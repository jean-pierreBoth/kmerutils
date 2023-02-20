//! This module describe structure used for Quality description, storing and loading
//!
//! As we want to spare memory we remap quality to a scale 0..8 (8 excluded) of 8 values.  
//! So waveletMatrix can encode values on 3 bits reducing to 0.56 memory load while
//! keeping all functionalies (access, rank and so on). Loading times get a factor 2 increase...
//!
//! qual 0x25 corresponds to p = 0.67. We do not expect quality under that threshold
//! if it occurs it will be remmapped to 0.  
//! qual 0x37 corresponds to p= 0.006. 
//! From 0x37 above it will be remmapped to 7
//!remap qualities to 1:6 qualities in  (0x25, 0x37)



use ::std::path::Path;
use ::std::cmp;

use ::wavelet_matrix::WaveletMatrix;

/// convert quality byte to a probablility with given threshold
pub fn quality_to_proba(q:u8, qmin:u8) -> f64 {
    10_f64.powf((qmin-q) as f64/10.0_f64)
}


/// As we want to spare memory we remap quality to a scale 0..8 (8 excluded) of 8 values
/// So waveletMatrix can encode values on 3 bits reducing to 0.56 memory load while
/// keeping all functionalies (access, rank and so on). Loading times get a factor 2 increase...
///
/// qual 0x25 corresponds to p = 0.67. We do not expect quality under that threshold
/// if it occurs it will be remmapped to 0
/// qual 0x37 corresponds to p= 0.006. 
/// From 0x37 above it will be remmapped to 7
/// remap qualities to 1:6 qualities in  (0x25, 0x37)
///

#[inline]
fn remap_quality8(q:u8) -> u8 {
    if q > 0x37 {
        return 7
    }
    else if q < 0x25 {
        return 0
    }
    else {
        let nqf = (cmp::min(q, 0x37) - 0x25) as f32 * 6.0_f32 / 18.0_f32;
        1 + nqf.floor() as u8
    }
}
/// a quality record can be either in raw data Vec\<u8\> or in WaveletMatrix

pub enum QualityMode {
    /// represent u8 encoding of quality
    Raw,
    /// represent WaveletMatrix encoding
    WM,
}


/// a sequence provides 2 methods :
///    - get the quality of a base given its rank.
///    - returns the representation mode

pub trait QSequence {
//    fn get_quality(&self, usize) -> Result<u8, u32>;
    fn get_mode(&self) -> QualityMode;
    fn get_read_num(&self) -> usize;
    fn len(&self) -> usize;
}




//  ========================================================================================
//        QSequenceWM
//
// As we will not always store all Qsequences it is useful to keep track
// associated read number
// We will possibly need a hash from num read to Qsequence
//
//  ========================================================================================


/// Wavelet matrix representation of a sequence (after 8bit mapping)
pub struct QSequenceWM {
    read_num : usize,
    pub qseq : WaveletMatrix,
}

impl QSequence for QSequenceWM {
    /// returns representation mode
    fn get_mode(&self) -> QualityMode {
        return QualityMode::WM
    }
    ///
    fn get_read_num(&self) -> usize {
        return self.read_num
    }
    ///
    fn len(&self) -> usize {
        self.qseq.len()
    }
}  // end of impl for QSequenceWM




impl QSequenceWM {
    pub fn new(read_n : usize, qv : &[u8]) -> QSequenceWM {
        let remapped : Vec<u64> = qv.iter().map(|q| remap_quality8(*q) as u64).collect();
        let qseqt = WaveletMatrix::new(&remapped);
         QSequenceWM { read_num: read_n, qseq: qseqt }
    }  // end of new

    /// get back to a raw sequence ... (This will run on a quality server or at least in a thread)
    /// The qualities are still mapped to a byte!
    pub fn decompress(&self) -> QSequenceRaw {
        let len = self.qseq.len();
        let mut vq = Vec::<u8>::with_capacity(len);
        for i in 0..len {
            vq.push(self.qseq.lookup(i) as u8);
        }
        QSequenceRaw{ read_num : self.read_num, qseq : vq}
    }

    /// This function returns the number of bits used in storage quality. (log2 of the max value)
    pub fn bit_len(&self) -> u8 {
        self.qseq.bit_len()
    }
}  // end of impl for QSequenceWM


//  ========================================================================================
//        QSequenceRaw
//  ========================================================================================


/// a quality sequence under raw representation
pub struct QSequenceRaw {
    pub read_num : usize,
    /// the remapped qualities in 0..7
    pub qseq : Vec<u8>,
}



impl QSequence for QSequenceRaw {
    ///
    fn get_mode(&self) -> QualityMode {
        return QualityMode::Raw
    }
    ///
    fn get_read_num(&self) -> usize {
        return self.read_num
    }
    ///
    fn len(&self) -> usize {
        self.qseq.len()
    }

}  // end of impl QSequence for QSequenceWM


impl QSequenceRaw {
    /// compress in format WM.
    pub fn compress_wm(&self) -> QSequenceWM {
        QSequenceWM::new(self.read_num, &(self.qseq))
    }

} // end of implementation block QSequenceRaw




pub fn load_quality_wm(filename: &String) ->  Result<Vec<QSequenceWM>, &'static str > {
    //
    println!("quality loading with needletail file : {} ", filename);
    //
    let path = Path::new(&filename);
    let f_info_res = path.metadata();
    match f_info_res {
        Ok(_meta) => (),
        Err(_e) => {
            println!("file does not exist: {:?}", filename);
            return Err("file does not exist");
        },
    }
    //
    let default_len = 1_000_000;
    let mut seq_array : Vec<QSequenceWM > = Vec::with_capacity(default_len);
    let mut n_qual = 0;
    let mut n_read = 0;
    let mut nb_bad_read = 0;    // reads with quality missing

    let mut reader = needletail::parse_fastx_file(&path).expect("expecting valid filename");
    while let Some(record) = reader.next() {
        n_read = n_read+1;
        //
        let seqrec = record.expect("invalid record");
        if let Some(qual) = seqrec.qual() { // match a reference
            n_qual += qual.len();
            // remap quality seq
            let newseq = QSequenceWM::new(n_read, qual);
            seq_array.push(newseq);
        }
        else {
            nb_bad_read += 1;
        }
        if seq_array.capacity() <= seq_array.len() + 100  {
            let old_len = seq_array.len() as f64;
            seq_array.reserve((old_len * 1.5) as usize);
        }
        if seq_array.len() % 200000 == 0 {
            println!(" nb rec read = {} ", seq_array.len());
        }
        //
    }
    //
    println!(" shrinking");
    seq_array.shrink_to_fit();
    println!(" shrinked ");
    //
    println!(" nb rec loaded = {} ", seq_array.len());
    println!("nb_qual {:?}", n_qual);
    println!("nb_read without quality {:?}", nb_bad_read);
    //
    return Ok(seq_array);
}  // end of load_quality
