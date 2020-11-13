//! This module contains computation of some statistics on base distribution, return times
//! between A/T pairs and C/G pairs.
// It gathers some heuristics on how to pairs blocks of read before going to sketching


use ::histogram::*;


use ::ndarray::*;


// rayon usage
use ::rayon::prelude::*;


// general use
use ::std::result;
use ::std::fs::OpenOptions;
use ::std::io::prelude::*;
//
use std::io;
use std::fmt;
use std::cmp;
use std::fs;
use std::mem;


// our includes

use crate::sequence::*;


// should be allocated via ndarray::Dim([sz as usize, 4]) so a sz by 4 matrix
// rust use row order as opposed to fortran which is row column

/// We keep in ReadBaseDistribution main statistics :
///  - readsizehisto  
///    an histogram of read length distribution
///   
///  - acgt_distribution  
///  For each  percentage and each base we store the fraction of reads in acgt_distribution.  
///  So  acgt_distribution[percent, cbase] = f means there is a fraction f of reads where cbase occurs at level percent%.   
///  base columns are in order a,c,g,t 
/// 

pub struct ReadBaseDistribution {
    pub readsizehisto : histogram::Histogram,
    /// upper value of histo
    upper_histo : usize,
    /// count number of values greater than upper , histo trackable value
    pub histo_out : usize,
    pub non_acgt : usize,
    /// array of dimensions 101,4.  101 beccause percentage can be in 0..100 included!
    pub acgt_distribution : Array2<f64>,
}



impl ReadBaseDistribution {
    /// readmaxsize is the maximum length represented in the histogram  
    /// prec : is the precision for histogram representation1,2 or 3
    pub fn new(readmaxsize : usize, prec : usize, d: ndarray::Ix2) -> ReadBaseDistribution {
        // TODO check that prec < log10(readmaxsize)
        // record histogram length from 1 to readmaxsize with slot of size readmaxsize/10**prec
        // lowest value arg in init must be >= 1
        assert!(prec >= 1, "precision for histogram construction should range >= 1");
        let histo = Histogram::configure().max_value(readmaxsize as u64).precision(prec as u32).build().unwrap();
        let array : Array2<f64> = Array2::zeros(d);
        ReadBaseDistribution{readsizehisto : histo, upper_histo : readmaxsize, histo_out : 0,
                             non_acgt : 0 as usize , acgt_distribution : array}
    }

    
    /// dump result in file of name name.
    /// first line is number of non acgt
    /// Then for each slot from 0 to 100 include , we get number of reads
    /// where a base (a,c,g,t) occurs at a given percentage.
    /// format is 4 column for bases in order a,c,g,t
    /// This file can be reloaded by Julia module BaseDistribution
    
    pub fn ascii_dump_acgt_distribution(&self, name : &String) -> result::Result<(), io::Error> {
        // open filename
        let fileres = OpenOptions::new().write(true).create(true).truncate(true).open(name);
        match fileres {
            Ok(mut file) => {
                //
                let (nbrow,_) = self.acgt_distribution.dim();
                for i in 0..nbrow {
                    write!(file, "{} {}  {}  {} \n", self.acgt_distribution[[i,0]] , self.acgt_distribution[[i,1]] ,
                           self.acgt_distribution[[i,2]] , self.acgt_distribution[[i,3]])?;  // ? get directly to io::Error !
                }
                file.flush()?;
                return Ok(());
            },
            Err(e) => {
                println!("could not open file {}", name);
                return Err(e);
            }
        } // end match
    } // end of ascii_dump_acgt_distribution

    /// dump nbslot values giving percentage of read in [ upper_histo * (i / nbslot) , upper_histo * (i+1)/nbslot]
    /// for i = 0..nbslot
    
    pub fn ascii_dump_readlen_distribution(&self, name : &str) -> result::Result<(), io::Error> {
        //
        let nb_entries = self.readsizehisto.entries() as usize;
        let maxlen = self.readsizehisto.maximum().unwrap() as usize;
        println!("ascii_dump_readlen_distribution nb_entries {}  maxlen {}", nb_entries, maxlen);
        let nbslot = nb_entries / 100;
        //
        let mut readsize : Vec<u64> = (0..(nbslot+1)).map(|_| 0u64).collect();
        for i in 0..(nbslot+1) {
            readsize[i] =  self.readsizehisto.percentile(100. * (i as f64/ nbslot as f64) as f64).unwrap();
        }
        let nb_points = 1000usize;
        let mut nb_read_vec = Vec::<usize>::with_capacity(nb_points);
        let mut abscisse = Vec::<usize>::with_capacity(nb_points);
        //
        let mut first_i = 0;
        let mut current_i = 0;
        for j in 0..nb_points {
            // find first i such that readsize[i] >= maxlen * (j/nbslot) 
            let threshold = (maxlen * j)/nb_points;
            while readsize[current_i]  < (threshold as u64)  && current_i < nbslot {
                current_i +=1;
            }
            if current_i < nbslot && current_i > first_i {
                let nb_in_slot = ((current_i - first_i) * nb_entries as usize) / nbslot; 
                nb_read_vec.push(nb_in_slot);
                abscisse.push(readsize[current_i] as usize);
            }
            first_i = current_i;
        }
        //
        let fileres = OpenOptions::new().write(true).create(true).truncate(true).open(name);
        match fileres {
            Ok(mut file) => {
                //
                for i in 0..nb_read_vec.len() {
                    write!(file, "{}  {} \n",abscisse[i] , nb_read_vec[i])?;  // ? get directly to io::Error !
                }
                file.flush()?;
                return Ok(());
            },
            Err(e) => {
                println!("could not open file {}", name);
                return Err(e);
            }
        } // end match
    } // end of ascii_readlen_distribution

    /// record read len in histogram keeping track of number of values outside histogram
    fn record_read_len(&mut self, sz : usize) -> () {
        let res = self.readsizehisto.increment(sz as u64);
        if !res.is_ok() {
            self.histo_out += 1;
        }
    } // end of record_read_len

} // end of implementation for BaseDistribution



/// computes base distribution statistics(in one thread) and dump result in file of name "bases.histo-1thread"
/// maxreadlen is maximum read length in histogram of read length distribution.
#[allow(dead_code)]
fn get_base_count(seq_array : &Vec<Sequence > , maxreadlen:usize) {
    //
    println!(" in get_base_count");
    let start_t = time::Instant::now();
    //
    let prec = (maxreadlen as f64).log10() as usize;
    let v_ref = &seq_array;
    let low = 0;
    let up = v_ref.len();
    let mut base_distribution = get_base_count_by_slice(v_ref.get(low..up).unwrap(), maxreadlen, prec).unwrap();       
    (*base_distribution).acgt_distribution *= 1./(seq_array.len() as f64);
    //
    let elapsed_t = start_t.elapsed().whole_seconds();
    println!(" elapsed time (s) in get_base_count {} ", elapsed_t);
    //
    base_distribution.ascii_dump_acgt_distribution(&String::from("bases.histo-1thread")).unwrap();
} // end of get_base_count



/// takes a slice on sequences and sendback base distribution on a channel
/// maxreadlen is the maximum size we record in histogram. It must be adapted to read length distribution
/// to get a significant historgram


fn get_base_count_by_slice(seq_array : &[Sequence] , maxreadlen : usize, prec : usize) -> result::Result< Box<ReadBaseDistribution> , ()> {
    //
    if cfg!(feature = "verbose_1") {
        println!(" in get_base_count_by_slice len : {} ", seq_array.len());
    }
    let start_t = time::Instant::now();
    // allocate a structure to store percentages, we want to cover an interval [0..100] so sz = 101
    let sz : usize = 101;

    let mut base_distribution : Box<ReadBaseDistribution>  =
        Box::new(ReadBaseDistribution::new(maxreadlen , prec, ndarray::Dim([sz as usize, 4])) );
    //
    let mut histo : Vec<u64> = vec![0,0,0,0];
    let mut nb_bad = 0;
    for i in 0..seq_array.len() {
        let decompressed = seq_array[i].decompress();
        base_distribution.record_read_len(decompressed.len());
        nb_bad = nb_bad + seq_array[i].base_count(&mut histo);
        for j in 0..4 {
            // compute percentage of each base in sequence
            let percent : usize = ((100 * histo[j]) as f64 / decompressed.len() as f64).round() as usize;
            (*base_distribution).acgt_distribution[ [percent ,j ] ] += 1.;
        }     
    }
    (*base_distribution).non_acgt = nb_bad as usize;
    //
    let elapsed_t = start_t.elapsed().whole_seconds();
    println!(" elapsed time (s) in get_base_count {} ", elapsed_t);
    //
    if cfg!(feature = "verbose_1") {
        println!(" out  get_base_count_par_slice nb_bad : {} ", nb_bad);
    }
    //
    Ok(base_distribution)
}





/// takes a slice on sequences and sendback base distribution on a channel
/// maxreadlen is the maximum size we record in histogram. It must be adapted to read length distribution
/// to get a significant historgram.  
/// Result is dumped in a file named : "bases.histo"
pub fn get_base_count_par (seq_array : &Vec<Sequence> ,  maxreadlen :usize, prec : usize) -> Option<Box<ReadBaseDistribution>> {
    //
    let nbthreads : usize = 2;
    let nb_cpus = num_cpus::get_physical();
    println!(" in get_base_count_par,  number of cpus found : {} " , nb_cpus);
    println!(" in get_base_count_par");
    let start_t = time::Instant::now();
    //
    // allocate nbthreads BaseDistribution. Do we need Cell or RefCell
    
    let v_ref = &seq_array;   // Vec n implemente pas Copy mais Ref oui!

    let distrib_collector : Vec<result::Result< Box<ReadBaseDistribution> , () >  > = (0..nbthreads).into_par_iter().map(|i|  {
        let low = (v_ref.len() / nbthreads)  * i;
        let up = (v_ref.len() / nbthreads)  * (i as usize + 1);
        println!(" in get_base_count_par {} {} ", low , up);
        let res: result::Result< Box<ReadBaseDistribution> , () > = get_base_count_by_slice(v_ref.get(low..up).unwrap(), maxreadlen , prec);
        res
    }).collect();

   
    // now we can reduce result summing all result in first 
    let sz = 101;    
    let mut base_distribution : Box<ReadBaseDistribution>  = Box::new(ReadBaseDistribution::new(maxreadlen , prec, ndarray::Dim([sz as usize, 4])) );
    
    for distrib in &distrib_collector {
        let ref_distrib = distrib.as_ref().unwrap();
        base_distribution.readsizehisto.merge(&ref_distrib.readsizehisto);
        // CAVEAT should check that upper_histo are coherent (the same) they are all initialized to maxreadlen
        base_distribution.upper_histo = maxreadlen;
        base_distribution.histo_out += ref_distrib.histo_out;
        base_distribution.non_acgt += ref_distrib.non_acgt;
        //   documentation suggest something like that should be possible
        base_distribution.acgt_distribution += &(ref_distrib.acgt_distribution);
    } // end of for on distrib_collector
    base_distribution.acgt_distribution *= 1./(seq_array.len() as f64);
    //
    let elapsed_t = start_t.elapsed().whole_seconds();
    println!(" elapsed time (s) in get_base_count {} ", elapsed_t);
    // dump
    base_distribution.ascii_dump_acgt_distribution(&String::from("bases.histo")).unwrap();
    //
    return Some(base_distribution);
} // end of get_base_count_par


//
//     Structures and Traits  for estimating Return Times for pair of bases
//     ========================================================
//




pub trait ReturnTimes {
    fn insert(delta : usize) -> bool;
    fn first() -> usize;  
}




// struct in Histogram or HdrHistogram are really two heavy for us as
// we will have one histogram per read.
// We will not collect statistics on more than 2^8 (=255) consecutive bases
// so u16 is sufficient for us. The structure ReturnTimes occupies
// INTER_TIMES_LENGTH (=12) u8 

/// struct ReturnTimes
/// This struct stores return times statistics for some conjugate pair AT (or CG) of bases by sequence.
/// This statistics is reverse complement invariant as we keep delta between A or T occurrence.
/// We can keep up to 2 ReturnTimesLR by read , one at the beginning of read one at end.



const RETURN_TIME_REC_MAGIC : [u8;1] = [0b11111111];            // 0xFF in hexa
const RETURN_TIME_FILE_MAGIC : [u8;4] = [0x00, 0xce, 0xab, 0xaf];

/// We count times between return up to 8, as probability that we have more than 8 consecutive C/G or A/T
/// in a read is roughly 1/2^8 and we sample on block size of length 250.
/// The number of searched bases can be computed by summing inter_times and adding 1. (slot/interval pb)
/// So the fraction of searched/total can be computed by (1+sum(inter_times))/window_size

const INTER_TIMES_LENGTH : usize = 8;

// chosen so that with 15% overlap on each side total block length is less than 256.
const DEFAULT_RET_BLOCK_SIZE : u8 = 190;


#[derive(Clone, Copy, Default)]
pub struct ReadBlockLocation {
    read_num : u32,
    /// read_pos is the num of block inside the read for which we sample data.
    /// So real pos in read_pos * size of BLOCK
    read_pos : u16,
}


impl fmt::Display for ReadBlockLocation {
    fn fmt(&self, dest : &mut fmt::Formatter) -> fmt::Result {
        write!(dest, "\n read position  ")?;
        write!(dest, "\n read: {}, pos : {} ", self.read_num, self.read_pos)
    } // end of fmt
}  // end impl fmt::Display


// size of struct is is 8 + 6 + 3 = 17 bytes
// So for short read of 101 bsases get a reduction of factor 5
// Note: total number of bases seen : 1 + sum(inter_times) so nb_seen just spares us the cost of summing inter_times


#[derive(Clone, Copy)]
pub struct ReturnTimesLR {
    inter_times : [u8;INTER_TIMES_LENGTH],
    /// position of beginning of stats
    pos : ReadBlockLocation,
    /// position of first occurrence of base pair in read after pos.read_pos!
    first : u8,
    /// number of bases seen of the searched pair.
    nb_seen : u8,
    /// number of bases scanned , less than 256 so fields first , above and nb_seen can be u8.
    /// But can be less than what user asked for if read is small or at end of a read.
    /// This field enables to get fraction of AT (or CG). 
    window_size : u8,
}  // end of struct ReturnTimes




impl Default for ReturnTimesLR {
    fn default() -> ReturnTimesLR {
        ReturnTimesLR {
            pos : ReadBlockLocation::default(),
            inter_times : [0u8;INTER_TIMES_LENGTH],
            first : 0,
            nb_seen : 0,
            window_size : DEFAULT_RET_BLOCK_SIZE,
        }
    }
}  // end of impl Default

impl fmt::Display for ReturnTimesLR {
    fn fmt(&self, dest : &mut fmt::Formatter) -> fmt::Result {
        self.pos.fmt(dest)?;
        write!(dest, "\n delta  count  ")?;
        for i in 0..INTER_TIMES_LENGTH {
            write!(dest, "\n {}       {}  " , i , self.inter_times[i])?;
        }
        write!(dest, "\n first : {},  nb_seen : {}, window_size : {}", self.first, self.nb_seen, self.window_size)
    } // end of fmt
}  // end impl fmt::Display



impl ReturnTimesLR {
    pub fn new(pos_arg: ReadBlockLocation, init_t : u8) -> ReturnTimesLR {
        ReturnTimesLR {
            pos: pos_arg,
            inter_times : [0u8;INTER_TIMES_LENGTH],
            first : init_t,
            nb_seen : 0,
            window_size : 0,
        }        
    }
    /// return size of internal buffer
    pub fn size(&self) -> usize {
        return self.inter_times.len();
    }
    /// get position of first. returns -1 if not seen 
    pub fn get_first(&self) -> i64 {
        if self.nb_seen > 0 {
            return self.first as i64;
        }
        else {
            return -1;
        }
    }
    
    /// returns number of searched bases found.
    pub fn get_nb_seen(&self) -> u8 {
        return self.nb_seen;
    }

    /// returns window_size
    pub fn get_window_size(&self) -> u8 {
        return self.window_size;        
    }
    
    /// insert a new delta (bases) between a couple (a | t)  or (c | g)
    /// every thing greate than INTER_TIMES_LENGTH-1 is stored in last slot.
    fn insert(&mut self, delta : usize) -> bool {
        let slot = delta.min(INTER_TIMES_LENGTH-1);
        self.inter_times[slot as usize] += 1;
        return true;
    } // end of fn insert

    /// reset all values to 0 if we want to reuse the structure.
    pub fn reset(&mut self, pos_arg: &ReadBlockLocation) -> () {
        self.pos = *pos_arg;
        self.first = 0;
        self.nb_seen = 0;
        self.window_size = 0;
        for i in 0..self.inter_times.len() {
            self.inter_times[i] = 0;
        }
    }  // end of reset
    
    /// fills in stat for a slice , always decompressed !
    /// we search returns for base_ret or conjugated 
    fn analyze_seq(&mut self, seq : &[u8],  lower_searched : u8, max_size : u8) -> std::result::Result<(), ()> {
        let mut lower_b;
        let mut last_seen : usize = 0;
        // 
        let imax = cmp::min(seq.len(), max_size as usize);
        self.window_size = imax as u8;
        // 
        for i in 0..imax {
            lower_b = get_ac_from_tg(seq[i]);
            if lower_b == lower_searched {
                if self.nb_seen == 0 {
                    self.first = i as u8;
                }
                else if self.nb_seen > 0 {
                    // we count number of base to skip to get another hit. (We spare one slot, so)
                    let delta = i - last_seen - 1;
                    self.insert(delta);
                }
                last_seen = i;               
                self.nb_seen += 1;
            } // we saw sthing searched
        }  // end of for on i
        //
        return Ok(());
    } // end of analyze_seq



    ///
    /// dump record in a binary file BufWriter<fs::File> but could be a cursor
    /// One record has size in bytes : 1(magic) + inter_times.len * u8 + 8 + (first , nb_seen, window_size)  i.e 24  bytes for
    /// inter_times.len == 12.
    ///
    pub fn dump_binary(&self, out : &mut dyn Write ) -> Result<usize, io::Error> {
        let mut nb_written = 0usize;
        // write magic at record beginning. Just one byte
        nb_written += out.write(&RETURN_TIME_REC_MAGIC)?;
        //
        // dump read position
        nb_written += out.write(unsafe { &mem::transmute::<u32, [u8;4]>(self.pos.read_num) } )?;
        nb_written += out.write(unsafe { &mem::transmute::<u16, [u8;2]>(self.pos.read_pos) } )?;
        //
        let out_slice = unsafe { ::std::slice::from_raw_parts(&self.inter_times[0] as *const u8,
                                                              self.inter_times.len()) };                     
        nb_written += out.write(out_slice)?;
        // write now 3 u16 in the following order first  and nb_seen and window_size
        // we need to dump real window_size as at end of read we can have less than default block size 
        nb_written += out.write(unsafe { &mem::transmute::<u8, [u8;1]>(self.first) } )?;
        nb_written += out.write(unsafe { &mem::transmute::<u8, [u8;1]>(self.nb_seen) } )?;
        nb_written += out.write(unsafe { &mem::transmute::<u8, [u8;1]>(self.window_size) } )?;
        //
        return Ok(nb_written);
    }
}  // end of implementation ReturnTimesLR


///
///  an iterator over ReturnTimesLR
///  
pub struct IterLR<'a> {
    rt: &'a ReturnTimesLR,
    /// current pos
    index: u8,
}



impl<'a> IterLR<'a> {
    pub fn new(rt_arg: &'a ReturnTimesLR) -> IterLR<'a> {
        IterLR {
            rt: rt_arg,
            index:0,
        } 
    }  // end of new
}  // end of implementation block



/// next returns an option containing if some a pair (delta, count)
impl<'a> Iterator for IterLR<'a> {
    type Item = (u8,u8);
    //
    fn next(&mut self) -> Option<(u8,u8)> {
        let limit = self.rt.inter_times.len() as u8;
        if self.index >= limit {
            return None;
        }
        else {
            // store current index so, we can increment index for next usage.
            let current = self.index;
            self.index += 1;
            return Some((current, self.rt.inter_times[current as usize]));
        }
    } // end of next
}


//
///  Return Times Writer
///
///  Collect raw sequences on the fly (before possible) compression
///  Internally the constructor :
///            1. opens a dump file to dump analyzed seq
///            2. Allocate a ReturnTimesLR analyze a sequence upon arrival, dump and reset ReturnTimesLR
///
///  Usage is : allocation of ReturnTimesWriter via new and then iterativley calls analyze
///          
//
//

#[allow(dead_code)]
pub struct ReturnTimesWriter {
    filename : String,
    return_t : ReturnTimesLR,
    bufw : io::BufWriter<fs::File>,
    /// lower base of each pair (AC/GT)
    a_or_c : u8,
    /// maximum window size to scan as asked by user. Can be less in reality if small reads are encountered
    max_window : u8,
}


impl ReturnTimesWriter {
    pub fn new(fname: String, searched_base : u8, max_size :u8) ->  result::Result<ReturnTimesWriter, io::Error> {
        let file = OpenOptions::new().write(true).create(true).truncate(true).open(&fname)?;
        let mut bufwrite = io::BufWriter::with_capacity(1_000_000_000, file);
        let rt = ReturnTimesLR::default();
        let lower_searched = get_ac_from_tg(searched_base);
        if lower_searched != b'A' && lower_searched != b'C' {
            println!(" bad base for return times analyze : {} ", lower_searched as char);
            return Err(io::Error::new(io::ErrorKind::Other, "bad base"));
        }
        //
        println!("allocating return times analyzer for base {} with window max {}", lower_searched as char, max_size);
        // now we initialize head of file, must dump magic and searched base.
        bufwrite.write(&RETURN_TIME_FILE_MAGIC)?;
        bufwrite.write(unsafe { &mem::transmute::<u8, [u8;1]>(max_size) } )?;
        bufwrite.write(unsafe { &mem::transmute::<u8, [u8;1]>(INTER_TIMES_LENGTH as u8) } )?;
        bufwrite.write(unsafe { &mem::transmute::<u8, [u8;1]>(searched_base) } )?;
        //
        Ok(ReturnTimesWriter {
            filename: fname,
            return_t: rt,
            bufw: bufwrite,
            a_or_c: lower_searched,
            max_window: max_size,
        })
    } // end of new
    //
    pub fn analyze(&mut self, read_num: usize, seq: &[u8]) -> result::Result<usize, io::Error> {
        // first we reset our stat
        let mut pos = ReadBlockLocation { read_num: read_num as u32, read_pos: 0};
        // loop inside seq.
        // Do we sample at beginning of segment?
        let mut nb_written = 0;
        let delta:usize = self.max_window as usize;
        let mut nb_sample = seq.len()/delta as usize;
        // ranges are semi open in Rust
        if seq.len() % delta > 0 {
            nb_sample += 1
        }
        for i in 0..nb_sample {
            let begin = i * delta;
            pos.read_pos = i as u16;
            self.return_t.reset(&pos);
            self.return_t.analyze_seq(&seq[begin..], self.a_or_c, self.max_window).unwrap();
            // now we dump, the record (recall it includes a magic)
            nb_written += self.return_t.dump_binary(&mut self.bufw)?;
        }
        return Ok(nb_written);
    }
}  // end of ReturnTimesWriter



/////////////////////////////////////////////////////////////////////////:

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_return_times_1() {
        let long_read = "ACATTGTACTTCGTTCAGTTACGTATTGCTTATTTTTCATGACCTCAGGGGTCGACGACATGTACCTACTGCACTCCAGC
CAGAACCTTAAAGACCCTGTATAAAAGCAAAAACCCATTGAGTCCAGCTAAGTGTATGGATATCAAAGCAAAATAATTGA
TTCAAACAGCTGCCGCATACCTGTTAACTCATTGTGTGTTATTCAGGTCTTCTTCAAACCTCTGCTTTTACAGCAATTGT
GGCCAATGTGAGACAGAATTTAGGCTAAAAACTTTGAGCAGTTGTCAATATTGTCATAGATAGAAGCCTATATGGCAATA
ATACGTCATAGCTGAAGCCAAACAAAGGCTATTAGTTCAATTTAAGACATGGTAGAACTCAAGTCTGGGGTGAATCAAGA
CTTTGAAGTACTTTATTCAACAGACAAGTTGCATTGCTTATCTCAGTATTATCTTTATAAAATCTTTCTTTTCTCTGGTT
TATTTCACATTTCTTGAGCTTTCATTATATCCACTACTTGGTTTGTGTATTAGATGAATTTTCCACTGAACTCATTGCTG
ATATCAAACTACAGCATTCATCATACCTTTCTGCCTGTTCTACGCACCACCTGCACAGAATTTCACTTCTTGACACCCTT
GTTTAAAGCCAGCAGCAACTGTGAAAACTTCAGCCAGGTAAATTTAAAGAGTTTAATTGAGCAATGAACAATTCACGAAT
CGGGCAGTTCCCAGAATCACAGCAGATTCACAGAGACTCCAGGGTGCCTCATGGTCAGAAACAATTATATGAACGAAACC
GTGATGTACAGGAACGGAACAGCTAGATTGATTGGTTACAGGTTGGTGTTGCCTTATTTGAACACAGTTTGAACGCTCAG
CAGTCTGATTGGTTAAGGAAGTATGACTCCTGGGTTTGGCCAACACTCAGCTGATACAGGCACATACTCTTTTTAGAGAC
GGTGTCTTGCTGTTGCCCAGGGTGGGTATTGGCATGATCTCCAGCTCACTGCAACCTCCGCCTCCAGGTTCAATGATTCT
CCTACCTCAACCTCCAGGTAGCTGGGATTACAGACTGGCGCCCTGCCTTAATAGCTGTACTTTTGATAGAGGCGAGGGTT
TCACCATGGTCAGGCTGATTCTCAAACTCCTGACTGGCCTCAGTCTGCCT";
        // we directly get a str
        // get rid of '\n' !!!! as we got the string from a (ascii fastq file)
        println!("long read length : {} ", long_read.len());
        let read_bytes = String::from(long_read).into_bytes();
        let mut filtered_bytes = Vec::new();
        for i in 0..read_bytes.len() {
            if read_bytes[i] != b'\n' {
                filtered_bytes.push(read_bytes[i]);
            }
        }
        // check we got rid of b'\n' get 1170 base string
        println!("filtered length : {} ", filtered_bytes.len());
        assert!(filtered_bytes.len() == read_bytes.len() - 14);

        let searched_base = b'A';
        let max_size = 255;
        let  mut rt_writer = ReturnTimesWriter::new(String::from("test_return_time_1.bin"), searched_base, max_size).unwrap();
        let res = rt_writer.analyze(0, &filtered_bytes);
        match res {
            Ok(nb_written) => {
                println!(" nb bytes written : {} ", nb_written);
            },
            Err(_) => {
                println!("an io error occurred");
            }
        }
        println!("rt_writer.return_t.get_first() : {}", rt_writer.return_t.get_first());
        assert!(rt_writer.return_t.get_first() == 2);
    }  // end of test_return_times_1


    
    #[test]
    fn test_return_times_2() {
        let long_read = "
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
TTCAAACAGCTGCCGCATACCTGTTAACTCATTGTGTGTTATTCAGGTCTTCTTCAAACCTCTGCTTTTACAGCAATTGT
GGCCAATGTGAGACAGAATTTAGGCTAAAAACTTTGAGCAGTTGTCAATATTGTCATAGATAGAAGCCTATATGGCAATA
ATTTCACATTTCTTGAGCCCCCATTATATCCACTACTTGGTTTGTGTATTAGATGAATTTTCCACTGAACTCATTGCTGA
CCTACCTCAACCTCCAGGTAGCTGGGATTACAGACTGGCGCCCTGCCTTAATAGCTGTACTTTTGATAGAGGCGAGGGTT
TCACCATGGTCAGGCTGATTCTCAAACTCCTGACTGGCCTCAGTCTGCCT";
        // we directly get a str
        // get rid of '\n' !!!! as we got the string from a (ascii fastq file)
        println!("long read length : {} ", long_read.len());
        let read_bytes = String::from(long_read).into_bytes();
        let mut filtered_bytes = Vec::new();
        for i in 0..read_bytes.len() {
            if read_bytes[i] != b'\n' {
                filtered_bytes.push(read_bytes[i]);
            }
        }
        // check we got rid of b'\n'
        println!("filtered length : {} ", filtered_bytes.len());
        assert!(filtered_bytes.len() == read_bytes.len() - 6);
        //    
        let searched_base = b'A';
        let max_size = 255;
        let  mut rt_writer = ReturnTimesWriter::new(String::from("test_return_time_2.bin"), searched_base, max_size as u8).unwrap();
        let res = rt_writer.analyze(0, &filtered_bytes);
        match res {
            Ok(nb_written) => {
                println!(" nb bytes written : {} ", nb_written);
            },
            Err(_) => {
                println!("an io error occurred");
            }
        }
        println!(" first occurence : {} ", rt_writer.return_t.get_first());
        // we analyze every block beginning at 200*i base. 3 blocks begins at 600 (0 indexation)
        // first occurence in last ( 3th)  block of 200 which begins at line 7 char 40. We get a T at 43
        assert!(rt_writer.return_t.get_first() == 0);
        // following assert depends on block_size fed up with it
//        assert!(rt_writer.return_t.inter_times[INTER_TIMES_LENGTH-1] == 0);
        //
        //
    }  // end of test_return_times_2

    
    //  Now we test iterator
    #[test]
    fn test_return_times_3() {
        let read = "GTTTAAAGCCAGCAGCAACTGTGAAAACTTCAGCCAGGTAAATTTAAAGAGTTTAATTGAGCAATGAACAATTCACGAA";
        //
        let searched_base = b'A';
        let max_size = 255;
        {   // we dump statistics to be able to check iterator in the following
            let  mut rt_writer = ReturnTimesWriter::new(String::from("test_return_time_3.bin"), searched_base, max_size).unwrap();
            let res = rt_writer.analyze(0, &String::from(read).as_bytes());
            match res {
                Ok(nb_written) => {
                    println!(" nb bytes written : {} ", nb_written);
                },
                Err(_) => {
                    println!("an io error occurred");
                }
            }
        }
        let mut rt = ReturnTimesLR::default();
        assert!(rt.inter_times.len() == INTER_TIMES_LENGTH);
        rt.analyze_seq(&String::from(read).as_bytes(), searched_base, max_size).unwrap();
        let rep_delta = vec![32, 11, 5, 2, 0];
        let mut slot = 0;
        //
        println!("\n test iterator");
        let mut iter = IterLR::new(&rt);
        while let Some(item) = iter.next() {
            if slot <  rep_delta.len() {
                println!("delta, count = {} {} ", item.0, item.1);
                assert!(item.1 == rep_delta[slot]);
                slot += 1;
            }
        }
        println!("slot  = {} ", slot);
        assert!(slot == 5);
        //
        println!("dump ReturnTimesLR : {}", rt);
        // the following assert is just to get all writes on stdout !! if needed
//        assert!(slot == 6);
    }  // end of test_return_times_3
        

} // end module test
