//! computes anchors in reads via inversible minhash and stores them in redis.




use ::std::rc::Rc;
use ::std::path::Path;
use ::std::fmt;

use log::debug;
use log::trace;
use log::error;

use redis::{Commands};

use std::str;
use std::cmp;

use fnv::FnvHasher;
use std::hash::{Hash};

pub use crate::redisbase::*;

use crate::base::{sequence::*, kmer::*, kmergenerator::*};

use crate::hashed::*;

use crate::minhash::*;


/// General parameters of anchor generation

#[derive(Default,Clone)]
pub struct AnchorsGeneratorParameters {
    ///
    fasta_name:String,
    /// window size for computing minhash
    window: u32,
    /// kmer numbers per windows
    nbkmer:u32,
    /// kmer size
    kmer_size:u16,
    /// overlap between windows
    overlap:u32,
}

impl AnchorsGeneratorParameters {
    pub fn new(fasta_name: String, window: u32, nbkmer:u32, kmer_size:u16, overlap:u32) -> AnchorsGeneratorParameters {
        AnchorsGeneratorParameters { fasta_name, window, nbkmer, kmer_size, overlap}
    }
    ///
    pub fn get_fasta_name(&self) -> String { self.fasta_name.clone()} 
    ///
    pub fn get_window(&self) -> u32 { self.window }
    ///
    pub fn get_nbkmer(&self) -> u32 { self.nbkmer }
    ///
    pub fn get_kmer_size(&self) -> u16 { self.kmer_size}
    ///
    pub fn get_overlap(&self) -> u32 { self.overlap}
} // end impl AnchorsGeneratorParameters


impl fmt::Display for AnchorsGeneratorParameters {
    fn fmt(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
        let mut _res = write!(formatter, "fastq file : {}", self.fasta_name);
        _res = write!(formatter, "slice size : {} ", self.window);
        _res = write!(formatter, "nb kmer to keep : {} ", self.nbkmer);
        write!(formatter, "overlap : {}", self.overlap)
    }
}




//===================================================================================
//   SliceAnchor
//
//  The first 3 fields gives key fields for redis, minhash is value.
//====================================================================================

    
/// slice anchor
#[derive(Clone)]
pub struct SliceAnchor<T:CompressedKmerT> {
    params: Rc<AnchorsGeneratorParameters>,
    //
    readnum:u32,
    // beginning of slice in read
    slicepos:u32,
    //
    minhash:Vec<InvHashCount<T>>
}


// readnum, slicepos will be keys in DB.
impl <T:CompressedKmerT> SliceAnchor<T> {
    pub fn new(params : &Rc<AnchorsGeneratorParameters> , readnum:u32, slicepos:u32, minhash:Vec<InvHashCount<T>>) -> Self {
        SliceAnchor {
            params:params.clone(),
            readnum:readnum,
            slicepos:slicepos,
            minhash: minhash,
       }
    } // end of new

    // get a slice identy key for redis
    fn get_key_for_redis(&self) -> SliceAnchorKeyRedis {
        SliceAnchorKeyRedis {
            filename: self.params.get_fasta_name(),
            process: String::from("invhash"),
            slice_size: self.params.get_window(),
            nb_bases: self.params.get_kmer_size(),
            readnum: self.readnum,
            slicepos: self.slicepos,            
        }        
    }  // end of get_key_for_redis


    // 
    fn get_value_for_redis(&self) -> SliceAnchorValueRedis {
        let nb_kmer = self.minhash.len();
        let mut anchor_value = Vec::<(u64,u8)>::with_capacity(nb_kmer);
        for km in &self.minhash {
            let couple = (km.hashed.get_hash() as u64, km.get_count());
            anchor_value.push(couple);
        }
        SliceAnchorValueRedis{hk_count:anchor_value}
    }  // end of get_value_for_redis


    
    // get MinhashKeyRedis for inverse indexing
    fn get_minhash_key_for_redis(&self) -> MinhashKeyRedis {
        let min_hash_val = self.minhash[0].hashed.get_hash();
        MinhashKeyRedis {
            filename: self.params.get_fasta_name(),
            process: String::from("invhash"),
            slice_size: self.params.get_window(),
            nb_bases: self.params.get_kmer_size(),
            minhash_val:min_hash_val as u64,
        }
    } // end of get_minhash_key_for_redis

    ///
    /// insert in redis base
    ///
    // CAVEAT We have a 2 step conversion
    // TO DO try to get some real information on what happened
    //
    fn redis_dump(&self, con: &mut redis::Connection) -> std::result::Result<u64, String> {
        // generate SliceAnchorKeyRedis        
        let key:SliceAnchorKeyRedis = self.get_key_for_redis();
        let inv_value = key.clone();
        // generate redis::Value with corresponding SliceAnchorKeyRedis obtained by into and then with Trait ToRedisArgs
        let value:SliceAnchorValueRedis = self.get_value_for_redis();
        // pass the hash command , key statisfies the  Trait ToRedisArgs, we can do either with redis::cmd
        // or use Command Trait on Connection
        // logging
        trace!("dumping key (read/slicepos {:?} {:?}", &key, &value);
        let res = con.hset::<&'static str , SliceAnchorKeyRedis, SliceAnchorValueRedis, redis::Value> (SLICE_ANCHOR_KEY, key, value);
        //         let res = redis::cmd("HSET").arg(key).arg(value).execute(con);
        if !res.is_ok() {
            return Err(String::from("error"));
        }
        //
        // inverse indexing from min of minhash to readnum:slicepos identity.
        //
        let inv_key = self.get_minhash_key_for_redis();
        // logging
        trace!("dumping inverse key (read/slicepos {:?} {:?}", &inv_key, &inv_value);
        let res = con.hset::<&'static str , MinhashKeyRedis , SliceAnchorKeyRedis , redis::Value> (MINHASH_1, inv_key, inv_value);        
        // decode result
        match res {
            Ok(redis::Value::Okay)        => Ok(1),
            _                             => Err(String::from("error")),
        }
    } // end of dump_in_redis_base

}  // end impl SliceAnchor<T>





// a function togenerate a SliceAnchor<T> from SliceAnchorKeyRedis sent to redis  and SliceAnchorValueRedis got from request.

pub fn create_sliceanchor_from_redis_key_value<T:CompressedKmerT+Hash>(params : &Rc<AnchorsGeneratorParameters>,
                                                                  key: &SliceAnchorKeyRedis,
                                                                  value: &SliceAnchorValueRedis) -> SliceAnchor<T> {
    //  allocate minhash
    let mut minhash:Vec<InvHashCount<T>> = Vec::with_capacity(value.hk_count.len());
    for i in 0..value.hk_count.len() {
        let invhash = InvHashedItem::<T>::new(value.hk_count[i].0);
        let v = InvHashCount::new(invhash, value.hk_count[i].1);
        minhash.push(v);
    }
    SliceAnchor::new(params, key.readnum, key.slicepos, minhash)
}

             
//    values to generate corresponds to
//   "prop:fn:process:pos:readnum:slicepos"         values : list of kmer_value, kmer_count spearated by :



//===============================================================================================

   
pub fn gen_anchor_mininvhash<T:CompressedKmerT>(params: &Rc<AnchorsGeneratorParameters>,
                                                    readnum:usize, slicepos:usize,
                                                    kmer_generator: & dyn KmerGenerationPattern<T> ,
                                                    seq: &Sequence) -> Option<SliceAnchor<T> >  {
    //  
    let nb_bases = params.get_kmer_size();
    let slice_size = params.get_window();
    let nbkmer = params.get_nbkmer() as usize;
    //
    if nb_bases <= 32 {
        let beg = slicepos as usize;
        let end = cmp::min(beg + slice_size as usize, seq.size()-1);
        let kmers : Vec<T>  = kmer_generator.generate_kmer_pattern_in_range(seq, beg, end);
        let mut minhash_a : MinInvHashCountKmer<T, FnvHasher>= MinInvHashCountKmer::new(nbkmer);
        minhash_a.sketch_kmer_slice(&kmers);
        // returns a Vec<InvHashCount<T> > ... type inference let us forget what we talk about.
        let hashcount = minhash_a.get_sketchcount();
        return Some(SliceAnchor::<T>::new(&params, readnum as u32, slicepos as u32, hashcount));
    }
    else {
        error!("gen_anchor_mininvhash > 32");
        return None;
    }
} // end of gen_anchor_mininvhash_kmer64bit



//=====================================================================================

///
/// gathers anchors for a read
///
pub struct ReadAnchors<T:CompressedKmerT> {
    // CAVEAT absolute num of read or number without non actg
    readnum:usize,
    slice_params:Rc<AnchorsGeneratorParameters>,
    anchors:Vec<SliceAnchor<T>>,
}



impl <T:CompressedKmerT> ReadAnchors<T> {
    pub fn new(slice_parameters: &Rc<AnchorsGeneratorParameters>) -> Self {
        ReadAnchors {
            readnum:0,
            slice_params:slice_parameters.clone(),
            anchors:Vec::new(),
        }
    }

    // careful to be coherent here between seq and readnum  CAVEAT
    fn initialize_from_sequence(&mut self, kmer_generator : & dyn KmerGenerationPattern<T>,  seq: &Sequence, readnum:usize) {
        //
        debug!("anchoring read {}", readnum);
        println!("using anchors parameters : {} {} ", &self.slice_params.get_window(), &self.slice_params.get_overlap());
        assert!(&self.slice_params.get_window() > &self.slice_params.get_overlap());
        assert!(self.slice_params.get_window() > 0);
        //
        self.readnum = readnum;
        //
        let seqlen = seq.size() as u32;
        debug!("anchoring sequence of length {}", seqlen);
        // get an idea of number of anchors we will have and allocate
        self.anchors = Vec::new();
        self.anchors.reserve(1 + ((seqlen/self.slice_params.get_window()) as usize));
        //
        let mut beg = 0;
        let mut anchor_opt;
        while beg < seqlen {
            if beg == 0 {
                debug!("anchoring beg {}", beg);
                anchor_opt = gen_anchor_mininvhash(&self.slice_params, self.readnum, beg as usize, kmer_generator, seq);
                beg = beg + self.slice_params.get_window() - self.slice_params.get_overlap();
            }
            else {
                debug!("anchoring beg {}", beg);
                anchor_opt = gen_anchor_mininvhash(&self.slice_params, self.readnum, beg as usize, kmer_generator, seq);            
                beg = beg + self.slice_params.get_window() - self.slice_params.get_overlap();
            }
            match anchor_opt {
                Some(anchor) =>  self.anchors.push(anchor),
                None => panic!("anchor creation failed"),
            }                
        }  // end while
        //
        self.anchors.shrink_to_fit();
        //
        trace!("nb anchors generated for read {}", self.anchors.len());
    } // end of initialize_from_sequence


    pub fn get_nb_slice(&self) -> usize {
        return self.anchors.len();
    }

    // keys is sequence of bytes consisting of readnum+slicepos
    // values is list of kmers values (without nb base encoding)
    // return the number of anchors dumped
    pub fn redis_dump(&self, db: &mut redis::Connection) -> usize {
        let mut nb_anchor_dumped = 0;
        //
        for anchor in &self.anchors {
            let res = anchor.redis_dump(db);
            match res {
                Ok(_) => { nb_anchor_dumped += 1;},
                Err(_) => {                    
                    error!("error occurred in redis dump");                   
                },
            }  // end match
        }  // end of for
        return nb_anchor_dumped;
    } // end of redis_dump
    
} // end of impl ReadAnchors<T>


//==============================================================================================



///
/// gathers anchors for a Fasta/q file
///


pub struct FastaAnchors<T:CompressedKmerT> {
    slice_params:Rc<AnchorsGeneratorParameters>,
    anchors:Vec<ReadAnchors<T>>,
    _redis_db:Option<String>,
    con:Option<redis::Connection>,
    store_anchor:bool,
}


impl  <T:CompressedKmerT> FastaAnchors<T> {
    /// redis_addr  : typically "redis://127.0.0.1/6379"
    pub fn new (slice_parameters:Rc<AnchorsGeneratorParameters>, store_anchor:bool, redis_addr:Option<String>) -> Self {
        let mut con = None;
        match  redis_addr {
            Some(ref name) => {
                let res_client = redis::Client::open(name.as_str());
                if res_client.is_ok() {
                    let res_con = res_client.unwrap().get_connection();
                    if res_con.is_ok() {
                        con = Some(res_con.unwrap());
                    }
                    else {
                        panic!("cannot get redis connection");
                    }
                }
            },
            None => (),
        }
        FastaAnchors {
            slice_params:Rc::clone(&slice_parameters),
            anchors:Vec::new(),
            _redis_db:redis_addr,
            con:con,
            store_anchor:store_anchor,
        }        
    } // end of new


    // on fly anchors computation, do not need any storage. Process sequence after sequence and go to redis.
    pub fn anchor_computation(& mut self) -> Result<usize, & 'static str> where KmerGenerator<T>: KmerGenerationPattern<T> {
        let fasta_name = self.slice_params.get_fasta_name();
        let path = Path::new(&fasta_name);
        let f_info_res = path.metadata();
        let filesize:u64;
        match f_info_res {
            Ok(meta) => {
                filesize = meta.len();           
                println!("file size: {:?}", filesize);
            },
            Err(_e) => {
                println!("file does not exist: {:?}", path);
                return Err("file does not exist");
            },
        }
        // Default : 2 bits / base
        let nb_bits:u8 = 2;
        let nb_bases = self.slice_params.get_kmer_size();
        let mut num_read = 0;
        let start_t = time::Instant::now();
        let mut reader = needletail::parse_fastx_file(&path).expect("expecting valid filename");
        while let Some(record) = reader.next() {
            let seqrec = record.expect("invalid record");
            let nb_bad = count_non_acgt(&seqrec.seq());
            if nb_bad == 0 {
                let newseq = Sequence::new(&seqrec.seq(), nb_bits);
                // get rid of non actg read
                num_read += 1;
                let kmer_generator = KmerGenerator::<T>::new(nb_bases as u8);
                let mut read_anchor = ReadAnchors::<T>::new(&self.slice_params);  // pass a ref to Rc which will be cloned.
                read_anchor.initialize_from_sequence(&kmer_generator, &newseq, num_read);
                if let Some(ref mut real_con) = self.con {
                    read_anchor.redis_dump(real_con);
                }
                if self.store_anchor {
                    self.anchors.push(read_anchor);
                }
            }
        }
        //
        let elapsed_t = start_t.elapsed().whole_seconds();
        println!(" elapsed time (s) in anchor computation/dump {} ", elapsed_t);
        //
        let _res = self.redis_bgrewriteaof();
        //
        Ok(1)
    }  // end of on_fly_computation

    pub fn get_nb_anchors(&self) -> usize {
        let mut nb_anchors = 0;
        for v in &(self.anchors) {
            nb_anchors += v.get_nb_slice();
        }
        nb_anchors
    }

    
    fn redis_bgrewriteaof(&mut self) -> usize {
        if self.con.is_some() {
            redis::cmd("BGREWRITEAOF").execute(self.con.as_mut().unwrap());
            return 1;
        }
        else {
            return 0;
        }
    }

    
}  // end of impl for FastaAnchors


 
