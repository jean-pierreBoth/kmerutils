//! key definitions and trait implementation for communicating with redis


//  TODO
//  SliceAnchorValueRedis and SliceAnchorKeyRedis MinhashKeyRedis should not have redis in their name
//  as the mechanism would be the same with any KV store RocksDB/TiKV
//  Only the trait implementation ToRedisArgs and FromRedisValue are specific to redis

extern crate redis;
use self::redis::{FromRedisValue,ToRedisArgs,RedisWrite};

use std::str;

use log::trace;


// REDIS on port 6379  dedicated to (my) genomic processing
// for key "properties"
//   "prop:user"                                    value string for user name
//   "prop:fn"                                      value string fastafile processed
//   "prop:fn:process"                              value string "anchor:minhash"
//   "prop:fn:process:nb_bases"                     value int
//   "prop:fn:process:slice_size"                   value int
//   "prop:fn:process:pos:readnum:slicepos"         values : list of (kmer_value, kmer_count):

// trait ToString implemented for u64 and  u8 enable easy generation of string for this last key value as it
// is implemented for T: fmt::Display + ?Sized and
// For the inverse parse uses FromStr (Cf rust documentation)
// String to u64   str.parse::<u64>()  


pub const FN_KEY          : &str  = "prop:fn";
pub const PROCESS_KEY     : &str  = "prop:fn:process";
pub const NB_BASES_KEY    : &str  = "prop:fn:process:bases";
pub const SLICE_SIZE_KEY  : &str  = "prop:fn:process:ssize";
pub const POS_KEY         : &str  = "prop:fn:process:readnum:slicepos";

// for inverse indexing


pub const MINHASH_1       : &str  = "prop:fn:process:minhash_1";
pub const MINHASH_2       : &str  = "prop:fn:process:minhash_2";

pub const SLICE_ANCHOR_KEY : &str = "prop:fn:process:ssize:bases:readnum:slicepos";





// internal structure we decode from Redis via FromRedisValue to get value part of SliceAnchor
// This enables getting the value for a given key.
//

#[derive(Debug)]
pub struct SliceAnchorValueRedis {
    pub hk_count: Vec<(u64,u8)>,
}

// beware each basic write_redis_args push a new vec<u8> coming fromm a call to into_bytes.
impl redis::ToRedisArgs for SliceAnchorValueRedis {
    fn write_redis_args<W:RedisWrite+?Sized> (&self, out: &mut W) {
        // minhash is a vector of ItemHash (u64) and count (u8)
        let mut key:Vec<u8> = Vec::new();
        let nb_kmer = self.hk_count.len();
        for i in 0..nb_kmer  {
            if i > 0 {
                key.push(b':');
            }
            key.extend_from_slice(self.hk_count[i].0.to_string().as_bytes());
            key.push(b',');
            key.extend_from_slice(self.hk_count[i].1.to_string().as_bytes());
            if i <  nb_kmer {
                key.push(b':');
            }
        }
        trace!("SliceAnchorValueRedis encoded sliceanchor for redis (kmer/count {:?}", &key);
        out.write_arg(&key); 
    } // end of write_redis_args
}  // end impl redis::ToRedisArgs for SliceAnchorValueRedis



// Given a redis Value, generate a vector of (invhashkmer, count)
// corresponding to (InvHashedKmer<T>, u8)
// We next have to go from u64 to InvHashedKmer<T> knowing the number of bases

impl FromRedisValue for SliceAnchorValueRedis {
    fn from_redis_value(v: &redis::Value) -> redis::RedisResult<Self> {
        let mut retvec = Vec::<(u64,u8)>::new();
        //
        match v {
            redis::Value::Data(ref bytes) => {
                // We must be in this case as we encode all values in one flat string in ToRedisArgs for SliceAnchorValueRedis.
                // We should not get a Value::Bulk
                // now we must decode vecu8 byte by byte. Inverse of write_redis_args
                let res = str::from_utf8(bytes);
                if res.is_err() {
                    Err(redis::RedisError::from((redis::ErrorKind::TypeError,"Not a String")))
                }
                else {
                    let str = res.unwrap();
                    let vcouple : Vec<&str> = str.split(':').collect(); // we get in v a vector of "kmer,count"
                    for couple in vcouple {
                        let vterms:Vec<&str> = couple.split(',').collect(); // we get in vterms[0] kmer, in vterm[1] count as strings
                        // convert from strings to u64 and u8
                        let resu64 = vterms[0].parse::<u64>();
                        let resu8 = vterms[1].parse::<u8>();
                        if resu64.is_err() || resu8.is_err() {
                            return Err(redis::RedisError::from((redis::ErrorKind::TypeError,"Cannot decode kmer or count")));
                        }
                        else { // everything is fine
                            let kmer_h = resu64.unwrap();
                            let count = resu8.unwrap();
                            retvec.push((kmer_h,count));
                            trace!("SliceAnchorValueRedis FromRedisValue pushing  (kmer/count {:?} {:?}", &kmer_h, &count);
                        }                   
                    } // end of for in couple
                    Ok(SliceAnchorValueRedis{hk_count: retvec})
                } // end else
            }, // case ref bytes
            // else return error 
            _ => Err(redis::RedisError::from((redis::ErrorKind::TypeError,"Not a Vec<u8>"))),        
        } // end match
    }  // end fn decode_minhash_from_redis_value(v: &redis::Value)
    
} // end impl FromRedisValue for SliceAnchorValueRedis 


// Structure for generating KEYS, implements ToRedisArgs

#[derive(Clone,Debug)]
pub struct SliceAnchorKeyRedis {
    pub filename:String,
    pub process:String,
    pub slice_size:u32,
    pub nb_bases:u16,
    pub readnum:u32,
    pub slicepos:u32
}




impl ToRedisArgs for SliceAnchorKeyRedis {
    // we concatenate all fields of SliceAnchorKeyRedis as string and return corresponding Vec<u8>
    fn write_redis_args<W:RedisWrite+?Sized>(&self, out: &mut W) {
        let mut key:Vec<u8> = Vec::new();
        // This suppose that filename is a valid utf-8 but filename are !?  CAVEAT
        key.extend_from_slice(self.filename.as_bytes());
        key.push(b':');
        key.extend_from_slice(self.process.as_bytes());
        key.push(b':');           
        // transfer other numerical values as [u8]
        key.extend_from_slice(self.slice_size.to_string().as_bytes());
        key.push(b':');           
        key.extend_from_slice(self.nb_bases.to_string().as_bytes());
        key.push(b':');           
        key.extend_from_slice(self.readnum.to_string().as_bytes());
        key.push(b':');           
        key.extend_from_slice(self.slicepos.to_string().as_bytes());
        //
        out.write_arg(&key);
    }  // end of write_redis_args
    
}  // end of impl ToRedisArgs for SliceAnchorKeyRedis


//=================================================================================
// For reverse request. Which read have a slice with a given minhash.
//=================================================================================


#[derive(Clone,Debug)]
pub struct MinhashKeyRedis {
    pub filename:String,
    pub process:String,
    pub slice_size:u32,
    pub nb_bases:u16,
    pub minhash_val:u64,
}



impl ToRedisArgs for  MinhashKeyRedis {
    // we concatenate all fields of SliceAnchorKeyRedis as string and return corresponding Vec<u8>
    fn write_redis_args<W:RedisWrite+?Sized>(&self, out: &mut W) {
        let mut key:Vec<u8> = Vec::new();
        // This suppose that filename is a valid utf-8 but filename are !?  CAVEAT
        key.extend_from_slice(self.filename.as_bytes());
        key.push(b':');
        key.extend_from_slice(self.process.as_bytes());
        key.push(b':');           
        // transfer other numerical values as [u8]
        key.extend_from_slice(self.slice_size.to_string().as_bytes());
        key.push(b':');           
        key.extend_from_slice(self.nb_bases.to_string().as_bytes());
        key.push(b':');           
        key.extend_from_slice(self.minhash_val.to_string().as_bytes());
        //
        out.write_arg(&key);
    }  // end of write_redis_args
    
}  // end of impl ToRedisArgs for MinhashKeyRedis
