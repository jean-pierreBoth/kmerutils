//! This module is mostly dedicated to doing simultanuously minhash and
//! related counting with compressed kmers.
//! 
//! This implementation of minhash is highly inspired by the finch module
//! but we need some genericity on types representing Kmer and Hash.  
//! We provide a generic module that relies on invertible hash
//! so we can avoid storing hashed kmer as hash value correspond bijectively to Kmers

//


// with the inspiration of the  finch module for NoHashHasher.
// perhaps use a bloomfilter instead of Hashmap. Especially if counts are not used in jaccard estimate.


// random generation numbers by decreasing speed XorShiftRng, ChaChaRng Isaac64Rng
// The individual items to store in the BinaryHeap




use crate::minhash::{MinHashDist,ItemHash};


use std::cmp::Ordering;

use std::collections::{BinaryHeap, HashMap};

use std::hash::{BuildHasherDefault, Hasher};

use std::mem;
use std::marker::PhantomData;


use log::trace;

use probminhash::invhash::*;

pub use crate::kmer::*;


/// Special Trait for hashed kmer with inversible hash.
// recall that PhantomData is zero sized

#[derive(Debug,Clone)]
pub struct InvHashedKmer<T:CompressedKmerT> {
    hash: ItemHash,
    t_marker: PhantomData<T>,
}



impl <T:CompressedKmerT> InvHashedKmer<T> {

    pub fn new(hash:ItemHash) -> Self {
        InvHashedKmer{hash:hash, t_marker:PhantomData,}
    }
    pub fn get_hash(&self) -> ItemHash {
        return self.hash;
    } // end of get impl
    
}  // end of impl InvHashedKmer





impl<T:CompressedKmerT> PartialEq for InvHashedKmer<T> {
    fn eq(&self, other: &InvHashedKmer<T>) -> bool {
        other.hash.eq(&self.hash)
    }
}



impl<T:CompressedKmerT> Eq for InvHashedKmer<T> {}


impl<T:CompressedKmerT> Ord for InvHashedKmer<T> {
    fn cmp(&self, other: &InvHashedKmer<T>) -> Ordering {
        self.hash.cmp(&other.hash)
    }
}


impl<T:CompressedKmerT> PartialOrd for InvHashedKmer<T> {
    fn partial_cmp(&self, other: &InvHashedKmer<T>) -> Option<Ordering> {
        Some(self.hash.cmp(&other.hash))
    }
}


//====================================================================================//



// size is 8 + 1 bytes !!
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct InvHashCountKmer<T:CompressedKmerT> {
    pub hashed: InvHashedKmer<T>,
    count: u8,
}

impl <T:CompressedKmerT> InvHashCountKmer<T> {
    pub fn new(hashed: InvHashedKmer<T>, count:u8) -> Self {
        InvHashCountKmer { hashed: hashed, count:count,}
    }
    pub fn get_count(&self) -> u8 {
        self.count
    }
}


//===================================================================================//





pub struct MinInvHashCountKmer<T:CompressedKmerT , H: Hasher+Default> {
    // heap to sort hashed values
    hashes: BinaryHeap<InvHashedKmer<T>>,
    // to keep a count how min hashed values encountered
    counts: HashMap<ItemHash, u8, BuildHasherDefault<H>>,
    total_count: u64,
    // number of hashed kmer to keep 
    size: usize,
}



impl <T:CompressedKmerT ,  H : Hasher+Default> MinInvHashCountKmer<T, H> {
    /// an allocator , size is capacity measured as  max number of hashed item
    /// store_kmer is to ask(or not) to keep the kmers hashed.
    /// as we use an invertible hasher for compressed kmers we do not need to keep track of kmers
    pub fn new(size: usize) -> Self {
        MinInvHashCountKmer {
            hashes: BinaryHeap::with_capacity(size + 1),
            counts: HashMap::with_capacity_and_hasher(size, BuildHasherDefault::<H>::default()),
            total_count: 0,
            size: size,
        }
    }  // end of new


    /// push an item in the sketching
    fn push(&mut self, item : &T) {
        //
        // hash with invertible hash 32 or 64 bit depending on size of T::Val
        let kmerval: T::Val = item.get_compressed_value();
        let new_hash = match mem::size_of::<T::Val>() {
            4 => {  let val_u = unsafe { mem::transmute_copy::<T::Val, u32>(&kmerval)};
                    int64_hash(val_u as u64)                                      
            },
            8 => {  let val_u = unsafe { mem::transmute_copy::<T::Val, u64>(&kmerval)};
                    int64_hash(val_u) as u64                                       
            },
            _ => panic!("bad size of kmer value"),
        };

        // do we insert
        let add_hash = match self.hashes.peek() {
            None => true,
            Some(old_max_hash) => (new_hash <= (*old_max_hash).hash) || (self.hashes.len() < self.size),
        };
        // if add_hash is true we must insert in hashes, 
        if add_hash {
            self.total_count += 1;
            if self.counts.contains_key(&new_hash) {
                // the item was already seen once.
                // let _lock = self.map_lock.lock().unwrap();
                let count = self.counts.entry(new_hash).or_insert(0u8);
                (*count) += 1;
            } else {
                // newhash is encountered for the first time
                self.hashes.push(InvHashedKmer {
                    hash: new_hash,
                    t_marker: PhantomData,
                });
                // 
                self.counts.insert(new_hash, 1u8);
                // just keep the number of minhash asked for
                if self.hashes.len() > self.size {
                    let hashitem = self.hashes.pop().unwrap();
                    let _old_count = self.counts.remove(&hashitem.hash).unwrap();
                }
            }
        } // end if add_hash        
    } // end push

    /// push a slice in the sketching
    pub fn sketch_kmer_slice(&mut self, to_sketch : &[T]) {
        trace!("sketching slice");
        to_sketch.into_iter().for_each(|x| self.push(x));
    } // end of sketch_slice


    /// returns a sorted vector of the sketch
    // We get a size 2 memory reduction with original minhash
    pub fn get_sketchcount(self) -> Vec<InvHashCountKmer<T> > {
        // this consumes the binary heap
        let mut vec = self.hashes.into_sorted_vec();
        
        let mut results = Vec::with_capacity(vec.len());
        for item in vec.drain(..) {
            let counts = *self.counts.get(&item.hash).unwrap();
            let counted_item = InvHashCountKmer {
                hashed: item,
                count: counts,
            };
            results.push(counted_item);
        }
        results
    }  // end of get_sketchcount

    
}  // end of impl MinInvHashCountKmer



/// compute different distances from sketch.
/// The arguments are supposed to come from get_sketchcount method that returns sorted (!!!) InvHashCountKmer
/// What do we do of counts? See ProbMinHash
pub fn minhash_distance<T:CompressedKmerT>(sketch1: &Vec<InvHashCountKmer<T> >, sketch2: &Vec<InvHashCountKmer<T> >) ->  MinHashDist {
    let mut i: usize = 0;
    let mut j: usize = 0;
    let mut common: u64 = 0;
    let mut total: u64 = 0;
    let sketch_size = sketch1.len();

    while i < sketch1.len() && j < sketch2.len() {
        if sketch1[i].hashed < sketch2[j].hashed {
            i += 1;
        } else if sketch2[j].hashed < sketch1[i].hashed {
            j += 1;
        } else {
            i += 1;
            j += 1;
            common += 1;
        }
        total += 1;
        if total >= sketch1.len() as u64 {
            break;
        }
    } // end while
    //
    // try to increase total up to asked sketch size
    //
    if total < sketch1.len() as u64 {
        // try to increase total.
        if i < sketch1.len() {
            total += (sketch1.len() - i) as u64;
        }
        if j < sketch1.len() {
            total += (sketch1.len() - j) as u64;
        }
        // now if ever total increase too much we truncate it
        if total > sketch_size as u64 {
            total = sketch_size as u64;
        }            
    }        
    //
    let containment: f64 = common as f64 / i as f64;
    let jaccard: f64 = common as f64 / total as f64;
    MinHashDist(containment, jaccard, common, total)
}  // end of minhash_distance


////////////////////////////////////////////////////////////////////////////////////////:


#[cfg(test)]
mod tests {
    use super::*;
    extern crate fnv;
    #[allow(unused_imports)]
    use self::fnv::FnvHasher; // extern fnv declared in test so we use self::fnv , if declared above we use super::fnv
    #[allow(unused_imports)]
    use crate::nohasher::NoHashHasher;
    

    use crate::kmergenerator::*;
    use crate::sequence::*;
    use env_logger;

    // initialize once log system for tests.
#[allow(dead_code)]
    fn init_log_test() {
        // we do not call unwrap as it is an error to init twice and possibly some one else initialized it
        let _ = env_logger::builder().is_test(true).try_init();
    }

    #[test]
    fn test_minhash_count_range_intersection_fnv() {
        init_log_test();
        // we construct 2 ranges [a..b] [c..d], with a<b, b < d, c<d sketch them and compute jaccard.
        // we should get something like max(b,c) - min(b,c)/ (b-a+d-c)
        //
        let str = String::from("TCAAAGGGAAACATTCAAAATCAGTATGCGCCCGTTCAGTTACGTATTGCTCTCGCTAATGAGATGGGCTGGGTACAGAG");
        let seq_bytes = str.as_bytes();
        let vkmer_a : Vec<Kmer16b32bit> = KmerGenerator::new(16).generate_kmer(&Sequence::new(&seq_bytes[0..80],2));
        let vkmer_b : Vec<Kmer16b32bit> = KmerGenerator::new(16).generate_kmer(&Sequence::new(&seq_bytes[60..],2));
        let mut minhash_a : MinInvHashCountKmer<Kmer16b32bit, FnvHasher>= MinInvHashCountKmer::new(5);
        let mut minhash_b : MinInvHashCountKmer<Kmer16b32bit, FnvHasher>= MinInvHashCountKmer::new(5);
        // now compute sketches
        println!("sketching a ");
        minhash_a.sketch_kmer_slice(&vkmer_a);
        println!("\n \n sketching b ");
        minhash_b.sketch_kmer_slice(&vkmer_b);
        let sketch_a = minhash_a.get_sketchcount();
        let sketch_b = minhash_b.get_sketchcount();
        // 
        let resdist = minhash_distance(&sketch_a, &sketch_b);
        trace!("distance minhash (contain, dist, common, total):  {}  {}   {}  {} ",
               resdist.0, resdist.1, resdist.2, resdist.3);
        assert!(resdist.2 > 0);
        //
    } // end of test_range_intersection


}  // end of mod test
