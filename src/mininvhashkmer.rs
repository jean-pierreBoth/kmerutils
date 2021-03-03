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




use crate::minhash::{MinHashDist};

use crate::hashed::*;


use std::collections::{BinaryHeap, HashMap};

use std::hash::{BuildHasherDefault, Hasher};

use std::marker::PhantomData;


use log::trace;

use probminhash::invhash::*;

pub use crate::base::kmer::*;



//===================================================================================//






////////////////////////////////////////////////////////////////////////////////////////:


#[cfg(test)]
mod tests {
    use super::*;
    extern crate fnv;
    #[allow(unused_imports)]
    use self::fnv::FnvHasher; // extern fnv declared in test so we use self::fnv , if declared above we use super::fnv
    #[allow(unused_imports)]
    use crate::nohasher::NoHashHasher;
    

    use crate::base::kmergenerator::*;
    use env_logger;

    // initialize once log system for tests.
    fn init_log_test() {
        // we do not call unwrap as it is an error to init twice and possibly some one else initialized it
        let _ = env_logger::builder().is_test(true).try_init();
    }

 

}  // end of mod test
