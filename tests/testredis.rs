extern crate kmerutils;

// to be launched by cargo test --test testredis [ -- --nocapture] 
// Cf Ch11-03 of rust book second edition

extern crate log;
use env_logger;

use std::rc::Rc;

pub use self::kmerutils::prelude::*;

use kmerutils::anchor::*;


//#[cfg(test)]

#[test]
fn test_mininvhash_redis() {
    //
    env_logger::Builder::from_default_env().init();
    //
    let anchor_parameters = Rc::new(AnchorsGeneratorParameters::new(String::from("./tests/2readEcoli.fastq"),
                                                                200, 3, 16,100));

    let mut f_anchors = FastaAnchors::<Kmer16b32bit>::new(Rc::clone(&anchor_parameters), true, None);
    let res = f_anchors.anchor_computation();
    println!("nb anchors computed {}", f_anchors.get_nb_anchors());
    //
    assert!(res.is_ok());
    
}   // end of test_mininvhash_redis
