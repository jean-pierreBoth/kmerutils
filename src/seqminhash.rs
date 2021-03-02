//! This module computes sketches of (part of) sequences via superminhash or minhash
//! of adjacent blocks in read sequences
//! See also jaccardweight module

use std::ops::Range;
use std::hash::{BuildHasherDefault};
use crate::base::{kmer::*, kmergenerator::*};

use crate::minhash::*;
use probminhash::invhash;
use probminhash::superminhasher::*;

use log::info;
use log::trace;

/// computes a sketch of kmer of size kmer_size generated from range in seq with Superminhash algorithm
pub fn sketch_seqrange_superminhash(seq: &Sequence, range:&Range<usize>, kmer_size:usize, sketch_size: usize) -> Vec<f64> {
    //
    info!("seqsketcher : entering superminhash_sketch_sequence");
    // default is invertible hash and then superminhash without any hashing
    let bh = BuildHasherDefault::<NoHashHasher>::default();
    let mut sminhash : SuperMinHash<u32, NoHashHasher>= SuperMinHash::new(sketch_size, &bh);
    //
    // generate all kmers include in range arg. dependance upon kmer_size 
    match kmer_size {
        16 => {
            let mut kmergen = KmerSeqIterator::<Kmer16b32bit>::new(16, &seq);
            kmergen.set_range(range.start, range.end).unwrap();
            loop {
                match kmergen.next() {
                    Some(kmer) => {
                        let canonical =  kmer.reverse_complement().min(kmer);
                        let hashval = invhash::int32_hash(canonical.0);
                        sminhash.sketch(&hashval).unwrap();
                    },
                    None => break,
                }
            } // end loop
            //
            sminhash.get_hsketch().clone()
        },
        9..=15 => {
            let mut kmergen = KmerSeqIterator::<Kmer32bit>::new(kmer_size as u8, &seq);
            kmergen.set_range(range.start, range.end).unwrap();
            loop {
                match kmergen.next() {
                    Some(kmer) => {
                        let canonical =  kmer.reverse_complement().min(kmer);
                        let hashval = invhash::int32_hash(canonical.0);
                        sminhash.sketch(&hashval).unwrap();
                    },
                    None => break,
                }
            } // end loop
            info!("seqsketcher : got a  superminhash_sketch");            
            sminhash.get_hsketch().clone()
        },
        _ => panic!("sketch_sequence_superminhash , unimplemented kmer_size {} {} {} ", kmer_size, file!(), line!()),
    }
} // end of sketch_sequence_superminhash


/// computes a sketch and count of kmers of size kmer_size generated from range in seq with Superminhash algorithm
pub fn sketch_seqrange_minhash(seq: &Sequence, range:&Range<usize>, kmer_size:usize, sketch_size: usize) -> Vec<HashCount<u32> > {
    //
    trace!("seqsketcher : entering sketch_seqrange_minhash");
    // default is invertible hash and then superminhash without any hashing
    let mut minhash : MinHashCount<u32, NoHashHasher> = MinHashCount::new(sketch_size, true);
    //
    // generate all kmers include in range arg. dependance upon kmer_size 
    match kmer_size {
        16 => {
            let mut kmergen = KmerSeqIterator::<Kmer16b32bit>::new(16, &seq);
            if let Err(_) = kmergen.set_range(range.start, range.end) {
                println!("sketch_seqrange_minhash: bad range, start = {} , end = {}", range.start, range.end);
                panic!("bad range");
            }
            loop {
                match kmergen.next() {
                    Some(kmer) => {
                        let canonical = kmer.reverse_complement().min(kmer);
                        let hashval = invhash::int32_hash(canonical.0);
//                        let hashval = kmer.0;
//                        trace!("pushed hasv : {} ",  hashval);
                        minhash.push(&hashval);
                    },
                    None => break,
                }
            } // end loop
            return minhash.get_sketchcount();
        },
        9..=15 => {
            let mut kmergen = KmerSeqIterator::<Kmer32bit>::new(kmer_size as u8, &seq);
            if let Err(_) = kmergen.set_range(range.start, range.end) {
                println!("sketch_seqrange_minhash: bad range, start = {} , end = {}", range.start, range.end);
                panic!("bad range");
            }
            
            loop {
                match kmergen.next() {
                    Some(kmer) => {
                        let canonical =  kmer.reverse_complement().min(kmer);
                        let hashval = invhash::int32_hash(canonical.0);
//                        debug!("pushed hasval : {} ",  hashval);
                        minhash.push(&hashval);
                    },
                    None => break,
                }
            } // end loop
            return minhash.get_sketchcount();
        },
        _ => panic!("sketch_sequence_minhash , unimplemented kmer_size {} {} {} ", kmer_size, file!(), line!()),
    }
}  // end of sketch_seqrange_minhash


////////////////////////////////////////////////////////////////////////////////////////////////////



#[cfg(test)]
mod tests {
    use super::*;
    //
    #[test]
    fn test_minhash_overlapping_ranges_kmer16b() {
        let seqstr = String::from("TCAAAGGGAAACATTCAAAATCAGTATGCGCCCGTTCAGTTACGTATTGCTCTCGCTAATGAGATGGGCTGGGTACAGAG");
        // transform to a Vec<u8> and then a Sequence
        let slu8 = seqstr.as_bytes();
        // get a sequence with 2 bits compression
        let seq = Sequence::new(&slu8,2);
        let range1 :  Range<usize> = Range{start:1, end:65};
        let kmer_size : usize = 16;
        let sketch_size: usize = 20;
        trace!(" sketching 1");
        // we get kmer from [1,16] to  [49,64] i.e 48 kmer
        let sk1 = sketch_seqrange_minhash(&seq, &range1, kmer_size, sketch_size);
        let range2 :  Range<usize> = Range{start:35, end:75};
        trace!("\n sketching 2");
         // we get kmer from [35,50] to  [59,74] i.e 25 kmer       
        let sk2 = sketch_seqrange_minhash(&seq, &range2, kmer_size, sketch_size);
        // we search  common 16-mers. they are the kmer from [35,50] to [49,64] ie 14 kmers
        // union of the 2 sets of kmer has cardinal 59
        let resdist:MinHashDist  = minhash_distance(&sk1, &sk2);
        //
        println!("distance minhash (contain, dist, common, total): {}  {}   {}  {} ",
               resdist.0, resdist.1, resdist.2, resdist.3);
        // sketch size id 20 i.e 1/3 of total 
        assert!(resdist.3 >= 3);
    } // end of  test_minhash_overlapping_ranges


    #[test]
    fn test_minhash_overlapping_ranges_kmer10b() {
        // a 80 character sequence
        let seqstr = String::from("TCAAAGGGAAACATTCAAAATCAGTATGCGCCCGTTCAGTTACGTATTGCTCTCGCTAATGAGATGGGCTGGGTACAGAG");
        // transform to a Vec<u8> and then a Sequence
        let slu8 = seqstr.as_bytes();
        // get a sequence with 2 bits compression
        let seq = Sequence::new(&slu8,2);
        let range1 :  Range<usize> = Range{start:1, end:65};
        let kmer_size : usize = 10;
        let sketch_size: usize = 20;
        trace!(" sketching 1");
        // we get kmer from [1,10] to  [55,64]] i.e 55 kmer
        let sk1 = sketch_seqrange_minhash(&seq, &range1, kmer_size, sketch_size);
        let range2 :  Range<usize> = Range{start:35, end:75};
        trace!("\n sketching 2");
         // we get kmer from [35,44] to  [65,74] i.e 30 kmer       
        let sk2 = sketch_seqrange_minhash(&seq, &range2, kmer_size, sketch_size);
        // we search  common 10-mers. they are the kmer from [35,44] to [55,64] ie 20 kmers
        // union of the 2 sets of kmer has cardinal 65
        let resdist:MinHashDist  = minhash_distance(&sk1, &sk2);
        //
        println!("distance super minhash (contain, dist, common, total): {} {} {} {}  ", resdist.0, resdist.1, resdist.2, resdist.3);
        //
        assert!(resdist.3 == 20);
    } // end of  test_superminhash_overlapping_ranges


    // Superminhash tests
    
    #[test]
    fn test_superminhash_overlapping_ranges_kmer16b() {
        let seqstr = String::from("TCAAAGGGAAACATTCAAAATCAGTATGCGCCCGTTCAGTTACGTATTGCTCTCGCTAATGAGATGGGCTGGGTACAGAG");
        // transform to a Vec<u8> and then a Sequence
        let slu8 = seqstr.as_bytes();
        // get a sequence with 2 bits compression
        let seq = Sequence::new(&slu8,2);
        let range1 :  Range<usize> = Range{start:1, end:65};
        let kmer_size : usize = 16;
        let sketch_size: usize = 50;
        trace!(" sketching 1");
        // we get kmer from [1,16] to  [49,64] i.e 48 kmer
        let sk1 = sketch_seqrange_superminhash(&seq, &range1, kmer_size, sketch_size);
        let range2 :  Range<usize> = Range{start:35, end:75};
        trace!("\n sketching 2");
         // we get kmer from [35,50] to  [59,74] i.e 25 kmer       
        let sk2 = sketch_seqrange_superminhash(&seq, &range2, kmer_size, sketch_size);
        // we search  common 16-mers. they are the kmer from [35,50] to [49,64] ie 14 kmers
        // union of the 2 sets of kmer has cardinal 59
        let resdist  = get_jaccard_index_estimate(&sk1, &sk2).unwrap();
        //
        println!("distance super minhash (contain, dist, common, total): {}  ", resdist);
        //
        assert!( resdist >= 0.15);
    } // end of  test_superminhash_overlapping_ranges


    // test sumper min hash with 16 base kmers.
    #[test]
    fn test_superminhash_overlapping_ranges_kmer10b() {
        let seqstr = String::from("TCAAAGGGAAACATTCAAAATCAGTATGCGCCCGTTCAGTTACGTATTGCTCTCGCTAATGAGATGGGCTGGGTACAGAG");
        // transform to a Vec<u8> and then a Sequence
        let slu8 = seqstr.as_bytes();
        // get a sequence with 2 bits compression
        let seq = Sequence::new(&slu8,2);
        let range1 :  Range<usize> = Range{start:1, end:65};
        let kmer_size : usize = 10;
        let sketch_size: usize = 20;
        trace!(" sketching 1");
        // we get kmer from [1,10] to  [55,64]] i.e 54 kmer
        let sk1 = sketch_seqrange_superminhash(&seq, &range1, kmer_size, sketch_size);
        let range2 :  Range<usize> = Range{start:35, end:75};
        trace!("\n sketching 2");
         // we get kmer from [35,44] to  [65,74] i.e 30 kmer       
        let sk2 = sketch_seqrange_superminhash(&seq, &range2, kmer_size, sketch_size);
        // we search  common 10-mers. they are the kmer from [35,44] to [55,64] ie 20 kmers
        // union of the 2 sets of kmer has cardinal 65 so we have an intersection of 20/65.
        let resdist  = get_jaccard_index_estimate(&sk1, &sk2).unwrap();
        //
        println!("distance super minhash : {}  ", resdist);
        //
        assert!(resdist >= 0.2);
    } // end of  test_superminhash_overlapping_ranges






    
}   // end of mod tests





