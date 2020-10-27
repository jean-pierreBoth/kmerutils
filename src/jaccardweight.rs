#[allow(unused_imports)]
use log::*;

#[allow(unused_imports)]
use std::hash::{BuildHasher, BuildHasherDefault, Hasher, Hash};
use indexmap::{IndexMap};
use fnv::{FnvBuildHasher};

use crate::nohasher::*;
use crate::kmer::*;
use crate::kmergenerator::*;


type FnvIndexMap<K, V> = IndexMap<K, V, FnvBuildHasher>;

use probminhash::probminhash::*;



/// compute jaccard probability index between a sequence and a vecor of sequences sequences for Kmer16b32bit
/// and returns a vector of Jaccard probability index
pub fn sketch_seqrange_probminhash3a_kmer16b32bit<F>(seqa: &Sequence, vseqb : Vec<&Sequence>, sketch_size: usize, fhash : F) -> Vec<f64> 
    where F : Fn(&Kmer16b32bit) -> u64 {
    //
    info!("seqsketcher : entering sketch_seqrange_probminhash3a");
    // a vector to return results
    let mut jaccard_vec = Vec::<f64>::with_capacity(vseqb.len());
    // default is invertible hash and then superminhash without any hashing
    let mut pminhasha = ProbMinHash3a::<usize,NoHashHasher>::new(sketch_size, 0);
    let mut wa : FnvIndexMap::<usize,f64> = FnvIndexMap::with_capacity_and_hasher(seqa.size(), FnvBuildHasher::default());
    //
    // generate all kmers include in range arg. dependance upon kmer_size 
    // seqa
    let mut kmergen = KmerSeqIterator::<Kmer16b32bit>::new(16, &seqa);
    kmergen.set_range(1, seqa.size()).unwrap();
    loop {
        match kmergen.next() {
            Some(kmer) => {
                let hashval = fhash(&kmer);
//                            let canonical =  kmer.reverse_complement().min(kmer);
//                            let hashval = probminhash::invhash::int32_hash(canonical.0);
                *wa.entry(hashval as usize).or_insert(0.) += 1.;
            },
            None => break,
        }
    }  // end loop
    pminhasha.hash_weigthed_idxmap(&wa);
    let siga = pminhasha.get_signature();
    // loop on vseqb to // with rayon
    for seqb in vseqb {
        let mut wb : FnvIndexMap::<usize,f64> = FnvIndexMap::with_capacity_and_hasher(seqb.size(), FnvBuildHasher::default());
        let mut kmergen = KmerSeqIterator::<Kmer16b32bit>::new(16, &seqb);
        kmergen.set_range(1, seqb.size()).unwrap();
        loop {
            match kmergen.next() {
                Some(kmer) => {
                    let hashval = fhash(&kmer);
    //                let canonical =  kmer.reverse_complement().min(kmer);
    //                let hashval = probminhash::invhash::int32_hash(canonical.0);
                    *wb.entry(hashval as usize).or_insert(0.) += 1.;
                },
                None => break,
            }
        }  // end loop 
        let mut pminhashb = ProbMinHash3a::<usize,NoHashHasher>::new(sketch_size, 0);
        pminhashb.hash_weigthed_idxmap(&wb);
        let sigb = pminhasha.get_signature();
        let jac64 = compute_probminhash_jaccard(siga, sigb);
        jaccard_vec.push(jac64);
    }
    return jaccard_vec;
    //
} // end of sketch_seqrange_probminhash3a_kmer16b32bit




pub fn sketch_seqrange_probminhash3a_kmer32bit<F>(seqa: &Sequence, vseqb : Vec<&Sequence>, sketch_size: usize, fhash : F) -> Vec<f64> 
    where F : Fn(&Kmer32bit) -> u32 {
    //
    info!("seqsketcher : entering sketch_seqrange_probminhash3a");
    // a vector to return results
    let mut jaccard_vec = Vec::<f64>::with_capacity(vseqb.len());
    // default is invertible hash and then superminhash without any hashing
    let mut pminhasha = ProbMinHash3a::<usize,NoHashHasher>::new(sketch_size, 0);
    let mut wa : FnvIndexMap::<usize,f64> = FnvIndexMap::with_capacity_and_hasher(seqa.size(), FnvBuildHasher::default());
    //
    // generate all kmers include in range arg. dependance upon kmer_size 
    // seqa
    let mut kmergen = KmerSeqIterator::<Kmer32bit>::new(16, &seqa);
    kmergen.set_range(1, seqa.size()).unwrap();
    loop {
        match kmergen.next() {
            Some(kmer) => {
                let hashval = fhash(&kmer);
//                            let canonical =  kmer.reverse_complement().min(kmer);
//                            let hashval = probminhash::invhash::int32_hash(canonical.0);
                *wa.entry(hashval as usize).or_insert(0.) += 1.;
            },
            None => break,
        }
    }  // end loop
    pminhasha.hash_weigthed_idxmap(&wa);
    let siga = pminhasha.get_signature();
    // loop on vseqb to // with rayon
    for seqb in vseqb {
        let mut wb : FnvIndexMap::<usize,f64> = FnvIndexMap::with_capacity_and_hasher(seqb.size(), FnvBuildHasher::default());
        let mut kmergen = KmerSeqIterator::<Kmer32bit>::new(16, &seqb);
        kmergen.set_range(1, seqb.size()).unwrap();
        loop {
            match kmergen.next() {
                Some(kmer) => {
                    let hashval = fhash(&kmer);
    //                let canonical =  kmer.reverse_complement().min(kmer);
    //                let hashval = probminhash::invhash::int32_hash(canonical.0);
                    *wb.entry(hashval as usize).or_insert(0.) += 1.;
                },
                None => break,
            }
        }  // end loop 
        let mut pminhashb = ProbMinHash3a::<usize,NoHashHasher>::new(sketch_size, 0);
        pminhashb.hash_weigthed_idxmap(&wb);
        let sigb = pminhasha.get_signature();
        let jac64 = compute_probminhash_jaccard(siga, sigb);
        jaccard_vec.push(jac64);
    }
    return jaccard_vec;
    //
} // end of sketch_seqrange_probminhash3a_kmer32bit





//===========================================================================================================

#[cfg(test)]
mod tests {
    use super::*;


    #[allow(dead_code)]
    fn log_init_test() {
        let _ = env_logger::builder().is_test(true).try_init();
    }

    #[test]
    // This function tests probability jaccard estimates on less than 16 bases 
    fn test_probminhash_kmer16b() {
        // initialize test logging
        log_init_test();
        //
        let seqstr = String::from("TCAAAGGGAAACATTCAAAATCAGTATGCGCCCGTTCAGTTACGTATTGCTCTCGCTAATGAGATGGGCTGGGTACAGAG");
        let seqabytes = seqstr.as_bytes();
        let seqa = Sequence::new(seqstr.as_bytes(),2);
        //
        let mut vecseqb = Vec::<&Sequence>::new();
        let seqb1 = Sequence::new(&seqabytes[1..40],2);   // half the length of seqa
        vecseqb.push(&seqb1);
        //
        let kmer_hash_fn = | kmer : &Kmer32bit | -> u32 {
            let canonical =  kmer.reverse_complement().min(*kmer);
            let hashval = probminhash::invhash::int32_hash(canonical.0);
            hashval
        };

        let vecsig = sketch_seqrange_probminhash3a_kmer32bit(&seqa, vecseqb, 50, kmer_hash_fn);
        debug!("vecsig {:?}", vecsig);
    }  // end of test_probminhash_kmer16b

}  // end of mod test

