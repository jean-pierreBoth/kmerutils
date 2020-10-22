#[allow(unused_imports)]
use log::{info,trace};

#[allow(unused_imports)]
use std::hash::{BuildHasher, BuildHasherDefault, Hasher, Hash};
use indexmap::{IndexMap};
use fnv::{FnvBuildHasher};

use crate::nohasher::*;
use crate::kmer::*;
use crate::kmergenerator::*;


type FnvIndexMap<K, V> = IndexMap<K, V, FnvBuildHasher>;

use probminhash::probminhash::*;

/// compute jaccard probability index between 2 sequences according to kmer size
pub fn sketch_seqrange_probminhash3a(seqa: &Sequence, seqb : &Sequence, kmer_size:usize, sketch_size: usize) -> f64 {
    //
    info!("seqsketcher : entering sketch_seqrange_probminhash3a");
    // default is invertible hash and then superminhash without any hashing
    let mut pminhasha = ProbMinHash3a::<usize,NoHashHasher>::new(sketch_size, 0);
    let mut pminhashb = ProbMinHash3a::<usize,NoHashHasher>::new(sketch_size, 0);
    let mut wa : FnvIndexMap::<usize,f64> = FnvIndexMap::with_capacity_and_hasher(70, FnvBuildHasher::default());
    let mut wb : FnvIndexMap::<usize,f64> = FnvIndexMap::with_capacity_and_hasher(70, FnvBuildHasher::default());
    //
    // generate all kmers include in range arg. dependance upon kmer_size 
    match kmer_size {
        16 => { // seqa
                let mut kmergen = KmerSeqIterator::<Kmer16b32bit>::new(16, &seqa);
                kmergen.set_range(1, seqa.size()).unwrap();
                loop {
                    match kmergen.next() {
                        Some(kmer) => {
                            let canonical =  kmer.reverse_complement().min(kmer);
                            let hashval = probminhash::invhash::int32_hash(canonical.0);
                            *wa.entry(hashval as usize).or_insert(0.) += 1.;
                        },
                         None => break,
                    }
                }  // end loop
                pminhasha.hash_weigthed_idxmap(&wa);
                let siga = pminhasha.get_signature();
                // seqb
                let mut kmergen = KmerSeqIterator::<Kmer16b32bit>::new(16, &seqb);
                kmergen.set_range(1, seqb.size()).unwrap();
                loop {
                    match kmergen.next() {
                        Some(kmer) => {
                            let canonical =  kmer.reverse_complement().min(kmer);
                            let hashval = probminhash::invhash::int32_hash(canonical.0);
                            *wb.entry(hashval as usize).or_insert(0.) += 1.;
                        },
                         None => break,
                    }
                }  // end loop 
                pminhashb.hash_weigthed_idxmap(&wb);
                let sigb = pminhasha.get_signature();
                let jac64 = compute_probminhash_jaccard(siga, sigb);
                return jac64;
        },  // end case 16
        _ => panic!("sketch_sequence_superminhash , unimplemented kmer_size {} {} {} ", kmer_size, file!(), line!()),
    }
    //
} // end of sketch_seqrange_probminhash3a_kmer1b32bit