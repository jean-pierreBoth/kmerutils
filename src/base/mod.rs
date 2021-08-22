//! Contains basic tools to describe sequences, kmers , generate and count Kmers  
//! It implements multithreaded counters with Bloom filters and Cuckoo filters


// textual scope of macro (without macro_export in which case global crate visibility)
#[macro_use]
pub mod nthash;

pub mod alphabet;
pub mod sequence;
pub mod kmer;
pub mod kmercount;
pub mod kmergenerator;

