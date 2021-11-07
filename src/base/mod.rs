//! Contains basic tools to describe DNA-sequences, kmers , generate and count Kmers  
//! It implements multithreaded counters with Bloom filters and Cuckoo filters


// textual scope of macro (without macro_export in which case global crate visibility)
#[macro_use]
pub mod nthash;

pub use sequence::*;
pub use kmertraits::*;
pub use kmer::*;
pub use kmer32bit::*;
pub use kmer16b32bit::*;
pub use kmer64bit::*;


pub use alphabet::*;

pub mod alphabet;
pub mod sequence;
pub mod kmertraits;
pub mod kmer;
pub mod kmer32bit;
pub mod kmer16b32bit;
pub mod kmer64bit;

pub mod kmercount;
pub mod kmergenerator;

