//! defines kmer traits


use std::io;
use std::hash::Hash;
use std::cmp::Ord;


/// A Kmer is defined with 2 traits, one for standard operations and
/// one for taking into account that our Kmers are most often 2 r bits encoded
/// so we need a trait for compression characteristics.
///
/// A trait all Kmer should implement. compressed or not.
pub trait KmerT {
    /// returns number of bases of kmer
    fn get_nb_base(&self) -> u8;
    /// returns the reverse complement kmer
    fn reverse_complement(&self) -> Self;
    /// push a (compressed! or 2 bit encoded) base at right end of kmer
    fn push(&self, base : u8) -> Self;
    /// each kmer type must know how to dump itself.
    fn dump(&self,  bufw : &mut dyn io::Write) -> io::Result<usize>;
}


/// all our Kmer implements this trait as they use 2/4 bit base encoding.
/// So basically our Kmer are copy, orderable, representable as hashable orderable compressed value.
/// which is a minimum
pub trait CompressedKmerT : KmerT+Ord+Copy  where Self::Val : Hash+Ord
{
    /// type of compressed value , u16, u32, u64.
    type Val;
    /// returns the max number of base supported by compressing strategy and size of Val
    fn get_nb_base_max() -> usize;
    /// return encoded value in type Val. In Val we have encoded number of bases and
    /// and kmer value. We do not get rid of number of base, and return raw value
    fn get_compressed_value(&self) -> Self::Val;
    /// get Kmer as a Vec<u8>
    fn get_uncompressed_kmer(&self) -> Vec<u8>;
    /// returns the size in bits of word supporting compressed kmer.
    /// not just number of base * size of base
    fn get_bitsize(&self) -> usize;
}  // end of trait CompressedKmer

