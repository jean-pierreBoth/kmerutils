//! This file implements Kmer for Amino Acid.
//! We use uncompressed kmers as generic slices
//! In fact KmerAA although stored in uncompress format trivially implements the trait CompressedKmerT
//! and so KmerSeqIterator and KmerGenerationPattern (see module Kmer) applies to KmerAA 

#![allow(unused)]

use std::io;
use std::io::{ErrorKind};

use std::cmp::Ordering;
use std::ops::{Range};

use indexmap::{IndexMap};
use fnv::FnvBuildHasher;
type FnvIndexMap<K, V> = IndexMap<K, V, FnvBuildHasher>;


#[allow(unused)]
use log::{debug,info,error};


use crate::base::kmertraits::*;

/// alphabet of RNA. 
pub struct Alphabet {
    pub bases: String,
}


impl Alphabet {
    pub fn new() -> Alphabet {
        Alphabet { bases : String::from("ACDEFGHIKLMNPQRSTVWY")}
    }
    //
    pub fn len(&self) -> u8 {
        return self.bases.len() as u8;
    }

    #[inline(always)]
    fn is_valid_base(&self, c: u8) -> bool {
        self.bases.find(c as char).is_some() 
    } // end is_valid_base
    
}  // end of impl Alphabet



#[derive(Copy,Clone,Hash)]
pub struct KmerAA<const N:usize> {
    aa : [u8; N],
} // end of struct KmerAA

impl <const N: usize> KmerAA<N> {

    pub fn new(init : [u8; N]) -> Self {
        KmerAA{aa: init}
    }

    fn get_uncompressed_value(&self) -> [u8;N] {
        self.aa
    }
}  // end of impl KmerAA



impl <const N:usize> KmerT for KmerAA<N> {

    fn get_nb_base(&self) -> u8 {
        N as u8
    } // end of get_nb_base

    // 
    fn push(&self, c : u8) -> Self {
        let mut nkmer = KmerAA::new(self.aa);
        nkmer.aa.rotate_left(1);
        nkmer.aa[N-1] = c;
        nkmer
    }  // end of push

    // TODO
    fn reverse_complement(&self) -> Self {
        panic!("KmerAA.KmerAA not yet implemented");
    } // end of reverse_complement


    // 
    fn dump(&self, bufw: &mut dyn io::Write) -> io::Result<usize> {
        // we can transform to a string.
        let str = String::from_utf8(self.aa.to_vec());
        if !str.is_ok() {
            log::error!("KmerAA conversion to a String failed");
            let err = std::io::Error::new(ErrorKind::Other, "KmerAA conversion to a String failed");
            return io::Result::Err(err);
        }
        let _ = bufw.write_fmt(format_args!("{:?}",str));
        // return length of string
        io::Result::Ok(str.unwrap().len())
    }
     
} // end of impl KmerT block for KmerAA


impl <const N : usize> PartialEq for KmerAA<N> {
    // we must check equality of field
    fn eq(&self, other: &KmerAA<N>) -> bool {
        self.aa.iter().eq(other.aa.iter())
    }
}  // end of impl PartialEq for KmerAA

impl <const N:usize> Eq for KmerAA<N> {}



/// We define ordering as a kind of "lexicographic" order by taking into account first number of base.
/// The more the number of base the greater. Then we have integer comparison between lower kmer part
/// which corresponds to lexicographic order as we have A < C < G < T in 2bit and 4bit encoding
/// 




impl <const N : usize> Ord for KmerAA<N> {

    fn cmp(&self, other: &KmerAA<N>) -> Ordering {
        if self.aa.iter().ge(other.aa.iter()) {
            Ordering::Greater
        }
        else {
            Ordering::Less
        }
    } // end cmp
} // end impl Ord for KmerAA 



impl <const N : usize> PartialOrd for KmerAA<N> {
    fn partial_cmp(&self, other: &KmerAA<N>) -> Option<Ordering> {
        Some(self.cmp(other))
    } // end partial_cmp
} // end impl Ord for KmerAA<N>



//=======================================================================

/// our sequence of Amino Acid is encoded on a byte (even if 5 bits are enough but we do not store sequences yet)
type SequenceAA = Vec<u8>;


pub struct KmerSeqIterator<'a, const N : usize> {
    /// size of kmer
    nb_base: usize,
    /// an iterator for base calling
    sequence: &'a SequenceAA,
    /// last position of last kmer returned. At the beginning its None
    previous: Option<usize>,
    ///
    range : Range<usize>,
} // end of KmerSeqIterator


impl<'a, const N : usize> KmerSeqIterator<'a, N> {

    pub fn new(seq : &'a SequenceAA) -> Self {
        let range = std::ops::Range{start : 0, end : seq.len() -1};
        KmerSeqIterator{nb_base : N, sequence : seq, previous : None, range}
    }

    /// iterates...
    pub fn next(&mut self) -> Option<KmerAA<N>> {
        None
    }

    /// defines the range of kmer generation.  
    /// All bases in kmer generated must be between in first..last last excluded!
    fn set_range(&mut self, first: usize, last:usize) -> std::result::Result<(),()> { 
        if last <= first || last > self.sequence.len() {
            return Err(());
        }
        else {
            self.range = Range{start:first, end:last};
            return Ok(());
        }
    } // end of set_range

} // end of impl block for KmerSeqIterator

//============================================================================


pub trait KmerGenerationPattern<T:KmerT> {
    /// generate all kmers included in 0..
    fn generate_kmer_pattern(&self, seq : & SequenceAA) -> Vec<T>;
    /// generate all kmers inclused in begin..end with end excluded as in rust conventions.
    fn generate_kmer_pattern_in_range(&self, seq : & SequenceAA, begin:usize, end:usize) -> Vec<T>;   
    /// generate kmers with their multiplicities
    fn generate_kmer_distribution(&self, seq : & SequenceAA) -> Vec<(T,usize)>;
}



//==========================================================================================


use std::marker::PhantomData;

pub struct KmerGenerator<T:KmerT> {
    /// size of kmer we generate
    pub kmer_size : u8,
    t_marker: PhantomData<T>,
}


impl  <T:KmerT> KmerGenerator<T> {
    pub fn new(ksize:u8) -> Self {
        KmerGenerator{kmer_size: ksize, t_marker : PhantomData}
    }
    /// generic driver for kmer generation
    pub fn generate_kmer (&self, seq : &SequenceAA) -> Vec<T> where Self: KmerGenerationPattern<T> {
        self.generate_kmer_pattern(seq)
    }
    /// generic driver for kmer generation
    pub fn generate_kmer_in_range(&self, seq : & SequenceAA, begin:usize, end:usize) -> Vec<T>
    where Self: KmerGenerationPattern<T> {
        self.generate_kmer_pattern_in_range(seq, begin, end)
    }
    /// generic driver for kmer distribution pattern
    pub fn generate_weighted_kmer(&self, seq : &SequenceAA) -> Vec<(T,usize)>  where Self : KmerGenerationPattern<T> {
        self.generate_kmer_distribution(seq)
    }
    ///
    pub fn get_kmer_size(&self) -> usize { self.kmer_size as usize}
}  // end of impl KmerGenerator



/*
    Now we have the basics of Kmer Traits we implement KmerSeqIterator and KmerGenerationPattern
 */




/// implementation of kmer generation pattern for KmerAA<N>
impl <const N : usize> KmerGenerationPattern<KmerAA<N>> for KmerGenerator<KmerAA<N>> {
    fn generate_kmer_pattern(&self, seq : &SequenceAA) -> Vec<KmerAA<N>> {
        if self.kmer_size as usize != N {
            panic!("KmerAA<N> has not the correct size!!");   // cannot happen !
        }
        // For a sequence of size the number of kmer is seq.size - kmer.size + 1  !!!
        // But it happens that "long reads" are really short 
        let nb_kmer = if seq.len() >= N { seq.len()-N+1} else {0};
        let mut kmer_vect = Vec::<KmerAA<N>>::with_capacity(nb_kmer);
        let mut kmeriter  = KmerSeqIterator::<N>::new(seq);
        loop {
            match kmeriter.next() {
                Some(kmer) => kmer_vect.push(kmer),
                None => break,
            }
        }
        //
        return kmer_vect;
    }  // end of generate_kmer_pattern

    /// generate all kmers associated to their multiplicity
    /// This is useful in the context of Jaccard Probability Index estimated with ProbminHash 
    fn generate_kmer_distribution(&self, seq : &SequenceAA) -> Vec<(KmerAA<N>,usize)> {
        if self.kmer_size as usize != N {
            panic!("KmerAA<N> has not the correct size!!");  // cannot happen
        }
        // For a sequence of size the number of kmer is seq.size - kmer.size + 1  !!!
        // But it happens that "long reads" are really short 
        let nb_kmer = if seq.len() >= N { seq.len()-N+1} else {0};
        let mut kmer_distribution : FnvIndexMap::<KmerAA<N>,usize> = FnvIndexMap::with_capacity_and_hasher(nb_kmer, FnvBuildHasher::default());
        let mut kmeriter = KmerSeqIterator::<N>::new(seq);
        loop {
            match kmeriter.next(){
                Some(kmer) => {
                    // do we store the kmer in the FnvIndexMap or a already hashed value aka nthash?
                    *kmer_distribution.entry(kmer).or_insert(0) += 1;
                },
                None => break,
            }
        }
        // convert to a Vec
        let mut hashed_kmers = kmer_distribution.keys();
        let mut weighted_kmer = Vec::<(KmerAA<N>,usize)>::with_capacity(kmer_distribution.len());
        loop {
            match hashed_kmers.next() {
                Some(key) => {
                    if let Some(weight) = kmer_distribution.get(key) {
                        // get back to Kmer16b32bit from 
                        weighted_kmer.push((*key,*weight));
                    };
                },
                None => break,
            }
        }
        //
        return weighted_kmer;
    }  // end of generate_kmer_pattern



    fn generate_kmer_pattern_in_range(&self, seq : &SequenceAA, begin:usize, end:usize) -> Vec<KmerAA<N>> {
        if self.kmer_size as usize != N {
            panic!("KmerAA<N> has not the correct size!!");   // cannot happen
        }
        if begin >= end {
            panic!("KmerGenerationPattern<'a, Kmer16b32bit>  bad range for kmer iteration");
        }
        // For a sequence of size the number of kmer is seq.size - kmer.size + 1  !!!
        // But it happens that "long reads" are really short 
        let nb_kmer = if seq.len() >= 16 { seq.len()-N+1} else {0};
        let mut kmer_vect = Vec::<KmerAA<N>>::with_capacity(nb_kmer);
        let mut kmeriter = KmerSeqIterator::<N>::new(seq);
        kmeriter.set_range(begin, end).unwrap();
        loop {
            match kmeriter.next() {
                Some(kmer) => kmer_vect.push(kmer),
                None => break,
            }
        }
        //
        return kmer_vect;
    }  // end of generate_kmer_pattern

}  // end of impl KmerGenerationPattern<'a, KmerAA<N>>




//===========================================================






#[cfg(test)]
mod tests {



    fn log_init_test() {
        let mut builder = env_logger::Builder::from_default_env();
        //    builder.filter_level(LevelFilter::Trace);
        let _ = builder.is_test(true).try_init();
    }

    // test iterator

    // test iterator with range
}  // end of mod tests