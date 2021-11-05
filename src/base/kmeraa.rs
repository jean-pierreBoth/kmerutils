//! This file implements Kmer for Amino Acid.
//! We use uncompressed kmers as generic slices
//! In fact KmerAA although stored in uncompress format trivially implements the trait CompressedKmerT
//! and so KmerSeqIterator and KmerGenerationPattern (see module Kmer) applies to KmerAA 

use std::io;
use std::io::{ErrorKind};

use std::cmp::Ordering;

#[allow(unused)]
use log::{debug,info,error};


use super::kmer::*;

#[derive(Copy,Clone)]
pub struct KmerAA<const N:usize> {
    aa : [u8; N],
} // end of struct KmerAA

impl <const N: usize> KmerAA<N> {

    pub fn new(init : [u8; N]) -> Self {
        KmerAA{aa: init}
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



impl <const N : usize> CompressedKmerT for  KmerAA<N> {
    type Val = [u8; N];

    fn get_nb_base_max() -> usize {
        N
    }

    /// return a clone of itself
    fn get_compressed_value(&self) ->  Self::Val {
        self.aa.clone()
    }

    fn get_uncompressed_kmer(&self) ->  Vec<u8> {
        self.aa.clone().to_vec()
    }

    fn get_bitsize(&self) -> usize {
        N*4
    }

}  // end of impl CompressedKmerT


/*
    Now we have the basics of Kmer Traits we implement KmerSeqIterator and KmerGenerationPattern
 */






//===========================================================


////////////////////////////////////////////////////////////////////////////////////////////////////



#[cfg(test)]
mod tests {

    use super::*;


    fn log_init_test() {
        let mut builder = env_logger::Builder::from_default_env();
        //    builder.filter_level(LevelFilter::Trace);
        let _ = builder.is_test(true).try_init();
    }


}