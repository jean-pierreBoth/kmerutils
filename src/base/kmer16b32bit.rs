//! module for Kmer16b32bit, encoding 16 base kmer on a u32
//! 


use std::mem;
use std::io;


use std::cmp::Ord;
use std::str::FromStr;

#[allow(unused)]
use log::{debug,trace};



// our includes

pub use super::{kmertraits::*};
pub use super::{nthash::*, alphabet::*};

/// type supporting 16 bases kmer as u32.

#[derive(Clone,Copy,PartialEq,Eq,PartialOrd,Ord,Hash,Debug)]
pub struct Kmer16b32bit(pub u32);


impl Kmer16b32bit {
    /// allocate a new kmer. initilized to 0
    pub fn new() -> Kmer16b32bit {
        Kmer16b32bit(0u32)
    } // end of new
}


impl KmerT for Kmer16b32bit {
    /// returns 12 as this struct is just for 16 base kmer
    #[inline(always)]
    fn get_nb_base(&self) -> u8 {
        16
    }
    /// just returns the reversed complement 16 bases kmer in 2bit encoding. Ch Hacker's delight.
    fn reverse_complement(&self) ->  Kmer16b32bit {
        // This depends on our choice for encoding ACGT as respectiveley  00  01 10 11 !!!
        // we use the reverse instruction which simplifies the previous impl based on  Hacker's delight 
        // we keep Hacker's delight trick to permut groups of 2bits to restore bases after reverse
        let mut revcomp = self.0 as u32;
        revcomp = !revcomp;
        //
        revcomp = revcomp.reverse_bits();
        revcomp = (revcomp & 0x55555555) << 1 | (revcomp & 0xAAAAAAAA) >> 1;
        // we complement
        Kmer16b32bit(revcomp)
    }
    /// push a base (2bits) at right end of kmer producing a new Kmer
    /// So arg base must be 2bit encoded !!! and there is no sure way to ensure arg is coherent
    fn push(&self, base : u8) -> Kmer16b32bit {
        // check base is encode ?
        let new_kmer = (self.0 << 2) | (base as u32 & 0b11);
        Kmer16b32bit(new_kmer)     
    }

    ///
    fn dump(&self, bufw: &mut dyn io::Write) -> io::Result<usize> {
        bufw.write(unsafe { &mem::transmute::<u32, [u8;4]>(self.0) } )
    }
    
} // end implementation 





impl CompressedKmerT for Kmer16b32bit {
    type Val = u32;
    ///
    fn get_nb_base_max() -> usize { 16 }
    /// a decompressing function mainly for test and debugging purpose
    fn get_uncompressed_kmer(&self) -> Vec<u8> {
        let alphabet = Alphabet2b::new();
        // we treat each block of 2 bis as u8 end call decoder of Alphabet2b
        let mut decompressed_kmer = Vec::<u8>::with_capacity(16);
        let mut base:u8;
        //
        let mut buf = self.0;
        for _ in 0..16 {
            buf = buf.rotate_left(2);
            base = (buf & 0b11) as u8; 
            decompressed_kmer.push(alphabet.decode(base));
        }
        return decompressed_kmer;
    }
    ///
    #[inline(always)]    
    fn get_compressed_value(&self) -> u32 {
        return self.0;
    }
    ///
    #[inline(always)]    
    fn get_bitsize(&self) -> usize { 32 }
}  // end of impl block for CompressedKmerT



impl KmerBuilder<Kmer16b32bit> for Kmer16b32bit {
    /// we know the number of bases by definition of the Kmer
    fn build(val: u32, _nb_base : u8) -> Kmer16b32bit {
        Kmer16b32bit(val)
    }
}

 
impl FromStr for Kmer16b32bit {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        // check that length is 16.
        if s.len() != 16 {
            return Err(String::from("length of kmer should be 16"));
        }
        let nb_bases=s.len();
        let mut kmer = Kmer16b32bit::new();
        let alphabet = Alphabet2b::new();
        //
        let sbytes = s.as_bytes();
        //
        for i in 0..nb_bases {
            if !alphabet.is_valid_base(sbytes[i]) {
                return Err(String::from("char not in ACGT"));
            }
            kmer = kmer.push(alphabet.encode(sbytes[i]));
        }
        Ok(kmer)
    } // end of from
} // end impl FromStr


//========================================================

#[cfg(test)]
mod tests {

#[allow(unused)]
use super::*;

#[test]
    fn test_reverse_complement_16b32bit() {
        //
        let to_reverse : Vec<Kmer16b32bit> = vec![
            // TACG_AGTA_GGAT_ACTA
            Kmer16b32bit(0b11000110_00101100_10100011_00011100),
            // ACTT_GGAA_CGTT_AATG
            Kmer16b32bit(0b00011111_10100000_01101111_00001110)
        ];
        //
        let reversedcomp : Vec<Kmer16b32bit> = vec! [
            // TAGT_ATCC_TACT_CGTA
            Kmer16b32bit(0b11001011_00110101_11000111_01101100),
            // CATT_AACG_TTCC_AAGT
            Kmer16b32bit(0b01001111_00000110_11110101_00001011)
        ];
        
        for i in 0..to_reverse.len() {
            let revcomp = to_reverse[i].reverse_complement();
            if revcomp != reversedcomp[i] {
                println!(" i  kmer reversed complement  : {} {:b}", i, revcomp.0);
            }
            assert!(revcomp == reversedcomp[i]);
        }
    } // end of test_reverse_complement_16b32bit()
}