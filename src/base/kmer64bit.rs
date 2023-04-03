//! implementation of Kmer64bit representing up to 28 bases encoded in u64
//! 
//! 

use std::mem;
use std::io;
use std::cmp::Ordering;
use std::cmp::Ord;
use std::str::FromStr;

#[allow(unused)]
use log::{debug,trace};

pub use super::kmertraits::*;
pub use super::{nthash::*, alphabet::*};


// For more than 16 bases the coding of number of bases in a separate field is more efficient
// This representation is consistent with Kmer32bit as self.0 gives the word supporting value of kmer

/// The type supporting Kmer for number of bases between 17 and 32.
#[derive(Clone,Copy,Debug, Hash)]
pub struct  Kmer64bit(pub u64, pub u8);


impl Kmer64bit {
    pub fn new(nb_base: u8) -> Kmer64bit {
        Kmer64bit(0u64, nb_base)
    }
}



impl PartialEq for Kmer64bit {
    // we must check number of bases and equality of field
    fn eq(&self, other: &Kmer64bit) -> bool {
        if (self.0 == other.0) & (self.1==other.1) { true } else {false}
    } // end of eq
}  // end of PartialEq implementation


// read the doc for Eq ....
impl Eq for Kmer64bit{}


/// We define ordering as a kind of "lexicographic" order by taking into account first number of base.
/// The more the number of base the greater. Then we have integer comparison between lower kmer part
/// which corresponds to lexicographic order as we have A < C < G < T in 2bit

impl Ord for Kmer64bit {
    fn cmp(&self, other: &Kmer64bit) -> Ordering {
        if self.1 != other.1 {
            return (self.1).cmp(&(other.1));
        }
        else {
            return (self.0).cmp(&(other.0));
        }
    } // end cmp
} // end impl Ord for Kmer32bit


impl PartialOrd for Kmer64bit {
    fn partial_cmp(&self, other: &Kmer64bit) -> Option<Ordering> {
        Some(self.cmp(other))
    } // end partial_cmp
} // end impl Ord for Kmer32bit




impl KmerT for Kmer64bit {
    /// retrieves the number of bases this kmer stores
    #[inline(always)]
    fn get_nb_base(&self) -> u8 {
        self.1
    }

    fn push(&self, base : u8) -> Kmer64bit {
        // shift left 2 bits, insert new base and enforce saved number of bases.
        // Note that the left shift pushes garbage bit between 4 upper bits and lower bits coding value.
        // We are cautious to clean them by using value_mask which reset those bit to 0!
        // It is useful when implementing PartialEq and compressed value 
        // We could use as a mask for value field : (0b1 << (2*self.get_nb_bases())) - 1 which enforce 0 bit between 4 upper bits
        // and lower bits coding value.
        let value_mask :u64 = (0b1 << (2*self.get_nb_base())) - 1;
        // shift left 2 bits, insert new base and enforce saved number of bases.
        let new_kmer = ((self.0 << 2) & value_mask) | (base as u64 & 0b11);
        trace!("in push new_kmer = {:b}",  new_kmer);
        Kmer64bit(new_kmer, self.1)
    }

    /// just returns the reversed complement 16 bases kmer in 2bit encoding. Ch Hacker's delight.
    fn reverse_complement(&self) ->  Kmer64bit {
        // we use the reverse instruction and 
        // we do a symetry as explained in Hacker's delight and complement.
        // This depends on our choice for encoding ACGT as respectiveley  00  01 10 11 !!!
        //
        let mut revcomp = self.0 as u64;
        revcomp = !revcomp;
        // then now we have to swap groups of 2 bits
        revcomp = revcomp.reverse_bits();
        revcomp = (revcomp & 0x55555555_55555555) << 1 | (revcomp & 0xAAAAAAAA_AAAAAAAA) >> 1;
        // We have to shift to the right 64-2*nb_bases as useful values are aligned right
        revcomp = revcomp >> (64 - 2 * self.1);
        Kmer64bit(revcomp, self.1)
    }


    fn dump(&self, bufw: &mut dyn io::Write) -> io::Result<usize> {
        bufw.write(unsafe { &mem::transmute::<u8, [u8;1]>(self.1) }).unwrap();
        bufw.write(unsafe { &mem::transmute::<u64, [u8;8]>(self.0) } )
    }    
}  // end of impl KmerT for Kmer64bit





impl CompressedKmerT for Kmer64bit {
    type Val = u64;
    /// This type can store 14 base at max
    fn get_nb_base_max() -> usize { 32 }
    /// a decompressing function mainly for test and debugging purpose
    fn get_uncompressed_kmer(&self) -> Vec<u8> {
        let nb_bases = self.1;
        let alphabet = Alphabet2b::new();
        // we treat each block of 2 bis as u8 end call decoder of Alphabet2b
        let mut decompressed_kmer = Vec::<u8>::with_capacity(nb_bases as usize);
        let mut base:u8;
        //
        let mut buf = self.0;
        // get the base coding part at left end of u32
        buf = buf.rotate_left((64 - 2 * nb_bases) as u32);
        for _ in 0..nb_bases {
            buf = buf.rotate_left(2);
            base = (buf & 0b11) as u8; 
            decompressed_kmer.push(alphabet.decode(base));
        }
        return decompressed_kmer;
    }
    /// return the pure value with part coding number of bases reset to 0.
    #[inline(always)]    
    fn get_compressed_value(&self) -> u64 {
        return self.0;
    }
    ///
    #[inline(always)]    
    fn get_bitsize(&self) -> usize { 64 }
}  // end of impl block of CompressedKmerT for Kmer64bit





impl KmerBuilder<Kmer64bit> for Kmer64bit {
    /// For Kmer64bit the number of bases is encoded in separate field 
    fn build(val: u64, nb_base : u8) -> Kmer64bit {
        Kmer64bit(val, nb_base)
    }
}



impl FromStr for Kmer64bit {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        // check that length is less than 32.
        if s.len() > 32 {
            return Err(String::from("kmer should be less than 32 long"));
        }
        let nb_bases=s.len();
        let mut kmer = Kmer64bit::new(nb_bases as u8);
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
    
}  // end of impl FromStr for Kmer64bit


//==================================================

#[cfg(test)]
mod tests {

#[allow(unused)]
use super::*;

    #[test]
    fn test_reverse_complement_12b_kmer64bit() {
        //
        let to_reverse : Vec<&'static str> = vec![
            // TACG_AGTA_GGAT
            "TACGAGTAGGAT",
            // ACTT_GGAA_CGTT
            "ACTTGGAACGTT"
        ];
        //
        let reversedcomp : Vec<&'static str> = vec! [
            // ATCC_TACT_CGTA
            "ATCCTACTCGTA",
            // AACG_TTCC_AAGT
            "AACGTTCCAAGT"
        ];
        
        for i in 0..to_reverse.len() {
            let revcomp: Kmer64bit = Kmer64bit::from_str(to_reverse[i]).unwrap().reverse_complement();
            let should_be: Kmer64bit = Kmer64bit::from_str(reversedcomp[i]).unwrap();
            if revcomp.0  !=  should_be.0 {
                println!(" i  kmer reversed complement  : {} {:b}", i, revcomp.0);
            }
            assert!(revcomp.0 == should_be.0 );
        }
    } // end of test_reverse_complement_kmer64bit  


#[test]
    fn test_reverse_complement_11b_kmer64bit() {
        //
        let to_reverse : Vec<&'static str> = vec![
            // TACG_AGTA_GGA
            "TACGAGTAGGA",
            // ACTT_GGAA_CGT
            "ACTTGGAACGT"
        ];
        //
        let reversedcomp : Vec<&'static str> = vec! [
            // TCC_TACT_CGTA
            "TCCTACTCGTA",
            // ACG_TTCC_AAGT
            "ACGTTCCAAGT"
        ];
        
        for i in 0..to_reverse.len() {
            let revcomp: Kmer64bit = Kmer64bit::from_str(to_reverse[i]).unwrap().reverse_complement();
            let should_be: Kmer64bit = Kmer64bit::from_str(reversedcomp[i]).unwrap();
            if revcomp.0  !=  should_be.0 {
                println!(" i  kmer reversed complement  : {} {:b}", i, revcomp.0);
            }
            assert_eq!(should_be.get_nb_base(), Kmer64bit::from_str(to_reverse[i]).unwrap().get_nb_base());
            assert!(revcomp.0 == should_be.0 );
        }
    } // end of test_reverse_complement_kmer32bit
}