//!  Kmer32bit 
//! 

use std::mem;
use std::io;
use std::cmp::Ordering;
use std::cmp::Ord;
use std::str::FromStr;

#[allow(unused)]
use log::{debug,trace};



// our includes

pub use super::kmertraits::*;
pub use super::{nthash::*, alphabet::*};

/// This is a type for Kmer less than 14 bases 2 bit encoded
/// Real number of base is encoded in upper 4 bits!!!!

#[derive(Clone,Copy,Debug, Hash)]
pub struct Kmer32bit(pub u32);

// we need to implement PartialEq, Eq PartialOrd and Ord for Kmer32bit.
// Note that as we have been very cautious in push method to not let garbage bits between upper fields used for coding
// number of bases and lower field for kmer coding, we can be brutal now and just test value
// with a mask equal to 0x0FFFFFFF Cf commentary in push method

impl PartialEq for Kmer32bit {
    // we must check number of bases and equality of field
    fn eq(&self, other: &Kmer32bit) -> bool {
        if self.0 & 0xF0000000 == other.0 & 0xF0000000 {
            // now we can test equality of kmer value part.
            // we have to check the lower 2*nb_bases bits of u32 word.
            if self.0 & 0x0FFFFFFF == other.0 & 0x0FFFFFFF {
                return true;
            }
            else {
                return false;
            }
        }
        return false;
    } // end of eq
}  // end of PartialEq implementation


// read the doc for Eq ....
impl Eq for Kmer32bit{}


/// We define ordering as a kind of "lexicographic" order by taking into account first number of base.
/// The more the number of base the greater. Then we have integer comparison between lower kmer part
/// which corresponds to lexicographic order as we have A < C < G < T in 2bit and 4bit encoding
impl Ord for Kmer32bit {
    fn cmp(&self, other: &Kmer32bit) -> Ordering {
        if self.0 & 0xF0000000 != other.0 & 0xF0000000 {
            return (self.0 & 0xF0000000).cmp(&(other.0 & 0xF0000000));
        }
        else {
            return (self.0 & 0x0FFFFFFF).cmp(&(other.0 & 0x0FFFFFFF));
        }
    } // end cmp
} // end impl Ord for Kmer32bit


impl PartialOrd for Kmer32bit {
    fn partial_cmp(&self, other: &Kmer32bit) -> Option<Ordering> {
        Some(self.cmp(other))
    } // end partial_cmp
} // end impl Ord for Kmer32bit



// now we can do the real job on Kmer32bit
//========================================


impl Kmer32bit {
    /// allocate a new kmer. Takes the number of base it will store as argument.
    pub fn new(nb_bases: u8) -> Kmer32bit {
        if nb_bases >= 15 {
            panic!("Kmer32bit cannot store more than 14 bases");
        }
        // insert 4 lower bits of nb_bases in 4 upper bytes of a u32
        let mut kmer = 0u32;
        kmer = kmer | ((nb_bases as u32) << 28);
        Kmer32bit(kmer)
    } // end of new

    /// fills in upper bits the number of bases this kmer will store
    pub fn set_nb_base(&mut self, nb_bases : u8) {
        if nb_bases >= 15 {
            panic!("Kmer32bit cannot store more than 14 bases");
        }
        // clean upper 4 bits
        self.0 = self.0 & 0x0FFFFFFF;
        // set upper bits
        self.0 = self.0 | ((nb_bases as u32) << 28)               
    } // end of  set_nb_base
}  // end of impl for 



impl KmerT for Kmer32bit {
    /// retrieves the number of bases this kmer stores
    #[inline(always)]
    fn get_nb_base(&self) -> u8 {
        (self.0.rotate_left(4) & 0b1111) as u8
    }
    /// push a base (2bits) at right end of kmer producing a new Kmer
    /// So arg base must be 2bit encoded !!! and we do not check that so be careful
    fn push(&self, base : u8) -> Kmer32bit {
        // store upper 4 bits
        let nb_bases_mask: u32 = self.0 & 0xF0000000;
        // shift left 2 bits, insert new base and enforce saved number of bases.
        // Note that the left shift pushes garbage bit between 4 upper bits and lower bits coding value.
        // We are cautious to clean them by using value_mask which reset those bit to 0!
        // It is useful when implementing PartialEq and compressed value 
        // We use a mask for value field : (0b1 << (2*self.get_nb_bases())) - 1 which enforce 0 bit between 4 upper bits
        // and lower bits coding value.
        let value_mask :u32 = (0b1 << (2*self.get_nb_base())) - 1;
        let mut new_kmer = ((self.0 << 2) & value_mask) | (base as u32 & 0b11);
        new_kmer = new_kmer | nb_bases_mask;
        // and enforce saved number of bases.
        trace!("in push new_kmer = {:b}",  new_kmer);
        Kmer32bit(new_kmer)
    }  // end of push

    /// just returns the reversed complement 16 bases kmer in 2bit encoding

    //  we can build upon the method used for Kmer16b32bit and then do the correct shift
    //  to get bases in the correct right part of a u32 and reset nb_bases in 4 upper bits.
    // 
    fn reverse_complement(&self) ->  Kmer32bit {
        let nb_bases_mask : u32 = self.0 & 0xF0000000;
        //
        // we do a symetry as explained in hacker's delight and complement.
        // Note that we skip the swap between 2 adjacent bits as base are encoded in blocks of 2 bits
        // This depends on our choice for encoding ACGT as respectiveley  00  01 10 11 !!!
        //
        let mut revcomp = self.0;
        revcomp = !revcomp;
        // the now we have to swap groups of 2 bits
        revcomp = (revcomp & 0x33333333)  <<   2 | (revcomp & 0xCCCCCCCC) >>  2;
        revcomp = (revcomp & 0x0F0F0F0F)  <<   4 | (revcomp & 0xF0F0F0F0) >>  4;
        revcomp = (revcomp & 0x00FF00FF)  <<   8 | (revcomp & 0xFF00FF00) >>  8;
        revcomp = (revcomp & 0x0000FFFF)  <<  16 | (revcomp & 0xFFFF0000) >> 16;
        // We have to shift to the right 32-2*nb_bases
        revcomp = revcomp >> (32 - 2 * nb_bases_mask.rotate_left(4));
        // and enforce nb_bases in upper 4 bits
        revcomp = (revcomp & 0x0FFFFFFF) | nb_bases_mask;
        //        
        Kmer32bit(revcomp)
    }
    /// we just do a raw write. Error prone when reloading. Any dump file must have a header
    /// describing number of bases! to distinguish from Kmer16b32bit
    fn dump(&self, bufw: &mut dyn io::Write) -> io::Result<usize> {
        bufw.write(unsafe { &mem::transmute::<u32, [u8;4]>((*self).0) } )
    }
} // end of impl KmerT for Kmer32bit



impl CompressedKmerT for Kmer32bit {
    type Val = u32;
    /// This type can store 14 base at max
    fn get_nb_base_max() -> usize { 14 }
    /// a decompressing function mainly for test and debugging purpose
    fn get_uncompressed_kmer(&self) -> Vec<u8> {
        let nb_bases = (self.0.rotate_left(4) & 0b1111) as usize;
        let alphabet = Alphabet2b::new();
        // we treat each block of 2 bis as u8 end call decoder of Alphabet2b
        let mut decompressed_kmer = Vec::<u8>::with_capacity(nb_bases);
        let mut base:u8;
        //
        let mut buf = self.0;
        // get the base coding part at left end of u32
        buf = buf.rotate_left((self.get_bitsize() - 2 * nb_bases) as u32);
        for _ in 0..nb_bases {
            buf = buf.rotate_left(2);
            base = (buf & 0b11) as u8; 
            decompressed_kmer.push(alphabet.decode(base));
        }
        return decompressed_kmer;
    }
    /// return the pure value with part coding number of bases reset to 0.
    #[inline(always)]    
    fn get_compressed_value(&self) -> u32 {
        // possibly be careful for garbage bits between 4 upper bits of word coding for nb_bases
        // and lower part coding value. We could compute a real mask by using nb_bases * 2 bits lower bits set 1.
        // BUT kmer is initialized to 0 and we have been careful in push method (See comments)
        return self.0 & 0x0FFFFFFF;
    }
    ///
    #[inline(always)]    
    fn get_bitsize(&self) -> usize { 32 }
}  // end of impl block for CompressedKmerT




 
impl FromStr for Kmer32bit {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        // check that length is less than 14.
        if s.len() > 14 {
            return Err(String::from("too long kmer"));
        }
        let nb_bases=s.len();
        let mut kmer = Kmer32bit::new(nb_bases as u8);
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
    
}  // end of impl FromStr for Kmer32bit


mod tests {

    #[allow(unused)]
    use super::*;

    #[test]
    fn test_reverse_complement_kmer32bit() {
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
            let revcomp: Kmer32bit = Kmer32bit::from_str(to_reverse[i]).unwrap().reverse_complement();
            let should_be: Kmer32bit = Kmer32bit::from_str(reversedcomp[i]).unwrap();
            if revcomp.0  !=  should_be.0 {
                println!(" i  kmer reversed complement  : {} {:b}", i, revcomp.0);
            }
            assert!(revcomp.0 == should_be.0 );
        }
    } // end of test_reverse_complement_kmer32bit


    #[test]
    fn test_ordandeq_kmer32() {
        let to_compare : Vec<&'static str> = vec![
            // TACG_AGTA_GGAT
            "TACGAGTAGGAT",
            // ACTT_GGAA_CGTT
            "ACTTGGAACGTT",
            // TACG_AGTA_GGAT,
            "TACGAGTAGGAT"
        ];
        let mut kmer_vect:Vec<Kmer32bit> = Vec::with_capacity(to_compare.len());
        for i in 0..to_compare.len() {
            let kmer: Kmer32bit = Kmer32bit::from_str(to_compare[i]).unwrap();
            trace!("  kmer i  : {} {:b}", i, kmer.0);
            kmer_vect.push(kmer);
        }
        // now some comparisons
        assert!(kmer_vect[0] == kmer_vect[2]);
        assert!(kmer_vect[0] > kmer_vect[1]);
    } // end of test_OrdandEq_kmer32


}