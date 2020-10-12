//! module kmer provides Kmer representation in 2bit encoding in 32, 64 bit unsigned integers
//!
//! In a 32 bit word we can encode a 16-mer. This type is Kmer16b32bit.
//! But we want also smaller Kmer, so we use also a 32 bit word but keep 4 (upper) bits to encode size of
//! kmer, it leaves us 28 bits for kmer representation so we can have kmers up to 14 bases. This is type Kmer32bit.
//! For larger kmer we use a Tuple(u64,u8) for kmer and number of base.
//! Rust has now u128 bits so we can go up to 64 bases-kmer with this strategy
//! A bit cumbersome but some space is spared and reverse complement
//! can be efficiently computed with bit symetry. (and soon use reverse instruction!)

use std::mem;
use std::io;
use std::cmp::Ordering;
use std::hash::Hash;
use std::cmp::Ord;
use std::str::FromStr;

use log::debug;
use log::trace;

use crate::nthash::*;

// our includes


pub use crate::sequence::*;

/// A Kmer is defined with 2 traits, one for standard operations and
/// one for taking into account that our Kmers are 2 bits encoded
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


/// all our Kmer implements this trait as they use 2 bit base encoding.
/// So basically our Kmer are copy, orderable, representable as hashable orderable compressed value.
/// which is a minimum
pub trait CompressedKmerT : KmerT+Ord+Copy  where Self::Val : Hash+Ord
{
    /// type of compressed value , u16, u32, u64.
    type Val;
    /// returns the max number of base supported by compressing strategy and size of Val
    fn get_nb_base_max() -> usize;
    /// return encoded value in type Val. Possibly in Val we have encoded number of bases and
    /// and kmer value. We get rid of number of base, and return pure encoded value.
    fn get_compressed_value(&self) -> Self::Val;
    /// get Kmer as a Vec<u8>
    fn get_uncompressed_kmer(&self) -> Vec<u8>;
    /// returns the size in bits of word supporting compressed kmer.
    /// not just number of base * size of base
    fn get_bitsize(&self) -> usize;
}  // end of trait CompressedKmer



/// 9 bytes structure. Cannot implement Trait CompressedKmer
//#[allow(dead_code)]
#[derive(Clone, Copy)]
pub struct KmerCoord {
    // num of read 
    pub read_num : u32,
    // pos in lower
    pub pos : u32,
}



//==============================================================================================


/// type supporting 16 bases kmer as u32.

#[derive(Clone,Copy,PartialEq,Eq,PartialOrd,Ord,Debug)]
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
        // we do a symetry as explained in Hacker's delight and complement.
        // Note that we skip the swap between 2 adjacent bits as base are encoded in blocks of 2 bits
        // This depends on our choice for encoding ACGT as respectiveley  00  01 10 11 !!!
        //
        let mut revcomp = self.0 as u32;
        revcomp = !revcomp;
        // the now we have to swap groups of 2 bits
        revcomp = (revcomp & 0x33333333)  <<   2 | (revcomp & 0xCCCCCCCC) >>  2;
        revcomp = (revcomp & 0x0F0F0F0F)  <<   4 | (revcomp & 0xF0F0F0F0) >>  4;
        revcomp = (revcomp & 0x00FF00FF)  <<   8 | (revcomp & 0xFF00FF00) >>  8;
        revcomp = (revcomp & 0x0000FFFF)  <<  16 | (revcomp & 0xFFFF0000) >> 16;
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


    
//===================================================================================

//  mask = 0xF0000000

// possiby we could use a u16 for smaller kmers.


/// This is a type for Kmer less than 14 bases 2 bit encoded
/// Real number of base is encoded in upper 4 bits!!!!

#[derive(Clone,Copy,Debug)]
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
        // We could use as a mask for value field : (0b1 << (2*self.get_nb_bases())) - 1 which enforce 0 bit between 4 upper bits
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

//===================================================================================


// For more than 16 bases the coding of number of bases in a separate field is more efficient
// This representation is consistent with Kmer32bit as self.0 gives the word supporting value of kmer

/// The type supporting Kmer for number of bases between 17 and 32.
#[derive(Clone,Copy,Debug)]
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
        // we do a symetry as explained in Hacker's delight and complement.
        // Note that we skip the swap between 2 adjacent bits as base are encoded in blocks of 2 bits
        // This depends on our choice for encoding ACGT as respectiveley  00  01 10 11 !!!
        //
        let mut revcomp = self.0 as u64;
        revcomp = !revcomp;
        // then now we have to swap groups of 2 bits
        revcomp = (revcomp & 0x33333333_33333333)  <<   2 | (revcomp & 0xCCCCCCCC_CCCCCCCC) >>  2;
        revcomp = (revcomp & 0x0F0F0F0F_0F0F0F0F)  <<   4 | (revcomp & 0xF0F0F0F0_F0F0F0F0) >>  4;
        revcomp = (revcomp & 0x00FF00FF_00FF00FF)  <<   8 | (revcomp & 0xFF00FF00_FF00FF00) >>  8;
        revcomp = (revcomp & 0x0000FFFF_0000FFFF)  <<  16 | (revcomp & 0xFFFF0000_FFFF0000) >> 16;
        revcomp = (revcomp & 0x00000000_FFFFFFFF)  <<  32 | (revcomp & 0xFFFFFFFF_00000000) >> 32;
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
}  // end of impl block for CompressedKmerT







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



//  implementation block of Trait NtHash
//========================================================================================


macro_rules! implement_nthash_for(
    ($ty:ty) => (
        impl  NtHash for  $ty {
            fn nthash_init(&self) -> u64 {
                let mut hval:u64 = 0;
                //
                let ksize = self.get_nb_base() as u32;
                let mut base:u8;
                let mut buf = (*self).0;
                // get the base coding part at left end of u32
                buf = buf.rotate_left(self.get_bitsize() as u32 - 2 * ksize);
                for i in 0..ksize {
                    buf = buf.rotate_left(2);
                    base = (buf & 0b11) as u8;
                    hval ^= base_map!(base,2).rotate_left(ksize-i-1);
                }
                return hval;          
            }
            //
            fn nthash_cycle(&mut self, hashval:u64, new_base:u8) -> u64 {
                let ksize = self.get_nb_base() as u32;
                let buf = (*self).0;
                // get the base coding part at left end of u32
                let old_base = (buf.rotate_left((self.get_bitsize() as u32 - 2 * ksize) + 2) & 0b11) as u8;
                self.push(new_base);
                //
                hashval.rotate_left(1) ^ base_map!(old_base,2).rotate_left(ksize as u32) ^ base_map!(new_base,2)
            }
            //    
            fn nthash_canonical_init(&self, fhash : &mut u64, rhash : &mut u64) -> (u64, u8) {
                *fhash = 0;
                *rhash = 0;
                //
                let ksize = self.get_nb_base() as u32;
                let mut base;
                let mut buf = (*self).0;
                buf = buf.rotate_left(self.get_bitsize() as u32 - 2 * ksize);
                for i in 0..ksize {
                    buf = buf.rotate_left(2);
                    base = (buf & 0b11) as u8;
                    (*fhash) = (*fhash) ^ (base_map!(base as usize, 2).rotate_left((ksize-i-1) as u32));
                    (*rhash) = (*rhash) ^ (base_map_complement!(base as usize, 2).rotate_left(i as u32));                            
                }
                if fhash <= rhash {
                    return (*fhash,0);
                }
                else {
                    return (*rhash,1);
                }
            }
            //
            fn nthash_canonical_cycle(&mut self, new_base:u8, fhash : &mut u64, rhash : &mut u64) -> (u64, u8) {
                *fhash = 0;
                *rhash = 0;
                //
                let ksize = self.get_nb_base() as u32;
                let buf = (*self).0;
                // get the base coding part at left end of u32
                let old_base = (buf.rotate_left((self.get_bitsize() as u32 - 2 * ksize) + 2) & 0b11) as u8;
                //
                *fhash = (*fhash).rotate_left(1) ^ base_map!(old_base,2).rotate_left(ksize as u32) ^ base_map!(new_base,2);
                *rhash = (*rhash).rotate_right(1) ^ (base_map_complement!(old_base as usize,2).rotate_left(ksize)) ^
                    (base_map_complement!(new_base as usize,2).rotate_left(ksize-1));
                // we must update to get old_base OK in future calls 
                self.push(new_base);
                //
                if fhash <= rhash {
                    return (*fhash,0);
                }
                else {
                    return (*rhash,1);
                }
            } // nthash_canonical_cycle

           
            fn nthash_mult_canonical_init(&self, fhash : &mut u64, rhash : &mut u64, hashed : &mut [u64]) -> u8 {
                *fhash = 0;
                *rhash = 0;
                //
                let res_hash = self.nthash_canonical_init(fhash, rhash);
                hashed[0] = res_hash.0;
                from_one_hash_val_to_mult_hash(self.get_nb_base() as u64 , hashed);
                // we return strand of minimum
                return res_hash.1;
            }  // end of nthash_mult_canonical_init
           
            fn nthash_mult_canonical_cycle(&mut self, new_base:u8, fhash : &mut u64, rhash : &mut u64, hashed : &mut [u64]) -> u8 {
                //
                let res_hash = self.nthash_canonical_cycle(new_base, fhash, rhash);
                hashed[0] = res_hash.0;
                from_one_hash_val_to_mult_hash(self.get_nb_base() as u64, hashed);
                //
                return res_hash.1;    
            } // end of nthash_mult_canonical_cycle
            
        }  // end of impl NtHash for Kmer32bit
    ) // end of match ty
);  // end of macro implement_nthash_for

implement_nthash_for!(Kmer32bit);
implement_nthash_for!(Kmer16b32bit);

//====================================================================================


// here we are suppose to know which read we refer to
// 4 bytes structure
#[allow(dead_code)]
pub struct KmerMinimizer {
    /// first bit is strand, then 7 bits for length of kmer, last 3 bytes for pos in read
    id : u32,  
}

/// 9 bytes structure
#[allow(dead_code)]
pub struct KmerHashed {
    /// hashed value created by invertible hash hash_64
    hashed : u64,
    /// nb base
    size : u8,
}

//===========================================================


////////////////////////////////////////////////////////////////////////////////////////////////////



#[cfg(test)]
mod tests {
    use super::*;
    //
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


    #[test]
    fn test_reverse_complement_kmer64bit() {
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
} // end of mod test
