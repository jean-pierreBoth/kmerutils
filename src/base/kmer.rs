//! module kmer provides Kmer representation using 2bit base encoding, and are stored in 32 or 64 bit unsigned integers
//!
//! In a 32 bit word we can encode a 16-mer. This type is Kmer16b32bit.  
//!   
//! For smaller Kmer, so we use also a 32 bit word but keep 4 (upper) bits to encode size of
//! kmer, it leaves us 28 bits for kmer representation so we can have kmers up to 14 bases. This is type Kmer32bit.  
//!    
//! For larger kmer we use a 2uple(u64,u8) for kmer and number of base.  
//! Rust has now u128 bits so we can go up to 64 bases-kmer with this strategy
//! A bit cumbersome but some space is spared and reverse complement
//! can be efficiently computed with bit symetry. (and soon use reverse instruction!)


#[allow(unused)]
use log::{debug,trace};



// our includes


pub use super::{nthash::*};

pub use super::{kmer32bit::*, kmer16b32bit::*, kmer64bit::*};


/// 9 bytes structure. Cannot implement Trait CompressedKmer
//#[allow(dead_code)]
#[derive(Clone, Copy)]
pub struct KmerCoord {
    // num of read 
    pub read_num : u32,
    // pos in lower
    pub pos : u32,
}


  


//  implementation block of Trait NtHash for Kmers using 2 bits encoded bases
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





#[cfg(test)]
mod tests {

#[allow(unused)]
use super::*;

//


    
  
} // end of mod test
