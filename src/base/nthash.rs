
////////////////////////////////////////////////////////////////////////////////////////////////
//! ntHash provides a recursive DNA kmer hashing:  adapted from Mohamadi Chu Birol BioInformatics 2016
//!
//! Kmer will be compressed most of the time, but possibly not...
//! so we must adapt .hpp code. We have seed initialization for 2/4/8 bit encoding for bases and complement base.


// We keep calls to rotate as it should be at least as fast as access to tables that
// are not guaranteed to be in cache, and rust rotate should call intel intel instruction

// shift for gerenerating multiple hash values
pub const MULTISHIFT : usize = 27;

// seed for gerenerating multiple hash values
pub const MULTISEED : u64 = 0x90b45d39fb6da1fa;

// 64-bit random SEED_ corresponding to bases and their complements

const SEED_A : u64 = 0x3c8bfbb395c60474;
const SEED_C : u64 = 0x3193c18562a02b4c;
const SEED_G : u64 = 0x20323ed082572324;
const SEED_T : u64 = 0x295549f54be24456;


// For 2bit encoded base we could get complement of base b by ^b&0b11 .... a bit longer than an array access ?
#[allow(dead_code)]
pub const OFFSET_COMP_2B : usize = 4;

// with an offset of OFFSET_COMP_2B  we get complementary base for a 2 bit encoded base b by accessing to BASE_MAPPING_2B[b+OFFSET_COMP_2B]
#[allow(dead_code)]
pub const BASE_MAPPING_2B: [u64;8] = [SEED_A, SEED_C, SEED_G, SEED_T,
   /* complement base */          SEED_T, SEED_G, SEED_C, SEED_A,
];


// with an offset of 8 we get complementary base for 4 bit encoded bases.
#[allow(dead_code)]
pub const OFFSET_COMP_4B : usize = 8;

// with an offset of OFFSET_COMP_4B  we get complementary base for a 4 bit encoded base b by accessing to BASE_MAPPING_4B[b+OFFSET_COMP_4B]
#[allow(dead_code)]
pub const BASE_MAPPING_4B: [u64;18] = [0, SEED_A, SEED_C, 0, SEED_G, 0, 0, 0, SEED_T, 
          /* complement base */    0, SEED_T, SEED_G, 0, SEED_T, 0, 0, 0, SEED_A
];




// for index less than 65 , 67, 81 and 84 we get complement base by adding offset 24

pub const OFFSET_COMP_8B : usize = 24;

// with an offset of OFFSET_COMP_8B  we get complementary base for a 4 bit encoded base b by accessing to BASE_MAPPING_8B[b+OFFSET_COMP_8B]
pub const BASE_MAPPING_8B: [u64;112] = [
    0, 0 ,0, 0,          0, 0, 0, 0,
    0, 0, 0, 0,          0, 0, 0, 0,
    0, 0, 0, 0,          0, 0, 0, 0,
    0, 0, 0, 0,          0, 0, 0, 0,

    0, 0, 0, 0,          0, 0, 0, 0,
    0, 0, 0, 0,          0, 0, 0, 0,
    0, 0, 0, 0,          0, 0, 0, 0,
    0, 0, 0, 0,          0, 0, 0, 0,

    0, SEED_A, 0, SEED_C,  0, 0, 0, 0,              /*   64..72  : A is 65, C is 67, G is 81 , and T is 84*/
    0,  0, 0, 0,           0, 0, 0, 0,              /*   72..80 */
    SEED_G, 0, 0, SEED_T,  0, 0, 0, 0,              /*   81..88 */

    0, SEED_T, 0, SEED_G,  0, 0, 0, 0,              /*  same block as the precedent one with bases complemented +24 translated */
    0,  0, 0, 0,           0, 0, 0, 0,          
    SEED_C, 0, 0, SEED_A,  0, 0, 0, 0,     
];



/// This function takes as arg a kmer size, an initial hval of kmer in hashed[0]
/// and computes more hash val in hashes[1..hashed_len]
/// it is part of nthash but could be used with inversible hash

#[inline(always)]
pub fn from_one_hash_val_to_mult_hash(ksize: u64, hashed: &mut [u64]) {
    let mut tmp_h : u64;
    //
    for i in 1..hashed.len() {
        // ^ comes after * but it is better showing it
        tmp_h = hashed[0] * (i as u64 ^ (ksize * MULTISEED)) as u64;
        tmp_h ^= tmp_h >> MULTISHIFT;
        hashed[i] = tmp_h;
    }  
}




/// trait describing nthash cyclic hashing.
// Will be implemented by str, and all kmer types

pub trait NtHash {
    /// first initialization hash value of a kmer type to hash
    fn nthash_init(&self) -> u64;
    /// cyclic hash, takes as arg previous hash value as hashval and new_base inserted at right of the kmer. 
    fn nthash_cycle(&mut self, hashval:u64, new_base:u8) -> u64;
    //
    /// function for initializing  canonical hash of a kmer (the minimum of hash du kmer et du kmer complement)
    /// It takes as arg the kmer to hash , returns forward and reverse hash value by reference
    /// the 2 uple((min(fhash, rhash), strand) with side = 0 fhash is the minimum 1 if rhash is the min
    fn nthash_canonical_init(&self, fhash : &mut u64, rhash : &mut u64) -> (u64, u8);
    //
    /// function for computing hash value of subsequent kmers by getting out old_base (leftmost) and inserting
    /// new_base at right end. It takes as arg the kmer size, old and new , forward and reverse hash
    /// previously computed. 
    /// It returns minimal value of hash and strand as init method and updated values of fhash (forward)
    /// and rhash (reverse complement)
    fn nthash_canonical_cycle(&mut self, new_base:u8, fhash : &mut u64, rhash : &mut u64) -> (u64, u8);
    //
    /// a method to compute initial multiple canonical sketches of a kmer.
    /// returns by reference in slice hashed the hashed.len() sketches asked and as value the strand where the
    /// sketches are taken from. 0 for forward and 1 for reverse.
    fn nthash_mult_canonical_init(&self, fhash : &mut u64, rhash : &mut u64, hashed : &mut [u64]) -> u8;
    //
    /// a method to cylcle on canonical multiple sketches of a initial kmer
    /// args are size of kmer , old_base getting out (left) of the kmer and new_base getting in at right of kmer.
    /// returns by reference in slice hashed the hashed.len() sketches asked and as value the strand where the
    /// sketches are taken from. 0 for forward and 1 for reverse.
    fn nthash_mult_canonical_cycle(&mut self, new_base:u8, fhash : &mut u64, rhash : &mut u64, hashed : &mut [u64]) -> u8;
}




///////////////////////////////////////////////////////////////////////////////////////////////


// macros to map a base to its seed
macro_rules! base_map(
    ($c:expr,$nb_bits:expr) => {
        match $nb_bits {
            8 => (BASE_MAPPING_8B[$c as usize] as u64),
            2 => (BASE_MAPPING_2B[$c as usize] as u64),
            4 => (BASE_MAPPING_4B[$c as usize] as u64),
            _ => panic!("cannot map base in nthash"),
        }
    }
);

// macro to map a base to its complement seed
macro_rules! base_map_complement(
    ($c:expr,$nb_bits:expr) => {
        match $nb_bits {
            8 => (BASE_MAPPING_8B[OFFSET_COMP_8B + $c as usize] as u64),
            2 => (BASE_MAPPING_2B[OFFSET_COMP_2B + $c as usize] as u64),
            4 => (BASE_MAPPING_4B[OFFSET_COMP_4B + $c as usize] as u64),
            _ => panic!("cannot map base in nthash"),
        }
    }
);


///////////////////////////////////////////////////////////////////////////////////////////////

/// nthash for uncompressed kmer, 8bit/base
/// This function does the first initialisation of hash value for a kmer passed as argument
//  It corresponds to getFhval in Chu-Birol

pub fn nthash_init_8b(kmer : &[u8]) -> u64 {
    let mut hval:u64 = 0;
    let ksize = kmer.len();
    for i in 0..ksize {
        //        hval ^= (BASE_MAPPING_8B[kmer[i] as usize] as u64).rotate_left((ksize- i-1) as u32 %64);
        hval ^= base_map!(kmer[i],8).rotate_left((ksize -i -1) as u32 %64);
    }
    return hval;   
} // end of nthash_u8


/// cyclic hashing of kmers.
/// We give as arguments : hash value of current kmer ,
///                        size of kmer the old_base and
///                        new base in 8bit encoding to be inserted
// at right end of kmer to get the new corresponding hash value. The old base (at leftmost byte is at byte ksize)
// hence the rotate arg
// It corresponds to NT64(const uint64_t fhVal, const unsigned char charOut, const unsigned char charIn, const unsigned k)
// in original Chu-Birol



#[inline(always)]
pub fn nthash_cycle_8b(hashval: u64, ksize: usize, old_base:u8, new_base:u8) -> u64 {
    hashval.rotate_left(1) ^ BASE_MAPPING_8B[old_base as usize].rotate_left((ksize%64) as u32) ^  BASE_MAPPING_8B[new_base as usize]
}



// now to hash the reverse complemnt


/// function to initialize reverse complement hash
//  It corresponds to getRhval in original ntHash Chu-Birol

pub fn nthash_rcomp_init_8b(kmer : &[u8]) -> u64 {
    let mut hval:u64 = 0;
    let ksize = kmer.len();
    for i in 0..ksize {
        hval ^= (BASE_MAPPING_8B[OFFSET_COMP_8B + kmer[i] as usize] as u64).rotate_left(i as u32 %64);
    }
    return hval;   
} // end of nthash_rcomp_init_8b


/// function to do cyclic hashing of reverse complement kmer knowing preceding hashval
/// We give as arguments :
///                     - hash value of reverse complement of current kmer ,
///                     - size of kmer
///                     - the old_base and
///                     - new base in 8bit encoding to be inserted

#[inline(always)]
pub fn nthash_rcomp_cycle_8b(hashval: u64, ksize: usize, old_base:u8, new_base:u8) -> u64 {
    hashval.rotate_right(1) ^ BASE_MAPPING_8B[OFFSET_COMP_8B + old_base as usize].rotate_right(1) ^
        BASE_MAPPING_8B[OFFSET_COMP_8B + new_base as usize].rotate_left((ksize -1) as u32 %64)
}




// canonical ntHash for one hash function
// it corresponds to  NTC64(const char * kmerSeq, const unsigned k, uint64_t& fhVal, uint64_t& rhVal)

/// computes initial canonical value for a kmer given in 8bit representation
/// args :
///        . kmer to hash
///        . ref to forward hash value
///        . ref to reverse hash value
///        . returns a 2-uple giving first the minimum of forward and reverse hash value and in second
///          value 0 if the minimum was in forward mode 1 in reverse.


pub fn nthash_canonical_init_8b(kmer : &[u8], fhash : &mut u64, rhash : &mut u64) -> (u64, u8) {
    *fhash = 0;
    *rhash = 0;
    //
    let ksize = kmer.len();
    for i in 0..ksize {
        (*fhash) = (*fhash) ^ (BASE_MAPPING_8B[kmer[i] as usize] as u64).rotate_left((ksize-i-1) as u32 %64);
        (*rhash) = (*rhash) ^ (BASE_MAPPING_8B[OFFSET_COMP_8B + kmer[i] as usize] as u64).rotate_left(i as u32 %64);
    }
    if fhash <= rhash {
        return (*fhash,0);
    }
    else {
        return (*rhash,1);
    }
} // end of nthash_canonical_init_8b



/// function for canonical hashing with 8bit base encoding

//  NTC64(const unsigned char charOut, const unsigned char charIn, const unsigned k, uint64_t& fhVal, uint64_t& rhVal)

pub fn nthash_canonical_cycle_8b(ksize: usize, old_base:u8, new_base:u8, fhash: &mut u64, rhash: &mut u64)  -> (u64, u8) {
    *fhash = nthash_cycle_8b(*fhash, ksize, old_base, new_base);
    *rhash = nthash_rcomp_cycle_8b(*rhash, ksize, old_base, new_base);
    if fhash <= rhash {
        return (*fhash,0);
    }
    else {
        return (*rhash,1);
    }
} // end of nthash_canonical_cycle_8b



//        function for canonical multiple hashing
//       =========================================



// Case where we need multiple canonical hash values. For example for (super)minhash with many hash function
// It is the most useful method, in fact for use in Superminhash. 
// Cf function NTMC64 in nthash.hpp  Mohamadi Chu Birol BioInformatics 2016



pub fn nthash_mult_canonical_init_8b(kmer : &[u8], fhash : &mut u64, rhash : &mut u64, hashed : &mut [u64]) -> u8 {
    *fhash = 0;
    *rhash = 0;
    let res_hash = nthash_canonical_init_8b(kmer, fhash, rhash);
    hashed[0] = res_hash.0;
    //
    from_one_hash_val_to_mult_hash(kmer.len() as u64, hashed);
    // we return strand of minimum
    return res_hash.1;
} // end of nthash_canonical_init_8b




pub fn nthash_mult_canonical_cycle_8b(ksize : usize, old_base:u8, new_base:u8, fhash : &mut u64, rhash : &mut u64, hashed : &mut [u64]) -> u8 {
    *fhash = 0;
    *rhash = 0;
    let res_hash = nthash_canonical_cycle_8b(ksize, old_base, new_base, fhash, rhash);
    hashed[0] = res_hash.0;
    //
    from_one_hash_val_to_mult_hash(ksize as u64, hashed);
    //
    return res_hash.1;
} // end of nthash_canonical_init_8b

///////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_nthash_simple_16bases() {
       // got a string 
        let seqstr = String::from("TCAAAGGGAAACATTCAAAATCAGTATGCGCCCGTTCAGTTACGTATTGCTCTCGCTAATGAGATGGGCTGGGTACAGAG");
        // transform to a Vec<u8> and then a Sequence
        let slu8 = seqstr.as_bytes();
        let kmer_size:usize = 16;
        let kmer_begin = 0;
        let mut kmer = &slu8[kmer_begin..kmer_size];
        println!("kmer initial {}" , String::from_utf8_lossy(kmer));
        let mut hval = nthash_init_8b(&kmer);
        for i in 1..seqstr.len() - kmer_size {
            println!("removed base {}, inserted base {}", slu8[i-1], slu8[i-1+kmer_size]);
            hval = nthash_cycle_8b(hval, kmer_size, slu8[i-1], slu8[i-1+kmer_size]);
            kmer = &slu8[i..i+kmer_size];
            let hash_cycled = nthash_init_8b(&kmer);
            println!(" i kmer hasval hashval_cycle {}  {}  {}  {}", i, String::from_utf8_lossy(kmer) , hval, hash_cycled);
            assert_eq!(hval, hash_cycled);
        }       
    }  // end of test_nthash_simple_16bases


    #[test]   
    fn test_nthash_canonical_16bases() {
       // got a string 
        let seqstr = String::from("TCAAAGGGAAACATTCAAAATCAGTATGCGCCCGTTCAGTTACGTATTGCTCTCGCTAATGAGATGGGCTGGGTACAGAG");
        // transform to a Vec<u8> and then a Sequence
        let slu8 = seqstr.as_bytes();
        let kmer_size:usize = 16;
        let kmer_begin = 0;
        let mut kmer = &slu8[kmer_begin..kmer_size];
        let mut fhash:u64 = 0;
        let mut rhash:u64 = 0;
        println!("kmer initial {}" , String::from_utf8_lossy(kmer));
        #[allow(unused_assignments)]
        let mut hval = nthash_canonical_init_8b(&kmer, &mut fhash, &mut rhash);
        for i in 1..seqstr.len() - kmer_size {
            println!("removed base {}, inserted base {}", slu8[i-1], slu8[i-1+kmer_size]);
            hval = nthash_canonical_cycle_8b(kmer_size, slu8[i-1], slu8[i-1+kmer_size], &mut fhash, &mut rhash);
            kmer = &slu8[i..i+kmer_size];
            let mut fhash_check:u64 = 0;
            let mut rhash_check:u64 = 0;
            let hash_cycled = nthash_canonical_init_8b(&kmer, &mut fhash_check, &mut rhash_check);
            println!(" i kmer hasval hashval_cycle {}  {}  {}  {}", i, String::from_utf8_lossy(kmer) , hval.0, hash_cycled.0);
            assert_eq!(hval, hash_cycled);
        }       
    }  // end of test_nthash_canonical_16bases


} // end of mod tests
