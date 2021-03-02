//! This module describes encoding of bases in 4 or 2 bits. 
//! Alphabet2b is used for compressing see [`super::kmer::CompressedKmerT`]
//! It provides compressing/decompressing utilities and iterator over a compressed sequence
//! of base




/// check if a base is ACGT
pub fn is_acgt(c : u8) -> bool {
    match c {
        b'A' | b'C' | b'G' | b'T' => true,
        _    => false,
    }
}


/// return lower conjugate from CG
pub fn get_ac_from_tg(c : u8) -> u8 {
    match c {
        b'T' => b'A',
        b'G' => b'C',
        _    => c,
    }
}

pub fn count_non_acgt(seq : &[u8]) -> usize {
    let non_acgt = seq.into_iter().fold(0u32, |acc, &b| acc+ !is_acgt(b) as u32);
    return non_acgt as usize;
}


/// We provide here small (and fast) utilities for encoding restructed DNA (A,T,C,G)
/// alphabet to 2 or 4 bits (or 8 bits) and utilities to manage corresponding sequences.
///
/// For a full generality need fall back to bio::alphabets (which use vec_map ...)
/// but our implementation is 2 or 3 times faster (and is not so well abstracted)

/// this trait provides basic encode/decode methods. methods panics if patterns are not in the alphabet.
pub trait BaseCompress {
    /// get a code on a reduced number of bits depending on implementation . Puts 1 in higher bit if fail
    fn encode(&self, c:u8) -> u8;
    /// given a bit pattern returns u8 giving char. Puts 1 in higher bit if fail
    fn decode(&self, c:u8) -> u8;
    /// get number bits used in encoding
    fn get_nb_bits(&self) -> u8;
    /// is char on the alphabet we can compress
    fn is_valid_base(&self, c: u8) -> bool;
    /// complement base
    fn complement(&self, c:u8) -> u8;
    /// pack a slice of 2 or 4 bases (depending on get_nb_bits) in a u8
    fn base_pack(&self, to_pack:&[u8]) -> u8;
    // pack a slice of 2 or 4 bases (depending on get_nb_bits) in a u8
//    fn base_unpack(&self, to_pack:&[u8]) -> u8;
}


//  Alphabet_2b


/// this structure compress to 2 bits the 4 bases ACGT
/// - A maps to 0b00
/// - C maps to 0b01
/// - G maps to 0b10
/// - T maps to 0b11
/// note : the lexicographic order is preserved and bases are conjugated

pub struct Alphabet2b {
    pub bases: String,
}


impl Alphabet2b {
    pub fn new() -> Alphabet2b {
        Alphabet2b { bases : String::from("ACGT")}
    }
    //
    pub fn len(&self) -> u8 {
        return 2;
    }
    //
    /// This function fills decompressed bases in slice unpacked.
    /// The slice must be larger than number of bases to be returned i.e 4 for Alphabet2b
    /// This mode to get the result is not the niciest but the fastest for repetitive call
    /// The caller must be aware that the byte packed must be fully initialized with consistent
    /// For non full bytes some bits can be garbage decode will return anything.
    /// It is the reason why in Sequence::new we fill the byte with a decodable pattern.
    //
    // Morally only some one that has consistently packed data can call unpack. Ugly but faster.
    // We could return a u32 but the caller would have to do decoding job
    //
    pub fn base_unpack(&self, packed : u8, unpacked : &mut [u8])  {
        // get each nibble and decode. just symetric as pack
        unpacked[3] = packed & 0b0011;
        unpacked[2] = (packed >> 2) & 0b0011;
        let packed2 = packed >> 4;
        unpacked[1] = packed2 & 0b0011;
        unpacked[0] = packed2 >> 2 & 0b11;

        unpacked[0] = self.decode(unpacked[0]);
        unpacked[1] = self.decode(unpacked[1]);
        unpacked[2] = self.decode(unpacked[2]);
        unpacked[3] = self.decode(unpacked[3]);
    }

    #[inline]
    pub fn nb_invalid_bases(&self, seq : &[u8]) -> u32 {
        let sum = seq.iter().fold(0u32, |acc, &b| acc+ !self.is_valid_base(b) as u32);
        return sum;
    }

}  // end impl Alphabet2b



impl BaseCompress for Alphabet2b {

    #[inline(always)]
    fn encode(&self, c:u8) -> u8 {
        match c {
            b'A' => 0b00,
            b'C' => 0b01,
            b'G' => 0b10,
            b'T' => 0b11,
            _    => panic!("pattern not a code in alpahabet_2b"),
        }
    } // end of function encode

    #[inline(always)]
    fn decode(&self, c:u8) -> u8 {
        match c {
            0b00 => b'A',
            0b01 => b'C',
            0b10 => b'G',
            0b11 => b'T',
            _    => panic!("pattern not a code in alpahabet_2b"),
        }
   }  // end of decode

    /// return base complement
    fn complement(&self, c:u8) -> u8 {
        match c {
            0b00 => 0b11,
            0b01 => 0b10,
            0b10 => 0b01,
            0b11 => 0b00,
            _    => panic!("pattern not a code in alpahabet_2b"),
        }
    } // end of complement

    #[inline(always)]
    fn  get_nb_bits(&self) -> u8 {
        return 2;
    }

    #[inline(always)]
    fn is_valid_base(&self, c: u8) -> bool {
        match c {
            b'A' | b'C' | b'G' | b'T' => true,
            _    => false,
        }
    } // end is_valid_base

    // we expect a slice of size at least 4 bytes
    fn base_pack(&self, to_pack: &[u8]) -> u8 {
        let mut packed = 0u8;
        for i in 0..4 {
            packed = packed | (self.encode(to_pack[i]) << 6 - 2 * i);
        }
        return packed;
    }

} // end implement section for Alphabet2b



//  Alphabet4b


/// this structure compress to 4 bits the 4 bases ACGT  
/// A maps to 0b0001 = 1 = 0x01  
/// C maps to 0b0010 = 2 = 0x02  
/// G maps to 0b0100 = 4 = 0x04  
/// T maps to 0b1000 = 8 = 0x08  

// note : the lexicographic order is preserved and bases are NOT conjugated
//         we can convert from Alphabet4b to Alphabet2b by counting (calling)  trailing_zeros
//         and converting form  Alphabet2b to Alphabet_4b by shifting
//         possibly we could encode A | C and so on

pub struct Alphabet4b {
    pub bases: String,
}


impl Alphabet4b {
    pub fn new() -> Alphabet4b {
        Alphabet4b { bases : String::from("ACGT")}
    }
    //
    pub fn len(&self) -> usize {
        return 4;
    }
    //
    /// This function fills decompressed bases in slice unpacked.
    /// The slice must be larger than number of bases to be returned i.e 2 for Alphabet4b
    /// This mode to get the result is not the niciest but the fastest for repetitive call
    /// The caller must be aware that the byte packed must be fully initialized with consistent pattern
    /// For non full bytes some bits can be garbage and decode will fail. Ugly but faster.
    pub fn base_unpack(&self, packed : u8, unpacked : &mut [u8])  {
        // get each nibble
//        println!(" received = {} " , packed);
        unpacked[1] = packed & 0x0F;
        unpacked[0] = (packed >> 4) & 0x0F;
        // as pack method does the encoding , decoding is done here also
        unpacked[0] = self.decode(unpacked[0]);
        unpacked[1] = self.decode(unpacked[1]);
        //
//        println!(" unpacked[0] = {} " , unpacked[0]);
//        println!(" unpacked[1] = {} " , unpacked[1]);
        return;
    }
    //
    #[inline]
    pub fn nb_invalid_bases(&self, seq : &[u8]) -> u32 {
        let sum = seq.iter().fold(0u32, |acc, &b| acc+ !self.is_valid_base(b) as u32);
        return sum;
    }
}  // end impl Alphabet4b


impl BaseCompress for Alphabet4b {

    #[inline(always)]
    fn encode(&self, c:u8) -> u8 {
        match c {
            b'A' => 0b0001,
            b'C' => 0b0010,
            b'G' => 0b0100,
            b'T' => 0b1000,
            b'N' => 0b1111,
            _    => panic!("char not in alpahabet4b"),
        }
    } // end of function encode

    #[inline(always)]
   fn decode(&self, c:u8) -> u8 {
        match c {
            0b0001 => b'A',
            0b0010 => b'C',
            0b0100 => b'G',
            0b1000 => b'T',
            0b1111 => b'N',   // just to fill up slices
            _    => panic!("pattern not a code in alpahabet_4b"),
        }
    }

    /// return base complement
    fn complement(&self, c:u8) -> u8 {
        match c {
            0b0001 => 0b1000,   // A -> T
            0b0010 => 0b0100,   // C -> G
            0b0100 => 0b0010,   // G -> C
            0b1000 => 0b0001,   // T -> A
            0b1111 => 0b1111,   // N -> N
            _        => panic!("pattern not a code in alpahabet_2b"),
        }
    } // end of complement


    #[inline(always)]
    fn  get_nb_bits(&self) -> u8 {
        return 4;
    }

    #[inline(always)]
    fn is_valid_base(&self, c: u8) -> bool {
        match c {
            b'A' | b'C' | b'G' | b'T' => true,
            _    => false,
        }
    } // end is_valid_base


    // we expect a slice of at least 2 bytes
    fn base_pack(&self, to_pack : &[u8]) -> u8 {
        debug_assert!(to_pack.len() >= 2);
        //
        let mut packed;
//        println!("packing = {}", to_pack[0]);
        packed = self.encode(to_pack[0]);
//        println!("packed = {}", packed);
        packed = packed << 4;
//        println!("packed = {}", packed);
//        println!("packing = {}", to_pack[1]);
        packed |= self.encode(to_pack[1]);
//        println!("packed = {}", packed);
        //
        return packed;
    }
} // end implement section Alphabet4b



///////////////////////////////////////////////////////////////////////////////////////
//         Alphabet8b
///////////////////////////////////////////////////////////////////////////////////////


/// uncompressed representation of sequence
pub struct Alphabet8b {
    pub bases: String,
}

impl Alphabet8b {
    pub fn new() -> Alphabet8b{
        Alphabet8b { bases : String::from("ACGT")}
    }
    //
    pub fn len(&self) -> usize {
        return 8;
    }
}



impl BaseCompress for Alphabet8b {

    #[inline(always)]
    fn encode(&self, c:u8) -> u8 {
        return c;
    } // end of function encode

    #[inline(always)]
   fn decode(&self, c:u8) -> u8 {
       return c;
    }

    /// return base complement
    fn complement(&self, c:u8) -> u8 {
        match c {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _    => panic!("pattern not a code in alpahabet_2b"),
        }
    } // end of complement


    #[inline(always)]
    fn  get_nb_bits(&self) -> u8 {
        return 8;
    }

    #[inline(always)]
    fn is_valid_base(&self, c: u8) -> bool {
        match c {
            b'A' | b'C' | b'G' | b'T' => true,
            _    => false,
        }
    } // end is_valid_base

    // we expect a slice of at least 2 bytes
    fn base_pack(&self, to_pack : &[u8]) -> u8 {
        //
        let packed = to_pack[0];
        //
        return packed;
    }
} // end implement section Alphabet8b


