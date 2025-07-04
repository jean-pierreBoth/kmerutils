//! This file implements KmerAA32bit and KmerAA64bit representing Kmer for Amino Acid.  
//!
//! We implement compression of bases on 5 bits stored in a u32 or a u64.
//! So KmerAA64bit can store up to 12 AA. For less than 6 AA a u32 is sufficient, use KmerAA32bit.  
//!
//! The module provides Kmer generation tools KmerSeqIterator and KmerGenerationPattern
//! as in module base.

#![allow(clippy::transmute_num_to_bytes)]

use std::mem;
use std::mem::size_of;

use std::io;

use std::str::FromStr;

use anyhow::*;
use fnv::{FnvBuildHasher, FnvHashMap};
use std::cmp::Ordering;
use std::ops::Range;

#[allow(unused)]
use log::{debug, error, info};

use crate::base::kmertraits::*;

/// alphabet of RNA is encoded from 1 to 20 according to lexicographic order.
pub struct Alphabet {
    pub bases: String,
}

/*
 We keep the 0 bit field

A = 00001
C = 00010
D = 00011
E = 00100
F = 00101
G = 00110
H = 00111
I = 01000
K = 01001
L = 01010
M = 01011
N = 01100
P = 01101
Q = 01111
R = 10000
S = 10001
T = 10010
V = 10011
W = 10100
Y = 10101
*/

impl Default for Alphabet {
    fn default() -> Self {
        Self::new()
    }
}

impl Alphabet {
    pub fn new() -> Alphabet {
        Alphabet {
            bases: String::from("ACDEFGHIKLMNPQRSTVWY"),
        }
    }
    //
    pub fn len(&self) -> u8 {
        self.bases.len() as u8
    }

    #[inline(always)]
    pub fn is_valid_base(&self, c: u8) -> bool {
        self.bases.find(c as char).is_some()
    } // end is_valid_base

    pub fn get_nb_bits(&self) -> u8 {
        5
    }

    // encode a base into its bit pattern and returns it in a u8
    fn encode(&self, c: u8) -> u8 {
        match c {
            b'A' => 0b00001,
            b'C' => 0b00010,
            b'D' => 0b00011,
            b'E' => 0b00100,
            b'F' => 0b00101,
            b'G' => 0b00110,
            b'H' => 0b00111,
            b'I' => 0b01000,
            b'K' => 0b01001,
            b'L' => 0b01010,
            b'M' => 0b01011,
            b'N' => 0b01100,
            b'P' => 0b01101,
            b'Q' => 0b01111,
            b'R' => 0b10000,
            b'S' => 0b10001,
            b'T' => 0b10010,
            b'V' => 0b10011,
            b'W' => 0b10100,
            b'Y' => 0b10101,
            _ => panic!("encode: not a code in alpahabet for amino acid: {:x}", c),
        } // end of match
    } // end of encode

    fn decode(&self, c: u8) -> u8 {
        match c {
            0b00001 => b'A',
            0b00010 => b'C',
            0b00011 => b'D',
            0b00100 => b'E',
            0b00101 => b'F',
            0b00110 => b'G',
            0b00111 => b'H',
            0b01000 => b'I',
            0b01001 => b'K',
            0b01010 => b'L',
            0b01011 => b'M',
            0b01100 => b'N',
            0b01101 => b'P',
            0b01111 => b'Q',
            0b10000 => b'R',
            0b10001 => b'S',
            0b10010 => b'T',
            0b10011 => b'V',
            0b10100 => b'W',
            0b10101 => b'Y',
            _ => panic!(
                "decode : pattern not a code in alpahabet for Amino Acid got : {:#b}",
                c & 0b11111
            ),
        }
    } // end of decode
} // end of impl Alphabet

//=======================================================================================
/// A Kmer of amino acids represented on 32 bits, it can store up to 6 AA.  
/// See also KmerAA64bit for less than 12 AA.
/// We implement Amino Acid Kmer as packed in a u32 using 5bits by base, so we can go up to 6 bases.

#[derive(Copy, Clone, Hash)]
pub struct KmerAA32bit {
    aa: u32,
    nb_base: u8,
} // end of struct KmerAA128bit

impl KmerAA32bit {
    pub fn new(nb_base: u8) -> Self {
        let nb_base_max = u32::BITS as usize / 5;
        if nb_base as usize >= nb_base_max {
            panic!(
                "For KmerAA32bit nb_base must be less or equal to {}",
                nb_base_max
            )
        }
        KmerAA32bit { aa: 0, nb_base }
    }
} // end of impl KmerAA128bit

impl KmerT for KmerAA32bit {
    fn get_nb_base(&self) -> u8 {
        self.nb_base
    } // end of get_nb_base

    //
    fn push(&self, c: u8) -> Self {
        // shift left 5 bits, insert new base and enforce 0 at upper bits
        let value_mask: u32 = (0b1 << (5 * self.get_nb_base())) - 1;
        // contrary to dna sequence base in seq is not encoded, we must encode it!!
        let encoded_base = Alphabet::new().encode(c);
        let new_kmer = ((self.aa << 5) & value_mask) | (encoded_base as u32 & 0b11111);
        log::trace!("after push {:#b}", new_kmer);
        KmerAA32bit {
            aa: new_kmer,
            nb_base: self.nb_base,
        }
    } // end of push

    // TODO
    fn reverse_complement(&self) -> Self {
        panic!("KmerAA32bit reverse_complement not yet implemented");
    } // end of reverse_complement

    fn dump(&self, bufw: &mut dyn io::Write) -> io::Result<usize> {
        let val = unsafe { mem::transmute::<u8, [u8; 1]>(self.nb_base) };
        bufw.write_all(&val).unwrap();
        let val = unsafe { mem::transmute::<u32, [u8; 4]>(self.aa) };
        bufw.write(&val)
    }
} // end of impl KmerT block for KmerAA128bit

impl PartialEq for KmerAA32bit {
    // we must check equality of field
    fn eq(&self, other: &KmerAA32bit) -> bool {
        (self.aa == other.aa) & (self.nb_base == other.nb_base)
    }
} // end of impl PartialEq for KmerAA128bit

impl Eq for KmerAA32bit {}

/// We define ordering as a kind of "lexicographic" order by taking into account first number of base.
/// The more the number of base the greater. Then we have integer comparison between aa parts
///
impl CompressedKmerT for KmerAA32bit {
    type Val = u32;

    fn get_nb_base_max() -> usize {
        u32::BITS as usize / 5
    }

    /// a decompressing function mainly for test and debugging purpose
    fn get_uncompressed_kmer(&self) -> Vec<u8> {
        let nb_bases = self.nb_base;
        let alphabet = Alphabet::new();
        // we treat each block of 5 bits as u8 end call decoder of Alphabet
        let mut decompressed_kmer = Vec::<u8>::with_capacity(nb_bases as usize);
        let mut base: u8;
        //
        let mut buf = self.aa;
        // get the base coding part at left end of u32
        log::debug!(
            "rotating left {}",
            Self::Val::BITS as usize - 5 * nb_bases as usize
        );
        buf = buf.rotate_left(8 * size_of::<Self::Val>() as u32 - 5 * nb_bases as u32);
        for _ in 0..nb_bases {
            buf = buf.rotate_left(5);
            base = (buf & 0b11111) as u8;
            decompressed_kmer.push(alphabet.decode(base));
        }
        decompressed_kmer
    }

    /// return the pure value with part coding number of bases reset to 0.
    #[inline(always)]
    fn get_compressed_value(&self) -> Self::Val {
        self.aa
    }

    #[inline(always)]
    fn get_bitsize(&self) -> usize {
        Self::Val::BITS as usize
    }
} // end of impl CompressedKmerT for KmerAA32bit

//===================================================================

impl Ord for KmerAA32bit {
    fn cmp(&self, other: &KmerAA32bit) -> Ordering {
        if self.nb_base != other.nb_base {
            (self.nb_base).cmp(&(other.nb_base))
        } else {
            (self.aa).cmp(&(other.aa))
        }
    } // end cmp
} // end impl Ord for KmerAA128bit

impl PartialOrd for KmerAA32bit {
    fn partial_cmp(&self, other: &KmerAA32bit) -> Option<Ordering> {
        Some(self.cmp(other))
    } // end partial_cmp
} // end impl Ord for KmerAA128bit

impl KmerBuilder<KmerAA32bit> for KmerAA32bit {
    /// for Kmer32bit we encode the number of bases in the 4 upper bits
    fn build(val: u32, nb_base: u8) -> KmerAA32bit {
        KmerAA32bit { aa: val, nb_base }
    }
}

//======================================================================

/// A Kmer of amino acids for less than 12 Amino Acid, stored on a u64.

#[derive(Copy, Clone, Hash)]
pub struct KmerAA64bit {
    aa: u64,
    nb_base: u8,
} // end of struct KmerAA64bit

impl KmerAA64bit {
    pub fn new(nb_base: u8) -> Self {
        if nb_base >= 12 {
            panic!("For KmerAA64bit nb_base must be less or equal to 12")
        }
        KmerAA64bit { aa: 0, nb_base }
    }
} // end of impl KmerAA64bit

impl KmerT for KmerAA64bit {
    fn get_nb_base(&self) -> u8 {
        self.nb_base
    } // end of get_nb_base

    //
    fn push(&self, c: u8) -> Self {
        // shift left 5 bits, insert new base and enforce 0 at upper bits
        let value_mask: u64 = (0b1 << (5 * self.get_nb_base())) - 1;
        // contrary to dna sequence base in seq is not encoded, we must encode it!!
        let encoded_base = Alphabet::new().encode(c);
        let new_kmer = ((self.aa << 5) & value_mask) | (encoded_base as u64 & 0b11111);
        log::trace!("after push {:#b}", new_kmer);
        KmerAA64bit {
            aa: new_kmer,
            nb_base: self.nb_base,
        }
    } // end of push

    // TODO
    fn reverse_complement(&self) -> Self {
        panic!("KmerAA64bit reverse_complement not yet implemented");
    } // end of reverse_complement

    fn dump(&self, bufw: &mut dyn io::Write) -> io::Result<usize> {
        let val = unsafe { mem::transmute::<u8, [u8; 1]>(self.nb_base) };
        bufw.write_all(&val).unwrap();
        //
        let val = unsafe { mem::transmute::<u64, [u8; 8]>(self.aa) };
        bufw.write(&val)
    }
} // end of impl KmerT block for KmerAA64bit

impl PartialEq for KmerAA64bit {
    // we must check equality of field
    fn eq(&self, other: &KmerAA64bit) -> bool {
        (self.aa == other.aa) & (self.nb_base == other.nb_base)
    }
} // end of impl PartialEq for KmerAA128bit

impl Eq for KmerAA64bit {}

impl Ord for KmerAA64bit {
    fn cmp(&self, other: &KmerAA64bit) -> Ordering {
        if self.nb_base != other.nb_base {
            (self.nb_base).cmp(&(other.nb_base))
        } else {
            (self.aa).cmp(&(other.aa))
        }
    } // end cmp
} // end impl Ord for KmerAA64bit

impl PartialOrd for KmerAA64bit {
    fn partial_cmp(&self, other: &KmerAA64bit) -> Option<Ordering> {
        Some(self.cmp(other))
    } // end partial_cmp
} // end impl Ord for KmerAA128bit

impl CompressedKmerT for KmerAA64bit {
    type Val = u64;

    fn get_nb_base_max() -> usize {
        12
    }

    /// a decompressing function mainly for test and debugging purpose
    fn get_uncompressed_kmer(&self) -> Vec<u8> {
        let nb_bases = self.nb_base;
        let alphabet = Alphabet::new();
        //
        let mut decompressed_kmer = Vec::<u8>::with_capacity(nb_bases as usize);
        let mut base: u8;
        //
        let mut buf = self.aa;
        // get the base coding part at left end of u32
        log::debug!("rotating left {}", 64 - 5 * nb_bases);
        buf = buf.rotate_left((64 - 5 * nb_bases) as u32);
        for _ in 0..nb_bases {
            buf = buf.rotate_left(5);
            base = (buf & 0b11111) as u8;
            decompressed_kmer.push(alphabet.decode(base));
        }
        decompressed_kmer
    }

    /// return the pure value with part coding number of bases reset to 0.
    #[inline(always)]
    fn get_compressed_value(&self) -> u64 {
        self.aa
    }

    #[inline(always)]
    fn get_bitsize(&self) -> usize {
        Self::Val::BITS as usize
    }
} // end of impl CompressedKmerT for KmerAA64bit

impl KmerBuilder<KmerAA64bit> for KmerAA64bit {
    /// for Kmer32bit we encode the number of bases in the 4 upper bits
    fn build(val: u64, nb_base: u8) -> KmerAA64bit {
        KmerAA64bit { aa: val, nb_base }
    }
} // end of  KmerBuilder<KmerAA64bit>

//=======================================================================

/// our sequence of Amino Acid is encoded on a byte (even if 5 bits are enough but we do not store sequences yet)
/// If necessary an implementation on bitvec could be used using a struct SeqIterator as a bridge from
/// Sequence to KmerSeqIterator
pub struct SequenceAA {
    seq: Vec<u8>,
}

impl SequenceAA {
    /// allocates and check for compatibility with alphabet
    pub fn new(str: &[u8]) -> Self {
        let alphabet = Alphabet::new();
        let _res = str.iter().map(|c| {
            if !alphabet.is_valid_base(*c) {
                log::error!("SequenceAA str is : {:?}", str);
                log::error!("SequenceAA character not in alphabet {}", c);
                std::process::abort();
            }
        });
        SequenceAA { seq: str.to_vec() }
    } // end of new

    pub fn len(&self) -> usize {
        self.seq.len()
    }

    pub fn is_empty(&self) -> bool {
        self.seq.len() == 0
    }
    /// return the the uncompressed lenght (maintained by analogy with DNA case)
    pub fn size(&self) -> usize {
        self.seq.len()
    }

    pub fn get_base(&self, pos: usize) -> u8 {
        if pos >= self.seq.len() {
            log::error!("\n base position after end of sequence , pos : {}", pos);
            log::error!("\n sequence(bytes) is : {:?} \n", self.seq);
            let str = String::from_utf8(self.seq.clone()).unwrap();
            log::error!("\n sequence(ascii) is  : {:?} \n", str);
            std::process::exit(1);
        } else {
            self.seq[pos]
        }
    } // end of get_base

    /// allocated a new sequence filtering out non valid bases
    pub fn new_filtered(buf: &[u8], alphabet: &Alphabet) -> Self {
        let nb_bases = buf.len();
        let mut seq = Vec::<u8>::with_capacity(nb_bases);
        for b in buf {
            if alphabet.is_valid_base(*b) {
                seq.push(*b);
            }
        }
        SequenceAA { seq: seq.to_vec() }
    }
} // end of SequenceAA

impl FromStr for SequenceAA {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        //
        let sbytes = s.as_bytes();
        let alphabet = Alphabet::new();
        //
        let _res = sbytes.iter().map(|c| {
            if !alphabet.is_valid_base(*c) {
                panic!("character not in alphabet {}", c);
            }
        });
        Ok(SequenceAA {
            seq: sbytes.to_vec(),
        })
    }
} // end of FromStr

impl ToString for SequenceAA {
    fn to_string(&self) -> String {
        //
        let stru8: Vec<u8> = self.seq.to_vec();
        std::string::String::from_utf8(stru8).unwrap()
    } // end of to_string
} //  end of ToString

//=========================================================================

// TODO factorize this trait?
/// A trait for generating Kmer along a sequence
pub trait KmerSeqIteratorT {
    /// KmerAA128bit
    type KmerVal;
    /// get next kmer or None
    fn next(&mut self) -> Option<Self::KmerVal>;
}

/// An iterator that produce Kmer along a sequence
pub struct KmerSeqIterator<'a, T>
where
    T: CompressedKmerT,
{
    /// size of kmer
    nb_base: usize,
    /// an iterator for base calling
    sequence: &'a SequenceAA,
    /// The alphabet needed to encode base from sequence into a (encoded!!) Kmer
    alphabet_aa: Alphabet,
    /// last position of last kmer returned. At the beginning its None
    previous: Option<T>,

    //
    range: Range<usize>,
    /// at present time, sequence for Amino Acid are not compressed, only Kmer so we do not need IterSequence as in mode base
    base_position: usize,
} // end of KmerSeqIterator

impl<'a, T> KmerSeqIterator<'a, T>
where
    T: CompressedKmerT,
{
    pub fn new(kmer_size: usize, seq: &'a SequenceAA) -> Self {
        let alphabet_aa = Alphabet::new();
        let range = std::ops::Range {
            start: 0,
            end: seq.len(),
        };
        let base_position = 0;
        KmerSeqIterator {
            nb_base: kmer_size,
            sequence: seq,
            alphabet_aa,
            previous: None,
            range,
            base_position,
        }
    }

    /// defines the range of kmer generation.  
    /// All bases in kmer generated must be between in first..last last excluded!
    pub fn set_range(&mut self, first: usize, last: usize) -> Result<(), anyhow::Error> {
        if last <= first || last > self.sequence.len() {
            log::error!(
                "KmerSeqIterator set_range first : {}, last : {}, length : {}, seq is : {:?}",
                first,
                last,
                self.sequence.len(),
                self.sequence.to_string()
            );
            Err(anyhow!("bad range for iterator"))
        } else {
            self.range = Range {
                start: first,
                end: last,
            };
            self.base_position = first;
            Ok(())
        }
    } // end of set_range
} // end of impl block for KmerSeqIterator

impl<'a, Kmer> KmerSeqIteratorT for KmerSeqIterator<'a, Kmer>
where
    Kmer: CompressedKmerT + KmerBuilder<Kmer>,
{
    type KmerVal = Kmer;

    /// iterates...
    fn next(&mut self) -> Option<Self::KmerVal> {
        // check for end of iterator
        if self.base_position >= self.sequence.len().min(self.range.end) {
            log::trace!(
                "iterator exiting at base pos {} range.end {} ",
                self.base_position,
                self.range.end
            );
            return None;
        }
        // check for too small seq
        if self.sequence.len() < self.nb_base {
            return None;
        }
        // now we know we are not at end of iterator
        // if we do not have a previous we have to contruct first kmer
        // we have to push a base.
        //
        if let Some(kmer) = self.previous {
            // in fact we have the base to push
            let next_base = self.sequence.get_base(self.base_position);
            log::trace!(
                " next pushing base : {}",
                char::from_u32(next_base as u32).unwrap()
            );
            self.previous = Some(kmer.push(next_base));
            self.base_position += 1;
            self.previous
        } else {
            // we are at beginning of kmer construction sequence, we must push kmer_size bases
            let kmer_size = self.nb_base;
            let nb_base_bits = 5;
            let pos = nb_base_bits * (kmer_size - 1);
            let next_base = self.sequence.get_base(self.base_position);
            let encoded_base = self.alphabet_aa.encode(next_base);
            self.base_position += 1;
            log::trace!(
                " init kmerAA base : {}",
                char::from_u32(next_base as u32).unwrap()
            );
            let mut new_kmer_val = <Kmer as CompressedKmerT>::Val::from(encoded_base) << pos;
            for i in 0..(kmer_size - 1) {
                let next_base = self.sequence.get_base(self.base_position);
                log::trace!(
                    " next kmerAA base : {}",
                    char::from_u32(next_base as u32).unwrap()
                );
                let encoded_base = self.alphabet_aa.encode(next_base);
                let base_val = <Kmer as CompressedKmerT>::Val::from(encoded_base)
                    << (pos - nb_base_bits - nb_base_bits * i);
                new_kmer_val = new_kmer_val | base_val;
                self.base_position += 1;
            }
            // we know by assert that self.nb_base is less than 12 so can be seen as a u8
            let new_kmer: Kmer =
                <Kmer as KmerBuilder<Kmer>>::build(new_kmer_val, self.nb_base as u8);
            self.previous = Some(new_kmer);
            self.previous
        }
    } // end of next
} // end of impl  KmerSeqIteratorT for KmerSeqIterator<'a, Kmer>

//============================================================================

/// A trait describing different mode of kmer production
pub trait KmerGenerationPattern<T: KmerT> {
    /// generate all kmers included in 0..
    fn generate_kmer_pattern(&self, seq: &SequenceAA) -> Vec<T>;
    /// generate all kmers inclused in begin..end with end excluded as in rust conventions.
    fn generate_kmer_pattern_in_range(&self, seq: &SequenceAA, begin: usize, end: usize) -> Vec<T>;
    /// generate kmers with their multiplicities
    fn generate_kmer_distribution(&self, seq: &SequenceAA) -> FnvHashMap<T, usize>;
}

//==========================================================================================

use std::marker::PhantomData;

pub struct KmerGenerator<T: KmerT> {
    /// size of kmer we generate
    pub kmer_size: u8,
    t_marker: PhantomData<T>,
}

impl<T: KmerT> KmerGenerator<T> {
    pub fn new(ksize: u8) -> Self {
        KmerGenerator {
            kmer_size: ksize,
            t_marker: PhantomData,
        }
    }
    /// generic driver for kmer generation
    pub fn generate_kmer(&self, seq: &SequenceAA) -> Vec<T>
    where
        Self: KmerGenerationPattern<T>,
    {
        self.generate_kmer_pattern(seq)
    }
    /// generic driver for kmer generation
    pub fn generate_kmer_in_range(&self, seq: &SequenceAA, begin: usize, end: usize) -> Vec<T>
    where
        Self: KmerGenerationPattern<T>,
    {
        self.generate_kmer_pattern_in_range(seq, begin, end)
    }
    /// generic driver for kmer distribution pattern
    pub fn generate_weighted_kmer(&self, seq: &SequenceAA) -> FnvHashMap<T, usize>
    where
        Self: KmerGenerationPattern<T>,
    {
        self.generate_kmer_distribution(seq)
    }
    //
    pub fn get_kmer_size(&self) -> usize {
        self.kmer_size as usize
    }
} // end of impl KmerGenerator

/// A utility to convert FnvHashMap<T,usize> to Vec<(T, usize)>
pub fn hashmap_count_to_vec_count<T: CompressedKmerT + std::hash::Hash + Eq>(
    kmer_distribution: &FnvHashMap<T, usize>,
) -> Vec<(T, usize)> {
    // convert to a Vec
    let hashed_kmers = kmer_distribution.keys();
    let mut weighted_kmer = Vec::<(T, usize)>::with_capacity(kmer_distribution.len());
    for key in hashed_kmers {
        if let Some(weight) = kmer_distribution.get(key) {
            weighted_kmer.push((*key, *weight));
        };
    }
    //
    weighted_kmer
} // end of hashmap_count_to_vec_count

/*
   Now we have the basics of Kmer Traits we implement KmerSeqIterator and KmerGenerationPattern
   Implementation  of Kmer Generation for KmerAA128bit
   ===================================================
*/

// We need a guess to allocate HashMap used with Kmer Generation
// for very long sequence we must avoid nb_kmer to sequence length! Find a  good heuristic
pub(super) fn get_nbkmer_guess(seq: &SequenceAA) -> usize {
    let nb = 10_000_000 * (1usize + seq.len().ilog2() as usize);

    seq.len().min(nb)
} // end of get_nbkmer_guess

// We need a guess to allocate HashMap used with Kmer Generation
// for vector of sequenc coming from a non concatnated file, we must avoid nb_kmer to sequence length! Find a  good heuristic
pub(super) fn get_nbkmer_guess_seqs(vseq: &[&SequenceAA]) -> usize {
    let total_nb_base = vseq.iter().fold(0, |acc, seq| acc + seq.size());
    total_nb_base.min(10_000_000 * (1usize + total_nb_base.ilog2() as usize))
} // end of get_nbkmer_guess_seqs

/// implementation of kmer generation pattern for KmerAA32bit\<N\>
impl KmerGenerationPattern<KmerAA32bit> for KmerGenerator<KmerAA32bit> {
    fn generate_kmer_pattern(&self, seq: &SequenceAA) -> Vec<KmerAA32bit> {
        if self.kmer_size as usize > KmerAA32bit::get_nb_base_max() {
            panic!(
                "KmerAA32bit cannot have size greater than {} !!",
                KmerAA32bit::get_nb_base_max()
            ); // cannot happen !
        }
        let kmer_size = self.kmer_size as usize;
        // For a sequence of size the number of kmer is seq.size - kmer.size + 1  !!!
        // But it happens that "long reads" are really short
        let nb_kmer = if seq.len() >= kmer_size {
            seq.len() - kmer_size + 1
        } else {
            0
        };
        let nb_kmer = nb_kmer.min(get_nbkmer_guess(seq));
        let mut kmer_vect = Vec::<KmerAA32bit>::with_capacity(nb_kmer);
        let mut kmeriter = KmerSeqIterator::<KmerAA32bit>::new(kmer_size, seq);
        while let Some(kmer) = kmeriter.next() {
            kmer_vect.push(kmer);
        }
        //
        kmer_vect
    } // end of generate_kmer_pattern

    /// generate all kmers associated to their multiplicity
    /// This is useful in the context of Jaccard Probability Index estimated with ProbminHash
    fn generate_kmer_distribution(&self, seq: &SequenceAA) -> FnvHashMap<KmerAA32bit, usize> {
        if self.kmer_size as usize > KmerAA32bit::get_nb_base_max() {
            panic!(
                "KmerAA32bit cannot have size greater than {} !!",
                KmerAA32bit::get_nb_base_max()
            ); // cannot happen !
        }
        // For a sequence of size the number of kmer is seq.size - kmer.size + 1  !!!
        // But it happens that "long reads" are really short
        let kmer_size = self.kmer_size as usize;
        //
        let nb_kmer = if seq.len() >= kmer_size {
            seq.len() - kmer_size + 1
        } else {
            0
        };
        let nb_kmer = nb_kmer.min(get_nbkmer_guess(seq));
        let mut kmer_distribution: FnvHashMap<KmerAA32bit, usize> =
            FnvHashMap::with_capacity_and_hasher(nb_kmer, FnvBuildHasher::default());
        let mut kmeriter = KmerSeqIterator::<KmerAA32bit>::new(kmer_size, seq);
        while let Some(kmer) = kmeriter.next() {
            // do we store the kmer in the FnvHashMap or a already hashed value aka nthash?
            *kmer_distribution.entry(kmer).or_insert(0) += 1;
        }
        //
        kmer_distribution
    } // end of generate_kmer_pattern

    fn generate_kmer_pattern_in_range(
        &self,
        seq: &SequenceAA,
        begin: usize,
        end: usize,
    ) -> Vec<KmerAA32bit> {
        if self.kmer_size as usize > KmerAA32bit::get_nb_base_max() {
            panic!(
                "KmerAA32bit cannot have size greater than {} !!",
                KmerAA32bit::get_nb_base_max()
            ); // cannot happen !
        }
        if begin >= end {
            panic!("KmerGenerationPattern<'a, KmerAA32bit> bad range for kmer iteration");
        }
        // For a sequence of size the number of kmer is seq.size - kmer.size + 1  !!!
        // But it happens that "long reads" are really short
        let kmer_size = self.kmer_size as usize;
        let nb_kmer = if seq.len() >= kmer_size {
            seq.len() - kmer_size + 1
        } else {
            0
        };
        let nb_kmer = nb_kmer.min(get_nbkmer_guess(seq));
        let mut kmer_vect = Vec::<KmerAA32bit>::with_capacity(nb_kmer);
        let mut kmeriter = KmerSeqIterator::<KmerAA32bit>::new(kmer_size, seq);
        kmeriter.set_range(begin, end).unwrap();
        while let Some(kmer) = kmeriter.next() {
            kmer_vect.push(kmer);
        }
        //
        kmer_vect
    } // end of generate_kmer_pattern
} // end of impl KmerGenerationPattern<'a, KmerAA32bit<N>>

/*
Implementation for Kmer64 bit
*/

/// implementation of kmer generation pattern for KmerAA64bit\<N\>
impl KmerGenerationPattern<KmerAA64bit> for KmerGenerator<KmerAA64bit> {
    fn generate_kmer_pattern(&self, seq: &SequenceAA) -> Vec<KmerAA64bit> {
        if self.kmer_size > 12 {
            panic!("KmerAA64bit cannot have size greater than 12!!"); // cannot happen !
        }
        let kmer_size = self.kmer_size as usize;
        // For a sequence of size the number of kmer is seq.size - kmer.size + 1  !!!
        // But it happens that "long reads" are really short
        let nb_kmer = if seq.len() >= kmer_size {
            seq.len() - kmer_size + 1
        } else {
            0
        };
        let nb_kmer = nb_kmer.min(get_nbkmer_guess(seq));
        let mut kmer_vect = Vec::<KmerAA64bit>::with_capacity(nb_kmer);
        let mut kmeriter = KmerSeqIterator::<KmerAA64bit>::new(kmer_size, seq);
        while let Some(kmer) = kmeriter.next() {
            kmer_vect.push(kmer);
        }
        //
        kmer_vect
    } // end of generate_kmer_pattern

    /// generate all kmers associated to their multiplicity
    /// This is useful in the context of Jaccard Probability Index estimated with ProbminHash
    fn generate_kmer_distribution(&self, seq: &SequenceAA) -> FnvHashMap<KmerAA64bit, usize> {
        if self.kmer_size as usize > 12 {
            panic!("KmerAA128bit cannot be greater than 12!!"); // cannot happen
        }
        // For a sequence of size the number of kmer is seq.size - kmer.size + 1  !!!
        // But it happens that "long reads" are really short
        let kmer_size = self.kmer_size as usize;
        //
        let nb_kmer = if seq.len() >= kmer_size {
            seq.len() - kmer_size + 1
        } else {
            0
        };
        let nb_kmer = nb_kmer.min(get_nbkmer_guess(seq));
        let mut kmer_distribution: FnvHashMap<KmerAA64bit, usize> =
            FnvHashMap::with_capacity_and_hasher(nb_kmer, FnvBuildHasher::default());
        let mut kmeriter = KmerSeqIterator::<KmerAA64bit>::new(kmer_size, seq);
        while let Some(kmer) = kmeriter.next() {
            // do we store the kmer in the FnvHashMap or a already hashed value aka nthash?
            *kmer_distribution.entry(kmer).or_insert(0) += 1;
        }
        //
        kmer_distribution
    } // end of generate_kmer_pattern

    fn generate_kmer_pattern_in_range(
        &self,
        seq: &SequenceAA,
        begin: usize,
        end: usize,
    ) -> Vec<KmerAA64bit> {
        if self.kmer_size as usize > 12 {
            panic!("KmerAA64bit cannot have size greater than 12"); // cannot happen
        }
        if begin >= end {
            panic!("KmerGenerationPattern<'a, KmerAA64bit>  bad range for kmer iteration");
        }
        // For a sequence of size the number of kmer is seq.size - kmer.size + 1  !!!
        // But it happens that "long reads" are really short
        let kmer_size = self.kmer_size as usize;
        let nb_kmer = if seq.len() >= kmer_size {
            seq.len() - kmer_size + 1
        } else {
            0
        };
        let nb_kmer = nb_kmer.min(get_nbkmer_guess(seq));
        let mut kmer_vect = Vec::<KmerAA64bit>::with_capacity(nb_kmer);
        let mut kmeriter = KmerSeqIterator::<KmerAA64bit>::new(kmer_size, seq);
        kmeriter.set_range(begin, end).unwrap();
        while let Some(kmer) = kmeriter.next() {
            kmer_vect.push(kmer);
        }
        //
        kmer_vect
    } // end of generate_kmer_pattern
} // end of impl KmerGenerationPattern<'a, KmerAA64bit<N>>

//===========================================================

#[cfg(test)]
mod tests {

    // to run with  cargo test -- --nocapture KmerAA128bit
    //  possibly with export RUST_LOG=INFO,kmerutils::rnautils=debug

    use super::*;

    fn log_init_test() {
        let mut builder = env_logger::Builder::from_default_env();
        //    builder.filter_level(LevelFilter::Trace);
        let _ = builder.is_test(true).try_init();
    }

    // test iterator
    #[test]
    fn test_seqaa_32bit_iterator_range() {
        log_init_test();
        //
        log::debug!("in test_seqaa_128bit_iterator_range");
        //
        let str = "MTEQIELIKLYSTRILALAAQMPHVGSLDNPDASAMKRSPLCGSKVTVDVIMQNGKITEFAQNVKACALGQAAASVAAQNIIGRTAEEVVRARDELAAMLKSGGPPPGPPFDGFEVLAPASEYKNRHASILLSLDATAEACASIAAQNSA";

        let seqaa = SequenceAA::from_str(str).unwrap();
        // ask for Kmer of size 4
        let mut seq_iterator = KmerSeqIterator::<KmerAA32bit>::new(4, &seqaa);
        // set a range
        seq_iterator.set_range(3, 10).unwrap(); // so that we have 4 4-Kmer  (4 = 10-1-kmer_size-3)
        // So we must havr from "QIEL"
        let mut kmer_num = 0;
        let kmer_res = ["QIEL", "IELI", "ELIK", "LIKL"];
        while let Some(kmer) = seq_iterator.next() {
            let k_uncompressed = kmer.get_uncompressed_kmer();
            let kmer_str = std::str::from_utf8(&k_uncompressed).unwrap();
            //            log::info!(" kmer {} = {:?}", kmer_num, kmer_str);
            if kmer_str != kmer_res[kmer_num] {
                log::error!(" kmer {} = {:?}", kmer_num, kmer_str);
                panic!(
                    "error in KmerAA32bit test::test_seq_aa_iterator \n at kmer num {}, got {:?} instead of {:?}",
                    kmer_num, kmer_str, kmer_res[kmer_num]
                );
            }
            kmer_num += 1;
        }
        // check iterator sees the end
        match seq_iterator.next() {
            Some(_kmer) => {
                panic!("iterator do not see end");
            }
            None => (),
        } // end match
    } // end of test_seqaa_iterator_range

    #[test]
    fn test_seqaa_64bit_iterator_range() {
        log_init_test();
        //
        log::debug!("in test_seqaa_64bit_iterator_range");
        //
        let str = "MTEQIELIKLYSTRILALAAQMPHVGSLDNPDASAMKRSPLCGSKVTVDVIMQNGKITEFAQNVKACALGQAAASVAAQNIIGRTAEEVVRARDELAAMLKSGGPPPGPPFDGFEVLAPASEYKNRHASILLSLDATAEACASIAAQNSA";

        let seqaa = SequenceAA::from_str(str).unwrap();
        // ask for Kmer of size 4
        let mut seq_iterator = KmerSeqIterator::<KmerAA64bit>::new(4, &seqaa);
        // set a range
        seq_iterator.set_range(3, 10).unwrap(); // so that we have 4 4-Kmer  (4 = 10-1-kmer_size-3)
        // So we must havr from "QIEL"
        let mut kmer_num = 0;
        let kmer_res = ["QIEL", "IELI", "ELIK", "LIKL"];
        while let Some(kmer) = seq_iterator.next() {
            let k_uncompressed = kmer.get_uncompressed_kmer();
            let kmer_str = std::str::from_utf8(&k_uncompressed).unwrap();
            //            log::info!(" kmer {} = {:?}", kmer_num, kmer_str);
            if kmer_str != kmer_res[kmer_num] {
                log::error!(" kmer {} = {:?}", kmer_num, kmer_str);
                panic!(
                    "error in KmerAA64bit test::test_seq_aa_iterator \n at kmer num {}, got {:?} instead of {:?}",
                    kmer_num, kmer_str, kmer_res[kmer_num]
                );
            }
            kmer_num += 1;
        }
        // check iterator sees the end
        match seq_iterator.next() {
            Some(_kmer) => {
                panic!("iterator do not see end");
            }
            None => (),
        } // end match
    } // end of test_seqaa_64bit_iterator_range

    // test we arrive at end correctly
    #[test]
    fn test_seqaa_iterator_end() {
        //
        log_init_test();
        log::debug!("in test_seqaa_iterator_end");
        //
        let str = "MTEQIELIKLYSTRILALAAQMPHVGSLDNPD";
        let seqaa = SequenceAA::from_str(str).unwrap();
        // ask for Kmer of size 8
        let mut last_kmer: &str = "toto";
        let mut seq_iterator = KmerSeqIterator::<KmerAA64bit>::new(8, &seqaa);
        let mut kmer_num = 0;
        let mut k_uncompressed;

        while let Some(kmer) = seq_iterator.next() {
            log::debug!("in test_seqaa_iterator_end iteration {}", kmer_num);
            k_uncompressed = kmer.get_uncompressed_kmer();
            last_kmer = std::str::from_utf8(&k_uncompressed).unwrap();
            log::debug!(" kmer {} = {:?}", kmer_num, last_kmer);
            kmer_num += 1;
        }
        if last_kmer != "VGSLDNPD" {
            log::info!(" last kmer seen {} = {:?}", kmer_num, last_kmer);
            panic!("test_seqaa_iterator_end did not get the correct last_kmer");
        }
    } // end of test_seqaa_iterator_end

    #[test]
    fn test_aastr_conversion() {
        //
        log_init_test();
        log::debug!("in test_seqaa_iterator_end");
        //
        let str = "MTEQIELIKLYSTRILALAAQMPHVGSLDNPD";
        let seqaa = SequenceAA::from_str(str).unwrap();
        //
        let orig = seqaa.to_string();
        assert_eq!(orig, String::from(str));
    }
} // end of mod tests
