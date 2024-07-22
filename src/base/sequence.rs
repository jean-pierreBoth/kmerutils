//! Describes a sequence and its associated compression/decompression schemes, iteration
//! and basic operations as reverse_complement.

pub use super::alphabet::*;

//
//===================================================================================
// Sequence
//===================================================================================
//

/// a sequence as a vector of byte (instead of a bitvec). Each byte contains 1,2 or 4 bases depending on compression
#[derive(Clone)]
pub struct Sequence {
    seq: Vec<u8>,
    /// first byte in description is nb_bits by base (2 or 4 or 8 depending on compression used, 8 means no compression) ,
    /// second byte in description is number of bases in last byte if last byte is incomplete. So description\[1\] == 0
    /// means all bytes are complete!!
    description: [u8; 2],
}

impl Sequence {
    /// This constructor is adapted to small sequences. For large sequence that concatenated the whole file
    /// prefer the function [`Self::with_capacity()`] and [Self::encode_and_add()]
    pub fn new(raw: &[u8], nb_bits: u8) -> Sequence {
        let nb_bases = raw.len();
        let nb_bases_by_byte = 8 / nb_bits as usize;
        let nb_full_bytes = nb_bases / nb_bases_by_byte;
        let nb_bases_in_last_byte = nb_bases - nb_full_bytes * nb_bases_by_byte;
        // if is an expression see doc if let
        let nb_bytes_needed = if nb_bases_in_last_byte > 0 {
            nb_full_bytes + 1
        } else {
            nb_full_bytes
        };
        //
        let mut seq: Vec<u8> = Vec::new();
        // now we must fill data
        // for bytes from 0 to nb_bytes_needed-1 we fill nb_base_by_byte
        // and for last byte we fill nb_bases_in_last_byte
        // we chose to use a match so we can have simple loops, instead of using virtual methods via traits.
        match nb_bits {
            8 => {
                // not much to do , but it must be done.
                seq.extend_from_slice(raw);
            }
            // 2 bits by base
            2 => {
                seq.reserve(nb_bytes_needed);
                // 4 bases per byte, we treat indexes i*4 , i*4+1 , i*4+2, i*4+3  in raw
                let alfa2b = Alphabet2b::new();
                // loop from 0 to nb_bytes_needed-1 include, here nb_bases_by_byte = 4
                for i in 0..nb_full_bytes {
                    let encoded4b = alfa2b.base_pack(&raw[4 * i..4 * (i + 1)]);
                    seq.push(encoded4b);
                }
                //
                // do not forget last byte, we have up to 3 bases to add so one more byte is  sufficient
                //
                if nb_bases_in_last_byte > 0 {
                    // a simple solution is to get a new slice of 4 b'A' (that will be encoded as 0)
                    // then fill the head of the slice with the last bytes of argument raw
                    // and call alfa2b.base_pack once more.
                    // When unpacking and decoding as we know nb_bases_in_last_byte we
                    // can garbage supplementary A added. Yeah real shit!
                    let mut last_bases = [b'A'; 4];
                    for i in 0..nb_bases_in_last_byte {
                        last_bases[i] = raw[4 * nb_full_bytes + i];
                    }
                    let encoded4b = alfa2b.base_pack(&last_bases);
                    seq.push(encoded4b);
                }
            } // end case 2
            // 4 bits by base
            4 => {
                seq.reserve(nb_bytes_needed);
                // 2 bases per byte , we treat indexes i*2 and i*2+1 in raw
                let alfa4b = Alphabet4b::new();
                for i in 0..nb_full_bytes {
                    let encoded2b = alfa4b.base_pack(&raw[2 * i..2 * (i + 1)]);
                    seq.push(encoded2b);
                }
                // do not forget last byte from 2*nb_bytes_needed to 2*nb_bytes_needed +
                if nb_bases_in_last_byte > 0 {
                    // we have exactly one base to add, but as decoding use always 2 bases we add a b'N'
                    // to match the last part of byte.
                    let mut last_to_encode = [0u8; 2];
                    last_to_encode[0] = raw[2 * nb_full_bytes];
                    last_to_encode[1] = b'Z';
                    let encoded2b = alfa4b.base_pack(&last_to_encode);
                    seq.push(encoded2b);
                }
            }

            _ => panic!(
                " error in Sequence::new , bad number of bits by base must be 2 4 or 8, got {} ",
                nb_bits
            ),
        }
        //
        Sequence {
            // struct field : value
            seq,
            description: [nb_bits, nb_bases_in_last_byte as u8],
        }
    } // end new

    #[inline(always)]
    pub fn nb_bits_by_base(&self) -> u8 {
        self.description[0]
    }

    #[inline(always)]
    fn nb_bases_in_last_byte(&self) -> u8 {
        self.description[1]
    }

    /// return a u8 encoded base. It is too expensive to allocate the corresponding decoding alphabet
    /// at each call. Use an iterator instead which allocates a decoder for the whole sequence.

    pub fn get_base(&self, pos: usize) -> u8 {
        let nb_bits = self.description[0] as usize;
        match nb_bits {
            8 => self.seq[pos],
            //
            4 | 2 => {
                let mask = (1 << nb_bits) - 1;
                let nb_base_by_byte = 8 / nb_bits;
                let byte = pos / nb_base_by_byte;
                let bit_offset = nb_bits * (pos % nb_base_by_byte);
                // as base are stacked from left to right, we must shift
                mask & (self.seq[byte] >> (8 - bit_offset - nb_bits))
            }
            //
            _ => panic!(" error in Sequence::new , bad number of bits by base must be 2 4 or 8"),
        }
    } // end of get_base

    /// decompress the whole sequence, returns a whole decoded sequence in a Vec.
    /// This function should not be called intensively as it is cost a reallocation
    //   could possibly return a Sequence structure with nb_bits_by_base=8 !!
    //
    pub fn decompress(&self) -> Vec<u8> {
        let mut seqvec: Vec<u8> = Vec::new();
        match self.nb_bits_by_base() {
            // case 4 bit by base
            4 => {
                // compression rate 2
                let alfa4b = Alphabet4b::new();
                let nb_full_bytes = if self.nb_bases_in_last_byte() > 0 {
                    self.seq.len() - 1
                } else {
                    self.seq.len()
                };
                let seqlen = nb_full_bytes * 2_usize + self.nb_bases_in_last_byte() as usize;
                log::trace!("seqlen = {} ", seqlen);
                seqvec = (0..seqlen).map(|_| 0).collect();
                // get a slice of 2 u8
                let small_slice_ref: &mut [u8] = &mut [0; 2]; //  slice of length initilized by 0
                let mut pos = 0;
                // unpack full bytes
                for i in 0..nb_full_bytes {
                    alfa4b.base_unpack(self.seq[i], small_slice_ref);
                    //
                    seqvec[pos..(pos + 2)].copy_from_slice(small_slice_ref);
                    pos += 2;
                }
                // unpack last byte if necessary
                if self.nb_bases_in_last_byte() > 0 {
                    alfa4b.base_unpack(self.seq[self.seq.len() - 1], small_slice_ref);
                    // exactly one base to push
                    seqvec[pos..pos + 1].copy_from_slice(&small_slice_ref[0..1]);
                }
            }
            // case 2 bit by base
            2 => {
                // compression rate 4
                let alfa2b = Alphabet2b::new();
                let nb_full_bytes = if self.nb_bases_in_last_byte() > 0 {
                    self.seq.len() - 1
                } else {
                    self.seq.len()
                };
                let seqlen = nb_full_bytes * 4_usize + self.nb_bases_in_last_byte() as usize;
                seqvec = (0..seqlen).map(|_| 0).collect();
                // get a slice of 4 u8
                let small_slice_ref: &mut [u8] = &mut [0; 4];
                let mut pos: usize = 0;
                // unpack full bytes
                let nb_full_bytes = if self.nb_bases_in_last_byte() > 0 {
                    self.seq.len() - 1
                } else {
                    self.seq.len()
                };
                for i in 0..nb_full_bytes {
                    alfa2b.base_unpack(self.seq[i], small_slice_ref);
                    //
                    seqvec[pos..(pos + 4)].copy_from_slice(small_slice_ref);
                    pos += 4;
                }
                // unpack last byte if necessary
                if self.nb_bases_in_last_byte() > 0 {
                    alfa2b.base_unpack(self.seq[self.seq.len() - 1], small_slice_ref);
                    // we have 1 2 or 3 bases more to push
                    let nb_bases = self.nb_bases_in_last_byte() as usize;
                    seqvec[pos..pos + nb_bases].copy_from_slice(&small_slice_ref[0..nb_bases]);
                }
            }
            // case no compression.
            8 => {
                // nothing to do, except a reallocation. Use a Rc ?
                seqvec.extend_from_slice(&self.seq);
            }
            // default case
            _ => panic!(" error in Sequence::new , bad number of bits by base must be 2 4 or 8"),
        }
        //
        seqvec
    }

    /// WARNING : return length of decompressed(!!!) i.e number ob bases sequence and not self.seq.len()
    #[inline]
    pub fn size(&self) -> usize {
        if self.nb_bases_in_last_byte() > 0 {
            (self.seq.len() - 1) * (8 / self.nb_bits_by_base() as usize)
                + self.nb_bases_in_last_byte() as usize
        } else {
            self.seq.len() * (8 / self.nb_bits_by_base() as usize)
        }
    }
    /// return length of compressed seq up to last block possibly not full. So it is an upper bound of useful size
    #[inline]
    pub fn compressed_length(&self) -> usize {
        self.seq.len()
    }

    /// Get a pointer to an object implementing alphabet.
    pub fn get_alphabet(&self) -> Box<dyn BaseCompress> {
        match self.nb_bits_by_base() {
            2 => Box::new(Alphabet2b::new()),
            4 => Box::new(Alphabet4b::new()),
            8 => Box::new(Alphabet8b::new()),
            _ => panic!("incoherent coding of sequence"),
        }
    } // end get_alphabet

    /// fast version of reverse complement for 2bit encoded sequence
    // if we had a 2bit encoded sequence (as we use for Kmer16b32bit) we could have done bit manipulation
    // with shift to the right of all bytes if last byte is incomplete
    // and then bit-symetry in byte followed by xor to complement bit and transposition
    // with have also 4bit and 8bit encoded sequences!!

    fn get_reverse_complement_2bitseq(&self) -> Sequence {
        let len = self.seq.len();
        let mut rev_vec: Vec<u8> = Vec::with_capacity(len);
        let nb_bases_in_last_byte = self.nb_bases_in_last_byte();
        let shift_amount;
        let shift_mask;
        let do_shift;
        //
        if nb_bases_in_last_byte > 0 {
            do_shift = true;
            shift_amount = 8 - 2 * nb_bases_in_last_byte;
            shift_mask = (0b1 << shift_amount) - 1;
        } else {
            do_shift = false;
            shift_amount = 0;
            shift_mask = 0;
        };
        // first transpose from vec to rev_vec with bit reversing and possible shift due
        for i in 0..len {
            let byte = self.seq[len - 1 - i];
            let mut rev_byte = byte;
            if do_shift {
                // take upper part of rev_byte and get it in lower part of byte
                rev_byte >>= shift_amount;
                // we have left space in rev_byte , take the lower shift_amount bits of
                // next byte and transfer it into upper part of current
                if i < len - 1 {
                    rev_byte |= (self.seq[len - 1 - i - 1] & shift_mask) << (8 - shift_amount);
                }
            }
            // swap adjacent bloc of 2 bits
            rev_byte = (rev_byte & 0x33) << 2 | (rev_byte & 0xCC) >> 2;
            // swap adjacent bloc of 4 bits
            rev_byte = (rev_byte & 0x0F) << 4 | (rev_byte & 0xF0) >> 4;
            //
            rev_byte = !rev_byte;
            rev_vec.push(rev_byte);
        }
        //
        Sequence {
            seq: rev_vec,
            description: self.description,
        }
    } // end of get_reverse_complement_2bitseq

    /// get reverse complement of sequence.

    pub fn get_reverse_complement(&self) -> Sequence {
        //
        if self.nb_bits_by_base() == 2 {
            return self.get_reverse_complement_2bitseq();
        }
        let mut rev_vec: Vec<u8> = Vec::with_capacity(self.size());
        // get a reverse iterator
        let mut iter = IterSequence::new(self, false);
        // get alphabet used in sequence, complement base
        let alphabet = self.get_alphabet();
        // transpose
        while let Some(base) = iter.next_back() {
            let cbase = alphabet.complement(base);
            rev_vec.push(alphabet.decode(cbase));
        }
        //
        Sequence::new(&rev_vec, self.nb_bits_by_base())
    } // end of get_reverse_complement

    /// get a function giving the good decoder for our sequence
    pub fn get_decoder(&self) -> Box<dyn Fn(u8) -> u8> {
        match self.nb_bits_by_base() {
            8 => Box::new(move |b| Alphabet8b::new().decode(b)),
            4 => Box::new(move |b| Alphabet4b::new().decode(b)),
            2 => Box::new(move |b| Alphabet2b::new().decode(b)),
            _ => panic!("should not occur"),
        }
    } // end of get_decoder

    // This method avoids a copy if sequence is not compressed.
    /// count number of each bases ACGT
    /// simple utility to count base (ACGT) proportions in a sequence
    /// Non ACGT bases are counted in nb_bad
    /// fill_histogram just to avoid copying a histo on stack
    pub fn base_count(&self, histo: &mut [u64]) -> usize {
        let mut nb_bad = 0;
        for i in 0..4 {
            histo[i] = 0;
        }
        // do we have to decompress
        let decompressed;
        let sliced_seq: &[u8] = if self.nb_bits_by_base() < 8 {
            decompressed = self.decompress();
            decompressed.as_slice()
        } else {
            self.seq.as_slice()
        };
        for c in sliced_seq {
            match c {
                b'A' => {
                    histo[0] += 1;
                } // A
                b'C' => {
                    histo[1] += 1;
                } // C
                b'G' => {
                    histo[2] += 1;
                } // G
                b'T' => {
                    histo[3] += 1;
                } // T
                _ => nb_bad += 1,
            };
        }
        nb_bad
    } // end of base_count

    /// allocate a Seq with capacity to store  nb_base encoded in nb_bits by base
    pub fn with_capacity(nb_bits: u8, nb_base: usize) -> Self {
        //
        assert!(nb_bits == 2 || nb_bits == 4 || nb_bits == 8);
        //
        let alloc_len = 2 + (nb_base * nb_bits as usize) / 8;
        log::trace!("Sequence with_capacity, nb bytes allocated : {}", alloc_len);
        let vec: Vec<u8> = Vec::with_capacity(alloc_len);
        let description = [nb_bits, 0];
        Sequence {
            seq: vec,
            description,
        }
    } // end of with_capacity

    /// shrink in case the call to with_capacity was too generous
    pub fn shrink_to_fit(&mut self) {
        self.seq.shrink_to_fit();
    }

    /// This function parse buffer to_add , filters out non ACGT, encode in alpabet associated to sequence and store in sequence
    /// The argument alphabet must correspond to the number of bits / base declared in Sequence initialization
    pub fn encode_and_add(&mut self, to_add: &[u8], alphabet: &dyn BaseCompress) {
        //
        log::trace!("encode_and_add, self.vec (vec in bytes) len : {}, capacity : {}, to_add length (bases): {:?}", self.seq.len(), self.seq.capacity(), to_add.len());
        log::trace!(" to_add : {:?}", to_add);
        //
        let nb_bits: usize = self.nb_bits_by_base() as usize;
        // do we need to grow self.vec ?
        if to_add.len() >= 8 * (self.seq.capacity() - self.seq.len()) / nb_bits {
            // we allocate nb_bits times more what we need
            let grow = to_add.len() - 8 * (self.seq.capacity() - self.seq.len()) / nb_bits;
            log::trace!("allocating nb new bytes : {:?}", grow);
            self.seq.try_reserve(grow).unwrap();
            log::trace!(
                "encode_and_add  after growing,  len {}, allocated : {} ",
                self.seq.len(),
                self.seq.capacity()
            );
        }
        let seqlen = self.seq.len();
        // do we have an incomplete byte in sequence ?
        let mut nb_scanned: usize = 0;
        let mut already = self.description[1] as usize;
        let mut to_encode: u8; // the byte we work on
                               // if we had a last incomplete byte we fill it
        if already > 0 {
            to_encode = self.seq[self.seq.len() - 1];
            nb_scanned += update_byte(
                &mut to_encode,
                &mut already,
                alphabet,
                &to_add[nb_scanned..],
            );
            self.seq[seqlen - 1] = to_encode;
        }
        // now we loop storing new encoded bytes in self.seq, we store 4 bases at each update until we encounter end of to_add
        // and we possibly fill an incomplete byte at end.
        while nb_scanned < to_add.len() {
            to_encode = 0;
            already = 0;
            nb_scanned += update_byte(
                &mut to_encode,
                &mut already,
                alphabet,
                &to_add[nb_scanned..],
            );
            if already > 0 {
                // we push only if sthing was inserted.
                // (ex : if a bad base arrive at beginning of byte and at end of record we must not push the byte
                self.seq.push(to_encode);
                // and if byte is full we reset already
                if already == (8 / nb_bits) {
                    // if byte is full we reset to 0
                    already = 0;
                }
            }
        }
        // now we must reset correct info in self.description to know how many bases we have in possibly incomplete last byte
        self.description[1] = already as u8;
    } // end of encode_and_add
} // end impl Sequence

// update a byte which has inside already base in it, encoding is done by alphabet
// return (number of base in examined in to_add, nb_inserted in byte)
#[inline]
fn update_byte(
    byte: &mut u8,
    already: &mut usize,
    alphabet: &dyn BaseCompress,
    to_add: &[u8],
) -> usize {
    let nb_bits = alphabet.get_nb_bits() as usize;
    // the number of base we can add is the number of bases in to_add limited by left space in byte
    // nb_max will return the number of bases examined in to_add and so the amount of advance the caller will have to use in subsequent calls in to_add
    let nb_max = (8 / nb_bits - *already).min(to_add.len());
    let mut shift = 0usize;
    let mut inserted = 0;
    // special care for N and non ACTG. We must fill byte so we loop until byte is filled or we are at end of to_add! Nasty bug
    while shift < to_add.len() && inserted < nb_max {
        match to_add[shift] as char {
            'A' | 'C' | 'T' | 'G' => {
                let encoded = alphabet.encode(to_add[shift]);
                log::trace!(
                    "encoding : {:?}, encoded : {}, byte before  : 0x{:x}",
                    to_add[shift] as char,
                    encoded,
                    byte
                );
                *byte |= encoded << (8 - nb_bits - *already * nb_bits);
                log::trace!(
                    "encoding : {:?}, encoded : {}, byte after : 0x{:x}",
                    to_add[shift] as char,
                    encoded,
                    byte
                );
                *already += 1;
                inserted += 1;
            }
            _ => {}
        };
        shift += 1;
    }
    shift
} // end of update_byte

//=======================================
//   Ierators
//========================================

/// an iterator sequence diving into byte and bits offset

pub struct IterSequence<'a> {
    /// the sequence i am an iterator of
    myseq: &'a Sequence,
    /// mask to extract a base
    mask: u8,
    /// specifies if iterator must send decoded base if sequence is compressed or send base
    /// as it is encoded in sequence
    must_decode: bool,
    /// The function to use for decoding a a bit group to a u8 base
    decoder: Option<Box<dyn Fn(u8) -> u8>>,
    /// The current byte numbered from 0 to seq.len() - 1
    byte: usize,
    /// The bit inside current byte numbered from 0 to 7
    bit: u8,
    /// last_byte to decode in compressed seq. So that self.byte must always be <= self.last_byte
    /// and self.byte and self.bit are initialized so that they correspond to first base
    last_byte: usize,
    /// position we are in after decoding of last byte. So it is 2, 4, 6, 8.
    last_bit: u8,
}

impl<'a> IterSequence<'a> {
    /// constructor of an iterator over a sequence. We choose at construction if
    /// we want next method returning a base compressed or uncompressed.
    pub fn new(seqarg: &'a Sequence, must_decode_arg: bool) -> IterSequence<'a> {
        //
        let nb_bases_by_byte = 8 / (seqarg.nb_bits_by_base() as usize);
        // case if sequence is small than nb_bases_by_byte
        let mut last_byte = if seqarg.size() >= nb_bases_by_byte {
            (seqarg.size() / nb_bases_by_byte) - 1
        } else {
            0
        };
        let mut last_bit = 8;
        //
        let nb_bits = seqarg.nb_bits_by_base();
        let mask: u8 = if nb_bits < 8 {
            (1 << nb_bits) - 1
        } else {
            0xFF
        };
        // could use seqarg.description[1] instead of seqarg.size() % nb_bases_by_byte
        if seqarg.size() % nb_bases_by_byte > 0 {
            last_byte += 1;
            last_bit = seqarg.nb_bits_by_base() * (seqarg.size() % nb_bases_by_byte) as u8;
        }
        log::trace!(
            "IterSequence constructor seq size last byte, last bit =  {} {} {} ",
            seqarg.size(),
            last_byte,
            last_bit
        );
        //
        let decoder = match must_decode_arg {
            true => Some(seqarg.get_decoder()),
            _ => None,
        };
        //
        IterSequence {
            myseq: seqarg,
            mask,
            must_decode: must_decode_arg,
            decoder,
            byte: 0,
            bit: 0,
            last_byte,
            last_bit,
        }
    } // end of new
    /// set range of iterator. By default it is from 0 to end of sequence. But this
    /// can be modified by this method. Careful , it is a range with end excluded as in rust usage!
    pub fn set_range(&mut self, begin: usize, end: usize) -> Result<(), ()> {
        if end <= begin || end > self.myseq.size() {
            return Err(());
        }
        // we must initialize self.byte and self.bit to a position in compressed sequence corresponding begin in
        // uncompressed sequence
        let nb_bases_by_byte = 8 / (self.myseq.nb_bits_by_base() as usize);
        // set begin byte/bit
        self.byte = begin / nb_bases_by_byte;
        // recall bits are numbered from 0 to 7 left to right.
        self.bit = self.myseq.nb_bits_by_base() * (begin % nb_bases_by_byte) as u8;
        // set end byte/bit, recall indexation begins at 0
        self.last_byte = (end / nb_bases_by_byte) - 1;
        self.last_bit = 8;
        if end % nb_bases_by_byte > 0 {
            self.last_byte += 1;
            self.last_bit = self.myseq.nb_bits_by_base() * (end % nb_bases_by_byte) as u8;
        };
        //
        //        log::debug!("IterSequence set range seq size  byte,  bit =  {} {} {} ", self.myseq.size(), self.byte , self.bit);
        //        log::debug!("IterSequence set range seq size last byte, last bit =  {} {} {} ", self.myseq.size(), self.last_byte , self.last_bit);
        //
        Ok(())
    }
    /// mostly used in Kmer generation to avoid decoding/encoding of each base
    ///
    pub fn set_decoding_state(&mut self, decode: bool) {
        self.must_decode = decode;
    }
}

// remind that self.description[0] is nb_bits / base
// description[1] is number of bases in last byte.
// The state is maintained so that we know when we are at end or we can do a read.
// After each read we must adjust state to self.byte to next byte.
// We are at end if self.byte >= self.myseq.seq.len() or if self.byte == self.myseq.seq.len()-1 &&
// Consequently if we enter with bit (==7) we are at end of sequence
// otherwise we should have bit = 0 i.e we are at beginning of new byte.

/// implementation of trait Iterator for IterSequence<'a>
impl<'a> Iterator for IterSequence<'a> {
    type Item = u8;
    //
    fn next(&mut self) -> Option<u8> {
        //        log::trace!("IterSequence entering next state : byte bit  = {} {} ", self.byte , self.bit);
        //
        if self.byte > self.last_byte || (self.byte == self.last_byte && self.bit >= self.last_bit)
        {
            //            log::trace!("IterSequence end of iterator {} {} ", self.byte , self.myseq.description[1]);
            return None;
        }
        //
        let endbit = if self.byte == self.last_byte && self.myseq.description[1] > 0 {
            // case where endbit is different from 8
            self.last_bit
        } else {
            8
        };
        //
        if self.bit >= endbit {
            //            log::trace!("IterSequence end of iterator byte, bit, tail : {} {} {} ", self.byte , self.bit, self.myseq.description[1]);
            return None;
        }
        //
        // now we know we have a next !!!
        //
        // log::trace!("IterSequence next : byte bit endbit  = {} {} {} ", self.byte , self.bit, endbit);
        // we must check consistency : endbit - bit >= nb_bits. moreover we could check parity of bit
        let nb_bits = self.myseq.nb_bits_by_base();
        assert!(endbit - self.bit >= nb_bits);
        let base = self.mask & (self.myseq.seq[self.byte] >> (8 - self.bit - nb_bits));
        // adjust state
        self.bit += nb_bits;
        if self.bit == endbit {
            // we are now at end of byte, increment byte that will cause eventually next call to return None
            self.byte += 1;
            if self.byte <= self.last_byte {
                // if not above last byte reset self.bit (else we stay at end , could put it at 8
                self.bit = 0;
            }
        }
        // return decoded or not
        match self.must_decode {
            true => {
                return Some((*(self.decoder.as_ref().unwrap()))(base));
            }
            false => Some(base),
        }
    } // end of next
} // end of impl<'a> Iterator

/// implementation of trait DoubleEndedIterator  for IterSequence<'a>

// read the description for DoubleEndedIterator!! It makes the doubleendediterator algorithms fully symetric

// remind that self.description[0] is nb_bits / base
// description[1] is number of bases in last byte.

// The state is maintained so that we know when we are at end or we can do a read.
// After each read we must adjust state to self.last_byte to next byte.
// We are at end if self.last_byte < self.byte  or if
// self.byte == self.last_byte && self.last_bit <= self.bit
// Consequently if we enter with bit (==7) we are at end of sequence
// otherwise we should have bit = 0 i.e we are at beginning of new byte.

impl<'a> DoubleEndedIterator for IterSequence<'a> {
    //
    #[allow(clippy::comparison_chain)]
    fn next_back(&mut self) -> Option<u8> {
        //        println!("DoubleEndedIterator IterSequence entering next_back state : last byte,  last bit  = {} {} ", self.last_byte , self.last_bit);
        // endbit is last bit that the reverse iterator can use , by default enbit is 0
        let mut endbit: u8 = 0;
        //
        if self.byte > self.last_byte {
            // println!("DoubleEndedIterator IterSequence end of iterator {} {}  return None", self.last_byte , self.myseq.description[1]);
            return None;
        }
        //
        else if self.byte == self.last_byte {
            if self.bit >= self.last_bit {
                return None;
            }
            // if last_byte equals self.byte, end_bit is self.bit as we cannot go furster than forward iterator position
            endbit = self.bit;
        }
        //
        if self.last_bit <= endbit {
            //            println!("DoubleEndedIterator IterSequence end of iterator byte, bit, tail : {} {} {} ", self.last_byte , self.last_bit, self.myseq.description[1]);
            return None;
        }
        // we must check consistency we know that last_bit > endbit
        // we must have :  last_bit - endbit >= nb_bits. moreover we could check parity of bit
        let nb_bits = self.myseq.nb_bits_by_base();
        assert!(self.last_bit - endbit >= nb_bits);
        //
        // now we know we have a next !!!
        //        println!("IterSequence next_back : last_byte last_bit endbit  = {} {} {} ", self.last_byte , self.last_bit, endbit);
        //
        let mask = if nb_bits < 8 {
            (1 << nb_bits) - 1
        } else {
            0xFF
        };
        let base = mask & (self.myseq.seq[self.last_byte] >> (8 - self.last_bit));
        // adjust state
        self.last_bit -= nb_bits;
        if self.last_bit == endbit {
            // we are now at end of byte, decrement last_byte that will cause eventually next call to return None
            if self.last_byte > self.byte {
                self.last_byte -= 1;
                if self.byte <= self.last_byte {
                    // if not under self.byte we reset self.last_bit at beginning (in backward path) of new byte.
                    self.last_bit = 8;
                }
            }
        }
        // return decoded or not
        //        println!("next_back returning {}",  (*(self.decoder))(base) );
        match self.must_decode {
            true => {
                return Some((*(self.decoder.as_ref().unwrap()))(base));
            }
            false => Some(base),
        }
    } // end of next
} // end of impl<'a> Iterator

/// implementation of IntoIterator generates an iterator with decoding flag on !

impl<'a> IntoIterator for &'a Sequence {
    type Item = u8;
    type IntoIter = IterSequence<'a>;
    // self here is a  &'a Sequence !!
    fn into_iter(self) -> Self::IntoIter {
        println!("IntoIterator");
        IterSequence::new(self, true)
    }
} // end of impl<'a> IntoIterator<'a>

/////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;
    //
    fn log_init_test() {
        let mut builder = env_logger::Builder::from_default_env();
        //    builder.filter_level(LevelFilter::Trace);
        let _ = builder.is_test(true).try_init();
    }
    //
    #[test]
    fn test_alphabet2b() {
        log_init_test();
        let alphabet2 = Alphabet2b::new();
        assert_eq!(alphabet2.encode(b'G'), 0b10);
        assert_eq!(alphabet2.decode(0b10), b'G');
    }

    #[test]
    fn encode4b_4bases_new() {
        log_init_test();
        //
        let v = vec![b'A', b'C', b'G', b'T'];
        let seq = Sequence::new(&v, 4);
        assert!(seq.seq.len() == 2);
        // test encoding
        println!("seq[0] = {:x}", seq.seq[0]);
        assert!(seq.seq[0] == 0x12);
        println!("seq[1] = {:x}", seq.seq[1]);
        assert!(seq.seq[1] == 0x48);
        assert!(seq.seq.len() == 2);
        assert!(seq.nb_bases_in_last_byte() == 0);
        // decoding test
        let decompressed = seq.decompress();
        println!("decompressed.len() = {} ", decompressed.len());
        assert!(decompressed.len() == 4);
        //
        println!("seq[0] = {:x}", seq.seq[0]);
        println!("seq[1] = {:x}", seq.seq[1]);
        //
        assert!(decompressed[0] == b'A');
        //
        assert!(decompressed[1] == b'C');
        //
        assert!(decompressed[2] == b'G');
        //
        assert!(decompressed[3] == b'T');
        //
    }

    #[test]
    fn encode4b_5bases_new() {
        log_init_test();
        //
        let v = vec![b'A', b'C', b'G', b'T', b'A'];
        let seq = Sequence::new(&v, 4);
        // test encoding
        println!("seq[0] = 0x{:x}", seq.seq[0]);
        assert!(seq.seq[0] == 0x12);
        println!("seq[1] = 0x{:x}", seq.seq[1]);
        assert!(seq.seq[1] == 0x48);
        // only upper part of byte is encoded
        println!("seq[2] = 0x{:x}", seq.seq[2]);
        assert!(seq.seq[2] & 0xF0 == 0x10);
        assert!(seq.seq.len() == 3);
        assert!(seq.nb_bases_in_last_byte() == 1);
        // decoding test
        let decompressed = seq.decompress();
        //
        assert!(decompressed.len() == 5);
        println!("seq[0] = 0x{:x}", seq.seq[0]);
        println!("seq[1] = 0x{:x}", seq.seq[1]);
        println!("seq[2] = 0x{:x}", seq.seq[2]);
        //
        assert!(decompressed[0] == b'A');
        //
        assert!(decompressed[1] == b'C');
        //
        assert!(decompressed[2] == b'G');
        //
        assert!(decompressed[3] == b'T');
        //
        assert!(decompressed[4] == b'A');
        //
        let mut ibase = 0;
        for b in &seq {
            log::info!(" ibase base = {} {} ", ibase, b);
            assert!(b == decompressed[ibase]);
            ibase += 1;
        }
        assert!(ibase == seq.size());
    } // end of encode4b_5bases

    #[test]
    fn encode2b_5bases() {
        log_init_test();
        //
        let v = vec![b'A', b'C', b'G', b'T', b'C'];
        let seq = Sequence::new(&v, 2);
        assert!(seq.seq.len() == 2);
        assert!(seq.nb_bases_in_last_byte() == 1);
        //
        // we get encoded one byte and 2 upper bits of another byte
        //
        // first byte encode is 0b0001_1011
        //
        println!("seq[0] = {}", seq.seq[0]);
        assert!(seq.seq[0] == 0x1B);
        println!("seq[1] = {}", seq.seq[1]);
        // second byte encode must begin with 0b01    corresponding to C
        // test the first 2 bits of seq.seq[1]
        assert!((seq.seq[1] >> 6) & 0b11 == 0b01);
        //
        println!(" decompressing");
        //
        let decompressed = seq.decompress();
        //
        println!("seq[0] = {}", seq.seq[0]);
        println!("seq[1] = {}", seq.seq[1]);
        //
        assert!(seq.nb_bases_in_last_byte() == 1);
        assert!(decompressed[0] == b'A');
        //
        assert!(decompressed[1] == b'C');
        //
        assert!(decompressed[2] == b'G');
        //
        assert!(decompressed[3] == b'T');
        //
        assert!(decompressed[4] == b'C');
        //  test get_base
        let alfa2b = Alphabet2b::new();
        println!("seq.get_base(0) = {}", alfa2b.decode(seq.get_base(0)));
        println!("seq.get_base(1) = {}", alfa2b.decode(seq.get_base(1)));
        println!("seq.get_base(4) = {}", alfa2b.decode(seq.get_base(4)));
        //
        assert!(seq.get_base(0) == alfa2b.encode(b'A'));
        assert!(seq.get_base(2) == alfa2b.encode(b'G'));
        assert!(seq.get_base(4) == alfa2b.encode(b'C'));
        // testing iterator
        let mut iterseq = IterSequence::new(&seq, true);
        let mut nb_ok = 0;
        loop {
            match iterseq.next() {
                Some(b1) => {
                    println!(" ibase base = {} {} ", nb_ok, b1);
                    assert!(b1 == decompressed[nb_ok]);
                    nb_ok += 1;
                }
                None => break,
            }
        }
        assert!(nb_ok == seq.size());
        // test into iterator
        let mut ibase = 0;
        for b in &seq {
            println!(" ibase base = {} {} ", ibase, b);
            assert!(b == decompressed[ibase]);
            ibase += 1;
        }
        assert!(ibase == seq.size());
    } // end of encode2b_5bases

    #[test]
    fn encode2b_4bases() {
        log_init_test();
        //
        let v = vec![b'A', b'C', b'G', b'T'];
        let seq = Sequence::new(&v, 2);
        assert!(seq.seq.len() == 1);
        assert!(seq.nb_bases_in_last_byte() == 0);
        // we get encoded one byte and 2 upper bits of another byte
        //
        // first byte encode is 0b0001_1011
        //
        println!("seq[0] = {}", seq.seq[0]);
        assert!(seq.seq[0] == 0x1B);
        //
        println!(" decompressing");
        //
        let decompressed = seq.decompress();
        //
        println!("seq[0] = {}", seq.seq[0]);
        //
        assert!(decompressed[0] == b'A');
        //
        assert!(decompressed[1] == b'C');
        //
        assert!(decompressed[2] == b'G');
        //
        assert!(decompressed[3] == b'T');
    } // end of encode2b_4bases

    use std::string::*;
    #[test]
    // this test also reverse iteration
    fn test_reverse_complement_sequence_8() {
        log_init_test();
        //
        let seqstr = String::from("TACGAGTAGGAT");
        let slu8 = seqstr.as_bytes();
        // get a sequence with 8 bits compression
        let seq = Sequence::new(slu8, 8);
        //
        let reverse_4 = seq.get_reverse_complement();
        let reverse = reverse_4.decompress();
        let reverse_str = String::from_utf8(reverse).unwrap();
        assert_eq!(reverse_str, String::from("ATCCTACTCGTA"));
    } // end of test_reverse_complement_sequence_8

    #[test]
    // this test also reverse iteration
    fn test_reverse_complement_sequence_2_12bases() {
        log_init_test();
        //
        let seqstr = String::from("TACGAGTAGGAT");
        let slu8 = seqstr.as_bytes();
        // get a sequence with 8 bits compression
        let seq = Sequence::new(slu8, 2);
        //
        let reverse_2 = seq.get_reverse_complement();
        let reverse = reverse_2.decompress();
        let reverse_str = String::from_utf8(reverse).unwrap();
        assert_eq!(reverse_str, String::from("ATCCTACTCGTA"));
    } // end of test_reverse_complement_sequence_2_12bases

    #[test]
    // this test also reverse iteration
    fn test_reverse_complement_sequence_2_14bases() {
        log_init_test();
        //
        let seqstr = String::from("TACGAGTAGGATCC");
        let slu8 = seqstr.as_bytes();
        // get a sequence with 8 bits compression
        let seq = Sequence::new(slu8, 2);
        //
        let reverse_2 = seq.get_reverse_complement();
        let reverse = reverse_2.decompress();
        let reverse_str = String::from_utf8(reverse).unwrap();
        assert_eq!(reverse_str, String::from("GGATCCTACTCGTA"));
    } // end of test_reverse_complement_sequence_2_14bases

    #[test]
    fn test_incremental_alpha2_14bases_seq_init() {
        log_init_test();
        // a 14 base sequence
        let mut seqstr = String::from("TACGAGTAGGATCC");
        let to_add = seqstr.as_bytes();
        // initialize with not enough!
        let mut seq_tocheck = Sequence::with_capacity(2, 4);
        //
        let alpha2b = Alphabet2b::new();
        //
        seq_tocheck.encode_and_add(to_add, &alpha2b);
        // we compare with normal initialization, check also for last byte
        let restored_str = String::from_utf8(seq_tocheck.decompress()).unwrap();
        assert_eq!(restored_str, seqstr);
        assert_eq!(seq_tocheck.nb_bases_in_last_byte(), 2);
        // now we test adding once more ...
        let seqstr2 = String::from("AAAGG");
        let to_add = seqstr2.as_bytes();
        seq_tocheck.encode_and_add(to_add, &alpha2b);
        // we add 5 bases
        assert_eq!(seq_tocheck.nb_bases_in_last_byte(), 3);
        let restored_str = String::from_utf8(seq_tocheck.decompress()).unwrap();
        seqstr.push_str(&seqstr2);
        assert_eq!(restored_str, seqstr);
    } // end of test_incremental_14b_seq_init

    #[test]
    fn test_encode_and_add_with_n() {
        log_init_test();
        // a 14 base sequence
        let seqstr = String::from("TCNGCAGTTGGATCCC");
        let to_add = seqstr.as_bytes();
        let mut seq_tocheck = Sequence::with_capacity(2, 4);
        //
        let alpha2b = Alphabet2b::new();
        //
        seq_tocheck.encode_and_add(to_add, &alpha2b);
        // we compare with normal initialization, check also for last byte, last byte contains 3
        let restored_str = String::from_utf8(seq_tocheck.decompress()).unwrap();
        assert_eq!(restored_str, String::from("TCGCAGTTGGATCCC"));
    } // end of test_encode_and_add_with_n

    #[test]
    fn test_encode_and_add_very_small_seq() {
        log_init_test();
        // a 14 base sequence
        let seqstr = String::from("TC");
        let to_add = seqstr.as_bytes();
        let mut seq_tocheck = Sequence::with_capacity(2, 4);
        //
        let alpha2b = Alphabet2b::new();
        //
        seq_tocheck.encode_and_add(to_add, &alpha2b);
        // we compare with normal initialization, check also for last byte, last byte contains 3
        let restored_str = String::from_utf8(seq_tocheck.decompress()).unwrap();
        assert_eq!(restored_str, String::from("TC"));
    } // end of test_encode_and_add_with_n

    #[test]
    fn test_incremental_alpha2_15bases_seq_init() {
        log_init_test();
        // a 14 base sequence
        let mut seqstr = String::from("TACGAGTAGGATCCC");
        let to_add = seqstr.as_bytes();
        // initialize with not enough!
        let mut seq_tocheck = Sequence::with_capacity(2, 4);
        //
        let alpha2b = Alphabet2b::new();
        //
        seq_tocheck.encode_and_add(to_add, &alpha2b);
        // we compare with normal initialization, check also for last byte, last byte contains 3
        let restored_str = String::from_utf8(seq_tocheck.decompress()).unwrap();
        assert_eq!(restored_str, seqstr);
        assert_eq!(seq_tocheck.nb_bases_in_last_byte(), 3);
        assert_eq!(seq_tocheck.size(), seqstr.len());
        // now we test adding once more ...
        let seqstr2 = String::from("AAAGG");
        let to_add = seqstr2.as_bytes();
        seq_tocheck.encode_and_add(to_add, &alpha2b);
        // we add 5 bases we get a 20 base sequence, the last byte is full
        log::info!("adding {:?}", seqstr2);
        assert_eq!(seq_tocheck.nb_bases_in_last_byte(), 0);
        assert_eq!(seq_tocheck.size(), seqstr.len() + seqstr2.len());
        let restored_str = String::from_utf8(seq_tocheck.decompress()).unwrap();
        seqstr.push_str(&seqstr2);
        assert_eq!(restored_str, seqstr);
    } // end of test_incremental_15b_seq_init

    #[test]
    fn encode4b_5bases_incr() {
        log_init_test();
        //
        let v = vec![b'A', b'C', b'G', b'T', b'A'];
        let mut seq = Sequence::with_capacity(4, 4);
        let alpha4b = Alphabet4b::new();
        seq.encode_and_add(&v, &alpha4b);
        // test encoding
        println!("seq[0] = 0x{:x}", seq.seq[0]);
        assert!(seq.seq[0] == 0x12);
        println!("seq[1] = 0x{:x}", seq.seq[1]);
        assert!(seq.seq[1] == 0x48);
        // only upper part of byte is encoded
        println!("seq[2] = 0x{:x}", seq.seq[2]);
        assert!(seq.seq[2] & 0xF0 == 0x10);
        assert!(seq.seq.len() == 3);
        assert!(seq.nb_bases_in_last_byte() == 1);
        // decoding test
        let decompressed = seq.decompress();
        //
        assert!(decompressed.len() == 5);
        println!("seq[0] = {:0x}", seq.seq[0]);
        println!("seq[1] = 0x{:x}", seq.seq[1]);
        println!("seq[2] = 0x{:x}", seq.seq[2]);
        //
        assert!(decompressed[0] == b'A');
        //
        assert!(decompressed[1] == b'C');
        //
        assert!(decompressed[2] == b'G');
        //
        assert!(decompressed[3] == b'T');
        //
        assert!(decompressed[4] == b'A');
        //
        let mut ibase = 0;
        for b in &seq {
            log::info!(" ibase base = {} {} ", ibase, b);
            assert!(b == decompressed[ibase]);
            ibase += 1;
        }
        assert!(ibase == seq.size());
    } // end of encode4b_5bases

    // a test for Alphabet4b although we do not use it
    #[test]
    fn test_incremental_alpha4_15bases_seq_init() {
        log_init_test();
        // a 15 base sequence
        let mut seqstr = String::from("TACGAGTAGGATCCC");
        let to_add = seqstr.as_bytes();
        //
        let mut seq_tocheck = Sequence::with_capacity(4, 4);
        //
        let alpha4b = Alphabet4b::new();
        //
        seq_tocheck.encode_and_add(to_add, &alpha4b);
        // we compare with normal initialization, check also for last byte, last byte contains 3
        let restored_str = String::from_utf8(seq_tocheck.decompress()).unwrap();
        assert_eq!(restored_str, seqstr);
        assert_eq!(seq_tocheck.nb_bases_in_last_byte(), 1);
        // now we test adding once more ...
        let seqstr2 = String::from("AAAGG");
        let to_add = seqstr2.as_bytes();
        seq_tocheck.encode_and_add(to_add, &alpha4b);
        // we add 5 bases we get a 20 base sequence, the last byte is full
        log::info!("adding {:?}", seqstr2);
        assert_eq!(seq_tocheck.nb_bases_in_last_byte(), 0);
        let restored_str = String::from_utf8(seq_tocheck.decompress()).unwrap();
        seqstr.push_str(&seqstr2);
        assert_eq!(restored_str, seqstr);
    } // end of test_incremental_15b_seq_init
} // end module test
