//! This file contains structure for Kmer generation from a compressed sequence
//! Essentially we provide for iterators which extract kmer from a range in a sequence
//
use log::{Level, trace};


use fnv::{FnvHashMap, FnvBuildHasher};

/// an IndexMap used for Kmer counting


pub use super::{kmertraits::*, kmer::*, sequence::*, kmer32bit::Kmer32bit, kmer16b32bit::Kmer16b32bit, kmer64bit::Kmer64bit};


pub trait KmerSeqIteratorT {
    /// Kmer32bit, Kmer16b32bit, Kmer64bit
    type KmerVal;
    /// get next kmer or None
    fn next(&mut self) -> Option<Self::KmerVal>;
}


// This structure must take care of encoding conversions
// as all Kmer must be 2 bits encoded
// The generator takes as argument size of kmer and a sequence to iterate over.
// The structure KmerSeqIterator is dependant upon implementation of Kmer representation.

/// a Kmer iterator over a sequence.
///
/// usage: see examples in tests
///
/// The structure [KmerGenerationPattern] is there to provide absraction over it.

pub struct KmerSeqIterator<'a , T > where T:CompressedKmerT + KmerBuilder<T> {
    /// size of kmer
    nb_base: u8,
    /// an iterator for base calling
    seqiter: IterSequence<'a>,
    /// last kmer returned. At the beginning its None
    previous: Option<T>,
}



impl<'a, T> KmerSeqIterator<'a, T>  where T:CompressedKmerT + KmerBuilder<T> {
    /// Constructor for a given sequence and kmersize
    pub fn new(ksize: u8, sequence: &'a Sequence) -> KmerSeqIterator<'a, T> {
        if ksize as usize > T::get_nb_base_max() {
            panic!("\n KmerSeqIterator cannot support so many bases for given kmer type, kmer size  {}", ksize);
        }
        let seqiter_arg = IterSequence::new(sequence, false);
        KmerSeqIterator{nb_base: ksize, seqiter:seqiter_arg, previous:None}
    } // end of new
    /// Set the range from which all kmer of a given size are to be extracted from the sequence associated to the iterator.
    pub fn set_range(&mut self, begin: usize, end: usize) -> std::result::Result<(),()> {
        self.seqiter.set_range(begin, end)
    }
} // end of impl for KmerSeqIterator



// KmerSeqIterator impl over generic Kmer

impl <'a, Kmer>  KmerSeqIteratorT for KmerSeqIterator<'a, Kmer> 
        where Kmer : CompressedKmerT + KmerBuilder<Kmer> {
    
    type KmerVal = Kmer;

    fn next(&mut self) -> Option<Kmer> {
       // check for end of iterator
       let next_base;
       match self.seqiter.next() {
           Some(b) => next_base = b,
           None => return None,
       }
        // now we know we are not at end of iterator
        // if we do not have a previous we have to contruct first kmer
        // we have to push a base.
        //
        if let Some(kmer) = self.previous {
            // in fact we have the base to push
            self.previous = Some(kmer.push(next_base));
            return self.previous;
        }
        else {
            // we are at beginning of kmer construction sequence we have first base
            // we need to place first base at the correct place.
            let kmer_size = self.nb_base as usize;            
            let pos = 2*(kmer_size -1);
            let mut new_kmer_val = <Kmer as CompressedKmerT>::Val::from(next_base) << pos;
            for i in 0..(kmer_size-1) {
                if let Some(next_base) = self.seqiter.next()  {
                    let base_val = <Kmer as CompressedKmerT>::Val::from(next_base) << (pos - 2 - 2*i);
                    new_kmer_val = new_kmer_val | base_val;
                }
                else {
                    return None;
                }
            } // end of for
            let new_kmer: Kmer = <Kmer as KmerBuilder<Kmer>>::build(new_kmer_val, self.nb_base);
            self.previous = Some(new_kmer);
            return Some(new_kmer);
        }
    }  // end of next

}



//=================== trait for kmer generation pattern ========================//



/// overload Kmer generation for various Kmer types. (see Rust Week 225).  
/// We define a generic trait defining a pattern for target type of kmer generation.  
/// NOTE: The sequence must be encode with 2 bit alphabet to generate compressed Kmer on 2 bits.  
/// The Kmer must have an encoding larger than the one used in the sequence encoding!  
/// **A panic is provoked in other cases!!**
/// 
pub trait KmerGenerationPattern<T:KmerT> {
    /// generate all kmers included in 0..
    /// Be cautious that all kmers are generated along the sequence. 
    /// So the length of the returned Vec is nearly the len of the sequence, 
    /// and for a long sequence this can be billions of Kmers. 
    /// For long sequence it can be useful to use generate_kmer_distribution
    fn generate_kmer_pattern(&self, seq : & Sequence) -> Vec<T>;
    /// generate all kmers inclused in begin..end with end excluded as in rust conventions.
    /// Length of returned vector corresponds to range length. 
    fn generate_kmer_pattern_in_range(&self, seq : & Sequence, begin:usize, end:usize) -> Vec<T>;   
    /// generate kmers with their multiplicities. 
    /// The length of returned Vec is n is the number of different Kmers in the sequence.
    /// **Note : We do not expect a sequence with a Kmer occuring more than 2^32 times!. This would cause an overflow!**
    fn generate_kmer_distribution(&self, seq : & Sequence) -> FnvHashMap<T,u32>;
}



use std::marker::PhantomData;

// a  structure just providing reference to type of generation pattern to use Kmer to generate
// The structure could be empty. It is the method generate_kmer that really counts.
// Via type inference we go from T to KmerGenerationPattern<'a, T> and so to the generation of correct T kmer
// Finally we put back T in KmerGenerator ??

/// This structure stores the size of kmer and provides overload of kmer generation for different types of Kmer.
///
///  usage: see an example in tests::test_generate_kmer16b32bit_pattern
///
/// 
/// NOTE: The sequence must be encoded with 2 bit alphabet to generate compressed Kmer on 2 bits.
///       The Kmer must have an encoding larger than the one used in the sequence encoding! 
///       A panic is provoked in other cases!!

pub struct KmerGenerator<T:KmerT> {
    /// size of kmer we generate
    pub kmer_size : u8,
    t_marker: PhantomData<T>,
}


impl  <T:KmerT> KmerGenerator<T> {
    pub fn new(ksize:u8) -> Self {
        KmerGenerator{kmer_size: ksize, t_marker : PhantomData}
    }
    /// generic driver for kmer generation
    pub fn generate_kmer (&self, seq : &Sequence) -> Vec<T> where Self: KmerGenerationPattern<T> {
        self.generate_kmer_pattern(seq)
    }
    /// generic driver for kmer generation
    pub fn generate_kmer_in_range(&self, seq : & Sequence, begin:usize, end:usize) -> Vec<T>
    where Self: KmerGenerationPattern<T> {
        self.generate_kmer_pattern_in_range(seq, begin, end)
    }
    /// generic driver for kmer distribution pattern
    pub fn generate_weighted_kmer(&self, seq : &Sequence) -> FnvHashMap<T,u32>  where Self : KmerGenerationPattern<T> {
        self.generate_kmer_distribution(seq)
    }
    ///
    pub fn get_kmer_size(&self) -> usize { self.kmer_size as usize}
}  // end of impl KmerGenerator


/// A utility to convert FnvHashMap<T,usize> to Vec<(T, usize)> 
pub fn hashmap_count_to_vec_count<T:CompressedKmerT+ std::hash::Hash>(kmer_distribution: &FnvHashMap<T,u32>) -> Vec<(T, u32)> {
        // convert to a Vec
        let mut hashed_kmers = kmer_distribution.keys();
        let mut weighted_kmer = Vec::<(T,u32)>::with_capacity(kmer_distribution.len());
        loop {
            match hashed_kmers.next() {
                Some(key) => {
                    if let Some(weight) = kmer_distribution.get(key) {
                        // get back to Kmer16b32bit from 
                        weighted_kmer.push((*key,*weight));
                    };
                },
                None => break,
            }
        }
        //
    return weighted_kmer;
}   // end of hashmap_count_to_vec_count



// We need a guess to allocate HashMap used with Kmer Generation
// for very long sequence we must avoid nb_kmer to sequence length! Find a  good heuristic
pub(super) fn get_nbkmer_guess(seq : &Sequence) -> usize {
    let nb = 10_000_000 * (1usize + seq.size().ilog2() as usize);
    let nb_kmer = seq.size().min(nb);
    return nb_kmer;
} // end of get_nbkmer_guess


// ======================= implementation for Kmer16b32bit ====================


/// implementation of kmer generation pattern for Kmer16b32bit
impl KmerGenerationPattern<Kmer16b32bit> for KmerGenerator<Kmer16b32bit> {
    fn generate_kmer_pattern(&self, seq : &Sequence) -> Vec<Kmer16b32bit> {
        if self.kmer_size != 16u8 {
            panic!("Kmer16b32bit has 16 bases!!");
        }
        if seq.nb_bits_by_base() != 2 {
            panic!("Sequence must be 2-bit encoded for 2-bit compressed kmer generation!!");
        }
        // For a sequence of size the number of kmer is seq.size - kmer.size + 1  !!!
        // But it happens that "long reads" are really short 
        let nb_kmer = if seq.size() >= 16 { seq.size()-16+1} else {0};
        let nb_kmer = nb_kmer.min(get_nbkmer_guess(seq));
        let mut kmer_vect = Vec::<Kmer16b32bit>::with_capacity(nb_kmer);
        let mut kmeriter  = KmerSeqIterator::<Kmer16b32bit>::new(self.kmer_size, seq);
        loop {
            match kmeriter.next() {
                Some(kmer) => kmer_vect.push(kmer),
                None => break,
            }
        }
        //
        return kmer_vect;
    }  // end of generate_kmer_pattern


    /// generate all kmers associated to their multiplicity
    /// This is useful in the context of Jaccard Probability Index estimated with ProbminHash 
    /// We are here with Kmers less than 15 bases, so we have at most 4**15 = 2**30 different kmers
    /// and so the multiplicity of a kmer cannot exceed the size of u32 used in FnvHashMap
    fn generate_kmer_distribution(&self, seq : &Sequence) -> FnvHashMap<Kmer16b32bit,u32> {
        if self.kmer_size != 16u8 {
            panic!("Kmer16b32bit has 16 bases!!");
        }
        if seq.nb_bits_by_base() != 2 {
            panic!("Sequence must be 2-bit encoded for 2-bit compressed kmer generation!!");
        }
        // For a sequence of size the number of kmer is seq.size - kmer.size + 1  !!!
        // But it happens that "long reads" are really short 
        let nb_kmer = if seq.size() >= 16 { seq.size()-16+1} else {0};
        let nb_kmer = nb_kmer.min(get_nbkmer_guess(seq));
        let mut kmer_distribution : FnvHashMap::<Kmer16b32bit,u32> = FnvHashMap::with_capacity_and_hasher(nb_kmer, FnvBuildHasher::default());
        let mut kmeriter = KmerSeqIterator::<Kmer16b32bit>::new(self.kmer_size, seq);
        loop {
            match kmeriter.next(){
                Some(kmer) => {
                    // we must convert Kmer16b32bit to u64 and be able to retrieve the original Kmer16b32bit
                    *kmer_distribution.entry(kmer).or_insert(0) += 1;
                },
                None => break,
            }
        }
        //
        return kmer_distribution;
    }  // end of generate_kmer_pattern



    fn generate_kmer_pattern_in_range(&self, seq : &Sequence, begin:usize, end:usize) -> Vec<Kmer16b32bit> {
        if self.kmer_size != 16u8 {
            panic!("Kmer16b32bit has 16 bases!!");
        }
        if seq.nb_bits_by_base() != 2 {
            panic!("Sequence must be 2-bit encoded for 2-bit compressed kmer generation!!");
        }
        if begin >= end {
            panic!("KmerGenerationPattern<'a, Kmer16b32bit>  bad range for kmer iteration");
        }
        // For a sequence of size the number of kmer is seq.size - kmer.size + 1  !!!
        // But it happens that "long reads" are really short 
        let nb_kmer = if seq.size() >= 16 { seq.size()-16+1} else {0};
        // TODO for very long sequence we must avoid nb_kmer to sequence length! Find a  good heuristic
        let nb_kmer = nb_kmer.min(get_nbkmer_guess(seq));
        let mut kmer_vect = Vec::<Kmer16b32bit>::with_capacity(nb_kmer);
        let mut kmeriter = KmerSeqIterator::<Kmer16b32bit>::new(self.kmer_size, seq);
        kmeriter.set_range(begin, end).unwrap();
        loop {
            match kmeriter.next() {
                Some(kmer) => kmer_vect.push(kmer),
                None => break,
            }
        }
        //
        return kmer_vect;
    }  // end of generate_kmer_pattern

}  // end of impl KmerGenerationPattern<'a, Kmer16b32bit>



// ======================= implementation for Kmer32bit ====================

/// implementation of kmer generation pattern for Kmer32bit

impl<'a> KmerGenerationPattern<Kmer32bit> for KmerGenerator<Kmer32bit> {
    fn generate_kmer_pattern(&self, seq : &Sequence) -> Vec<Kmer32bit> {
        if self.kmer_size > 14u8 {
            panic!("Kmer16b32bit has less than 14 bases!!");
        }
        if seq.nb_bits_by_base() != 2 {
            panic!("Sequence must be 2-bit encoded for 2-bit compressed kmer generation!!");
        }
        // For a sequence of size the number of kmer is seq.size - kmer.size + 1  !!!
        // But it happens that "long reads" are really short
        let kmer_size = self.kmer_size as usize;
        let nb_kmer = if seq.size() >= kmer_size { seq.size()- kmer_size+1} else {0};
        let nb_kmer = nb_kmer.min(get_nbkmer_guess(seq));
        let mut kmer_vect = Vec::<Kmer32bit>::with_capacity(nb_kmer);
        let mut kmeriter = KmerSeqIterator::<Kmer32bit>::new(self.kmer_size, seq);
        loop {
            match kmeriter.next(){
                Some(kmer) => kmer_vect.push(kmer),
                None => break,
            }
        }
        //
        return kmer_vect;
    }  // end of generate_kmer_pattern


    /// generate all kmers associated to their multiplicity
    /// This is useful in the context of Jaccard Probability Index estimated with ProbminHash
    /// We are here with Kmers less than 15 bases, so we have at most 4**15 = 2**30 different kmers
    /// and so the multiplicity of a kmer cannot exceed the size of u32 used in FnvHashMap
    fn generate_kmer_distribution(&self, seq : &Sequence) -> FnvHashMap<Kmer32bit,u32> {
        if self.kmer_size > 14u8 {
            panic!("Kmer32bit has more than 14 bases!!");
        }
        if seq.nb_bits_by_base() != 2 {
            panic!("Sequence must be 2-bit encoded for 2-bit compressed kmer generation!!");
        }
        // For a sequence of size the number of kmer is seq.size - kmer.size + 1  !!!
        // But it happens that "long reads" are really short 
        let nb_kmer = if seq.size() >= 16 { seq.size()-16+1} else {0};
        let nb_kmer = nb_kmer.min(get_nbkmer_guess(seq));
        let mut kmer_distribution : FnvHashMap::<Kmer32bit,u32> = FnvHashMap::with_capacity_and_hasher(nb_kmer, FnvBuildHasher::default());
        let mut kmeriter = KmerSeqIterator::<Kmer32bit>::new(self.kmer_size, seq);
        loop {
            match kmeriter.next(){
                Some(kmer) => {
                    trace!(" storing kmer {} ", String::from_utf8_lossy(kmer.get_uncompressed_kmer().as_slice()));                    // we must convert Kmer64bit to u64 and be able to retrieve the original Kmer64bit
                    *kmer_distribution.entry(kmer).or_insert(0) += 1;
                },
                None => break,
            }
        }

        //
        return kmer_distribution;
    }  // end of generate_kmer_pattern



    fn generate_kmer_pattern_in_range(&self, seq : &Sequence, begin:usize, end:usize) -> Vec<Kmer32bit> {
        if self.kmer_size > 14u8 {
            panic!("Kmer16b32bit has less than 14 bases!!");
        }
        if seq.nb_bits_by_base() != 2 {
            panic!("Sequence must be 2-bit encoded for 2-bit compressed kmer generation!!");
        }
        //
        if end <= begin {
            panic!("KmerGenerationPattern<'a, Kmer32bit>:generate_kmer_pattern_in_range incoherent range");
        }
        // For a sequence of size the number of kmer is seq.size - kmer.size + 1  !!!
        // But it happens that "long reads" are really short
        let kmer_size = self.kmer_size as usize;
        let nb_kmer = if seq.size() >= kmer_size { seq.size()- kmer_size+1} else {0};
        let nb_kmer = nb_kmer.min(get_nbkmer_guess(seq));
        let mut kmer_vect = Vec::<Kmer32bit>::with_capacity(nb_kmer);
        let mut kmeriter = KmerSeqIterator::<Kmer32bit>::new(self.kmer_size, seq);
        kmeriter.set_range(begin, end).unwrap();
        // as we have set range in kmeriter we iter as long as we get a kmer
        loop {
            match kmeriter.next(){
                Some(kmer) => {
                    trace!(" storing kmer {} ", String::from_utf8_lossy(kmer.get_uncompressed_kmer().as_slice()));                    // we must convert Kmer64bit to u64 and be able to retrieve the original Kmer64bit
                    kmer_vect.push(kmer);
                },
                None => break,
            }
        }
        //
        return kmer_vect;
    }  // end of generate_kmer_pattern_in_range

    
}  // end of impl<'a> KmerGenerationPattern<'a, Kmer32bit>

// ======================= implementation for Kmer64bit ====================


impl KmerGenerationPattern<Kmer64bit> for KmerGenerator<Kmer64bit> {
    
    fn generate_kmer_pattern(&self, seq : & Sequence) -> Vec<Kmer64bit> {
        if self.kmer_size > 32u8 {
            panic!("Kmer64bit has less than 32 bases!!");
        }
        if seq.nb_bits_by_base() != 2 {
            panic!("Sequence must be 2-bit encoded for 2-bit compressed kmer generation!!");
        }
        // For a sequence of size the number of kmer is seq.size - kmer.size + 1  !!!
        // But it happens that "long reads" are really short
        let mut nb_generated : u64 = 0;
        let kmer_size = self.kmer_size as usize;
        let nb_kmer = if seq.size() >= kmer_size { seq.size()- kmer_size+1} else {0};
        let nb_kmer = nb_kmer.min(get_nbkmer_guess(seq));
        let mut kmer_vect = Vec::<Kmer64bit>::with_capacity(nb_kmer);
        let mut kmeriter = KmerSeqIterator::<Kmer64bit>::new(self.kmer_size, seq);
        loop {
            match kmeriter.next(){
                Some(kmer) =>  { 
                    kmer_vect.push(kmer);
                    nb_generated += 1;
                    if log::log_enabled!(Level::Debug) {
                        if nb_generated % 1_000_000 == 0 {
                            log::debug!("nb kmer generated  : {}", nb_generated);
                        }
                    }
                }
                None => break,
            }
        }
        //
        return kmer_vect;
    }  // end of generate_kmer_pattern


    /// generate all kmers associated to their multiplicity
    /// This is useful in the context of Jaccard Probability Index estimated with ProbminHash 
    /// Note : The multiplicity of a given kmer is limited by the size of u32. (must be less than 2**32)
    fn generate_kmer_distribution(&self, seq : &Sequence) -> FnvHashMap<Kmer64bit,u32> {
        if self.kmer_size > 32u8 {
            panic!("Kmer64bit has less than 32 bases!!");
        }
        if seq.nb_bits_by_base() != 2 {
            panic!("Sequence must be 2-bit encoded for 2-bit compressed kmer generation!!");
        }
        // For a sequence of size the number of kmer is seq.size - kmer.size + 1  !!!
        // But it happens that "long reads" are really short 
        let mut nb_generated : u64 = 0;
        let kmer_size = self.kmer_size as usize;
        let nb_kmer = if seq.size() >= kmer_size { seq.size()- kmer_size+1} else {0};
        let nb_kmer = nb_kmer.min(get_nbkmer_guess(seq));
        let mut kmer_distribution : FnvHashMap::<Kmer64bit,u32> = FnvHashMap::with_capacity_and_hasher(nb_kmer, FnvBuildHasher::default());
        let mut kmeriter = KmerSeqIterator::<Kmer64bit>::new(self.kmer_size, seq);
        let mut nb_base = 0;
        loop {
            match kmeriter.next(){
                Some(kmer) => {
                    trace!(" storing {} ", String::from_utf8_lossy(kmer.get_uncompressed_kmer().as_slice()));                    // we must convert Kmer64bit to u64 and be able to retrieve the original Kmer64bit
                    *kmer_distribution.entry(kmer).or_insert(0) += 1;
                    if nb_base == 0 {
                        nb_base = kmer.1;
                    }
                    nb_generated += 1;
                    if log::log_enabled!(Level::Debug) {
                        if nb_generated % 1_000_000 == 0 {
                            log::debug!("nb kmer generated  : {}, nb_different kmers {}", nb_generated, kmer_distribution.len());
                        }
                    }
                },
                None => break,
            }
        }   
        //
        return kmer_distribution;
    }  // end of generate_kmer_pattern



    /// Note : The multiplicity of a given kmer is limited by the size of u32. (must be less than 2**32)
    fn generate_kmer_pattern_in_range(&self, seq : &Sequence, begin:usize, end:usize) -> Vec<Kmer64bit> {
        if self.kmer_size > 32u8 {
            panic!("Kmer64bit has less than 32 bases!!");
        }
        if seq.nb_bits_by_base() != 2 {
            panic!("Sequence must be 2-bit encoded for 2-bit compressed kmer generation!!");
        }
        if end <= begin {
            panic!("KmerGenerationPattern<'a, Kmer64bit>:generate_kmer_pattern_in_range incoherent range");
        }
        // For a sequence of size the number of kmer is seq.size - kmer.size + 1  !!!
        // But it happens that "long reads" are really short
        let kmer_size = self.kmer_size as usize;
        let nb_kmer = if seq.size() >= kmer_size { seq.size()- kmer_size+1} else {0};
        let nb_kmer = nb_kmer.min(get_nbkmer_guess(seq));
        let mut kmer_vect = Vec::<Kmer64bit>::with_capacity(nb_kmer);
        let mut kmeriter = KmerSeqIterator::<Kmer64bit>::new(self.kmer_size, seq);
        kmeriter.set_range(begin, end).unwrap();
        // as we have set range in kmeriter we iter as long as we get a kmer
        loop {
            match kmeriter.next(){
                Some(kmer) => kmer_vect.push(kmer),
                None => break,
            }
        }
        //
        return kmer_vect;
    }  // end of generate_kmer_pattern_in_range

    
} // end of impl KmerGenerationPattern


//===========================================================================

// function to be called from any where
pub fn generate_all_kmer16b32bit(seqvec : &Vec<Sequence>) {
    //
    println!(" in generating kmer ...");
    let start_t = std::time::Instant::now();

    let mut nb_kmer : u64 = 0;
    
    for seq in seqvec {
        let vkmer : Vec<Kmer16b32bit> = KmerGenerator::new(16 as u8).generate_kmer(&seq);
        nb_kmer += vkmer.len() as u64;
    }
    
    let elapsed_t = start_t.elapsed().as_secs();
    println!(" elapsed time (s) in kmer generation {} ", elapsed_t);
    println!(" nb kmer generated {}", nb_kmer);
}


pub fn generate_all_kmer32bit(kmer_size:u8, seqvec : &Vec<Sequence>) {
    println!(" in generating kmer ...");
    let start_t = std::time::Instant::now();

    let mut nb_kmer : u64 = 0;
    
    for seq in seqvec {
        let vkmer : Vec<Kmer32bit> = KmerGenerator::new(kmer_size).generate_kmer(&seq);
        nb_kmer += vkmer.len() as u64;
    }
    
    let elapsed_t = start_t.elapsed().as_secs();
    println!(" elapsed time (s) in kmer generation {} ", elapsed_t);
    println!(" nb kmer generated {}", nb_kmer);
}


pub fn generate_all_kmer64bit(kmer_size:u8, seqvec : &Vec<Sequence>) {
    println!(" in generating kmer ...");
    let start_t = std::time::Instant::now();

    let mut nb_kmer : u64 = 0;
    
    for seq in seqvec {
        let vkmer : Vec<Kmer64bit> = KmerGenerator::new(kmer_size).generate_kmer(&seq);
        nb_kmer += vkmer.len() as u64;
    }
    
    let elapsed_t = start_t.elapsed().as_secs();
    println!(" elapsed time (s) in kmer generation {} ", elapsed_t);
    println!(" nb kmer generated {}", nb_kmer);
}

//=======================================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[allow(dead_code)]
    fn log_init() {
        let _ = env_logger::builder().is_test(true).try_init();
    }
        //
    //   11010000001010100000000100111101
    //    T C A A A G G G A A A C A T T C
    //   In this test all bytes are complete
    #[test]
    fn test_gen_kmer16b32bit_80bases() {
        log_init();
        // got a string of 80 bases shoud have 66 Kmers OK
        let seqstr = String::from("TCAAAGGGAAACATTCAAAATCAGTATGCGCCCGTTCAGTTACGTATTGCTCTCGCTAATGAGATGGGCTGGGTACAGAG");
        log::info!("test_gen_kmer16b32bit_80bases generating 16b kmer from : {}", seqstr);
        // transform to a Vec<u8> and then a Sequence
        let slu8 = seqstr.as_bytes();
        // get a sequence with 2 bits compression
        let seq = Sequence::new(&slu8,2);
        let mut kmergen32 = KmerSeqIterator::<Kmer16b32bit>::new(16, &seq);
        //
        for i in 0..seqstr.len()-16+1 {
            let kmer = kmergen32.next().unwrap();
            println!("kmer i = {}  {:b}  " , i , kmer.0);
            let v = kmer.get_uncompressed_kmer();
            println!("kmer i = {}  {:?}  " , i , v);
            // recall from_utf8_lossy return a Cow('a,str)
            let str = String::from_utf8_lossy(v.as_slice());
            println!("kmer i = {}  {}  " , i , str);
            assert_eq!(str, seqstr[i..i+16]);
        }        
    }  // end of test_gen_kmer


    
    #[test]
    fn test_gen_kmer16b32bit_50bases() {
        log_init();
        // got a string of 80 bases shoud have 66 Kmers OK
        let seqstr = String::from("TCAAAGGGAAACATTCAAAATCAGTATGCGCCCGTTCAGTTACGTATTGC");
        log::info!("test_gen_kmer16b32bit_50bases generating 16b kmer from : {}", seqstr);
        // transform to a Vec<u8> and then a Sequence
        let slu8 = seqstr.as_bytes();
        // get a sequence with 2 bits compression
        let seq = Sequence::new(&slu8,2);
        let mut kmergen32 = KmerSeqIterator::<Kmer16b32bit>::new(16, &seq);
        //
        for i in 0..seqstr.len()-16+1 {
            let kmer = kmergen32.next().unwrap();
            println!("kmer i = {}  {:b}  " , i , kmer.0);
            let v = kmer.get_uncompressed_kmer();
            println!("kmer i = {}  {:?}  " , i , v);
            // recall from_utf8_lossy return a Cow('a,str)
            let str = String::from_utf8_lossy(v.as_slice());
            println!("kmer i = {}  {}  " , i , str);
            assert_eq!(str, seqstr[i..i+16]);
        }
        // check iterator sees the end
        match kmergen32.next() {
            Some(kmer) => {
                let v = kmer.get_uncompressed_kmer();
                println!("kmer =  {:?}  " , v);
                panic!("iterator do not see end");
            },
            None => (),
        } // end match
    }  // end of test_gen_kmer


    #[test]
    fn test_gen_kmer16b32bit_50bases_range_iterator() {
        log_init();
        // got a string of 50 bases shoud have 66 Kmers OK
        let seqstr = String::from("TCAAAGGGAAACATTCAAAATCAGTATGCGCCCGTTCAGTTACGTATTGC");
        // transform to a Vec<u8> and then a Sequence
        let slu8 = seqstr.as_bytes();
        // get a sequence with 2 bits compression
        let seq = Sequence::new(&slu8,2);
        let mut kmergen32 = KmerSeqIterator::<Kmer16b32bit>::new(16, &seq);
        //
        kmergen32.set_range(3,25).unwrap();
        //
        for i in 3..25-16+1 {
            match kmergen32.next() {
                Some(kmer) => {
                    println!("kmer i = {}" , i);
                    let v = kmer.get_uncompressed_kmer();
                    println!("          kmer =  {:?}  " , v);
                    // recall from_utf8_lossy return a Cow('a,str)
                    let str = String::from_utf8_lossy(v.as_slice());
                    println!("          kmer i = {}  {}  " , i , str);
                    assert_eq!(str, seqstr[i..i+16]);
                },
                None => {
                    println!(" i =  {}", i);
                    panic!("iterator do not see end");
                },
            } // end match
        }
        // check iterator sees the end
        match kmergen32.next() {
            Some(kmer) => {
                let v = kmer.get_uncompressed_kmer();
                println!("kmer =  {:?}  " , v);
                panic!("iterator do not see end");
            },
            None => (),
        } // end match
    }  // end of test_gen_kmer


    
    // test generation of kmer of size 11
    #[test]
    fn test_gen_kmer32bit_50bases() {
        // got a string of 50 bases , genrate 11 Kmers, we shoud have 40 Kmers OK
        let seqstr = String::from("TCAAAGGGAAACATTCAAAATCAGTATGCGCCCGTTCAGTTACGTATTGC");
        // transform to a Vec<u8> and then a Sequence
        let slu8 = seqstr.as_bytes();
        // get a sequence with 2 bits compression
        let seq = Sequence::new(&slu8,2);
        let kmer_size : usize = 11;
        let mut kmergen_11b = KmerSeqIterator::<Kmer32bit>::new(kmer_size as u8, &seq);
        //
        for i in 0..seqstr.len()-kmer_size+1 {
            let kmer = kmergen_11b.next().unwrap();
            println!("kmer i = {}  {:b}  " , i , kmer.0);
            let v = kmer.get_uncompressed_kmer();
            println!("kmer i = {}  {:?}  " , i , v);
            // recall from_utf8_lossy return a Cow('a,str)
            let str = String::from_utf8_lossy(v.as_slice());
            println!("kmer i = {}  {}  " , i , str);
            assert_eq!(str, seqstr[i..i+kmer_size]);
        }
        // check iterator sees the end
        match kmergen_11b.next() {
            Some(kmer) => {
                let v = kmer.get_uncompressed_kmer();
                println!("kmer =  {:?}  " , v);
                panic!("iterator do not see end");
            },
            None => (),
        } // end match
    }  // end of test_gen_kmer


/* test_generate_weighted_kmer32bit should output the following with --nocapture
  but not in this order as we switched from IndexMap to FnvHashMap in indexing Kmers!!!!

    kmer [84, 67, 65]  TCA ,   weight 4
    kmer [67, 65, 65]  CAA ,   weight 2
    kmer [65, 65, 65]  AAA ,   weight 4
    kmer [65, 65, 71]  AAG ,   weight 1
    kmer [65, 71, 71]  AGG ,   weight 1
    kmer [71, 71, 71]  GGG ,   weight 1
    kmer [71, 71, 65]  GGA ,   weight 1
    kmer [71, 65, 65]  GAA ,   weight 1
    kmer [65, 65, 67]  AAC ,   weight 1
    kmer [65, 67, 65]  ACA ,   weight 1
    kmer [67, 65, 84]  CAT ,   weight 1
    kmer [65, 84, 84]  ATT ,   weight 2
    kmer [84, 84, 67]  TTC ,   weight 2
    kmer [65, 65, 84]  AAT ,   weight 1
    kmer [65, 84, 67]  ATC ,   weight 1
    kmer [67, 65, 71]  CAG ,   weight 2
    kmer [65, 71, 84]  AGT ,   weight 2
    kmer [71, 84, 65]  GTA ,   weight 2
    kmer [84, 65, 84]  TAT ,   weight 2
    kmer [65, 84, 71]  ATG ,   weight 1
    kmer [84, 71, 67]  TGC ,   weight 1
    kmer [71, 67, 71]  GCG ,   weight 1
    kmer [67, 71, 67]  CGC ,   weight 1
    kmer [71, 67, 67]  GCC ,   weight 1
    kmer [67, 67, 67]  CCC ,   weight 1
    kmer [67, 67, 71]  CCG ,   weight 1
    kmer [67, 71, 84]  CGT ,   weight 2
    kmer [71, 84, 84]  GTT ,   weight 2
    kmer [84, 84, 65]  TTA ,   weight 1
    kmer [84, 65, 67]  TAC ,   weight 1
    kmer [65, 67, 71]  ACG ,   weight 1
*/

    #[test]
    fn test_generate_weighted_kmer32bit() {
        //
        log_init();
        println!("test_generate_weighted_kmer32bit");
        // got a string of 48 bases shoud have 66 Kmers OK
        let seqstr = String::from("TCAAAGGGAAACATTCAAAATCAGTATGCGCCCGTTCAGTTACGTATT");
        // transform to a Vec<u8> and then a Sequence
        let slu8 = seqstr.as_bytes();
        // get a sequence with 2 bits compression
        let seq = Sequence::new(&slu8,2);
        let weighted_kmer_h : FnvHashMap<Kmer32bit,u32> = KmerGenerator::new(3).generate_weighted_kmer(&seq);
        let weighted_kmer = hashmap_count_to_vec_count(&weighted_kmer_h);
        for x in weighted_kmer {
            let ukmer = (x.0).get_uncompressed_kmer();
            println!("kmer {:?}  {} ,   weight {}", ukmer , String::from_utf8_lossy(ukmer.as_slice()), x.1);
        }
    } // end of test_generate_weighted_kmer32bit
    

    #[test]
    fn test_generate_weighted_kmer64bit() {
        //
        log_init();
        // got a string with repetition of the beginning to have some long kmers with multiplicity >= 2
        let seqstr = String::from("TCAAAGGGAAACATTCAAAATCAGTATGCGCCCGTTCAGTTACGTATTTCAAAGGGAAACATTCAAAATCAG");
        // transform to a Vec<u8> and then a Sequence
        let slu8 = seqstr.as_bytes();
        // get a sequence with 2 bits compression
        let seq = Sequence::new(&slu8,2);
        let weighted_kmer_h : FnvHashMap<Kmer64bit, u32> = KmerGenerator::new(15).generate_weighted_kmer(&seq);
        let weighted_kmer = hashmap_count_to_vec_count(&weighted_kmer_h);
        // first 9 kmers have weight 2 else weight 1
        for x in weighted_kmer {
            let ukmer = (x.0).get_uncompressed_kmer();
            let searched = String::from_utf8_lossy(ukmer.as_slice()).into_owned();
            println!("kmer {:?} position : {},   weight {}", ukmer , searched, x.1);
            if x.1 == 2 {
                // check we have it
                let pos1 = seqstr.find(&searched).unwrap();
                // check we have it once more
                let pos2 = &seqstr[pos1+1..].find(&searched).unwrap();
                // and no more
                if let Some(pos3) = seqstr[pos1+1+pos2+1..].find(&searched) {
                    log::error!("found kmer : {} more than 2 times, last found at {}", searched, pos1+1+pos2+1+pos3);
                }       
                assert!(&seqstr[pos1+1+pos2+1..].find(&searched).is_none());
            }
            if x.1 == 1 {
                // check we have it
                let pos = &seqstr.find(&searched).unwrap();
                // check we have it once 
                assert!(&seqstr[pos+1..].find(&searched).is_none());
            }
        }
    } // end of test_generate_weighted_kmer64bit



    #[test]
    fn test_generate_kmer16b32bit_pattern() {
        // got a string of 80 bases shoud have 66 Kmers OK
        let seqstr = String::from("TCAAAGGGAAACATTCAAAATCAGTATGCGCCCGTTCAGTTACGTATT");
        // transform to a Vec<u8> and then a Sequence
        let slu8 = seqstr.as_bytes();
        // get a sequence with 2 bits compression
        let seq = Sequence::new(&slu8,2);
        //
        let vkmer : Vec<Kmer16b32bit> = KmerGenerator::new(16).generate_kmer(&seq);
        assert_eq!(vkmer.len(), seqstr.len() - 16 +1);
        // check first kmer
        let v = vkmer[0].get_uncompressed_kmer();
        let str = String::from_utf8_lossy(v.as_slice());
        assert_eq!(str, seqstr[0..16]);
        // check last kmer
        let lastpos = seqstr.len() - 16;
        let v = vkmer[vkmer.len()-1].get_uncompressed_kmer();
        let str = String::from_utf8_lossy(v.as_slice());
        assert_eq!(str, seqstr[lastpos..seqstr.len()]);
    }

        
    #[test]
    fn test_generate_kmer64bit_pattern() {
        // got a string of 48 bases shoud have 21 Kmers OK
        let seqstr = String::from("TCAAAGGGAAACATTCAAAATCAGTATGCGCCCGTTCAGTTACGTATT");
        // transform to a Vec<u8> and then a Sequence
        let slu8 = seqstr.as_bytes();
        // get a sequence with 2 bits compression
        let seq = Sequence::new(&slu8,2);
        //
        let kmer_size : usize = 21;
        let vkmer : Vec<Kmer64bit> = KmerGenerator::new(kmer_size as u8).generate_kmer(&seq);
        assert_eq!(vkmer.len(), seqstr.len() - kmer_size +1);
        // check first kmer
        let v = vkmer[0].get_uncompressed_kmer();
        let str = String::from_utf8_lossy(v.as_slice());
        assert_eq!(str, seqstr[0..kmer_size]);
        // check last kmer
        let lastpos = seqstr.len() - kmer_size;
        let v = vkmer[vkmer.len()-1].get_uncompressed_kmer();
        let str = String::from_utf8_lossy(v.as_slice());
        assert_eq!(str, seqstr[lastpos..seqstr.len()]);
    }

    // test generation of kmer of size 11
    #[test]
    fn test_gen_kmer64bit_50bases() {
        // got a string of 50 bases , genrate 21 Kmers, we shoud have 48 Kmers OK
        let seqstr = String::from("TCAAAGGGAAACATTCAAAATCAGTATGCGCCCGTTCAGTTACGTATTGC");
        // transform to a Vec<u8> and then a Sequence
        let slu8 = seqstr.as_bytes();
        // get a sequence with 2 bits compression
        let seq = Sequence::new(&slu8,2);
        let kmer_size : usize = 21;
        let mut kmergen_21b = KmerSeqIterator::<Kmer64bit>::new(kmer_size as u8, &seq);
        //
        for i in 0..seqstr.len()-kmer_size+1 {
            let kmer = kmergen_21b.next().unwrap();
            println!("kmer i = {}  {:b}  " , i , kmer.0);
            let v = kmer.get_uncompressed_kmer();
            println!("kmer i = {}  {:?}  " , i , v);
            // recall from_utf8_lossy return a Cow('a,str)
            let str = String::from_utf8_lossy(v.as_slice());
            println!("kmer i = {}  {}  " , i , str);
            assert_eq!(str, seqstr[i..i+kmer_size]);
        }
        // check iterator sees the end
        match kmergen_21b.next() {
            Some(kmer) => {
                let v = kmer.get_uncompressed_kmer();
                println!("kmer =  {:?}  " , v);
                panic!("iterator do not see end");
            },
            None => (),
        } // end match
    }  // end of test_gen_kmer64bit_50bases

    
}  // end of mod tests
