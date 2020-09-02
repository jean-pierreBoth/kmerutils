//! This file contains structure for Kmer generation from a compressed sequence
//! Essentially we provide for iterators which extract kmer from a range in a sequence
//



pub use crate::kmer::*;

// This structure must take care of encoding conversions
// as all Kmer must be 2 bits encoded
// The generator takes as argument an iterator over a sequence.
// The structure KmerSeqIterator is dependant upon implementation of Kmer representation.
// The structure KmerGenerationPattern is there to provide absraction over it.
/// a Kmer iterator over a sequence.
///
/// usage:
/// 
///    let mut kmeriter = KmerSeqIterator::<Kmer16b32bit>::new(16, &seq);
///            match kmeriter.next(){
///                Some(kmer) => ...,
///                None => break,
///            }
/// 
///
pub struct KmerSeqIterator<'a , T > where T:CompressedKmerT {
    /// size of kmer
    nb_base: u8,
    /// an iterator for base calling
    seqiter: IterSequence<'a>,
    /// last kmer returned. At the beginning its None
    previous: Option<T>,
}



impl<'a, T> KmerSeqIterator<'a, T>  where T:CompressedKmerT {
    /// Constructor with default mode: no decoding of base!! for faster kmer generation, and default length
    pub fn new(ksize: u8, sequence: &'a Sequence) -> KmerSeqIterator<'a, T> {
        if ksize as usize > T::get_nb_base_max() {
            panic!("\n KmerSeqIterator cannot support so many bases for given kmer type, kmer size  {}", ksize);
        }
        let seqiter_arg = IterSequence::new(sequence, false);
        KmerSeqIterator{nb_base: ksize, seqiter:seqiter_arg, previous:None}
    } // end of new
    /// Set the range from which all kmer of a given size are to be extracted.
    pub fn set_range(&mut self, begin: usize, end: usize) -> std::result::Result<(),()> {
        self.seqiter.set_range(begin, end)
    }
} // end of impl for KmerSeqIterator




impl<'a>  KmerSeqIterator<'a, Kmer16b32bit> {
    /// return next kmer if any.
    pub fn next(&mut self) -> Option<Kmer16b32bit> {
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
            // we need 15 more!!
            let mut new_kmer_val:u32 = (next_base as u32) << 30;
            for i in 0..15 {
                if let Some(next_base) = self.seqiter.next()  {
                    new_kmer_val = new_kmer_val | ((next_base as u32) << (28-2*i));
                }
                else {
                    return None;
                }
            } // end of for
            let new_kmer = Kmer16b32bit(new_kmer_val);
            self.previous = Some(new_kmer);
            return Some(new_kmer);   
        } // end else        
    }  // end of next
} // end of impl<'a,T> KmerSeqIterator<'a, Kmer16b32bit>




impl<'a>  KmerSeqIterator<'a, Kmer32bit> {
    pub fn next(&mut self) -> Option<Kmer32bit> {
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
//            info!("in Kmer32bit next nb_base, base = {} {:b} ", kmer_size, next_base);
            let pos = 2*(kmer_size -1);
            let mut new_kmer:u32 = (kmer_size as u32) << 28;
            new_kmer = new_kmer | (next_base as u32) << pos;
            for i in 0..(kmer_size-1) {
                if let Some(next_base) = self.seqiter.next()  {
                    new_kmer = new_kmer | ((next_base as u32) << (pos-2-2*i));
                }
                else {
                    return None;
                }
            } // end of for
//            info!("in Kmer32bit next , kmer  = {:b} ", new_kmer);
            self.previous = Some(Kmer32bit(new_kmer));
            return Some(Kmer32bit(new_kmer));   
        } // end else        
    }  // end of next
} // end of impl<'a,T> KmerSeqIterator<'a, Kmer16b32bit>






impl<'a>  KmerSeqIterator<'a, Kmer64bit> {
    pub fn next(&mut self) -> Option<Kmer64bit> {
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
//            info!("in Kmer64bit next nb_base, base = {} {:b} ", kmer_size, next_base);
            let pos = 2*(kmer_size -1);
            let mut new_kmer = 0u64;
            new_kmer = new_kmer| (next_base as u64) << pos;
            for i in 0..(kmer_size-1) {
                if let Some(next_base) = self.seqiter.next()  {
                    new_kmer = new_kmer | ((next_base as u64) << (pos-2-2*i));
                }
                else {
                    return None;
                }
            } // end of for
//            info!("in Kmer64bit next , kmer  = {:b} ", new_kmer);
            self.previous = Some(Kmer64bit(new_kmer, self.nb_base));
            return Some(Kmer64bit(new_kmer, self.nb_base));   
        } // end else        
    }  // end of next
} // end of impl<'a,T> KmerSeqIterator<'a, Kmer64bit>




//=================== trait for kmer generation pattern ========================//



/// overload Kmer generation as in Rust Week 225
/// We define a generic trait defining a pattern for target type of kmer generation

pub trait KmerGenerationPattern<T:KmerT> {
    /// generate all kmers included in 0..
    fn generate_kmer_pattern(&self, seq : & Sequence) -> Vec<T>;
    /// generate all kmers inclused in begin..end with end excluded as in rust conventions.
    fn generate_kmer_pattern_in_range(&self, seq : & Sequence, begin:usize, end:usize) -> Vec<T>;   
}


// a  structure just providing reference to type of generation pattern to use Kmer to generate
// The structure could be empty. It is the method generate_kmer that really counts.
// Via type inference we go from T to KmerGenerationPattern<'a, T> and so to the generation of correct T kmer
// Finally we put back T in KmerGenerator ??

/// This structure store the size of kmer and provides overload of kmer generation for different types of Kmer.
///
///  usage:
///
///  ``` let vkmer : Vec<Kmer16b32bit> = KmerGenerator::new::<Kmer16b32bit>(16).generate_kmer(&seq);```
///


use std::marker::PhantomData;

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
    pub fn get_kmer_size(&self) -> usize { self.kmer_size as usize}
}  // end of impl KmerGenerator


// ======================= implementation for Kmer16b32bit ====================


/// implementation of kmer generation pattern for Kmer16b32bit
impl KmerGenerationPattern<Kmer16b32bit> for KmerGenerator<Kmer16b32bit> {
    fn generate_kmer_pattern(&self, seq : &Sequence) -> Vec<Kmer16b32bit> {
        if self.kmer_size != 16u8 {
            panic!("Kmer16b32bit has 16 bases!!");
        }
        // For a sequence of size the number of kmer is seq.size - kmer.size + 1  !!!
        // But it happens that "long reads" are really short 
        let nb_kmer = if seq.size() >= 16 { seq.size()-16+1} else {0};
        let mut kmer_vect = Vec::<Kmer16b32bit>::with_capacity(nb_kmer);
        let mut kmeriter = KmerSeqIterator::<Kmer16b32bit>::new(self.kmer_size, seq);
        loop {
            match kmeriter.next(){
                Some(kmer) => kmer_vect.push(kmer),
                None => break,
            }
        }
        //
        return kmer_vect;
    }  // end of generate_kmer_pattern


    fn generate_kmer_pattern_in_range(&self, seq : &Sequence, begin:usize, end:usize) -> Vec<Kmer16b32bit> {
        if self.kmer_size != 16u8 {
            panic!("Kmer16b32bit has 16 bases!!");
        }
        if begin >= end {
            panic!("KmerGenerationPattern<'a, Kmer16b32bit>  bad range for kmer iteration");
        }
        // For a sequence of size the number of kmer is seq.size - kmer.size + 1  !!!
        // But it happens that "long reads" are really short 
        let nb_kmer = if seq.size() >= 16 { seq.size()-16+1} else {0};
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
        // For a sequence of size the number of kmer is seq.size - kmer.size + 1  !!!
        // But it happens that "long reads" are really short
        let kmer_size = self.kmer_size as usize;
        let nb_kmer = if seq.size() >= kmer_size { seq.size()- kmer_size+1} else {0};
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

    
    fn generate_kmer_pattern_in_range(&self, seq : &Sequence, begin:usize, end:usize) -> Vec<Kmer32bit> {
        if self.kmer_size > 14u8 {
            panic!("Kmer16b32bit has less than 14 bases!!");
        }
        //
        if end <= begin {
            panic!("KmerGenerationPattern<'a, Kmer32bit>:generate_kmer_pattern_in_range incoherent range");
        }
        // For a sequence of size the number of kmer is seq.size - kmer.size + 1  !!!
        // But it happens that "long reads" are really short
        let kmer_size = self.kmer_size as usize;
        let nb_kmer = if seq.size() >= kmer_size { seq.size()- kmer_size+1} else {0};
        let mut kmer_vect = Vec::<Kmer32bit>::with_capacity(nb_kmer);
        let mut kmeriter = KmerSeqIterator::<Kmer32bit>::new(self.kmer_size, seq);
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

    
}  // end of impl<'a> KmerGenerationPattern<'a, Kmer32bit>

// ======================= implementation for Kmer64bit ====================


impl KmerGenerationPattern<Kmer64bit> for KmerGenerator<Kmer64bit> {
    
    fn generate_kmer_pattern(&self, seq : & Sequence) -> Vec<Kmer64bit> {
        if self.kmer_size > 32u8 {
            panic!("Kmer64bit has less than 32 bases!!");
        }
        // For a sequence of size the number of kmer is seq.size - kmer.size + 1  !!!
        // But it happens that "long reads" are really short
        let kmer_size = self.kmer_size as usize;
        let nb_kmer = if seq.size() >= kmer_size { seq.size()- kmer_size+1} else {0};
        let mut kmer_vect = Vec::<Kmer64bit>::with_capacity(nb_kmer);
        let mut kmeriter = KmerSeqIterator::<Kmer64bit>::new(self.kmer_size, seq);
        loop {
            match kmeriter.next(){
                Some(kmer) => kmer_vect.push(kmer),
                None => break,
            }
        }
        //
        return kmer_vect;
    }  // end of generate_kmer_pattern

    
    fn generate_kmer_pattern_in_range(&self, seq : &Sequence, begin:usize, end:usize) -> Vec<Kmer64bit> {
        if self.kmer_size > 32u8 {
            panic!("Kmer64bit has less than 32 bases!!");
        }
        if end <= begin {
            panic!("KmerGenerationPattern<'a, Kmer64bit>:generate_kmer_pattern_in_range incoherent range");
        }
        // For a sequence of size the number of kmer is seq.size - kmer.size + 1  !!!
        // But it happens that "long reads" are really short
        let kmer_size = self.kmer_size as usize;
        let nb_kmer = if seq.size() >= kmer_size { seq.size()- kmer_size+1} else {0};
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
    println!(" in generating kmer ...");
    let start_t = time::Instant::now();

    let mut nb_kmer : u64 = 0;
    
    for seq in seqvec {
        let vkmer : Vec<Kmer16b32bit> = KmerGenerator::new(16 as u8).generate_kmer(&seq);
        nb_kmer += vkmer.len() as u64;
    }
    
    let elapsed_t = start_t.elapsed().whole_seconds();
    println!(" elapsed time (s) in kmer generation {} ", elapsed_t);
    println!(" nb kmer generated {}", nb_kmer);
}


pub fn generate_all_kmer32bit(kmer_size:u8, seqvec : &Vec<Sequence>) {
    println!(" in generating kmer ...");
    let start_t = time::Instant::now();

    let mut nb_kmer : u64 = 0;
    
    for seq in seqvec {
        let vkmer : Vec<Kmer32bit> = KmerGenerator::new(kmer_size).generate_kmer(&seq);
        nb_kmer += vkmer.len() as u64;
    }
    
    let elapsed_t = start_t.elapsed().whole_seconds();
    println!(" elapsed time (s) in kmer generation {} ", elapsed_t);
    println!(" nb kmer generated {}", nb_kmer);
}


pub fn generate_all_kmer64bit(kmer_size:u8, seqvec : &Vec<Sequence>) {
    println!(" in generating kmer ...");
    let start_t = time::Instant::now();

    let mut nb_kmer : u64 = 0;
    
    for seq in seqvec {
        let vkmer : Vec<Kmer64bit> = KmerGenerator::new(kmer_size).generate_kmer(&seq);
        nb_kmer += vkmer.len() as u64;
    }
    
    let elapsed_t = start_t.elapsed().whole_seconds();
    println!(" elapsed time (s) in kmer generation {} ", elapsed_t);
    println!(" nb kmer generated {}", nb_kmer);
}

//=======================================================================================

#[cfg(test)]
mod tests {
    use super::*;
    //
    //   11010000001010100000000100111101
    //    T C A A A G G G A A A C A T T C
    //   In this test all bytes are complete
    #[test]
    fn test_gen_kmer16b32bit_80bases() {
        // got a string of 80 bases shoud have 66 Kmers OK
        let seqstr = String::from("TCAAAGGGAAACATTCAAAATCAGTATGCGCCCGTTCAGTTACGTATTGCTCTCGCTAATGAGATGGGCTGGGTACAGAG");
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
        // got a string of 80 bases shoud have 66 Kmers OK
        let seqstr = String::from("TCAAAGGGAAACATTCAAAATCAGTATGCGCCCGTTCAGTTACGTATTGC");
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
        // got a string of 60 bases , genrate 11 Kmers, we shoud have 49 Kmers OK
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
