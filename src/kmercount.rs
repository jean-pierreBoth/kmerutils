//! This module contains trait/struct to count kmers
// We take a slice of kmer and count them


use log::info;


use ::std::mem;
use ::std::io;
use ::std::io::{Write,Read};
use ::std::fs;
use ::std::fs::OpenOptions;
use ::std::marker::PhantomData;

use crossbeam_utils::thread::*;

use time::*;
use ::cuckoofilter::*;
use ::bloom::*;

use std::collections::hash_map::DefaultHasher;
use metrohash::{MetroHash64};
use fnv::*;

pub use crate::kmergenerator::*;

// a magic for file dump of unique kmers
// this type of dump is useful if we just want to know if we can garbage a kmer.
// file contains position at present time as we expect a link with qualities
pub const COUNTER_UNIQUE:u32 = 0xcea2bbdd;



// a magic for multiple kmers.
// file will have counts and no position
pub const COUNTER_MULTIPLE:u32 = 0xcea2bbff;


// to be made generic

/// Trait kmer counter. absraction of basic  request over simple case
/// and theaded case where we maintain a pool of counters

pub trait KmerCountT {
    /// the type of kmer we count
    type Kmer;
    /// return estimated count (0 if not seen)
    fn get_count(&self, kmer: Self::Kmer) -> u32;
    /// insert a kmer. to be templated
    fn insert_kmer(&mut self, kmer: Self::Kmer);
    /// returns number of distinck kmers
    fn get_nb_distinct(&self) -> u64;
    /// returns the number of unique kmers
    fn get_nb_unique(&self) -> u64;
}



// for multi threading


// our modules

use crate::invhash::*;

/// The structure to count kmers. 
pub struct KmerCounter<Kmer>  {
    /// number of bits for element in bloom filter
    bloom_f_nb_bits: u8,
    /// The false positive rate required for. standard is 0.03
    _fpr: f32,
    /// a cuckoo filter to keep track of kmer encountered only once
    cuckoo_f: CuckooFilter<MetroHash64>,
    /// a counting bloom filter to keep track of kmer encountered at least twice
    cbloom_f: CountingBloomFilter,
    /// total number of disctint elements seen
    nb_distinct: u64,
    ///
    _kmertype: PhantomData<Kmer>,
}


impl <Kmer> KmerCounter<Kmer> where Kmer: CompressedKmerT
{
    pub fn new(fpr_arg: f32, capacity: usize, nb_bits: usize) -> KmerCounter<Kmer> {
        //
        KmerCounter {  bloom_f_nb_bits: nb_bits as u8 , _fpr: fpr_arg, 
                       cuckoo_f: CuckooFilter::with_capacity(capacity as usize),
                       cbloom_f: CountingBloomFilter::with_rate(nb_bits, fpr_arg, capacity as u32),
                       nb_distinct:0,
                       _kmertype: PhantomData,
        }
    } // end of new


    /// get count for a kmer. returns only something if kmer has been seen at least twice.
    pub fn get_above2_count(&self, kmer: Kmer) -> u32 {
        if self.cbloom_f.contains(& kmer.get_compressed_value()) {
            self.cbloom_f.estimate_count(& kmer.get_compressed_value())
        }
        else { 0 }
    }  // get_above2_count



    /// a method to free memory occupied by the counter of unique kmers
    pub fn eliminate_once_kmer(&mut self) {
        let occupied_size = self.cuckoo_f.memory_usage();
        let nbitems = self.cuckoo_f.len();
        println!(" once kmer counter : size = {} , nbitems = {}",  occupied_size, nbitems);
        // purge
        self.cuckoo_f = CuckooFilter::with_capacity(0 as usize);
        println!(" after purging size = {}", self.cuckoo_f.memory_usage());       
    }

    /// returns number of bits used for a count
    pub fn get_count_nb_bits(&self) -> u8 {
        self.bloom_f_nb_bits
    }
} // end impl KmerCounter



/// This function dumps to a file the content of a kmer counter.
/// We only dump kmers seen more than once: given a kmer if it is not multiple, then it is unique
/// It appears that generating kmer is relatively fast compared to counting.
/// So we take the counter, regenerate kmers and fetch them in counter before dumping.
///
/// Args are
///   - the counter to dump
///   - the file name to dump in
///   - the sequence list from which we count kmer
///   - the generator to use for to regenerate kmer and search them in counter. (Of course it must generates the same kmer as those
///     counted in counter.
///
/// format of file is :
///    - a u32 magic : COUNTER_MULTIPLE (0xcea2bbff)
///    - kmer size as u8
///    - number of bytes for count 1 or 2.
///    - number of kmers as a u64. Note this an approximate value (due to false positive for multiple kmer count)
///    - then list of (kmer , count as u8 or u16 depending on number of bytes for count)
//
pub fn dump_in_file_multiple_kmer<Kmer>(counter: &KmerCounter<Kmer>, fname: &String,
                                          seqvec : &Vec<Sequence>,
                                          kmer_generator : &KmerGenerator<Kmer>) -> std::result::Result<u64, io::Error>
     where Kmer: CompressedKmerT ,
           KmerGenerator<Kmer>: KmerGenerationPattern<Kmer>
{
    //
    println!(" dumping kmer in file : {} ", fname);
    //
    let magic: u32  = COUNTER_MULTIPLE;
    let nb_bytes_by_count;
    if counter.bloom_f_nb_bits  <= 8 {
        nb_bytes_by_count = 1;
    }
    else if counter.bloom_f_nb_bits <= 16 {
        nb_bytes_by_count = 2;
    }
    else {
        panic!("kmercount::dump_in_file_counted_kmer , can only dump kmer count using less than 16 bits for count"); 
    }
    // as kmer generation is fast we generate once more all kmers and check for those that are in once_f
    let file = OpenOptions::new().write(true).create(true).truncate(true).open(&fname)?;
    let mut bufw : io::BufWriter<fs::File> = io::BufWriter::with_capacity(1_000_000_000, file);
    let mut nb_kmer_dumped : u64 = 0;
    let kmer_size = kmer_generator.get_kmer_size();
    // write magic, kmer_size and number of kmers and size of in bytes used for counting in bloomfilter.
    bufw.write(unsafe { &mem::transmute::<u32, [u8;4]>(magic) } )?;
    bufw.write(unsafe { &mem::transmute::<u8, [u8;1]>(kmer_size as u8) })?;
    bufw.write(unsafe { &mem::transmute::<u8, [u8;1]>(nb_bytes_by_count as u8) })?;
    // total number of kmer 
    let nb_kmer_to_dump : u64 = counter.get_nb_distinct() - counter.get_nb_unique();
    println!("dumping nb kmer {}", nb_kmer_to_dump);
    bufw.write(unsafe { &mem::transmute::<u64, [u8;8]>(nb_kmer_to_dump) } )?;
    //
    // now we must generate kmer once more. But we must keep track of already dumped kmer.
    //
    let mut cuckoo_dumped = CuckooFilter::<MetroHash64>::with_capacity((2 * nb_kmer_to_dump) as usize);
    for seq in seqvec {
        let vkmer : Vec<Kmer> = kmer_generator.generate_kmer(&seq);
        for kmer in &vkmer {
            // get canonical kmer
            let kmin = kmer.reverse_complement().min(*kmer);
            let key = kmin.get_compressed_value();
            // the following requires a read access on filter
            let already = cuckoo_dumped.contains(& key);
            if !already {
                // then we can dump, and it will be dumped once.
                let count = counter.get_above2_count(kmin) as u8;
                // now we enter a section that requires a write access on filter. Much less frequent especially
                // as we go further in kmer sequence
                if count >= 2 {
                    // needs a write access on cuckoo
                    cuckoo_dumped.add(&key).unwrap();
                    kmin.dump(&mut bufw)?;  // dump kmer
                    match nb_bytes_by_count {
                        1 => {
                            bufw.write(unsafe { &mem::transmute::<u8, [u8;1]>(count as u8) })?;
                        },
                        2 => {
                            bufw.write(unsafe { &mem::transmute::<u16, [u8;2]>(count as u16) })?;
                        },
                        _ => panic!("kmercount::dump_in_file_counted_kmer incoherent number of bytes for count"),
                    };                
                    nb_kmer_dumped += 1;
                    if nb_kmer_dumped % 100_000_000 == 0 {
                        println!("dump_kmer_counter, number of kmer dumped : {} ", nb_kmer_dumped);
                    }
                } // end count >= 2
            } // end if new         
        } // end of for on kmer
    } // end of for on  seq
    //
    println!("dump_kmer_counter, number of kmer to dump, number of kmer dumped : {}  {} ", nb_kmer_to_dump, nb_kmer_dumped);
    //
    Ok(nb_kmer_dumped)
}  // end of dump_in_file


// =========================================================================
//           implementation of trait KmerCountT
// =========================================================================

impl <Kmer> KmerCountT for KmerCounter<Kmer>
    where Kmer:CompressedKmerT
{
    type Kmer=Kmer;
    /// insert a Kmer
    fn insert_kmer(&mut self, kmer: Kmer) {
        //
        // we should avoid hashing twice (once for each counter) the same kmer...
        // could implement a hasher that does nothing for a type representing a already hashed 32bit kmer
        // we store a kmer in cuckoo filter when we see it for the first time.
        // when seen for the second time we delete it from the cucckoo and insert it twice in cbloom_f
        //
        if self.cbloom_f.contains(& kmer.get_compressed_value()) {
            self.cbloom_f.insert(& kmer.get_compressed_value());    // insert once more
        }
        else {
            // if test_and_add succeded, i.e returns an Ok(true)!!!)  we have one more new kmer
            let inserted = match self.cuckoo_f.test_and_add(& kmer.get_compressed_value()) {
                Ok(true) => true,
                _ => false,
            };
            if inserted {
                self.nb_distinct += 1;
            }
            else {
                // remove from cuckoo and transfer to cbloom_f
                self.cuckoo_f.delete(& kmer.get_compressed_value());
                self.cbloom_f.insert(& kmer.get_compressed_value());                
                self.cbloom_f.insert(& kmer.get_compressed_value());                
            }
        }  // end not yet inserted in cbloom       
    } // end of insert_kmer

     /// get count for a kmer. returns also kmer with .
    fn get_count(&self, kmer: Kmer) -> u32 {
        if self.cbloom_f.contains(& kmer.get_compressed_value()) {
            return self.cbloom_f.estimate_count(& kmer.get_compressed_value());
        }
        else {
            self.cuckoo_f.contains(& kmer.get_compressed_value()) as u32
        }
    } // end of get_count

    /// returns number of different kmers.
    fn get_nb_distinct(&self) -> u64 {
        self.nb_distinct
    }

    /// returns the number of unique kmers
    fn get_nb_unique(&self) -> u64 {
        self.cuckoo_f.len() as u64
    }
} // end of impl KmerCountT for KmerCounter



// This function counts 16-mers from a vector of Sequence
fn count_kmer<Kmer>(seqvec : &Vec<Sequence>,  kmer_generator : &KmerGenerator<Kmer>) -> KmerCounter<Kmer>
     where Kmer: CompressedKmerT ,
           KmerGenerator<Kmer>: KmerGenerationPattern<Kmer>
{
    println!(" in counting kmer ...");
    let start_t = Instant::now();

    let mut nb_kmer : u64 = 0;
    // rule of thumb could be seqvec.len() * 1000 for variable length read
    let capacity = 1_200_000_000;
    let nb_bits = 8;
    let mut nbseq = 0u64;
    //
    let mut kmer_counter  = KmerCounter::new(0.03, capacity, nb_bits);
    println!("space occupied by cuckoo filter {}", kmer_counter.cuckoo_f.memory_usage());
    //
    for seq in seqvec {
        let vkmer : Vec<Kmer> = kmer_generator.generate_kmer(&seq);
        for kmer in &vkmer {
            let kmin = kmer.reverse_complement().min(*kmer);
            kmer_counter.insert_kmer(kmin);
        }
        nb_kmer += vkmer.len() as u64;
        nbseq += 1;
        if nbseq % 100_000 == 0 {
            println!(" nb  read treated = {} ", nbseq);
            println!(" nb once kmer = {} nb distinct = {} nb generated = {} ", kmer_counter.cuckoo_f.len(), kmer_counter.nb_distinct, nb_kmer);
        }        
    }
    let elapsed_t = start_t.elapsed().whole_seconds();
    println!(" elapsed time (s) in kmer32b counting {} ", elapsed_t);
    println!(" nb kmer generated {}", nb_kmer);
    //
    kmer_counter
}


/// counts 16-kmer
pub fn count_kmer16b32bit(seqvec : & Vec<Sequence>) -> KmerCounter<Kmer16b32bit> {
    let kmer_generator = KmerGenerator::new(16);
    let counter:KmerCounter<Kmer16b32bit> = count_kmer(seqvec, &kmer_generator);
    //
    counter
}

/// counts kmer of size <= 14
pub fn count_kmer32bit(seqvec : & Vec<Sequence>, kmer_size:usize) -> KmerCounter<Kmer32bit> {
    if kmer_size >= 15 {
        panic!("count_kmer32bit cannot count kmer of size greater than 14, size received : {}", kmer_size);
    }
    let kmer_generator = KmerGenerator::new(kmer_size as u8);
    let counter:KmerCounter<Kmer32bit> = count_kmer(seqvec, &kmer_generator);
    counter
}


/// counts kmer up to 32.
pub fn count_kmer64bit(seqvec : & Vec<Sequence>, kmer_size:usize) -> KmerCounter<Kmer64bit> {
    if kmer_size > 32 {
        panic!("count_kmer64bit cannot count kmer of size greater than 14, size received : {}", kmer_size);
    }
    else if kmer_size <= 16 {
        panic!("count_kmer64bit is too expensive, use count_kmer32bit for kmer of size  : {}", kmer_size);
        
    }
    let kmer_generator = KmerGenerator::new(kmer_size as u8);
    let counter:KmerCounter<Kmer64bit> = count_kmer(seqvec, &kmer_generator);
    counter
}






// ======================================================================
//      To Parallelize Kmer counting we use a pool of KmerCounter
// =====================================================================

/// To Parallelize Kmer counting we use a pool of KmerCounter and we need a function to dispatch each kmer
/// to its counter.
///
/// We get one counter by thread. kmer are dispatched via
/// a modulo number of threads after one more pass of hash
/// to ensure equidistribution of loads. Of course
/// The routine filling the data must use the same scheme
/// to dispatch kmer as the one used to answer request.





/// object counted by KmerCounterPool must be dispatchable
pub trait DispatchableT {
    type ToDispatch;
    ///
    fn dispatch(&self, nb_receiver:usize) -> usize;
}



impl DispatchableT for Kmer16b32bit {
    /// compressed value is a u32
    type ToDispatch=u32;
    ///
    #[inline(always)]
    fn dispatch(&self, nb_receiver:usize) -> usize {
        (int32_hash(self.get_compressed_value()) % nb_receiver as u32) as usize
    }
}


impl DispatchableT for Kmer32bit {
    /// compressed value is a u32
    type ToDispatch=u32;
    ///
    #[inline(always)]
    fn dispatch(&self, nb_receiver:usize) -> usize {
       (int32_hash(self.get_compressed_value()) % nb_receiver as u32) as usize
    }
}


impl DispatchableT for Kmer64bit {
    /// compressed value is a u32
    type ToDispatch=u64;
    ///
    #[inline(always)]
    fn dispatch(&self, nb_receiver:usize) -> usize {
        (int64_hash(self.get_compressed_value()) % nb_receiver as u64) as usize
    }
}


/// We maintain a pool of counters. One per thread.

pub struct KmerCounterPool<Kmer>  {
    // one counter per thread
    pub counters: Vec<Box<KmerCounter<Kmer> >>,
}


impl<Kmer> KmerCounterPool<Kmer> {
    pub fn new (counters:Vec<Box<KmerCounter<Kmer>>>) -> KmerCounterPool<Kmer> {
        // CAVEAT we should check that all counters have same number of bits in bloom filter
        KmerCounterPool{counters: counters}
    }

    /// get count for a kmer only if kmer has been seen at least twice.
    /// do not check for cuckoo filter as we want only count >= 2. So it is a little faster.
    // possibly put it trait KmerCountT
    pub fn get_above2_count(&self, kmer: Kmer) -> u32
    where Kmer: CompressedKmerT+DispatchableT
    {
        let loc = kmer.dispatch(self.counters.len());        
        self.counters[loc as usize].get_above2_count(kmer)
    }  // get_above2_count
    
    /// returns number of bits used for a count
    pub fn get_count_nb_bits(&self) -> u8 
        where Kmer: CompressedKmerT {
        if self.counters.len() > 0 {
            self.counters[0].as_ref().get_count_nb_bits()
        }
        else { 0}
    }

    /// This function dumps in fname the content of KmerCounterPool
    /// We dump in fact only kmer that are not unique.
    /// For unique Kmer it is their position that is interesting and can be correlated with Qualities and
    /// this is done by stucture KmerFilter1 for 16 bases Kmers
    ///
    /// format of file is :
    ///    - a u32 magic : COUNTER_MULTIPLE:u32 = 0xcea2bbff
    ///    - kmer size as u8
    ///    - nb_bytes_by_count as u8
    ///    - number of kmer dumped as u64. Note this an approximate value (due to false positive for multiple kmer count)
    ///    - then list (kmers,count)
    ///

    pub fn dump_kmer_counter (&self, fname: &String, seqvec : &Vec<Sequence>,
                             kmer_generator : &KmerGenerator<Kmer>) -> std::result::Result<usize, io::Error>
    where Kmer: CompressedKmerT+DispatchableT,
          KmerGenerator<Kmer>: KmerGenerationPattern<Kmer>,
    {
        //
        println!(" in KmerCounterPool::dump_kmer_counter, dump file : {} ", fname);
        let start_t = time::Instant::now();
        //
        let magic: u32  = COUNTER_MULTIPLE;
        // as kmer generation is fast we generate once more all kmers and check for those that are in once_f
        let file = OpenOptions::new().write(true).create(true).truncate(true).open(&fname)?;
        let mut bufw : io::BufWriter<fs::File> = io::BufWriter::with_capacity(1_000_000_000, file);
        let mut nb_kmer_dumped = 0;
        let kmer_size = kmer_generator.get_kmer_size();
        // as bloom f allocated have 8 bits per slot ... but we dump that to have the possibility to change.
        let nb_bytes_by_count = 1;
        // write magic, kmer_size and number of kmers
        bufw.write(unsafe { &mem::transmute::<u32, [u8;4]>(magic) } )?;
        bufw.write(unsafe { &mem::transmute::<u8, [u8;1]>(kmer_size as u8) })?;
        bufw.write(unsafe { &mem::transmute::<u8, [u8;1]>(nb_bytes_by_count as u8) })?;
        // how many kmer we have to dump?  
        let nb_kmer_to_dump = self.get_nb_distinct() - self.get_nb_unique();
        bufw.write(unsafe { &mem::transmute::<u64, [u8;8]>(nb_kmer_to_dump) } )?;
        // now we must generate kmer once more. But we must keep track of already dumped kmer.
        let mut cuckoo_dumped = CuckooFilter::<MetroHash64>::with_capacity((2.0 * nb_kmer_to_dump as f64) as usize);
        for seq in seqvec {
            let vkmer : Vec<Kmer> = kmer_generator.generate_kmer(&seq);
            for kmer in &vkmer {
                // get canonical kmer
                let kmin = kmer.reverse_complement().min(*kmer);
                let key = kmin.get_compressed_value();
                // the following requires a read access on filter
                let already = cuckoo_dumped.contains(& key);
                if !already {
                    let count = self.get_above2_count(kmin) as u8;
                    // now we enter a section that requires a write access on filter. Much less frequent especially
                    // as we go further in kmer sequence
                    if count >= 2 {
                        cuckoo_dumped.add(&key).unwrap();
                        kmin.dump(&mut bufw)?;
                        // dump count
                        bufw.write(unsafe { &mem::transmute::<u8, [u8;1]>(count) })?; 
                        nb_kmer_dumped += 1;
                        if nb_kmer_dumped % 100_000_000 == 0 {
                            println!("dump_kmer_counter, number of kmer dumped : {} ", nb_kmer_dumped);
                        }
                    } // end test on count
                } // end test on new
            } // end of for on kmer
        } // end of for on seq
        let elapsed_t = start_t.elapsed().whole_seconds();
        println!("dump_kmer_counter, number of kmer dumped : {} number to dump {} ", nb_kmer_dumped, nb_kmer_to_dump);
        println!(" elapsed time (s) {} ", elapsed_t);
        //
        Ok(nb_kmer_dumped as usize)
    }  // end of dump_in_file


}  // end of impl KmerCounterPool


// Now we can have a generic implementation of KmerCountT for KmerCounterPool
// once we can impose DispatchableT constraint on Kmer

impl <Kmer> KmerCountT for  KmerCounterPool<Kmer>
   where Kmer:CompressedKmerT+DispatchableT
{
    // in type equation, left Kmer is associated type of KmerCountT , right is template parameter
    type Kmer=Kmer;
    ///
    fn insert_kmer(&mut self, kmer: Kmer) {
        // dispatch
        let loc = kmer.dispatch(self.counters.len());        
        self.counters[loc as usize].insert_kmer(kmer);
    }  // end of insert_kmer
    
    fn get_count(&self, kmer: Kmer) -> u32 {
        let loc = kmer.dispatch(self.counters.len());        
        self.counters[loc as usize].get_count(kmer)
    }   

    fn get_nb_distinct(&self) -> u64 {
        self.counters.iter().map(|v| v.get_nb_distinct()).sum()
    }

    /// returns number of unique kmers (loops on internal counters)
    fn get_nb_unique(&self) -> u64 {
        self.counters.iter().map(|v| v.get_nb_unique()).sum()
    }
}



///////////////////////////////////////////////////////////////////////


    
/// This function dumps in fname the content of KmerCounterPool
/// We dump in fact only kmer that are not unique.
/// For unique Kmer it is their position that is interesting and can be correlated with Qualities and
/// this is done by stucture KmerFilter1 for 16 bases Kmers
///
/// format of file is :
///    - a u32 magic : COUNTER_MULTIPLE:u32 = 0xcea2bbff
///    - kmer size as u8
///    - number of kmers as a u64. Set to 0u64 if unknown else > 0 (to get fast allocation and check on file when reloading)
///    - then list of (kmers,numseq,numkmer)
///

pub fn dump_kmer_counter<Kmer>(counter_pool: &KmerCounterPool<Kmer>, fname: &String, seqvec : &Vec<Sequence>,
                               kmer_generator : &KmerGenerator<Kmer>) -> std::result::Result<usize, io::Error>
    where Kmer: CompressedKmerT+DispatchableT,
          KmerGenerator<Kmer>: KmerGenerationPattern<Kmer>,
{
    //
    println!(" dumping multiple kmer in file : {} ", fname);
    //
    let magic: u32  = COUNTER_MULTIPLE;
    // as kmer generation is fast we generate once more all kmers and check for those that are in once_f
    let file = OpenOptions::new().write(true).create(true).truncate(true).open(&fname)?;
    let mut bufw : io::BufWriter<fs::File> = io::BufWriter::with_capacity(1_000_000_000, file);
    let mut nb_kmer_dumped : usize = 0;
    let kmer_size = kmer_generator.get_kmer_size();
    // as bloom f allocated have 8 bits per slot ... but we dump that to have the possibility to change.
    let nb_bytes_by_count = 1;
    // write magic, kmer_size and number of kmers
    bufw.write(unsafe { &mem::transmute::<u32, [u8;4]>(magic) } )?;
    bufw.write(unsafe { &mem::transmute::<u8, [u8;1]>(kmer_size as u8) })?;
    bufw.write(unsafe { &mem::transmute::<u8, [u8;1]>(nb_bytes_by_count as u8) })?;
    // how many kmer we have to dump?  
    let nb_kmer_to_dump = counter_pool.get_nb_distinct() - counter_pool.get_nb_unique();
    println!(" threaded dump multiple kmer , nb kmer to dump : {} ", nb_kmer_to_dump);
    bufw.write(unsafe { &mem::transmute::<u64, [u8;8]>(nb_kmer_to_dump as u64) })?;
        // now we must generate kmer once more. But we must keep track of already dumped kmer.
    let mut cuckoo_dumped = CuckooFilter::<MetroHash64>::with_capacity((1.5 * nb_kmer_to_dump as f64) as usize);
    bufw.write(unsafe { &mem::transmute::<u64, [u8;8]>(nb_kmer_to_dump) } )?;
    // now we must generate kmer once more
    for seq in seqvec {
        let vkmer : Vec<Kmer> = kmer_generator.generate_kmer(&seq);
        for kmer in &vkmer {
            let kmin = kmer.reverse_complement().min(*kmer);
            let key = kmin.get_compressed_value();
            // the following requires a read access on filter
            let already = cuckoo_dumped.contains(& key);
            if !already {
                // we rely on KmerCounterPool being a partition. counter_pool.get_count does the dispatching job as Kmer:DispatchableT
                let count = counter_pool.get_count(kmin) as u8;
                if count >= 2 {
                    // needs a write access on cuckoo
                    cuckoo_dumped.add(&key).unwrap();
                    // needs a mutex on bufw and nb_kmer_dumped
                    kmin.dump(&mut bufw)?;
                    // dump count
                    bufw.write(unsafe { &mem::transmute::<u8, [u8;1]>(count as u8) })?; 
                    nb_kmer_dumped += 1;
                }
            } // end if !already
        } // end of for kmer
    } // end of for on seq
    println!("dump_kmer_counter, number of kmer to dump {} , dumped : {} ", nb_kmer_to_dump, nb_kmer_dumped);
    //
    Ok(nb_kmer_dumped as usize)
}  // end of dump_in_file





pub fn threaded_dump_kmer_counter<Kmer>(counter_pool: &KmerCounterPool<Kmer>, fname: &String, seqvec : &Vec<Sequence>,
                               nb_threads:usize, kmer_size:usize) -> std::result::Result<usize, io::Error>
    where Kmer: CompressedKmerT+DispatchableT+Send+Sync,
          Kmer::Val: Send+Sync,
          KmerGenerator<Kmer>: KmerGenerationPattern<Kmer>,
{
    //
    println!(" threaded dump multiple kmer in file : {} ", fname);
    let start_t = time::Instant::now();
    //
    let magic: u32  = COUNTER_MULTIPLE;
    // as kmer generation is fast we generate once more all kmers and check for those that are in once_f
    let file = OpenOptions::new().write(true).create(true).truncate(true).open(&fname)?;
    let mut bufw : io::BufWriter<fs::File> = io::BufWriter::with_capacity(1_000_000_000, file);
    // as bloom f allocated have 8 bits per slot ... but we dump so that to have the possibility to change.
    let nb_bytes_by_count;
    if counter_pool.counters[0].bloom_f_nb_bits <= 8 {
        nb_bytes_by_count = 1;
    }
    else if counter_pool.counters[0].bloom_f_nb_bits <= 16 {
        nb_bytes_by_count = 2;
    }
    else {    
        panic!("threaded_dump_kmer_counter bloom filter has more than 16 bits/item");
    }
    //
    // write magic, kmer_size and number of kmers
    bufw.write(unsafe { &mem::transmute::<u32, [u8;4]>(magic) } )?;
    bufw.write(unsafe { &mem::transmute::<u8, [u8;1]>(kmer_size as u8) })?;
    bufw.write(unsafe { &mem::transmute::<u8, [u8;1]>(nb_bytes_by_count as u8) })?;
    // how many kmer we have to dump?  
    let nb_kmer_to_dump = counter_pool.get_nb_distinct() - counter_pool.get_nb_unique();
    bufw.write(unsafe { &mem::transmute::<u64, [u8;8]>(nb_kmer_to_dump as u64) })?;
    println!(" threaded dump multiple kmer , nb kmer to dump : {} ", nb_kmer_to_dump);
    // we set count to u16 beccause of high coverage and high error rates in long reads that force us to reduce kmer size!!
    let (s, r) = crossbeam_channel::bounded::<(Kmer, u16)>(1_000_000);
    //
    crossbeam_utils::thread::scope(|scope| {
        //
        let mut join_handles = Vec::new();
        // channel creation, type of message is a tuple (compressed value of kmer, u8)
        // cuckoo filter creation
        // senders of kmer.
        for i in 0..nb_threads {
            let sender_i = s.clone();
            let sender_handle = scope.spawn(move |_|   {
                let nb_to_dump_i = (counter_pool.counters[i].get_nb_distinct() - counter_pool.counters[i].get_nb_unique()) as usize;
                let mut nb_sent : usize = 0;
                let filter_size: u32 = (1.25 * nb_kmer_to_dump as f32) as u32 / nb_threads as u32;
                let mut filter : CuckooFilter<MetroHash64> = CuckooFilter::with_capacity(filter_size as usize);
                let kmer_generator = KmerGenerator::new(kmer_size as u8);
                let mut numseq = 0;
                for seq in seqvec {
                    let vkmer : Vec<Kmer> = kmer_generator.generate_kmer(&seq);
                    for kmer in &vkmer {
                        let kmin = kmer.reverse_complement().min(*kmer);
                        if kmin.dispatch(nb_threads) == i {
                            let key = kmin.get_compressed_value();                        
                            // the following requires a read access on filter
                            let already;
                            already = filter.contains(& key);
                            if !already {
                                // we rely on KmerCounterPool being a partition. counter_pool.get_count does the dispatching job as Kmer:DispatchableT
                                let count = counter_pool.get_count(kmin) as u16;
                                if count >= 2 {
                                    // now we can dump
                                    filter.add(& key).unwrap();
                                    nb_sent += 1;
                                    let msg = (kmin, count);
                                    sender_i.send(msg);
                                } // end if counted >= 2
                            } // end if !already
                        } //  end if dispatch
                    } // end of for on kmer
                    numseq += 1;
                    // an exit test based on the fact we know how many each thread has to dump
                    if nb_sent >= nb_to_dump_i {
                        println!("thread i {} has done its job nb kmer sent {} , exiting at numseq {} ", i, nb_sent, numseq);
                        return Box::new(nb_sent);
                    }
                    //
                    if numseq % 100_000_000 == 0 {
                        println!("thread i {} processed seq {}", i, numseq);
                    }
                } // end of for on seq
                // return nb_sent
                Box::new(nb_sent)
            }); // end of scope.thread
            join_handles.push(sender_handle);
        } // end of for on sender threads
        // each sender_i is dropped by its thread, now we drop original sender, so receptor
        // will exit nicely
        drop(s);
        //
        // receptor thread
        //
        let receptor_handle = scope.spawn(move |_| {
            let mut nb_received : u64 = 0;
            let mut idump = 1;
            r.for_each(|msg| { let _res = msg.0.dump(&mut bufw);
                               bufw.write(unsafe { &mem::transmute::<u16, [u8;2]>(msg.1)} ).unwrap();
                               nb_received += 1;
                               if (nb_received * 10) >= idump * nb_kmer_to_dump {
                                   // an impression every 10%
                                   println!("receptor recieved {} msg, to dump {} ", nb_received, nb_kmer_to_dump);
                                   idump += 1;
                               }
            }); // end for_each
            Box::new(nb_received as usize)
        });// end of reciever thread.
        join_handles.push(receptor_handle);
        //
        // now we must merge
        //
        let mut nb_sent = 0;
        let mut nb_msg_vec : Vec<usize> = Vec::with_capacity(nb_threads+1);
        for handle in join_handles {
            nb_msg_vec.push(*(handle.join().unwrap()));
        }
        for _ in 0..nb_msg_vec[nb_msg_vec.len()-1] {
            nb_sent += 1;
        }
        println!("nb msg sent {}, nb msg received {} ", nb_sent, nb_msg_vec[nb_msg_vec.len()-1]);
    } // end of closure in scope
    ).unwrap();  // end of scope
    //
    let elapsed_t = start_t.elapsed().whole_seconds();
    println!(" elapsed time (s) {} ", elapsed_t);
    //
    Ok(1)
} // end of threaded_dump_kmer_counter




/// This function counts K-mers from a vector of Sequence

pub fn count_kmer_thread_independant<Kmer>(seqvec : &Vec<Sequence>, nb_threads:usize, kmer_size:usize) -> KmerCounterPool<Kmer>
where Kmer: CompressedKmerT+DispatchableT+Send,
      KmerGenerator<Kmer>: KmerGenerationPattern<Kmer>
{
    println!(" in counting kmer ... : count_kmer_thread_independant");
    //
    let capacity = 1_000_000_000;
    let nb_bits = 8;
    //
    // as kmer generation seems to be (very fast) and all time is spent in counting
    // each thread generates the whole set of kmers but we associate to each thread a filter to select
    // the kmers it takes. At the end of all threads we must gather results.
    // As bloom filters are union able its OK. All single kmer in cuckoo will stay unique as filters
    // realize a partition, but we will have to tranfer contents to one of cuckoo.
    //
        
    let poolcounters :KmerCounterPool<Kmer> = crossbeam_utils::thread::scope(|scope| {
        let start_t = time::Instant::now();
        let mut join_handles: Vec<Box<ScopedJoinHandle<KmerCounter<Kmer>>>> = Vec::with_capacity(nb_threads as usize);
        for i in 0..nb_threads {
            let mut nbseq_i = 0;
            // can generate thread
            let handle = scope.spawn(move |_| {
                let mut kmer_counter_i  = KmerCounter::new(0.03, capacity/(nb_threads as usize), nb_bits);
                let kmer_generator = KmerGenerator::new(kmer_size as u8);
                for seq in seqvec {
                    let vkmer : Vec<Kmer> = kmer_generator.generate_kmer(&seq);
                    for kmer in &vkmer {
                        let kmin = kmer.reverse_complement().min(*kmer);
                        // the idea to dispatch kmer after one more pass of hashing
                        // is taken from H.Li. Pb is that inequal repartition of AT/CG give a non uniform
                        // repartion of kmin.
                        if kmin.dispatch(nb_threads) == i  {
                            kmer_counter_i.insert_kmer(kmin);
                        }
                    }
                    nbseq_i += 1;
                    if nbseq_i % 5_000_000 == 0 {
                        println!(" nb  read treated by thread i = {} {} ", i, nbseq_i);
                    }
                } // end of for on seq
                kmer_counter_i
            });  // end of spawn
            println!(" thread  i= {} ", i);
            join_handles.push(Box::new(handle));
        } // end of for on threads
        //
        // now we must merge
        let mut kmer_counters : Vec<Box<KmerCounter<Kmer>>>  = Vec::with_capacity(nb_threads as usize);
        let mut nb_unique = 0u64;
        let mut nb_distinct = 0u64;
        for handle in join_handles {
            let thread_counter = handle.join().unwrap();
            nb_unique += thread_counter.cuckoo_f.len() as u64;
            nb_distinct += thread_counter.nb_distinct;
            println!(" nb once kmer = {} nb distinct = {}", thread_counter.cuckoo_f.len(), thread_counter.nb_distinct);           
            kmer_counters.push(Box::new(thread_counter));
            // transfer from thread_counter to kmer_counter
        }
        let elapsed_t = start_t.elapsed().whole_seconds();
        println!(" merging counters nb once kmer = {} nb distinct = {}", nb_unique, nb_distinct);           
        println!(" elapsed time (s) in threaded  counting scope {} ", elapsed_t);
        KmerCounterPool::new(kmer_counters)
    }).unwrap();  // end of scope
    // print some statistics
    println!(" merging counters nb distinct = {}", poolcounters.get_nb_distinct());
    //
    poolcounters
} // end of count_kmer_16b32bit_threaded



/// This function counts K-mers from a vector of Sequence
/// One thread generates all kmer and dispatch them to counter threads.
/// via messages.
// In fact as number of threads increase, it is more efficient than count_kmer_thread_independant

use std::cell::{RefCell};
use crossbeam_channel::{Sender};

pub fn count_kmer_threaded_one_to_many<Kmer>(seqvec : &Vec<Sequence>, nb_threads:usize, count_size:usize, kmer_size:usize) -> KmerCounterPool<Kmer>
where Kmer: CompressedKmerT+DispatchableT+Send,
      KmerGenerator<Kmer>: KmerGenerationPattern<Kmer>
{
    println!(" in counting kmer ... : count_kmer_threaded_one_to_many, kmer size : {}", kmer_size);
    //
    let capacity = 1_000_000_000;
    let nb_bits = count_size;
    //
    // In this counting function one thread generates all kmers and dispatch to a counter thread.
    // At the end of all threads we must gather results.
    // As bloom filters are union able its OK. All single kmer in cuckoo will stay unique as filters
    // realize a partition, but we will have to tranfer contents to one of cuckoo.
    // Filters must be a partition of the kmer set!!!
    //
        
    let poolcounters:KmerCounterPool<Kmer> = crossbeam_utils::thread::scope(|scope| {
        //
        let mut channels : Vec<RefCell< Sender<Kmer>  > >  = Vec::new();
        //
        let start_t = time::Instant::now();
        let mut join_handles: Vec<Box<ScopedJoinHandle<KmerCounter<Kmer>>>> = Vec::with_capacity(nb_threads as usize);
        // receptor threads
        for i in 0..nb_threads {
            let channel = crossbeam_channel::bounded::<Kmer>(1_000_000);
            channels.push(RefCell::new(channel.0));
            // can generate thread
            let receiver = channel.1;
            let receptor_handle = scope.spawn( move |_| {
                let mut kmer_counter_i  = KmerCounter::new(0.03, capacity/(nb_threads as usize), nb_bits);
                let mut nb_received = 0u64;
                receiver.for_each(|msg| { // just insert
                    kmer_counter_i.insert_kmer(msg);
                    nb_received += 1;
                    if nb_received % 500_000_000 == 0 {
                        println!(" nb  kmer treated  by thread i = {} {} ", i, nb_received);
                    }
                }); // end for_each
                println!(" nb  read kmer by thread i = {} {} ", i, nb_received);
                kmer_counter_i
            });// end of reciever thread.
            join_handles.push(Box::new(receptor_handle));
        }  // end of for i on receiver threads
        //        
        // sender thread
        //
        let kmer_generator = KmerGenerator::new(kmer_size as u8);
        let mut nbseq = 0;
        for seq in seqvec {
            let vkmer : Vec<Kmer> = kmer_generator.generate_kmer(&seq);
            for kmer in &vkmer {
                let kmin = kmer.reverse_complement().min(*kmer);
                // the idea to dispatch kmer after one more pass of hashing
                // is taken from H.Li. Pb is that inequal repartition of AT/CG give a non uniform
                // repartion of kmin.
                let i = kmin.dispatch(nb_threads);
                channels[i].borrow().send(kmin);
            }
            nbseq += 1;
            if nbseq % 10_000_000 == 0 {
                println!(" nb  read treated by sender  =  {} ", nbseq);
            }
        } // end of for on seq
        println!(" nb  read treated by sender  =  {} ", nbseq);
        for c in channels {
            drop(c.borrow());
        }
        //
        // now we must merge
        let mut kmer_counters : Vec<Box<KmerCounter<Kmer>>>  = Vec::with_capacity(nb_threads as usize);
        let mut nb_unique = 0u64;
        let mut nb_distinct = 0u64;
        for handle in join_handles {
            let thread_counter = handle.join().unwrap();
            nb_unique += thread_counter.cuckoo_f.len() as u64;
            nb_distinct += thread_counter.nb_distinct;
            println!(" nb once kmer = {} nb distinct = {}", thread_counter.cuckoo_f.len(), thread_counter.nb_distinct);           
            kmer_counters.push(Box::new(thread_counter));
            // transfer from thread_counter to kmer_counter
        }
        let elapsed_t = start_t.elapsed().whole_seconds();
        println!(" merging counters nb once kmer = {} nb distinct = {}", nb_unique, nb_distinct);           
        println!(" elapsed time (s) in count_kmer_threaded_one_to_many  {} ", elapsed_t);
        KmerCounterPool::new(kmer_counters)
    }).unwrap();  // end of scope
    //
    poolcounters
} // end of count_kmer_threaded_one_to_many




//==========================================================================================
//     Structure pour juste filter les Kmer a moindre cout en evitant le bloom filter
//
//     To be mage generic as KmerCounter is turns out to be really useful.
//==========================================================================================    

/// A Kmer filter to keep track of unique kmers without using bloom filter
/// Purpose  : try to find correlation between unique Kmer and quality.

/// It is up to the user to clean field all_f to free some memory by calling
/// It is 25% faster than using KmerCounter which use a bloom filter for keeping track of non unique kmers"

pub struct KmerFilter1 {
    /// kmer size in number of base
    kmer_size: u8,
    /// a cuckoo filter to keep track of kmer encountered only once
    once_f: CuckooFilter<DefaultHasher>,
    /// all kmer encountered
    all_f: CuckooFilter<DefaultHasher>,
}


impl KmerFilter1 {
    pub fn new(size: u8, capacity: u32) ->  KmerFilter1 {
        //
        KmerFilter1 {
            kmer_size: size,
            once_f: CuckooFilter::with_capacity(capacity as usize),
            all_f: CuckooFilter::with_capacity(2 * capacity as usize)
        }
    }
    /// insert a Kmer16b32bit
    pub fn insert_kmer16b32bit(&mut self, kmer: Kmer16b32bit) {
        //
        // we should avoid hashing twice (once for each counter) the same kmer...
        // could implement a hasher that does nothing for a type representing a already hashed 32bit kmer
        // we store a kmer in cuckoo filter when we see it for the first time.
        // when seen for the second time we delete it from the cf and insert it twice in cbloom_f
        // As we store hashed kmer by invertible hash in the bloom filter the hash is done twice...
        //
        if self.all_f.contains(& kmer.0) {
            if self.once_f.contains(& kmer.0) {
                self.once_f.delete(& kmer.0);  // remove beccause this kmer is not unique
            }
        }
        else {
            // it was seen 0 or once
            if !self.once_f.contains(& kmer.0) {
                // it was never seen,  we insert it
                self.all_f.add(&kmer.0).unwrap();                
                self.once_f.add(&kmer.0).unwrap();                
            }
            
        }  // end not yet inserted in cbloom       
    } // end of insert_kmer

    /// This function dumps in fname the only once kmers. We want information on position to correlate with Quality. 
    /// format of file is :
    ///    - a u32 magic : 0xcea2bbdd
    ///    - kmer size as u8
    ///    - number of kmers as a u64. Set to 0u64 if unknown else > 0 (to get fast allocation and check on file when reloading)
    ///    - then list of (kmers,numseq,numkmer)
    //
    pub fn dump_in_file_once_kmer16b32bit(&self, fname: &String, seqvec : &Vec<Sequence>) -> std::result::Result<usize, io::Error> {
        //
        println!(" dumping once_kmer in file : {} ", fname);
        //
        let magic: u32  = 0xcea2bbdd;
        let kmer_generator = KmerGenerator::<Kmer16b32bit>::new(16);
        // as kmer generation is fast we generate once more all kmers and check for those that are in once_f
        let file = OpenOptions::new().write(true).create(true).truncate(true).open(&fname)?;
        let mut bufw : io::BufWriter<fs::File> = io::BufWriter::with_capacity(1_000_000_000, file);
        let mut nb_kmer_dumped = 0;
        // write magic, kmer_size and number of kmers
        bufw.write(unsafe { &mem::transmute::<u32, [u8;4]>(magic) } )?;
        bufw.write(unsafe { &mem::transmute::<u8, [u8;1]>(self.kmer_size) })?;
        let nb_kmer : u64 = self.once_f.len() as u64;
        bufw.write(unsafe { &mem::transmute::<u64, [u8;8]>(nb_kmer) } )?;
        let mut numseq = 0u32;
        // now we must generate kmer once more
        for seq in seqvec {
            let vkmer : Vec<Kmer16b32bit> = kmer_generator.generate_kmer(&seq);
            let mut numkmer = 0u32;
            for kmer in &vkmer {
                let kmin = kmer.reverse_complement().min(*kmer);
                if self.once_f.contains(&kmin.0) {
                    (kmin as Kmer16b32bit).dump(&mut bufw)?;
                    // dump numseq and numkmer 4 bytes each
                    bufw.write(unsafe { &mem::transmute::<u32, [u8;4]>(numseq) })?; 
                    bufw.write(unsafe { &mem::transmute::<u32, [u8;4]>(numkmer)} )?;
                    nb_kmer_dumped += 1;
                }
                numkmer += 1;
            }
            numseq += 1;
        }
        println!("dump_in_file_once_kmer16b32bit, number of kmer dumped : {} ", nb_kmer_dumped);
        //
        Ok(nb_kmer_dumped)
    }  // end of dump_in_file

    /// to free memory occupied by kmer seen more than once
    pub fn eliminate_not_once_kmer(&mut self) {
        unimplemented!();
    }
    
} // end of implementation KmerFilter1


/// a function for counting identifying unique kmers.
pub fn filter1_kmer_16b32bit(seqvec : &Vec<Sequence>) -> KmerFilter1 {
   println!(" in filtering unique kmer ...");
    let start_t = time::Instant::now();

    let mut nb_kmer : u64 = 0;
    let capacity = 1000_000_000;
    let mut nbseq = 0u64;
    //
    let mut kmer_counter = KmerFilter1::new(16u8, capacity);
    println!("space occupied by KmerFilter1 filter {}", 2* kmer_counter.once_f.memory_usage());
    //
    for seq in seqvec {
        let vkmer : Vec<Kmer16b32bit> = KmerGenerator::<Kmer16b32bit>::new(16).generate_kmer(&seq);
        for kmer in &vkmer {
            // do not forget reverse complement manip
            let kmin = kmer.reverse_complement().min(*kmer);
            kmer_counter.insert_kmer16b32bit(kmin);
        }
        nb_kmer += vkmer.len() as u64;
        nbseq += 1;
       if nbseq % 100_000 == 0 {
           println!(" nb  read treated = {} ", nbseq);
           println!(" nb once kmer = {} nb distinct = {} nb generated = {} ", kmer_counter.once_f.len(), kmer_counter.all_f.len() , nb_kmer);
        }        
    }
    let elapsed_t = start_t.elapsed().whole_seconds();
    println!(" elapsed time (s) in filter1_kmer_32b counting {} ", elapsed_t);
    println!(" nb kmer generated {}", nb_kmer);
    //
    kmer_counter
}

//==============================================================================================
//  Structure to reload in a hashmap the dump of counting filter
//  Temporary implementation. What is the cost?
//                            Can we see a correlation between unique kmers and low qualities.
//==============================================================================================


#[derive(Copy,Clone,PartialEq)]
pub enum CounterType {
    /// just list unique kmers
    Unique,
    /// list only multiple kmer with count
    Multiple,
}


// as long as we store Kmer16b32bit we need only slots of size u32.
// from julia it is the field positions that is really useful.
// hmap could be useful to perturbate kmer and check status of new kmer.
/// Works only for up to 32 bit kmers Now. To be made generic.
/// This structure is used to reload a dump of counting filter for Kmer16b32bit
/// It associate a slot to Kmer16b32bit where we can find coordinate or count of kmer.


pub struct  KmerCountReload {
    ///
    _c_type: CounterType,
    /// kmer size
    kmer_size: u32,
    /// nb_kmer in file
    nb_kmer: usize,
    /// a hashmap mapping kmer to a slot in description array
    // The indirect acces to positions via hmap cost us some memory but enables Julia to ask positions of unique kmers via position
    h_pos: Option<fnv::FnvHashMap<u32,u32> >,
    ///
    positions: Vec<KmerCoord>,
    /// This is the field that corresponds to count of multiple kmers.
    h_count: Option<fnv::FnvHashMap<u32,u16 > >,
}  // 


impl KmerCountReload {
    /// a constructor that just does allocation once we know size
    pub fn new(c_type: CounterType, kmer_size: u32, size:usize) -> KmerCountReload {
        let mut h_pos : Option<fnv::FnvHashMap<u32,u32> > = None;
        let mut h_count : Option<fnv::FnvHashMap<u32,u16 > > = None;
        let mut positions : Vec<KmerCoord> = Vec::new();
        match c_type {
            CounterType::Unique   => {
                h_pos = Some(FnvHashMap::with_capacity_and_hasher(size, Default::default()));
                positions.reserve(size);
            } ,
            CounterType::Multiple => {
                h_count = Some(fnv::FnvHashMap::with_capacity_and_hasher(size, Default::default()));
            },
        };
        KmerCountReload{_c_type: c_type,  kmer_size:kmer_size, nb_kmer:size, h_pos: h_pos, positions: positions, h_count: h_count}
    }  // end of new
    //
    // insert a kmer data vectors positions and count and insert resulting insertion slot in hashmap 
    fn insert_kmer_pos(&mut self, kmer:u32, numseq:u32, kmer_pos:u32) {
        // get last free slot
        let slot = self.positions.len();
        let coord = KmerCoord{read_num:numseq, pos:kmer_pos};
        // insert in hmap (kmer,slot)
        self.h_pos.as_mut().unwrap().insert(kmer,slot as u32);
        self.positions.push(coord);
    }
    ///
    fn insert_kmer_count(&mut self, kmer:u32, count : u16) {
        self.h_count.as_mut().unwrap().insert(kmer, count);
    }
    //
    pub fn get_kmer_size(&self) -> u32 { self.kmer_size }
    /// returns the number of kmer declared in file
    pub fn get_nb_kmer(&self) -> usize { self.nb_kmer }
    //
    //
    // private methods for multiple kmer case , part of constructing strategy.
    // We cannot know exactly how many kmers we have to reload as we use
    // another cuckoo/bloom filter when dumping. The number of Kmer
    // in the dump file is an approximation subject to false positive error.
    // (There should some more in the file than announced in head of file.
    // So we must loop until EOF encounter.
    //
    pub fn load_multiple_kmers_from_file(fname: &String) -> Option<Box<KmerCountReload> >  {
        //
        info!(" in  KmerCountReload::load_multiple_kmers");
        let start_t = time::Instant::now();
        //
        let file_r = OpenOptions::new().read(true).open(&fname);
        let file;
        match file_r {
            Ok(sthing) => {file = sthing},
            Err(_e) => { println!("KmerCountReload::load_multiple_kmers_from_file cannot open file {} ", fname);
                         return None;
            }
        }
        let mut bufr: io::BufReader<fs::File> = io::BufReader::with_capacity(1_000_000_000, file);
        // we need some 4 bytes buffer
        let mut buf_4bytes_1 = [0u8;4];
        let mut buf_1bytes_1 = [0u8;1];
        let mut buf_2bytes_1 = [0u8;2];
        let mut io_res;
        let magic;
        let nb_kmer:u64;
        let kmer_size;
        let nb_bytes_by_count;
        //
        // now read header
        //
        // write magic (4 bytes), kmer_size, size of counter,  and number of kmers
        io_res = bufr.read_exact(&mut buf_4bytes_1);
        if io_res.is_err() {
            println!("KmerCountReload::load_multiple_kmers_from_file could no read magic");
            return None;
        }
        else { // to be replaced by from_bytes as soon as API goes from nightly to stable
            magic = unsafe { mem::transmute::<[u8;4], u32>(buf_4bytes_1)}
        }
        let c_type = match magic {
            COUNTER_MULTIPLE   => CounterType::Multiple,
            _ => {
                println!("KmerCountReload::load_multiple_kmers_from_file unknow magic {:x} ", magic);
                return None;
            }
        };
        // read kmer size 1 byte
        io_res = bufr.read_exact(&mut buf_1bytes_1);
        if io_res.is_err() {
            println!("KmerCountReload::load_from_file could no read kmer_size");
            return None;
        }
        else { // to be replaced by from_bytes as soon as API goes from nightly to stable
            kmer_size = unsafe { mem::transmute::<[u8;1], u8>(buf_1bytes_1)};
            println!("KmerCountReload::load_multiple_kmers_from_file got kmer size {}", kmer_size);
        }
        // read counter size (one byte)
        io_res = bufr.read_exact(&mut buf_1bytes_1);
        if io_res.is_err() {
            println!("KmerCountReload::load_multiple_kmers_from_file could no read kmer_size");
            return None;
        }
        else { // to be replaced by from_bytes as soon as API goes from nightly to stable
            nb_bytes_by_count = unsafe { mem::transmute::<[u8;1], u8>(buf_1bytes_1)};
            println!("KmerCountReload::load_multiple_kmers_from_file got counter size {}", nb_bytes_by_count);
        }
        if nb_bytes_by_count > 2 {
            println!("load_multiple_kmers, kmer count on more than 2 bytes not yet implemented");
            return None;
        }
        // read number of kmers , 8 bytes
        let mut buf_8bytes =  [0u8;8];
        io_res = bufr.read_exact(&mut buf_8bytes);
        if io_res.is_err() {
            println!("KmerCountReload::load_multiple_kmers_from_file could read number of kmers");
            return None;
        }
        else { // to be replaced by from_bytes as soon as API goes from nightly to stable
            nb_kmer = unsafe { mem::transmute::<[u8;8], u64>(buf_8bytes)};
            println!("KmerCountReload::load_multiple_kmers_from_file got nb kmer {}", nb_kmer);
        }
        // increase a little nb_kmer to take into account false positive in dumping
        let nb_kmer_plus = (1.03 * nb_kmer as f64) as u64;
        println!("KmerCountReload::load_multiple_kmers_from_file allocating with kmer_size nb_kmer {} {}", kmer_size, nb_kmer_plus); 
        let mut kmer_count = KmerCountReload::new(c_type, kmer_size as u32, nb_kmer_plus as usize);
        //
        // Now we loop and read
        // As the dump uses a count filter , the number of kmer dumped is exact at some %.
        // So we loop until EOF.
        //
        let mut count:u16;
        let mut kmer;
        let mut nb_kmer_read = 0;
        //
        loop {
            // read kmer , depends upon size
            io_res = bufr.read_exact(&mut buf_4bytes_1);
            // if there is an EOF it is just now
            match io_res {
                Err(e) => {
                    if e.kind() == io::ErrorKind::UnexpectedEof {
                        // can check order of magnitude of kmer read
                        println!("KmerCountReload::load_multiple_kmers_from_file , EOF encountered");
                        println!(" nb kmer declared in file {} , nb loaded {}", nb_kmer, nb_kmer_read);  
                        let elapsed_t = start_t.elapsed().whole_seconds();
                        println!(" elapsed time (s) in KmerCountReload:: {} ", elapsed_t);
                        return Some(Box::new(kmer_count));
                    }
                    else {
                        println!("KmerCountReload::load_from_file could no read kmer");
                        return None;
                    };
                },
                Ok(()) => {
                    kmer = unsafe { mem::transmute::<[u8;4], u32>(buf_4bytes_1)};
                },
            };
            // read count on 1 or 2 bytes.
            if nb_bytes_by_count == 1 {
                io_res = bufr.read_exact(&mut buf_1bytes_1);
                if io_res.is_err() {
                    println!("KmerCountReload::load_from_file could no seq num");
                    return None;
                }
                else { // to be replaced by from_bytes as soon as API goes from nightly to stable
                    count = unsafe { mem::transmute::<[u8;1], u8>(buf_1bytes_1)} as u16
                }
            }
            else {
                io_res = bufr.read_exact(&mut buf_2bytes_1);
                if io_res.is_err() {
                    println!("KmerCountReload::load_from_file could no seq num");
                    return None;
                }
                else { // to be replaced by from_bytes as soon as API goes from nightly to stable
                    count = unsafe { mem::transmute::<[u8;2], u16>(buf_2bytes_1)} as u16
                }
            }
            // now we can insert kmer.
            if nb_kmer_read % 100_000_000 == 0 {
                println!("KmerCountReload::load_multiple_kmers  nb_kmer kmer, count {} {}", nb_kmer_read, count);                
            }
            kmer_count.insert_kmer_count(kmer, count);
            nb_kmer_read += 1;
        } // end loop
        //
    }  // end of load_multiple_kmers_from_file

    
    /// reload the dump of a counting filter in KmerCountReload
    // Read characteristics of dump and then call specialized functions depending on dump type
    pub fn load_unique_kmer_from_file(fname: &String) -> Option<Box<KmerCountReload> > {
        //
        info!(" in  KmerCountReload::load_unique_kmer_from_file {} ", fname);
        println!(" in  KmerCountReload::load_unique_kmer_from_file");
        let start_t = time::Instant::now();
        //
        let file_r = OpenOptions::new().read(true).open(&fname);
        let file;
        match file_r {
            Ok(sthing) => {file = sthing},
            Err(_e) => { println!("KmerCountReload::load_from_file cannot open file {} ", fname);
                         return None;
            }
        }
        let mut bufr: io::BufReader<fs::File> = io::BufReader::with_capacity(1_000_000_000, file);
        //
        let mut io_res;
        let magic;
        let nb_kmer:u64;
        let kmer_size;
        // we need some 4 bytes buffer
        let mut buf_4bytes_1 = [0u8;4];
        // write magic (4 bytes), kmer_size and number of kmers
        io_res = bufr.read_exact(&mut buf_4bytes_1);
        if io_res.is_err() {
            println!("KmerCountReload::load_from_file could no read magic");
            return None;
        }
        else { // to be replaced by from_bytes as soon as API goes from nightly to stable
            magic = unsafe { mem::transmute::<[u8;4], u32>(buf_4bytes_1)}
        }
        let c_type = match magic {
            COUNTER_UNIQUE   => CounterType::Unique,
            _ => {
                println!("KmerCountReload::load_from_file unknow magic {:x} ", magic);
                return None;
            }
        };
        // read kmer size 1 byte
        let mut b=[0u8;1];
        io_res = bufr.read_exact(&mut b);
        if io_res.is_err() {
            println!("KmerCountReload::load_from_file could no read kmer_size");
            return None;
        }
        else { // to be replaced by from_bytes as soon as API goes from nightly to stable
            kmer_size = unsafe { mem::transmute::<[u8;1], u8>(b)};
            println!("KmerCountReload::load_from_file got kmer size {}", kmer_size);
        }
        // read number of kmers , 8 bytes
        let mut buf_8bytes =  [0u8;8];
        io_res = bufr.read_exact(&mut buf_8bytes);
        if io_res.is_err() {
            println!("KmerCountReload::load_from_file could read number of kmers");
            return None;
        }
        else { // to be replaced by from_bytes as soon as API goes from nightly to stable
            nb_kmer = unsafe { mem::transmute::<[u8;8], u64>(buf_8bytes)};
            println!("KmerCountReload::load_from_file got nb kmer {}", nb_kmer);
        }
        // increase a little nb_kmer to take into account false positive in dumping
        let nb_kmer_plus = (1.03 * nb_kmer as f64) as u64;
        println!("KmerCountReload::load_from_file allocating with kmer_size nb_kmer {} {}", kmer_size, nb_kmer_plus); 
        let mut kmer_count = KmerCountReload::new(c_type, kmer_size as u32, nb_kmer_plus as usize);
        //
       // Now we loop and read
        let mut numseq;
        let mut kmerpos;
        let mut kmer;
        //
        for ikmer in 0..nb_kmer {
            // read kmer , depends upon size
            io_res = bufr.read_exact(&mut buf_4bytes_1);
            if io_res.is_err() {
                println!("KmerCountReload::load_from_file could no read kmer");
                return None;
            }
            else { // to be replaced by from_bytes as soon as API goes from nightly to stable
                kmer = unsafe { mem::transmute::<[u8;4], u32>(buf_4bytes_1)}
            }
            // read numseq
            io_res = bufr.read_exact(&mut buf_4bytes_1);
            if io_res.is_err() {
                println!("KmerCountReload::load_from_file could no seq num");
                return None;
            }
            else { // to be replaced by from_bytes as soon as API goes from nightly to stable
                numseq = unsafe { mem::transmute::<[u8;4], u32>(buf_4bytes_1)}
            }
            // read kmer pos
            io_res = bufr.read_exact(&mut buf_4bytes_1);
            if io_res.is_err() {
                println!("KmerCountReload::load_from_file error reading kmer pos");
                return None;
            }
            else { // s we have kmer_size, we will be able to generalize.
                kmerpos = unsafe { mem::transmute::<[u8;4], u32>(buf_4bytes_1)};
            }
            // now we can insert kmer.
            if ikmer == 0 {
                println!("KmerCountReload::load_from_file first kmer, numseq, kpos {} {} {} ", kmer, numseq, kmerpos);                
            }
            else {
                if ikmer % 10_000_000 == 0 {
                    println!("KmerCountReload::load_unique_kmer_from_file  kmer rank  kmer, numseq, kpos {} {} {} {} ",
                             ikmer, kmer, numseq, kmerpos);                
                }
            }
           kmer_count.insert_kmer_pos(kmer, numseq, kmerpos);
        }
        //
        let elapsed_t = start_t.elapsed().whole_seconds();
        println!(" elapsed time (s) in KmerCountReload::load_from_file {} ", elapsed_t);
        println!(" nb kmer loaded {}", nb_kmer);
        //
        Some(Box::new(kmer_count))
    } // end of KmerCountReload::load_from_file
    //
    // Some utilities just to be accessed from Julia
    //
    /// returns the position of a (unique!!) kmer. To analyze quality/kmer correlation.
    pub fn get_unique_kmer_coord(&self, _kmer:Kmer16b32bit) -> Option<KmerCoord> {
        unimplemented!();        
    }
    /// get a kmer position from its rank. This is mostly for julia interface
    pub fn get_coord_from_rank(&self, rank:usize) -> Option<KmerCoord> {
        if rank < self.positions.len() {
            Some(self.positions[rank])
        }
        else {
            None
        }
    }  // end of get_coord_from_rank

    /// finally returns an option on counts
    pub fn get_multi_kmer_counts(&self) -> Option<Vec<u16>> {
        if self._c_type != CounterType::Multiple || self.h_count == None {
            None
        }
        // avoid the call to unwrap and use ref to avoid the move out! Rust opaque tricks...
        else if let Some(ref hmap) = self.h_count {
            let v = hmap.values().map(|x| *x).collect();
            Some(v)
        }
        else { None}
    }  // end of get_multi_kmer_counts


    
}  // end of KmerCountReload



///////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;
    use std::ops::Range;
    use rand::distributions::*;
    #[test]
    fn test_kmer_counter() {
        let fp_rate = 0.03;
        let nb_random = 1_000_000;
        //
        let mut kmer_counter = KmerCounter::new(fp_rate, 10_000_000, 8);
        //
        let seqstr = String::from("TCAAAGGGAAACATTCAAAATCAGTATGCGCCCGTTCAGTTACGTATTGCTCTCGCTAATGAGATGGGCTGGGTACAGAG");
        // transform to a Vec<u8> and then a Sequence
        let slu8 = seqstr.as_bytes();
        // get a sequence with 2 bits compression
        let seq = Sequence::new(&slu8,2);
        let vkmer : Vec<Kmer16b32bit> = KmerGenerator::new(16).generate_kmer(&seq);
        // now we will sample 1_000_000 kmers
        let between = Uniform::<usize>::new(2, vkmer.len());        
        let mut rng = rand::thread_rng();
        let mut xsi;
        // we insert once vkmer[0]
        kmer_counter.insert_kmer(vkmer[0]);
        // we insert once vkmer[1]
        kmer_counter.insert_kmer(vkmer[1]);
        for _ in 0..nb_random {
            xsi = between.sample(&mut rng);
            let kmer = vkmer[xsi];
            kmer_counter.insert_kmer(kmer);
        }
        // now we insert once more vkmer[1]
        kmer_counter.insert_kmer(vkmer[1]);
        //
        // now check counts for first 2 kmers
        //
        println!(" kmer 0 : {}", kmer_counter.get_count(vkmer[0]));
        println!(" kmer 1 : {}", kmer_counter.get_count(vkmer[1]));
        assert_eq!(kmer_counter.get_count(vkmer[0]), 1);
        assert_eq!(kmer_counter.get_count(vkmer[1]), 2);
        //
        //
        let countvec: Vec<u32> = (0..vkmer.len()).map(|i| kmer_counter.get_count(vkmer[i])).collect();
        let maxpos = (0..countvec.len()).max_by_key(|&x| countvec[x]).unwrap() as usize;
        let minpos = (0..countvec.len()).min_by_key(|&x| countvec[x]).unwrap() as usize;
        println!(" kmer maxpos {} minpos {} ", maxpos, minpos);
        // now we can check the error rate ?
        for i in 2..countvec.len() {
            println!(" pos count {}  {} ", i, countvec[i]);
        }
        println!(" kmer min {} count {} ", minpos, countvec[minpos]);
        println!(" kmer max {} count {} ", maxpos, countvec[maxpos]);
        println!(" nb unique {}, nb distinct {}", kmer_counter.get_nb_unique(), kmer_counter.get_nb_distinct());

//        assert_eq!(1,0);       
    } // end of test_kmer_counter_alloc


    
    #[test]
    fn test_false_positive() {
        let mut kmer_counter = KmerCounter::new(0.03, 10_000_000, 8);
        //
        let seqstr = String::from("TCAAAGGGAAACATTCAAAATCAGTATGCGCCCGTTCAGTTACGTATTGCTCTCGCTAATGAGATGGGCTGGGTACAGAG");
        // transform to a Vec<u8> and then a Sequence
        let slu8 = seqstr.as_bytes();
        // get a sequence with 2 bits compression
        let seq = Sequence::new(&slu8,2);
        // seq has 80 chars, we get 65 = 80-16+1 kmers of 16 bases
        let vkmer : Vec<Kmer16b32bit> = KmerGenerator::new(16).generate_kmer(&seq);
        println!(" got nb kmers : {} ", vkmer.len());
        // now we will sample 1_000_000 kmers in upper half of kmers so from kmer 31 to 
        let between = Uniform::from(Range{start: vkmer.len()/2, end: vkmer.len()});        
        let mut rng = rand::thread_rng();
        let nb_random = 1_000_000;
        let mut xsi;
        for _ in 0..nb_random {
            xsi = between.sample(&mut rng);
            let kmer = vkmer[xsi];
            kmer_counter.insert_kmer(kmer);
        }
        // now check if we see something between 0 and vkmer.len()/2 - 1
        let countvec: Vec<u32> = (0..vkmer.len()).map(|i| kmer_counter.get_count(vkmer[i])).collect();
        let maxpos = (0..countvec.len()).max_by_key(|&x| countvec[x]).unwrap() as usize;
        let minpos = (0..countvec.len()).min_by_key(|&x| countvec[x]).unwrap() as usize;
        println!(" kmer maxpos {} minpos {} ", maxpos, minpos);
        for i in 0..countvec.len() {
            if i < countvec.len() /2 {
                assert!(countvec[i] == 0);
            }
            else {
                assert!(countvec[i] == 255);
            }
            println!(" pos count {}  {} ", i, countvec[i]);
        }
        println!(" kmer min {} count {} ", minpos, countvec[minpos]);
        println!(" kmer max {} count {} ", maxpos, countvec[maxpos]);

    } // end of test_false_positive


    
} // end mod tests
