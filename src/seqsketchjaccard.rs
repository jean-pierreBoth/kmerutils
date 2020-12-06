//! This module provides sequence signature computation and Jaccard probability index using the probminhash crate.  
//! The Jaccard probability index is a Jaccard index between sequences taking into account multiplicity of kmers. 
//! It is adapted to short read sequences, for long read it is necessary to split sequences into blocks, see module seqblcoksketch.
//! 
//! The kmers of a given size are generated for each sequence, kmers lists are hashed by the probminhash algorithm 
//! and a jaccard index between sequences is computed.
//! 
//! It is also possible to ask for the common Kmers found in the signature of 2 sequences
//! 
//! 
//! 

#![allow(dead_code)]



#[allow(unused_imports)]
use log::*;

use std::io;
use std::io::{Write,Read};
use std::io::{ErrorKind};

use std::fs;
use std::fs::OpenOptions;
use std::fmt::{Debug};

#[allow(unused_imports)]
use std::hash::{BuildHasher, BuildHasherDefault, Hasher, Hash};
use indexmap::{IndexMap};
use fnv::{FnvBuildHasher};

use crate::nohasher::*;

use crate::kmergenerator::*;

use rayon::prelude::*;

type FnvIndexMap<K, V> = IndexMap<K, V, FnvBuildHasher>;

use probminhash::probminhasher::*;


/// given 2 weighted set given as IndexMap, compute weighted jaccard index and return common objects if any or None
pub fn compute_probminhash3a_jaccard<D,H, Hidx>(idxa : &IndexMap<D, f64, Hidx>, idxb: &IndexMap<D, f64, Hidx>, 
                                                        sketch_size : usize, return_object: bool)  -> (f64, Option<Vec<D>>)
            where   D : Copy+Eq+Hash+Debug+Default,
                    H:  Hasher+Default ,
                    Hidx : std::hash::BuildHasher {
        //
        let mut pminhasha = ProbMinHash3a::<D,H>::new(sketch_size, D::default());
        pminhasha.hash_weigthed_idxmap(idxa);
        let mut pminhashb = ProbMinHash3a::<D,H>::new(sketch_size, D::default());
        pminhashb.hash_weigthed_idxmap(idxb);
        let siga = pminhasha.get_signature();
        let sigb = pminhashb.get_signature();
        let jac : f64;
        if !return_object {
            jac = compute_probminhash_jaccard(&siga, &sigb);
            return (jac,None);
        }
        else {
            return probminhash_get_jaccard_objects(&siga, &sigb);
        }
}  // end of compute_probminhash3a_jaccard


/// given 2 signatures of objects sets, compute weighted jaccard index and return common objects if any or None
pub fn probminhash_get_jaccard_objects<D:Eq+Copy>(siga : &Vec<D>, sigb : &Vec<D>) -> (f64, Option<Vec<D>>) {
    let sig_size = siga.len();
    assert_eq!(sig_size, sigb.len());
    //
    let mut common_objects =  Vec::<D>::new();
    let mut inter = 0;
    for i in 0..siga.len() {
        if siga[i] == sigb[i] {
            inter += 1;
            common_objects.push(siga[i]);
        }
    }
    let jp = inter as f64/siga.len() as f64;
    //
    if jp > 0. {
        return (jp,Some(common_objects));
    }
    else {
        return (0., None);
    }
}  // end of compute_probminhash_jaccard





/// compute jaccard probability index between a sequence and a vector of sequences for Kmer16b32bit
/// and returns a vector of Jaccard probability index.  
/// This function is threaded (with Rayon) so it is best used with a vector of sequence of sufficient size
pub fn jaccard_index_probminhash3a_kmer16b32bit<F>(seqa: &Sequence, vseqb : &Vec<Sequence>, sketch_size: usize, fhash : F) -> Vec<f64> 
    where F : Fn(&Kmer16b32bit) -> u32 + Sync + Send{
    //
    info!("seqsketcher : entering sketch_seqrange_probminhash3a_kmer16b32bit");
    // a vector to return results
    let mut jaccard_vec = Vec::<f64>::with_capacity(vseqb.len());
    for _ in 0..vseqb.len() {
        jaccard_vec.push(0.);
    }
    // default is invertible hash and then superminhash without any hashing
    let mut pminhasha = ProbMinHash3a::<usize,NoHashHasher>::new(sketch_size, 0);
    let mut wa : FnvIndexMap::<usize,f64> = FnvIndexMap::with_capacity_and_hasher(seqa.size(), FnvBuildHasher::default());
    //
    // generate all kmers include in range arg. dependance upon kmer_size 
    // seqa
    let mut kmergen = KmerSeqIterator::<Kmer16b32bit>::new(16, &seqa);
    kmergen.set_range(0, seqa.size()).unwrap();
    loop {
        match kmergen.next() {
            Some(kmer) => {
                let hashval = fhash(&kmer);
                *wa.entry(hashval as usize).or_insert(0.) += 1.;
            },
            None => break,
        }
    }  // end loop
    pminhasha.hash_weigthed_idxmap(&wa);
    let siga = pminhasha.get_signature();
    // loop on vseqb to // with rayon
    let comput_closure = | seqb : &Sequence, i:usize | -> (usize,f64) {
        let mut wb : FnvIndexMap::<usize,f64> = FnvIndexMap::with_capacity_and_hasher(seqb.size(), FnvBuildHasher::default());
        let mut kmergen = KmerSeqIterator::<Kmer16b32bit>::new(16, &seqb);
        kmergen.set_range(0, seqb.size()).unwrap();
        loop {
            match kmergen.next() {
                Some(kmer) => {
                    let hashval = fhash(&kmer);
                    *wb.entry(hashval as usize).or_insert(0.) += 1.;
                },
                None => break,
            }
        }  // end loop 
        let mut pminhashb = ProbMinHash3a::<usize,NoHashHasher>::new(sketch_size, 0);
        pminhashb.hash_weigthed_idxmap(&wb);
        let sigb = pminhashb.get_signature();
        let jac = compute_probminhash_jaccard(siga, sigb);        
        return (i,jac);
    };
    let jac_with_rank : Vec::<(usize,f64)> = (0..vseqb.len()).into_par_iter().map(|i| comput_closure(&vseqb[i],i)).collect();
    // re-order from jac_with_rank to jaccard_vec as the order of return can be random!!
    for i in 0..jac_with_rank.len() {
        let slot = jac_with_rank[i].0;
        jaccard_vec[slot] = jac_with_rank[i].1;
    }
    return jaccard_vec;
    //
} // end of sketch_seqrange_probminhash3a_kmer16b32bit


/// This function computes and return signatures of a vector of sequences by generating kmers of size kmer_size. 
/// The sketch is done with probminhash3a algorithm.   
/// The size of signature of each sequence is sketch_size.  
/// fhash is any hash function, but usually it is identity, invhash on kmer or on min of kmer and reverse complement.
/// These are the hash function that make possible to get back to the original kmers (or at least partially in the case using the min)
/// 
pub fn sketch_probminhash3a_kmer32bit<F>(vseq : &Vec<Sequence>, sketch_size: usize, kmer_size : u8, fhash : F) -> Vec<Vec<u32> >
    where F : Fn(&Kmer32bit) -> u32 + Send + Sync {
    //
    let comput_closure = | seqb : &Sequence, i:usize | -> (usize,Vec<u32>) {
        let mut wb : FnvIndexMap::<u32,f64> = FnvIndexMap::with_capacity_and_hasher(seqb.size(), FnvBuildHasher::default());
        let mut kmergen = KmerSeqIterator::<Kmer32bit>::new(kmer_size, &seqb);
        kmergen.set_range(0, seqb.size()).unwrap();
        loop {
            match kmergen.next() {
                Some(kmer) => {
                    let hashval = fhash(&kmer);
                    *wb.entry(hashval as u32).or_insert(0.) += 1.;
                },
                None => break,
            }
        }  // end loop 
        let mut pminhashb = ProbMinHash3a::<u32,NoHashHasher>::new(sketch_size, 0);
        pminhashb.hash_weigthed_idxmap(&wb);
        let sigb = pminhashb.get_signature();
        // get back from usize to Kmer32bit ?. If fhash is inversible possible, else NO.
        return (i,sigb.clone());
    };
    //
    let sig_with_rank : Vec::<(usize,Vec<u32>)> = (0..vseq.len()).into_par_iter().map(|i| comput_closure(&vseq[i],i)).collect();
    // re-order from jac_with_rank to jaccard_vec as the order of return can be random!!
    let mut jaccard_vec = Vec::<Vec<u32>>::with_capacity(vseq.len());
    for _ in 0..vseq.len() {
        jaccard_vec.push(Vec::new());
    }
    // CAVEAT , boxing would avoid the clone?
    for i in 0..sig_with_rank.len() {
        let slot = sig_with_rank[i].0;
        jaccard_vec[slot] = sig_with_rank[i].1.clone();
    }
    jaccard_vec
}  // end of sketchprobminhash3a_kmer32bit




/// This function computes and return signatures of a vector of sequences by generating kmers of size kmer_size. 
/// The sketch is done with probminhash3a algorithm.   
/// The size of signature of each sequence is sketch_size.  
/// fhash is any hash function, but usually it is identity, invhash on kmer or on min of kmer and reverse complement.
/// These are the hash function that make possible to get back to the original kmers (or at least partially in the case using the min)
/// 
pub fn sketch_probminhash3_kmer32bit<F>(vseq : &Vec<Sequence>, sketch_size: usize, kmer_size : u8, fhash : F) -> Vec<Vec<u32> >
    where F : Fn(&Kmer32bit) -> u32 + Send + Sync {
    //
    let comput_closure = | seqb : &Sequence, i:usize | -> (usize,Vec<u32>) {
        let mut wb : FnvIndexMap::<u32,f64> = FnvIndexMap::with_capacity_and_hasher(seqb.size(), FnvBuildHasher::default());
        let mut kmergen = KmerSeqIterator::<Kmer32bit>::new(kmer_size, &seqb);
        kmergen.set_range(0, seqb.size()).unwrap();
        loop {
            match kmergen.next() {
                Some(kmer) => {
                    let hashval = fhash(&kmer);
                    *wb.entry(hashval as u32).or_insert(0.) += 1.;
                },
                None => break,
            }
        }  // end loop 
        let mut pminhashb = ProbMinHash3a::<u32,NoHashHasher>::new(sketch_size, 0);
        pminhashb.hash_weigthed_idxmap(&wb);
        let sigb = pminhashb.get_signature();
        // get back from usize to Kmer32bit ?. If fhash is inversible possible, else NO.
        return (i,sigb.clone());
    };
    //
    let sig_with_rank : Vec::<(usize,Vec<u32>)> = (0..vseq.len()).into_par_iter().map(|i| comput_closure(&vseq[i],i)).collect();
    // re-order from jac_with_rank to jaccard_vec as the order of return can be random!!
    let mut jaccard_vec = Vec::<Vec<u32>>::with_capacity(vseq.len());
    for _ in 0..vseq.len() {
        jaccard_vec.push(Vec::new());
    }
    // CAVEAT , boxing would avoid the clone?
    for i in 0..sig_with_rank.len() {
        let slot = sig_with_rank[i].0;
        jaccard_vec[slot] = sig_with_rank[i].1.clone();
    }
    jaccard_vec
}  // end of sketchprobminhash3_kmer32bit







/// Compute jaccard probability index between a sequence and a vector of sequences for Kmer32bit with probminhash3a.      
/// It returns a vector of Jaccard probability index.
/// the fhash function is a hash function.  
/// The function is threaded with the Rayon crate.
pub fn jaccard_index_probminhash3a_kmer32bit<F>(seqa: &Sequence, vseqb : &Vec<Sequence>, sketch_size: usize, 
                    kmer_size : u8, fhash : F) -> Vec<f64> 
                    where F : Fn(&Kmer32bit) -> u32 + Send + Sync {
    //
    info!("seqsketcher : entering compute_jaccard_index_probminhash3a_kmer32bit");
    // a vector to return results
    let mut jaccard_vec = Vec::<f64>::with_capacity(vseqb.len());
    for _ in 0..vseqb.len() {
        jaccard_vec.push(0.);
    }
    // default is invertible hash and then superminhash without any hashing
    let mut pminhasha = ProbMinHash3a::<usize,NoHashHasher>::new(sketch_size, 0);
    let mut wa : FnvIndexMap::<usize,f64> = FnvIndexMap::with_capacity_and_hasher(seqa.size(), FnvBuildHasher::default());
    //
    // generate all kmers include in range arg. dependance upon kmer_size 
    // seqa
    let mut kmergen = KmerSeqIterator::<Kmer32bit>::new(kmer_size, &seqa);
    kmergen.set_range(0, seqa.size()).unwrap();
    loop {
        match kmergen.next() {
            Some(kmer) => {
                let hashval = fhash(&kmer);
                trace!(" kmer in seqa {:?}, hvalval  {:?} ", kmer.get_uncompressed_kmer(), hashval);
                *wa.entry(hashval as usize).or_insert(0.) += 1.;
            },
            None => break,
        }
    }  // end loop
    pminhasha.hash_weigthed_idxmap(&wa);
    let siga = pminhasha.get_signature();
    trace!("siga = {:?}", siga);
    // loop on vseqb to // with rayon
    let comput_closure = | seqb : &Sequence, i:usize | -> (usize,f64) {
        let mut wb : FnvIndexMap::<usize,f64> = FnvIndexMap::with_capacity_and_hasher(seqb.size(), FnvBuildHasher::default());
        let mut kmergen = KmerSeqIterator::<Kmer32bit>::new(kmer_size, &seqb);
        kmergen.set_range(0, seqb.size()).unwrap();
        loop {
            match kmergen.next() {
                Some(kmer) => {
                    let hashval = fhash(&kmer);
                    *wb.entry(hashval as usize).or_insert(0.) += 1.;
                },
                None => break,
            }
        }  // end loop 
        let mut pminhashb = ProbMinHash3a::<usize,NoHashHasher>::new(sketch_size, 0);
        pminhashb.hash_weigthed_idxmap(&wb);
        let sigb = pminhashb.get_signature();
        let jac = compute_probminhash_jaccard(siga, sigb);        
        return (i,jac);
    };
    let jac_with_rank : Vec::<(usize,f64)> = (0..vseqb.len()).into_par_iter().map(|i| comput_closure(&vseqb[i],i)).collect();
    // re-order from jac_with_rank to jaccard_vec as the order of return can be random!!
    for i in 0..jac_with_rank.len() {
        let slot = jac_with_rank[i].0;
        jaccard_vec[slot] = jac_with_rank[i].1;
    }
    return jaccard_vec;
} // end of sketch_seqrange_probminhash3a_kmer32bit



/// Compute jaccard probability index between a sequence and a vector of sequences for Kmer32bit with probminhash3.      
/// It returns a vector of Jaccard probability index.
/// the fhash function is a hash function.  
/// The function is threaded with the Rayon crate.
pub fn jaccard_index_probminhash3_kmer32bit<F>(seqa: &Sequence, vseqb : &Vec<Sequence>, sketch_size: usize, 
                    kmer_size : u8, fhash : F) -> Vec<f64> 
                    where F : Fn(&Kmer32bit) -> u32 + Send + Sync {
    //
    info!("seqsketcher : entering compute_jaccard_index_probminhash3a_kmer32bit");
    // a vector to return results
    let mut jaccard_vec = Vec::<f64>::with_capacity(vseqb.len());
    for _ in 0..vseqb.len() {
        jaccard_vec.push(0.);
    }
    // default is invertible hash and then superminhash without any hashing
    let mut pminhasha = ProbMinHash3::<usize,NoHashHasher>::new(sketch_size, 0);
    let mut wa : FnvIndexMap::<usize,f64> = FnvIndexMap::with_capacity_and_hasher(seqa.size(), FnvBuildHasher::default());
    //
    // generate all kmers include in range arg. dependance upon kmer_size 
    // seqa
    let mut kmergen = KmerSeqIterator::<Kmer32bit>::new(kmer_size, &seqa);
    kmergen.set_range(0, seqa.size()).unwrap();
    loop {
        match kmergen.next() {
            Some(kmer) => {
                let hashval = fhash(&kmer);
                trace!(" kmer in seqa {:?}, hvalval  {:?} ", kmer.get_uncompressed_kmer(), hashval);
                *wa.entry(hashval as usize).or_insert(0.) += 1.;
            },
            None => break,
        }
    }  // end loop
    pminhasha.hash_weigthed_idxmap(&wa);
    let siga = pminhasha.get_signature();
    trace!("siga = {:?}", siga);
    // loop on vseqb to // with rayon
    let comput_closure = | seqb : &Sequence, i:usize | -> (usize,f64) {
        let mut wb : FnvIndexMap::<usize,f64> = FnvIndexMap::with_capacity_and_hasher(seqb.size(), FnvBuildHasher::default());
        let mut kmergen = KmerSeqIterator::<Kmer32bit>::new(kmer_size, &seqb);
        kmergen.set_range(0, seqb.size()).unwrap();
        loop {
            match kmergen.next() {
                Some(kmer) => {
                    let hashval = fhash(&kmer);
                    *wb.entry(hashval as usize).or_insert(0.) += 1.;
                },
                None => break,
            }
        }  // end loop 
        let mut pminhashb = ProbMinHash3::<usize,NoHashHasher>::new(sketch_size, 0);
        pminhashb.hash_weigthed_idxmap(&wb);
        let sigb = pminhashb.get_signature();
        let jac = compute_probminhash_jaccard(siga, sigb);        
        return (i,jac);
    };
    let jac_with_rank : Vec::<(usize,f64)> = (0..vseqb.len()).into_par_iter().map(|i| comput_closure(&vseqb[i],i)).collect();
    // re-order from jac_with_rank to jaccard_vec as the order of return can be random!!
    for i in 0..jac_with_rank.len() {
        let slot = jac_with_rank[i].0;
        jaccard_vec[slot] = jac_with_rank[i].1;
    }
    return jaccard_vec;
} // end of sketch_seqrange_probminhash3_kmer32bit


//============================================== 
//   Dump utilities
//==============================================


const MAGIC_SIG_DUMP : u32 = 0xceabeadd;

// CAVEAT should go to serde/bson

// dumps in an open write buffer a vector os signatures
pub fn dump_signatures_block_u32(signatures : &Vec<Vec<u32>>, out : &mut dyn Write) -> io::Result<()> {
    for i in 0..signatures.len() {
        for j in 0..signatures[i].len() {
            out.write(& signatures[i][j].to_le_bytes()).unwrap();
        }
    }  // end of for i
    //
    return Ok(());
} // end of dump_signatures


/// initialize dump file. Nota we intialize with size of key signature : 4 bytes.  
/// 
/// Format of file is :
/// -  MAGIC_SIG_DUMP as u32
/// -  sig_size 4 or 8 dumped as u32 according to type of signature Vec<u32> or Vec<u64>
/// -  sketch_size  : length of vecteur dumped as u32
/// -  kmer_size    : as u32
pub fn create_signature_dump(dumpfname:&String, kmer_size : u8, sketch_size: usize) -> io::BufWriter<fs::File> {
    let dumpfile_res = OpenOptions::new().write(true).create(true).truncate(true).open(&dumpfname);
    let dumpfile;
    if dumpfile_res.is_ok() {
        dumpfile = dumpfile_res.unwrap();
    } else {
        println!("cannot open {}", dumpfname);
        std::process::exit(1);
    }
    let sig_size : u32 = 4;
    let sketch_size_u32 = sketch_size as u32;
    let kmer_size_u32 = kmer_size as u32;
    let mut sigbuf : io::BufWriter<fs::File> = io::BufWriter::with_capacity(1_000_000_000, dumpfile);
    sigbuf.write(& MAGIC_SIG_DUMP.to_le_bytes()).unwrap();
    sigbuf.write(& sig_size.to_le_bytes()).unwrap();
    sigbuf.write(& sketch_size_u32.to_le_bytes()).unwrap();
    sigbuf.write(& kmer_size_u32.to_le_bytes()).unwrap();
    //
    return sigbuf;
}




#[allow(dead_code)]
/// structure to reload a file consisting of sketch
struct SigSketchFileReader {
    fname:String,
    /// signature size in bytes. 4 for u32, 8 for u64
    sig_size : u8,
    /// the number of sketch by object hashed
    sketch_size:usize,
    /// size of kmers used in sketching.
    kmer_size : u8,
    /// read buffer 
    signature_buf:io::BufReader<fs::File>
}



impl SigSketchFileReader {
    /// initialize the fields fname, sketch_size, kmer_size and allocates signature_buf but signatures will be read by next.
    #[allow(dead_code)]
    pub fn new(fname:&String) -> Result<SigSketchFileReader, String> {
        let dumpfile_res = OpenOptions::new().read(true).open(&fname);
        let dumpfile;
        if dumpfile_res.is_ok() {
            dumpfile = dumpfile_res.unwrap();
        } else {
            println!("cannot open {}", fname);
            std::process::exit(1);
        }
        let mut signature_buf : io::BufReader<fs::File> = io::BufReader::with_capacity(1_000_000_000, dumpfile);
        let mut buf_u32 = [0u8;4];
        let mut io_res;
        // check magic
        io_res = signature_buf.read_exact(&mut buf_u32);
        if io_res.is_err() {
            println!("SigSketchFileReader could no read magic");
            return Err(String::from("SigSketchFileReader could no read magic"));
        }
        let magic = u32::from_le_bytes(buf_u32);
        if magic != MAGIC_SIG_DUMP {
            println!("file {} is not a dump of signature", fname);
            return Err(String::from("file is not a dump of signature"));
        }
        //
        // read sig_size
        //
        io_res = signature_buf.read_exact(&mut buf_u32);
        if io_res.is_err() {
            println!("SigSketchFileReader could no read sketch_size");
            return Err(String::from("SigSketchFileReader could no read sketch_size"));
        }
        let sig_size = u32::from_le_bytes(buf_u32);
        if sig_size != 4 {
            println!("SigSketchFileReader could no read sketch_size");
            return Err(String::from("SigSketchFileReader , sig_size != 4 not yet implemented"));            
        }
        //
        // read sketch_size
        //
        io_res = signature_buf.read_exact(&mut buf_u32);
        if io_res.is_err() {
            println!("SigSketchFileReader could no read sketch_size");
            return Err(String::from("SigSketchFileReader could no read sketch_size"));
        }
        let sketch_size = u32::from_le_bytes(buf_u32);
        trace!("read sketch size {}", sketch_size);
        //
        // check kmer_size
        //
        io_res = signature_buf.read_exact(&mut buf_u32);
        if io_res.is_err() {
            println!("SigSketchFileReader could no read kmer_size");
            return Err(String::from("SigSketchFileReader could no read kmer_size"));
        }
        let kmer_size = u32::from_le_bytes(buf_u32);
        trace!("read kmer_size {}", kmer_size);
        //

        Ok(SigSketchFileReader{fname: fname.clone() , sig_size: sig_size as u8 , sketch_size: sketch_size as usize, kmer_size: kmer_size as u8, signature_buf})
    } // end of new

    /// return kmer_size used sketch dump
    #[allow(dead_code)]
    pub fn get_kmer_size(&self) -> u8 {
        self.kmer_size
    }

    /// returns number of base signature per object
    #[allow(dead_code)]
    pub fn get_signature_length(&self) -> usize {
        self.sketch_size
    }

    /// returns size in bytes of base sketch : 4 or 8
    #[allow(dead_code)]
    pub fn get_signature_size(&self) -> usize {
        self.sig_size as usize
    }
    /// emulates iterator API. Return next object's signature (a Vec<u32> ) if any, None otherwise.
    #[allow(dead_code)]
    pub fn next(&mut self) -> Option<Vec<u32> > {
        let nb_bytes = self.sketch_size * std::mem::size_of::<u32>();
        let mut buf : Vec<u8> = (0..nb_bytes).map(|_| 0u8).collect();

        let io_res = self.signature_buf.read_exact(buf.as_mut_slice());
        //
        if io_res.is_err() {
            // we check that we got EOF or rust ErrorKind::UnexpectedEof
            match io_res.err().unwrap().kind() {
                ErrorKind::UnexpectedEof => return None,
                        _            =>  { 
                                        println!("an unexpected error occurred reading signature buffer");
                                        std::process::exit(1);
                                    }
            }
        }
        else {
            let sig = Vec::<u32>::with_capacity(self.sketch_size);
            return Some(sig);        
        }
    } // end of next
     
} // end of impl SigSketchFile


// ====================================================================================================
//   Some tests
// ====================================================================================================


#[cfg(test)]
mod tests {
    
    use super::*;


#[allow(dead_code)]
    fn log_init_test() {
        let mut builder = env_logger::Builder::from_default_env();
        //    builder.filter_level(LevelFilter::Trace);
        let _ = builder.is_test(true).try_init();
    }


    #[test]
    // This function tests probability jaccard estimates on less than 16 bases 
    fn test_probminhash_kmer_smallb() {
        // initialize test logging
        log_init_test();
        let kmer_size = 5;
        let sketch_size = 50;
        // 80 bases
        let seqstr = String::from("TCAAAGGGAAACATTCAAAATCAGTATGCGCCCGTTCAGTTACGTATTGCTCTCGCTAATGAGATGGGCTGGGTACAGAG");
        let seqabytes = seqstr.as_bytes();
        let seqa = Sequence::new(seqstr.as_bytes(),2);
        //
        let mut vecseqb = Vec::<Sequence>::new();
        // seqb1 has 40-kmer_size in common with seqstr. Jaccard index should be (40-kmer)/(80-kmer_size)
        let seqb1 = Sequence::new(&seqabytes[0..40],2);   // half the length of seqa
        vecseqb.push(seqb1);
        //
        let seqarevcomp = seqa.get_reverse_complement();
        vecseqb.push(seqarevcomp.clone());
        let reverse_str = String::from_utf8(seqarevcomp.decompress()).unwrap();
        println!("\n reverse string : {}", reverse_str);
        // for this hash function seqarevcomp sjaccard signature  should be : [ (40-kmer_size)/(80-kmer_size) , 1.]
        let jac_theo_0 = (40-kmer_size) as f64 / (80-kmer_size) as f64;
        let kmer_revcomp_hash_fn = | kmer : &Kmer32bit | -> u32 {
            let canonical =  kmer.reverse_complement().min(*kmer);
            let hashval = probminhash::invhash::int32_hash(canonical.0);
            hashval
        };  
        // for this hash function (40-kmer)/(80-kmer_size) [ (40-kmer)/(80-kmer_size) , epsil]
        // epsil being the proba there is a common small kmer between seqstr and revers comp which can occur
        let kmer_identity = | kmer : &Kmer32bit | -> u32  {
            let hashval = kmer.0;
            hashval
        };
        let vecsig = jaccard_index_probminhash3a_kmer32bit(&seqa, &vecseqb, sketch_size, kmer_size, kmer_revcomp_hash_fn);
        debug!("vecsig with revcomp hash  {:?}", vecsig);
        assert!(vecsig[0] >= 0.75 * jac_theo_0);
        assert!(vecsig[1] >= 1.);
        // now we try with identity hash
        println!("calling with identity hash");
        let vecsig = jaccard_index_probminhash3a_kmer32bit(&seqa, &vecseqb, sketch_size, kmer_size, kmer_identity);
        debug!("vecsig with identity  {:?}", vecsig);
        assert!(vecsig[0] >= 0.75 * jac_theo_0);
        // get the kmer in intersection if any between seqa and its reverse complement.
        if vecsig[1] > 0. {
            // means we have a kmer in common in seqa and reverse complement of seqa. We check it
            println!("got intersection with reverse complement seq");
            let mut wa : FnvIndexMap::<usize,f64> = FnvIndexMap::with_capacity_and_hasher(seqa.size(), FnvBuildHasher::default());
            let mut pminhasha = ProbMinHash3a::<usize,NoHashHasher>::new(sketch_size, 0);
            // generate all kmers include in range arg. dependance upon kmer_size in seqa 
            let mut kmergen = KmerSeqIterator::<Kmer32bit>::new(kmer_size, &seqa);
            kmergen.set_range(0, seqa.size()).unwrap();
            loop {
                match kmergen.next() {
                    Some(kmer) => {
                        let hashval = kmer_identity(&kmer);
                        debug!(" kmer in seqa {:?}, hvalval  {:?} ", kmer.get_uncompressed_kmer(), hashval);
                        *wa.entry(hashval as usize).or_insert(0.) += 1.;
                    },
                    None => break,
                }
            }  // end loop
            pminhasha.hash_weigthed_idxmap(&wa);
            //
            let mut wb : FnvIndexMap::<usize,f64> = FnvIndexMap::with_capacity_and_hasher(seqarevcomp.size(), FnvBuildHasher::default());
            let mut pminhashb = ProbMinHash3a::<usize,NoHashHasher>::new(sketch_size, 0);
            // generate all kmers include in range arg. dependance upon kmer_size 
            let mut kmergen = KmerSeqIterator::<Kmer32bit>::new(kmer_size, &seqarevcomp);
            kmergen.set_range(0, seqarevcomp.size()).unwrap();
            loop {
                match kmergen.next() {
                    Some(kmer) => {
                        let hashval = kmer_identity(&kmer);
                        trace!(" kmer in seqrevcomp {:?}, hvalval  {:?} ", kmer.get_uncompressed_kmer(), hashval);
                        *wb.entry(hashval as usize).or_insert(0.) += 1.;
                    },
                    None => break,
                }
            }  // end loop    
            pminhashb.hash_weigthed_idxmap(&wb);
            let (jac, common) = probminhash_get_jaccard_objects(pminhasha.get_signature(), pminhashb.get_signature());
            debug!("jac for common objects = {}", jac);
            if jac > 0. {
                // with kmer size = 5 we have ACGTA and TACGT that are common!
                trace!("common kemrs {:?}", common.unwrap());
            }
        }  // end search of intersecting kmers
        //
        assert!(vecsig[1] <= 0.1);

    }  // end of test_probminhash_kmer_smallb


    #[test]
    fn test_probminhash_kmer_16b32bit_serial() {
        // initialize test logging
        log_init_test();
        // 80 bases
        let kmer_size = 16;
        let seqstr = String::from("TCAAAGGGAAACATTCAAAATCAGTATGCGCCCGTTCAGTTACGTATTGCTCTCGCTAATGAGATGGGCTGGGTACAGAG");
        let seqabytes = seqstr.as_bytes();
        let seqa = Sequence::new(seqstr.as_bytes(),2);
        //
        let mut vecseqb = Vec::<Sequence>::new();
        // seqb1 has 40-kmer_size in common with seqstr. Jaccard index should be (40-kmer)/(80-kmer_size)
        let seqb1 = Sequence::new(&seqabytes[0..40],2);   // half the length of seqa
        vecseqb.push(seqb1);
        //
        let seqarevcomp = seqa.get_reverse_complement();
        vecseqb.push(seqarevcomp.clone());
        let reverse_str = String::from_utf8(seqarevcomp.decompress()).unwrap();
        println!("\n reverse string : {}", reverse_str);
        // for this hash function seqarevcomp sjaccard signature  should be : [ (40-kmer_size)/(80-kmer_size) , 1.]
        let jac_theo_0 = (40-kmer_size) as f64 / (80-kmer_size) as f64;
        let kmer_revcomp_hash_fn = | kmer : &Kmer16b32bit | -> u32 {
            let canonical =  kmer.reverse_complement().min(*kmer);
            let hashval = probminhash::invhash::int32_hash(canonical.0);
            hashval
        };  
        // for this hash function (40-kmer)/(80-kmer_size) [ (40-kmer)/(80-kmer_size) , epsil]
        // epsil being the proba there is a common small kmer between seqstr and revers comp which can occur
        let kmer_identity = | kmer : &Kmer16b32bit | -> u32  {
            let hashval = kmer.0;
            hashval
        };
        let vec_0 = jaccard_index_probminhash3a_kmer16b32bit(&seqa, &vec![vecseqb[0].clone()], 50, kmer_revcomp_hash_fn);
        let vec_1 = jaccard_index_probminhash3a_kmer16b32bit(&seqa, &vec![vecseqb[1].clone()], 50, kmer_revcomp_hash_fn);
        let mut vecsig = Vec::<f64>::with_capacity(2);
        vecsig.push(vec_0[0]);
        vecsig.push(vec_1[0]);
        debug!("vecsig with revcomp hash  {:?}", vecsig);
        assert!(vecsig[0] >= 0.75 * jac_theo_0);
        assert!(vecsig[1] >= 1.);
        //
        println!("calling with identity hash");
        let vecsig = jaccard_index_probminhash3a_kmer16b32bit(&seqa, &vecseqb, 50, kmer_identity);
        debug!("vecsig with identity  {:?}", vecsig);
        assert!(vecsig[0] >= 0.75 * jac_theo_0);
        assert!(vecsig[1] <= 0.1);
    }  // end of test_probminhash_kmer_16b32bit


#[test]
// This test checks for parallel computation of signature with the same sequences as  test_probminhash_kmer_16b32bit
   fn test_probminhash_kmer_16b32bit_par() {
        log_init_test();
        // 80 bases
        let kmer_size = 16;
        let seqstr = String::from("TCAAAGGGAAACATTCAAAATCAGTATGCGCCCGTTCAGTTACGTATTGCTCTCGCTAATGAGATGGGCTGGGTACAGAG");
        let seqabytes = seqstr.as_bytes();
        let seqa = Sequence::new(seqstr.as_bytes(),2);
        //
        let mut vecseqb = Vec::<Sequence>::new();
        // seqb1 has 40-kmer_size in common with seqstr. Jaccard index should be (40-kmer)/(80-kmer_size)
        let seqb1 = Sequence::new(&seqabytes[0..40],2);   // half the length of seqa
        vecseqb.push(seqb1);
        //
        let seqarevcomp = seqa.get_reverse_complement();
        vecseqb.push(seqarevcomp);
        let kmer_revcomp_hash_fn = | kmer : &Kmer16b32bit | -> u32 {
            let canonical =  kmer.reverse_complement().min(*kmer);
            let hashval = probminhash::invhash::int32_hash(canonical.0);
            hashval
        }; 
        let vec_jac = jaccard_index_probminhash3a_kmer16b32bit(&seqa, &vecseqb, kmer_size, kmer_revcomp_hash_fn);
        debug!("vecsig with revcomp hash  {:?}", vec_jac);
        let jac_theo_0 = (40-kmer_size) as f64 / (80-kmer_size) as f64;
        assert!(vec_jac[0] >= 0.75 * jac_theo_0);
        assert!(vec_jac[1] >= 1.);
    }



//
//========================================================================================================
//   io tests
//==========================================================================================================


// This tests reload of a signature dump (if a test file is present) 

#[test]
fn test_reload_sketch_file() {
    log_init_test();
    //
    let fname = String::from("/home.1/jpboth/Rust/kmerutils/Runs/umpsigk8s200");
    let sketch_reader_res = SigSketchFileReader::new(&fname);
    // check result with a as_ref to avoid consuming value
    let sketch_reader_res_ref = sketch_reader_res.as_ref();
    if let Some(msg) = sketch_reader_res_ref.err()  {
        println!("test_reload_sketch_file, error with file : {} {} ", fname, msg);
        return;
    }
    // get a mut on result
    let mut sketch_reader_ref = sketch_reader_res.ok().unwrap();
    //
    println!("kmer size : {}", sketch_reader_ref.get_kmer_size());
    println!("sig length : {}", sketch_reader_ref.get_signature_length());
    println!("sig size : {}", sketch_reader_ref.get_signature_size());
    //
    let mut nbread = 0;
    while let Some(_sig) = sketch_reader_ref.next()  {
        nbread += 1;
        if nbread % 100000 == 0 {
            println!("loaded nb sig : {}", nbread);
        }
    }
    println!("loaded nb sig : {}", nbread);
} // end test_reload_sketch_file





}  // end of mod test

