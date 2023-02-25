//! This module provides sequence signature computation and Jaccard probability index using the probminhash crate.  
//! The Jaccard probability index is a Jaccard index between sequences taking into account multiplicity of kmers. 
//!    
//! For long (many Gbytes) sequences it may be useful to split sequences into blocks, see module seqbloocksketch.
//! 
//! The kmers of a given size are generated for each sequence, kmers lists are hashed by the probminhash algorithm 
//! and a jaccard weighted index between sequences is computed. See [Probminhash](https://crates.io/crates/probminhash)
//! 
//! It is also possible to ask for the common Kmers found in the signature of 2 sequences
//! 
//! 
//! 



use log::*;

use std::io;
use std::io::{Write,Read,ErrorKind, BufReader, BufWriter };

use std::fs;
use std::fs::OpenOptions;
use std::fmt::{Debug};
use std::path::{Path, PathBuf};

use std::hash::{BuildHasherDefault, Hasher, Hash};

use serde::{Deserialize, Serialize};
use serde_json::{to_writer};

use indexmap::{IndexMap};
use fnv::{FnvHashMap, FnvBuildHasher};

use num;

use crate::nohasher::*;

use crate::base::{kmer::*, kmergenerator::*, kmergenerator::KmerSeqIteratorT};

use rayon::prelude::*;


use probminhash::{probminhasher::*, superminhasher::SuperMinHash};
use probminhash::jaccard::compute_probminhash_jaccard;


// We need a guess to allocate HashMap used with Kmer Generation
// for very long sequence we must avoid nb_kmer to sequence length! Find a  good heuristic
fn get_nbkmer_guess(seq : &Sequence) -> usize {
    let nb = 100_000_000 * (1usize + seq.size().ilog2() as usize);
    let nb_kmer = seq.size().min(nb);
    return nb_kmer;
} // end of get_nbkmer_guess



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


/// This structure describes the kmer size used in computing sketches and the number of sketch we want.
#[derive(Serialize,Deserialize,Copy,Clone)]
pub struct SeqSketcher {
    kmer_size : usize,
    sketch_size : usize
}  // end of SeqSketcher


impl SeqSketcher {
    /// 
    pub fn new(kmer_size: usize, sketch_size : usize) -> Self {
        SeqSketcher{kmer_size, sketch_size}
    }

    /// returns kmer size
    pub fn get_kmer_size(&self) -> usize {
        self.kmer_size
    }

    /// return sketch size
    pub fn get_sketch_size(&self) -> usize {
        self.sketch_size
    }  
    
    /// serialized dump
    pub fn dump_json(&self, filename : &String) -> Result<(), String> {
        //
        let filepath = PathBuf::from(filename.clone());
        //
        log::info!("dumping sketching parameters in json file : {}", filename);
        //
        let fileres = OpenOptions::new().write(true).create(true).truncate(true).open(&filepath);
        if fileres.is_err() {
            log::error!("SeqSketcher dump : dump could not open file {:?}", filepath.as_os_str());
            println!("SeqSketcher dump: could not open file {:?}", filepath.as_os_str());
            return Err("SeqSketcher dump failed".to_string());
        }
        // 
        let mut writer = BufWriter::new(fileres.unwrap());
        let _ = to_writer(&mut writer, &self).unwrap();
        //
        Ok(())
    } // end of dump


    /// reload from a json dump
    pub fn reload_json(dirpath : &Path) -> Result<SeqSketcher, String> {
        log::info!("in reload_json");
        //
        let filepath = dirpath.join("sketchparams_dump.json");
        let fileres = OpenOptions::new().read(true).open(&filepath);
        if fileres.is_err() {
            log::error!("Sketcher reload_json : reload could not open file {:?}", filepath.as_os_str());
            println!("Sketcher reload_json: could not open file {:?}", filepath.as_os_str());
            return Err("Sketcher reload_json could not open file".to_string());            
        }
        //
        let loadfile = fileres.unwrap();
        let reader = BufReader::new(loadfile);
        let sketch_params:SeqSketcher = serde_json::from_reader(reader).unwrap();
        //
        log::info!("SeqSketcher reload, kmer_size : {}, sketch_size : {}", 
            sketch_params.get_kmer_size(), sketch_params.get_sketch_size());     
        //
        Ok(sketch_params)
    } // end of reload_json





    /// This function computes and return signatures of a vector of sequences by generating kmers of size kmer_size. 
    /// The sketch is done with probminhash3a algorithm. See the crate [probminhash](https://crates.io/crates/probminhash)  
    /// It is a generic implementation of probminhash3a  against our standard compressed Kmer types.  
    /// Kmer::Val is the base type u32, u64 on which compressed kmer representations relies.  
    ///    
    /// fhash is any hash function, but usually it is identity, invhash on kmer or on min of kmer and reverse complement.
    /// These are the hash functions that make possible to get back to the original kmers (or at least partially in the case using the min)..
    /// 
    /// The argument type of the hashing function F specify the type of Kmer to generate along the sequence.  
    pub fn sketch_probminhash3a_compressedkmer<'b, Kmer : CompressedKmerT + KmerBuilder<Kmer>, F>(&self, vseq : &'b Vec<&Sequence>, fhash : F) -> Vec<Vec<Kmer::Val> >
        where F : Fn(&Kmer) -> Kmer::Val + Send + Sync,
              Kmer::Val : num::PrimInt + Send + Sync + Debug,
              KmerGenerator<Kmer> :  KmerGenerationPattern<Kmer> {
        //
        log::debug!("entering sketch_probminhash3a_compressedkmer");
        //
        let comput_closure = | seqb : &Sequence, i:usize | -> (usize,Vec<Kmer::Val>) {
            // if we get very large sequence (many Gb length) we must be cautious on size of hashmap; i.e about number of different kmers!!! 
            let nb_kmer = get_nbkmer_guess(seqb);
            let mut wb : FnvHashMap::<Kmer::Val,u64> = FnvHashMap::with_capacity_and_hasher(nb_kmer, FnvBuildHasher::default());
            let mut kmergen = KmerSeqIterator::<Kmer>::new(self.kmer_size as u8, &seqb);
            kmergen.set_range(0, seqb.size()).unwrap();
            loop {
                match kmergen.next() {
                    Some(kmer) => {
                        let hashval = fhash(&kmer);
                        *wb.entry(hashval).or_insert(0) += 1;
                    },
                    None => break,
                }
            }  // end loop 
            let mut pminhashb = ProbMinHash3a::<Kmer::Val,NoHashHasher>::new(self.sketch_size, 
                <Kmer::Val>::default());
            pminhashb.hash_weigthed_hashmap(&wb);
            let sigb = pminhashb.get_signature();
            // get back from usize to Kmer32bit ?. If fhash is inversible possible, else NO.
            return (i,sigb.clone());
        };
        //
        let sig_with_rank : Vec::<(usize,Vec<Kmer::Val>)> = (0..vseq.len()).into_par_iter().map(|i| comput_closure(vseq[i],i)).collect();
        // re-order from jac_with_rank to jaccard_vec as the order of return can be random!!
        let mut jaccard_vec = Vec::<Vec<Kmer::Val>>::with_capacity(vseq.len());
        for _ in 0..vseq.len() {
            jaccard_vec.push(Vec::new());
        }
        // CAVEAT , boxing would avoid the clone?
        for i in 0..sig_with_rank.len() {
            let slot = sig_with_rank[i].0;
            jaccard_vec[slot] = sig_with_rank[i].1.clone();
        }
        jaccard_vec
    }  // end of sketch_probminhash3a_compressedkmer



    //   Probminhash3
    //  ==============

    /// This function computes and return signatures of a vector of sequences by generating kmers of size kmer_size. 
    /// The sketch is done with probminhash3 algorithm.   
    /// The size of signature of each sequence is sketch_size.  
    /// fhash is any hash function, but usually it is identity, invhash on kmer or on min of kmer and reverse complement.  
    /// These are the hash function that make possible to get back to the original kmers (or at least partially in the case using the min).  
    /// 
    /// The argument type of the hashing function F specify the type of Kmer to generate along the sequence.  

    pub fn sketch_probminhash3_kmer32bit<F>(&self, vseq : &Vec<&Sequence>, fhash : F) -> Vec<Vec<u32> >
        where F : Fn(&Kmer32bit) -> u32 + Send + Sync {
        //
        let comput_closure = | seqb : &Sequence, i:usize | -> (usize,Vec<u32>) {
            // if we get very large sequence (many Gb length) we must be cautious on size of hashmap; i.e about number of different kmers!!! 
            let nb_kmer = get_nbkmer_guess(seqb);
            let mut wb : FnvHashMap::<u32,f64> = FnvHashMap::with_capacity_and_hasher(nb_kmer, FnvBuildHasher::default());
            let mut kmergen = KmerSeqIterator::<Kmer32bit>::new(self.kmer_size as u8, &seqb);
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
            let mut pminhashb = ProbMinHash3a::<u32,NoHashHasher>::new(self.sketch_size, 0);
            pminhashb.hash_weigthed_hashmap(&wb);
            let sigb = pminhashb.get_signature();
            // get back from usize to Kmer32bit ?. If fhash is inversible possible, else NO.
            return (i,sigb.clone());
        };
        //
        let sig_with_rank : Vec::<(usize,Vec<u32>)> = (0..vseq.len()).into_par_iter().map(|i| comput_closure(vseq[i],i)).collect();
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



    //  Superminhash

    /// a generic implementation of superminhash  against our standard compressed Kmer types.  
    /// Kmer::Val is the base type u32, u64 on which compressed kmer representations relies.
    /// F is a hash function returning morally a u32, usize or u64.  
    /// The argument type of the hashing function F specify the type of Kmer to generate along the sequence.  
    pub fn sketch_superminhash_compressedkmer<'b, Kmer : CompressedKmerT + KmerBuilder<Kmer>, F, T>(&self, vseq : &'b Vec<&Sequence>, fhash : F) -> Vec<Vec<f64> >
        where F : Fn(&Kmer) -> T + Send + Sync,
              T : std::hash::Hash, 
              Kmer::Val : num::PrimInt + Send + Sync + Debug,
              KmerGenerator<Kmer> :  KmerGenerationPattern<Kmer> {
        //
        log::debug!("entering sketch_superminhash_compressedkmer");
        //
        let comput_closure = | seqb : &Sequence, i:usize | -> (usize,Vec<f64>) {
            //
            log::debug!(" in sketch_superminhash_compressedkmer, closure");
            //
            let bh = BuildHasherDefault::<fnv::FnvHasher>::default();
            let mut sminhash : SuperMinHash<T, fnv::FnvHasher>= SuperMinHash::new(self.sketch_size, &bh);

            let mut kmergen = KmerSeqIterator::<Kmer>::new(self.kmer_size as u8, &seqb);
            kmergen.set_range(0, seqb.size()).unwrap();
            loop {
                match kmergen.next() {
                    Some(kmer) => {
                        let hashval = fhash(&kmer);
                        if sminhash.sketch(&hashval).is_err() {
                            log::error!("could not hash kmer : {:?}", kmer.get_uncompressed_kmer());
                            std::panic!("could not hash kmer : {:?}", kmer.get_uncompressed_kmer());
                        }
                    },
                    None => break,
                }
            }  // end loop 
            let sigb = sminhash.get_hsketch();
            // get back from usize to Kmer32bit ?. If fhash is inversible possible, else NO.
            return (i,sigb.clone());
        };
        //
        let sig_with_rank : Vec::<(usize,Vec<f64>)> = (0..vseq.len()).into_par_iter().map(|i| comput_closure(vseq[i],i)).collect();
        // re-order from jac_with_rank to jaccard_vec as the order of return can be random!!
        let mut jaccard_vec = Vec::<Vec<f64>>::with_capacity(vseq.len());
        for _ in 0..vseq.len() {
            jaccard_vec.push(Vec::new());
        }
        // CAVEAT , boxing would avoid the clone?
        for i in 0..sig_with_rank.len() {
            let slot = sig_with_rank[i].0;
            jaccard_vec[slot] = sig_with_rank[i].1.clone();
        }
        jaccard_vec
    } // end of sketch_superminhash_compressedkmer




    /// initialize dump file. Nota we intialize with size of key signature : 4 bytes.  
    /// 
    /// Format of file is :
    /// -  MAGIC_SIG_DUMP as u32
    /// -  sig_size 4 or 8 dumped as u32 according to type of signature Vec\<u32\> or Vec\<u64\>
    /// -  sketch_size  : length of vecteur dumped as u32
    /// -  kmer_size    : as u32
    /// 
    pub fn create_signature_dump(&self, dumpfname:&String) -> io::BufWriter<fs::File> {
        let dumpfile_res = OpenOptions::new().write(true).create(true).truncate(true).open(&dumpfname);
        let dumpfile;
        if dumpfile_res.is_ok() {
            dumpfile = dumpfile_res.unwrap();
        } else {
            println!("cannot open {}", dumpfname);
            std::process::exit(1);
        }
        let sig_size : u32 = 4;
        let sketch_size_u32 = self.sketch_size as u32;
        let kmer_size_u32 = self.kmer_size as u32;
        let mut sigbuf : io::BufWriter<fs::File> = io::BufWriter::with_capacity(1_000_000_000, dumpfile);
        sigbuf.write(& MAGIC_SIG_DUMP.to_le_bytes()).unwrap();
        sigbuf.write(& sig_size.to_le_bytes()).unwrap();
        sigbuf.write(& sketch_size_u32.to_le_bytes()).unwrap();
        sigbuf.write(& kmer_size_u32.to_le_bytes()).unwrap();
        //
        return sigbuf;
    }  // end of create_signature_dump



}  // end of impl SeqSketcher




/// Compute jaccard probability index between a sequence and a vector of sequences for all CompressedKmer  with probminhash3a.      
/// It returns a vector of Jaccard probability index.
/// the fhash function is a hash function.  
/// The function is threaded with the Rayon crate.
pub fn jaccard_index_probminhash3a< Kmer : CompressedKmerT + KmerBuilder<Kmer>, F>(seqa: &Sequence, vseqb : &Vec<Sequence>, sketch_size: usize, kmer_size : u8, fhash : F) -> Vec<f64> 
    where F : Fn(&Kmer) -> Kmer::Val + Send + Sync,
              Kmer::Val : num::PrimInt + Send + Sync + Debug,
              KmerGenerator<Kmer> :  KmerGenerationPattern<Kmer>{
    //
    debug!("seqsketcher : entering compute_jaccard_index_probminhash3a");
    // a vector to return results
    let mut jaccard_vec = Vec::<f64>::with_capacity(vseqb.len());
    for _ in 0..vseqb.len() {
        jaccard_vec.push(0.);
    }
    // default is invertible hash and then superminhash without any hashing
    let mut pminhasha = ProbMinHash3a::<<Kmer as CompressedKmerT>::Val,NoHashHasher>::new(sketch_size, Kmer::Val::default());
    // if we get very large sequence (many Gb length) we must be cautious on size of hashmap; i.e about number of different kmers!!! 
    let nb_kmer = get_nbkmer_guess(seqa);
    let mut wa : FnvHashMap::<Kmer::Val,u64> = FnvHashMap::with_capacity_and_hasher(nb_kmer, FnvBuildHasher::default());
    //
    // generate all kmers include in range arg. dependance upon kmer_size 
    // seqa
    let mut kmergen = KmerSeqIterator::<Kmer>::new(kmer_size, &seqa);
    kmergen.set_range(0, seqa.size()).unwrap();
    loop {
        match kmergen.next() {
            Some(kmer) => {
                let hashval = fhash(&kmer);
                trace!(" kmer in seqa {:?}, hvalval  {:?} ", kmer.get_uncompressed_kmer(), hashval);
                *wa.entry(hashval).or_insert(0) += 1;
            },
            None => break,
        }
    }  // end loop
    pminhasha.hash_weigthed_hashmap(&wa);
    let siga = pminhasha.get_signature();
    trace!("siga = {:?}", siga);
    // loop on vseqb to // with rayon
    let comput_closure = | seqb : &Sequence, i:usize | -> (usize,f64) {
        // if we get very large sequence (many Gb length) we must be cautious on size of hashmap; i.e about number of different kmers!!! 
        let nb_kmer = get_nbkmer_guess(seqb);
        let mut wb : FnvHashMap::<Kmer::Val, u64> = FnvHashMap::with_capacity_and_hasher(nb_kmer, FnvBuildHasher::default());
        let mut kmergen = KmerSeqIterator::<Kmer>::new(kmer_size, &seqb);
        kmergen.set_range(0, seqb.size()).unwrap();
        loop {
            match kmergen.next() {
            Some(kmer) => {
                let hashval = fhash(&kmer);
                *wb.entry(hashval).or_insert(0) += 1;
            },
            None => break,
            }
        }  // end loop 
        let mut pminhashb = ProbMinHash3a::<Kmer::Val,NoHashHasher>::new(sketch_size, Kmer::Val::default());
        pminhashb.hash_weigthed_hashmap(&wb);
        let sigb = pminhashb.get_signature();
        let jac = compute_probminhash_jaccard(siga, sigb);        
        return (i,jac);
    };
    //
    let jac_with_rank : Vec::<(usize,f64)> = (0..vseqb.len()).into_par_iter().map(|i| comput_closure(&vseqb[i],i)).collect();
    // re-order from jac_with_rank to jaccard_vec as the order of return can be random!!
    for i in 0..jac_with_rank.len() {
        let slot = jac_with_rank[i].0;
        jaccard_vec[slot] = jac_with_rank[i].1;
    }
    return jaccard_vec;
} // end of sketch_seqrange_probminhash3a



/// Compute jaccard probability index between a sequence and a vector of sequences for Kmer32bit with probminhash3.      
/// It returns a vector of Jaccard probability index.
/// the fhash function is a hash function.  
/// The function is threaded with the Rayon crate.
pub fn jaccard_index_probminhash3_kmer32bit<F>(seqa: &Sequence, vseqb : &Vec<Sequence>, sketch_size: usize, 
                    kmer_size : u8, fhash : F) -> Vec<f64> 
                    where F : Fn(&Kmer32bit) -> u32 + Send + Sync {
    //
    debug!("seqsketcher : entering compute_jaccard_index_probminhash3a_kmer32bit");
    // a vector to return results
    let mut jaccard_vec = Vec::<f64>::with_capacity(vseqb.len());
    for _ in 0..vseqb.len() {
        jaccard_vec.push(0.);
    }
    // default is invertible hash and then superminhash without any hashing
    let mut pminhasha = ProbMinHash3::<usize,NoHashHasher>::new(sketch_size, 0);
    // if we get very large sequence (many Gb length) we must be cautious on size of hashmap; i.e about number of different kmers!!! 
    let nb_kmer = get_nbkmer_guess(seqa);
    let mut wa : FnvHashMap::<usize,f64> = FnvHashMap::with_capacity_and_hasher(nb_kmer, FnvBuildHasher::default());
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
    pminhasha.hash_weigthed_hashmap(&wa);
    let siga = pminhasha.get_signature();
    trace!("siga = {:?}", siga);
    // loop on vseqb to // with rayon
    let comput_closure = | seqb : &Sequence, i:usize | -> (usize,f64) {
        // if we get very large sequence (many Gb length) we must be cautious on size of hashmap; i.e about number of different kmers!!! 
        let nb_kmer = get_nbkmer_guess(seqb);
        let mut wb : FnvHashMap::<usize,f64> = FnvHashMap::with_capacity_and_hasher(nb_kmer, FnvBuildHasher::default());
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
        pminhashb.hash_weigthed_hashmap(&wb);
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

// dumps in an open write buffer a vector of signatures
pub fn dump_signatures_block_u32(signatures : &Vec<Vec<u32>>, out : &mut dyn Write) -> io::Result<()> {
    for i in 0..signatures.len() {
        for j in 0..signatures[i].len() {
            out.write(& signatures[i][j].to_le_bytes()).unwrap();
        }
    }  // end of for i
    //
    return Ok(());
} // end of dump_signatures




/// structure to reload a file consisting of sketch
pub struct SigSketchFileReader {
    _fname:String,
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
    pub fn new(fname:&String) -> Result<SigSketchFileReader, String> {
        let dumpfile_res = OpenOptions::new().read(true).open(&fname);
        let dumpfile;
        if dumpfile_res.is_ok() {
            dumpfile = dumpfile_res.unwrap();
        } else {
            println!("cannot open {}", fname);
            return Err(String::from("SigSketchFileReader : could not open dumpfile"))
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

        Ok(SigSketchFileReader{_fname: fname.clone() , sig_size: sig_size as u8 , sketch_size: sketch_size as usize, kmer_size: kmer_size as u8, signature_buf})
    } // end of new

    /// return kmer_size used sketch dump
    pub fn get_kmer_size(&self) -> u8 {
        self.kmer_size
    }

    /// returns number of base signature per object
    pub fn get_signature_length(&self) -> usize {
        self.sketch_size
    }

    /// returns size in bytes of base sketch : 4 or 8
    pub fn get_signature_size(&self) -> usize {
        self.sig_size as usize
    }
    /// emulates iterator API. Return next object's signature (a Vec\<u32\> ) if any, None otherwise.
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
     
} // end of impl SigSketchFileReader


// ====================================================================================================
//   Some tests
// ====================================================================================================


#[cfg(test)]
mod tests {
    
    use super::*;
//    use probminhash::superminhasher::compute_superminhash_jaccard;

// we define compute_superminhash_jaccard to avoid bumping version of probminhash now!
// TODO use probminhash::superminhasher::compute_superminhash_jaccard ASAP probminhash gets to 0.1.7
#[inline]
    fn compute_superminhash_jaccard(hsketch: &Vec<f64>  , other_sketch: &Vec<f64>)  -> Result<f64, ()>  {
        return probminhash::superminhasher::get_jaccard_index_estimate(hsketch, other_sketch);
    }


    fn log_init_test() {
        let mut builder = env_logger::Builder::from_default_env();
        //    builder.filter_level(LevelFilter::Trace);
        let _ = builder.is_test(true).try_init();
    }


    #[test]
    // This function tests probability jaccard estimates on kmers of size less than 16 bases 
    fn test_pminhasha_kmer_smallb() {
        // initialize test logging
        log_init_test();
        log::info!("test_probminhasha_kmer_smallb");
        //
        let kmer_size = 5;
        let sketch_size = 4000;
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
        log::info!("jac_theo_0 : {:.3e}", jac_theo_0);
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
        let vecsig = jaccard_index_probminhash3a(&seqa, &vecseqb, sketch_size, kmer_size, kmer_revcomp_hash_fn);
        log::info!("vecsig with revcomp hash  {:?}", vecsig);
        assert!(vecsig[0] >= 0.75 * jac_theo_0);
        assert!(vecsig[1] >= 1.);
        // now we try with identity hash
        println!("calling with identity hash");
        let vecsig = jaccard_index_probminhash3a(&seqa, &vecseqb, sketch_size, kmer_size, kmer_identity);
        debug!("vecsig with identity  {:?}", vecsig);
        assert!(vecsig[0] >= 0.75 * jac_theo_0);
        // get the kmer in intersection if any between seqa and its reverse complement.
        if vecsig[1] > 0. {
            // means we have a kmer in common in seqa and reverse complement of seqa. We check it
            println!("got intersection with reverse complement seq");
            let mut wa : FnvHashMap::<u32,f64> = FnvHashMap::with_capacity_and_hasher(seqa.size(), FnvBuildHasher::default());
            let mut pminhasha = ProbMinHash3a::<u32,NoHashHasher>::new(sketch_size, 0);
            // generate all kmers include in range arg. dependance upon kmer_size in seqa 
            let mut kmergen = KmerSeqIterator::<Kmer32bit>::new(kmer_size, &seqa);
            kmergen.set_range(0, seqa.size()).unwrap();
            loop {
                match kmergen.next() {
                    Some(kmer) => {
                        let hashval = kmer_identity(&kmer);
                        debug!(" kmer in seqa {:?}, hvalval  {:?} ", String::from_utf8(kmer.get_uncompressed_kmer()).unwrap(), hashval);
                        *wa.entry(hashval).or_insert(0.) += 1.;
                    },
                    None => break,
                }
            }  // end loop
            pminhasha.hash_weigthed_hashmap(&wa);
            //
            let mut wb : FnvHashMap::<u32,f64> = FnvHashMap::with_capacity_and_hasher(seqarevcomp.size(), FnvBuildHasher::default());
            let mut pminhashb = ProbMinHash3a::<u32,NoHashHasher>::new(sketch_size, 0);
            // generate all kmers include in range arg. dependance upon kmer_size 
            let mut kmergen = KmerSeqIterator::<Kmer32bit>::new(kmer_size, &seqarevcomp);
            kmergen.set_range(0, seqarevcomp.size()).unwrap();
            loop {
                match kmergen.next() {
                    Some(kmer) => {
                        let hashval = kmer_identity(&kmer);
                        trace!(" kmer in seqrevcomp {:?}, hvalval  {:?} ", kmer.get_uncompressed_kmer(), hashval);
                        *wb.entry(hashval).or_insert(0.) += 1.;
                    },
                    None => break,
                }
            }  // end loop    
            pminhashb.hash_weigthed_hashmap(&wb);
            let (jac, common) = probminhash_get_jaccard_objects(pminhasha.get_signature(), pminhashb.get_signature());
            debug!("jac for common objects = {}", jac);
            if jac > 0. {
                // with kmer size = 5 we have ACGTA and TACGT that are common!
                debug!("common kemrs {:?}", common.unwrap());
            }
        }  // end search of intersecting kmers
        //
        assert!(vecsig[1] <= 0.1);

    }  // end of test_probminhasha_kmer_smallb


    #[test]
    fn test_pminhasha_k16b32bit_serial() {
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
        debug!("\n reverse string : {}", reverse_str);
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
        let vec_0 = jaccard_index_probminhash3a(&seqa, &vec![vecseqb[0].clone()], 50, 16,  kmer_revcomp_hash_fn);
        let vec_1 = jaccard_index_probminhash3a(&seqa, &vec![vecseqb[1].clone()], 50, 16, kmer_revcomp_hash_fn);
        let mut vecsig = Vec::<f64>::with_capacity(2);
        vecsig.push(vec_0[0]);
        vecsig.push(vec_1[0]);
        info!("vecsig with revcomp hash  {:?}", vecsig);
        assert!(vecsig[0] >= 0.75 * jac_theo_0);
        assert!(vecsig[1] >= 1.);
        //
        info!("calling with identity hash");
        let vecsig = jaccard_index_probminhash3a(&seqa, &vecseqb, 50, 16, kmer_identity);
        info!("vecsig with identity  {:?}", vecsig);
        assert!(vecsig[0] >= 0.75 * jac_theo_0);
        assert!(vecsig[1] <= 0.1);
    }  // end of test_probminhash_kmer_16b32bit


#[test]
// This test checks for parallel computation of signature with the same sequences as  test_probminhash_kmer_16b32bit
   fn test_pminhash_kmer64bit_serial() {
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
        let kmer_revcomp_hash_fn = | kmer : &Kmer64bit | -> u64 {
            let canonical =  kmer.reverse_complement().min(*kmer);
            let hashval = probminhash::invhash::int64_hash(canonical.0);
            hashval
        }; 
        let vec_jac = jaccard_index_probminhash3a(&seqa, &vecseqb, 50, 16, kmer_revcomp_hash_fn);
        let jac_theo_0 = (40-kmer_size) as f64 / (80-kmer_size) as f64;
        info!("vecsig with revcomp hash  {:?} jaccard theo : {:.3e}", vec_jac, jac_theo_0);
        assert!(vec_jac[0] >= 0.75 * jac_theo_0);
        assert!(vec_jac[1] >= 1.);
    }



    #[test]
    fn test_superminhash_kmer_16b32bit_serial() {
        // initialize test logging
        log_init_test();
        // 80 bases
        let kmer_size = 16;
        let sketch_size = 100;
        //
        let mut vecseq = Vec::<&Sequence>::new();
        //
        let seqstr = String::from("TCAAAGGGAAACATTCAAAATCAGTATGCGCCCGTTCAGTTACGTATTGCTCTCGCTAATGAGATGGGCTGGGTACAGAG");
        let seqabytes = seqstr.as_bytes();
        let seqa = Sequence::new(seqstr.as_bytes(),2);
        vecseq.push(&seqa);
        //
        // seqb1 has 40-kmer_size in common with seqa. Jaccard index seqa should be (40-kmer)/(80-kmer_size)
        let seqb1 = Sequence::new(&seqabytes[0..40],2);   // half the length of seqa
        vecseq.push(&seqb1);
        //
        let seqarevcomp = seqa.get_reverse_complement();
        vecseq.push(&seqarevcomp);
        let reverse_str = String::from_utf8(seqarevcomp.decompress()).unwrap();
        log::debug!("\n reverse string : {}", reverse_str);
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
        // do a superminhash sketching of sequence with kmer_revcomp_hash_fn
        let sketcher = SeqSketcher::new(kmer_size, sketch_size);
        let sig_vec = sketcher.sketch_superminhash_compressedkmer(&vecseq,kmer_revcomp_hash_fn);
        // now we can compute jaccard index between sig_vec[0] and the 2 others i.e sig_vec[1] and sig_vec[2]
        let d_01 = compute_superminhash_jaccard(&sig_vec[0], &sig_vec[1]).unwrap();
        let d_02 = compute_superminhash_jaccard(&sig_vec[0], &sig_vec[2]).unwrap();
        debug!("seqa with revcomp hash  {:?}", sig_vec[0]);
        debug!("seqb with revcomp hash  {:?}", sig_vec[1]);
        debug!("seq rev comp with revcomp hash  {:?}", sig_vec[2]);
        debug!("ditances  with revcomp hash  {:.3e}  {:.3e}", d_01, d_02);
        assert!(d_01 >= 0.75 * jac_theo_0);
        assert!(d_02 >= 1.);
        //
        // do a superminhash sketching of sequence with kmer_identity
        println!("calling with identity hash");
        let sig_vec = sketcher.sketch_superminhash_compressedkmer(&vecseq,kmer_identity);
        let d_01 = compute_superminhash_jaccard(&sig_vec[0], &sig_vec[1]).unwrap();
        let d_02 = compute_superminhash_jaccard(&sig_vec[0], &sig_vec[2]).unwrap();
        debug!("vecsig with identity  {:?}", sig_vec);
        info!("ditances  with revcomp hash  {:.3e}  {:.3e}", d_01, d_02);
        assert!(d_01 >= 0.75 * jac_theo_0);
        assert!(d_02 <= 0.1);
    }  // end of test_superminhash_kmer_16b32bit_serial

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
    if !sketch_reader_res.is_ok() {
        return;
    }
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

