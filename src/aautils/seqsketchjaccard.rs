//! provide minimal tool to sketch RNA sequences by probminhash3a


//#![allow(unused)]

use std::io::{BufReader, BufWriter };


use std::fs::OpenOptions;
use std::fmt::Debug;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};
use serde_json::{to_writer};

use fnv::{FnvHashMap, FnvBuildHasher};

use num;

#[allow(unused)]
use crate::nohasher::*;

use crate::base::{kmertraits::*};
use crate::aautils::{kmeraa::*};

use rayon::prelude::*;


use probminhash::probminhasher::*;



// TODO this should be factorized with DNA case.

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





    /// A generic version of sketching with probminhash3a on compressed kmer for amino acids
    pub fn sketch_probminhash3a_compressedkmeraa<'b, Kmer : CompressedKmerT + KmerBuilder<Kmer>, F>(&self, vseq : &'b Vec<&SequenceAA>, fhash : F) -> Vec<Vec<Kmer::Val> >
        where F : Fn(&Kmer) -> Kmer::Val + Send + Sync,
              Kmer::Val : num::PrimInt + Send + Sync + Debug,
              KmerGenerator<Kmer> :  KmerGenerationPattern<Kmer> {

            let comput_closure = | seqb : &SequenceAA, i:usize | -> (usize,Vec<Kmer::Val>) {
                // if we get very large sequence (many Gb length) we must be cautious on size of hashmap; i.e about number of different kmers!!! 
                let nb_kmer = get_nbkmer_guess(seqb);
                let mut wb : FnvHashMap::<Kmer::Val,u64> = FnvHashMap::with_capacity_and_hasher(nb_kmer, FnvBuildHasher::default());
                let mut kmergen = KmerSeqIterator::<Kmer>::new(self.kmer_size, &seqb);
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
    }  // end of sketch_probminhash3a_compressedkmeraa

} // end of SeqSketcher (RNA case)



//=========================================================


#[cfg(test)]
mod tests {

use super::*;
use std::str::FromStr;

    fn log_init_test() {
        let mut builder = env_logger::Builder::from_default_env();
        //    builder.filter_level(LevelFilter::Trace);
        let _ = builder.is_test(true).try_init();
    }

    #[test]
    fn test_seqaa_probminhash_64bit() {
        log_init_test();
        //
        log::debug!("test_seqaa_probminhash");
        //
        let str1 = "MTEQIELIKLYSTRILALAAQMPHVGSLDNPDASAMKRSPLCGSKVTVDVIMQNGKITFDGFEVLAPASEYKNRHASILLSLDATAEACASIAAQNSA";
        // The second string is the first half of the first repeated
        let str2 = "MTEQIELIKLYSTRILALAAQMPHVGSLDNPDASAMKRSPLCGSKVMTEQIELIKLYSTRILALAAQMPHVGSLDNPDASAMKRSPLCGSKV";

        let seq1 = SequenceAA::from_str(str1).unwrap();
        let seq2 = SequenceAA::from_str(str2).unwrap();
        let vseq = vec![&seq1, &seq2];
        let kmer_size = 5;
        let sketch_size = 400;
        let sketcher = SeqSketcher::new(kmer_size, sketch_size);
        let nb_alphabet_bits = Alphabet::new().get_nb_bits();
        // we need a hash function from u128 to f64
        let kmer_hash_fn = | kmer : &KmerAA64bit | -> <KmerAA64bit as CompressedKmerT>::Val {
            let mask : <KmerAA64bit as CompressedKmerT>::Val = num::NumCast::from::<u64>((0b1 << nb_alphabet_bits*kmer.get_nb_base()) - 1).unwrap();
            let hashval = kmer.get_compressed_value() & mask;
            hashval
        };
        let mask : u64 = num::NumCast::from::<u64>((0b1 << nb_alphabet_bits*kmer_size as u8) - 1).unwrap();
        log::debug!("mask = {:b}", mask);
        //
        log::info!("calling sketch_probminhash3a_compressedKmerAA64bit");
        let signatures = sketcher.sketch_probminhash3a_compressedkmeraa(&vseq, kmer_hash_fn); 
        // get distance between the 2 strings  
        let sig1 = &signatures[0];
        let sig2 = &signatures[1];
        //
        let inter : u64 = sig1.iter().zip(sig2.iter()).map(|(a,b)| if a==b {1} else {0}).sum();
        let dist = inter as f64/sig1.len() as f64;
        log::info!("inter : {:?} length {:?} jaccard distance {:?}", inter, sig1.len(), dist );
        assert!( (dist-0.5).abs() < 1./10.);
    } // end of test_seqaa_probminhash_64bit



    #[test]
    fn test_seqaa_probminhash_32bit() {
        log_init_test();
        //
        log::debug!("test_seqaa_probminhash");
        //
        let str1 = "MTEQIELIKLYSTRILALAAQMPHVGSLDNPDASAMKRSPLCGSKVTVDVIMQNGKITFDGFEVLAPASEYKNRHASILLSLDATAEACASIAAQNSA";
        // The second string is the first half of the first repeated
        let str2 = "MTEQIELIKLYSTRILALAAQMPHVGSLDNPDASAMKRSPLCGSKVMTEQIELIKLYSTRILALAAQMPHVGSLDNPDASAMKRSPLCGSKV";

        let seq1 = SequenceAA::from_str(str1).unwrap();
        let seq2 = SequenceAA::from_str(str2).unwrap();
        let vseq = vec![&seq1, &seq2];
        let kmer_size = 5;
        let sketch_size = 400;
        let sketcher = SeqSketcher::new(kmer_size, sketch_size);
        let nb_alphabet_bits = Alphabet::new().get_nb_bits();
        // we need a hash function from u128 to f64
        let kmer_hash_fn = | kmer : &KmerAA32bit | -> <KmerAA32bit as CompressedKmerT>::Val {
            let mask : <KmerAA32bit as CompressedKmerT>::Val = num::NumCast::from::<u32>((0b1 << nb_alphabet_bits*kmer.get_nb_base()) - 1).unwrap();
            let hashval = kmer.get_compressed_value() & mask;
            hashval
        };
        let mask : u64 = num::NumCast::from::<u32>((0b1 << nb_alphabet_bits*kmer_size as u8) - 1).unwrap();
        log::debug!("mask = {:b}", mask);
        //
        log::info!("calling sketch_probminhash3a_compressedKmerAA32bit");
        let signatures = sketcher.sketch_probminhash3a_compressedkmeraa(&vseq, kmer_hash_fn); 
        // get distance between the 2 strings  
        let sig1 = &signatures[0];
        let sig2 = &signatures[1];
        //
        let inter : u64 = sig1.iter().zip(sig2.iter()).map(|(a,b)| if a==b {1} else {0}).sum();
        let dist = inter as f64/sig1.len() as f64;
        log::info!("inter : {:?} length {:?} jaccard distance {:?}", inter, sig1.len(), dist );
        assert!( (dist-0.5).abs() < 1./10.);
    } // end of test_seqaa_probminhash_32bit



    #[test]
    fn test_seqaa_probminhash_gen() {
        log_init_test();
        //
        log::debug!("test_seqaa_probminhash");
        //
        let str1 = "MTEQIELIKLYSTRILALAAQMPHVGSLDNPDASAMKRSPLCGSKVTVDVIMQNGKITFDGFEVLAPASEYKNRHASILLSLDATAEACASIAAQNSA";
        // The second string is the first half of the first repeated
        let str2 = "MTEQIELIKLYSTRILALAAQMPHVGSLDNPDASAMKRSPLCGSKVMTEQIELIKLYSTRILALAAQMPHVGSLDNPDASAMKRSPLCGSKV";

        let seq1 = SequenceAA::from_str(str1).unwrap();
        let seq2 = SequenceAA::from_str(str2).unwrap();
        let vseq = vec![&seq1, &seq2];
        let kmer_size = 5;
        let sketch_size = 400;
        let sketcher = SeqSketcher::new(kmer_size, sketch_size);
        let nb_alphabet_bits = Alphabet::new().get_nb_bits();
        // we need a hash function from u128 to f64
        let kmer_hash_fn = | kmer : &KmerAA32bit | -> <KmerAA32bit as CompressedKmerT>::Val {
            let mask : <KmerAA32bit as CompressedKmerT>::Val = num::NumCast::from::<u32>((0b1 << nb_alphabet_bits*kmer.get_nb_base()) - 1).unwrap();
            let hashval = kmer.get_compressed_value() & mask;
            hashval
        };
        let mask : u64 = num::NumCast::from::<u32>((0b1 << nb_alphabet_bits*kmer_size as u8) - 1).unwrap();
        log::debug!("mask = {:b}", mask);
        //
        log::info!("calling sketch_probminhash3a_compressed_kmeraa for KmerAA32bit");
        let signatures = sketcher.sketch_probminhash3a_compressedkmeraa(&vseq, kmer_hash_fn); 
        // get distance between the 2 strings  
        let sig1 = &signatures[0];
        let sig2 = &signatures[1];
        //
        let inter : u64 = sig1.iter().zip(sig2.iter()).map(|(a,b)| if a==b {1} else {0}).sum();
        let dist = inter as f64/sig1.len() as f64;
        log::info!("inter : {:?} length {:?} jaccard distance {:?}", inter, sig1.len(), dist );
        assert!( (dist-0.5).abs() < 1./10.);
    } // end of test_seqaa_probminhash_32bit


}  // end of mod tests in rnautils::seqsketchjaccard