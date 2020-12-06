//! This file split a sequence in blocks and sketch each block.I t is then possible to compute 
//! Jaccard Probability Distance between blocks.
//! This module is adaptef to long reads with variable length.
//! 

#![allow(dead_code)]


#[allow(unused_imports)]

use log::*;



#[allow(unused_imports)]
use std::hash::{BuildHasher, BuildHasherDefault, Hasher, Hash};
use indexmap::{IndexMap};
use fnv::{FnvBuildHasher};

use std::io;
use std::io::{Write,Read};

use std::fs;
use std::fs::OpenOptions;

use crate::nohasher::*;

use crate::kmergenerator::*;

use rayon::prelude::*;

type FnvIndexMap<K, V> = IndexMap<K, V, FnvBuildHasher>;

use probminhash::probminhasher::*;
use hnsw_rs::prelude::*;


const MAGIC_BLOCKSIG_DUMP : u32 = 0xceabbadd;

pub struct BlockSeqSketcher {
    block_size : usize,
    kmer_size : usize,
    sketch_size : usize
}



type BlockSketchedSeq = (usize, Vec<BlockSketched>);


impl BlockSeqSketcher {
    //
    fn new(block_size : usize, kmer_size: usize, sketch_size:usize) -> BlockSeqSketcher {
        BlockSeqSketcher{block_size, kmer_size, sketch_size}
    }
    /// skecth a sequence and returns a Vector of BlockSketched
    fn blocksketch_sequence<F>(& self, numseq: usize, seq : &Sequence, fhash : &F)  -> BlockSketchedSeq 
        where F : Fn(&Kmer32bit) -> u32 + Sync + Send {
        //
        let nb_blocks = if seq.size() % self.block_size == 0 { seq.size() / self.block_size} else  { 1 };
        // estimate number of block to preallocated result
        let mut sketch = Vec::<BlockSketched>::with_capacity(nb_blocks);
        //
        // get a kmer generator generate all kmers include in range arg. dependance upon kmer_size 
        // 
        let mut kmergen = KmerSeqIterator::<Kmer32bit>::new(self.kmer_size as u8, &seq);
        kmergen.set_range(0, seq.size()).unwrap();
        // loop on blocks and compute signature of each block
        for numblock in 0..nb_blocks {
            let mut pminhasha = ProbMinHash3a::<u32,NoHashHasher>::new(self.sketch_size, 0);
            let mut wa : FnvIndexMap::<u32,f64> = FnvIndexMap::with_capacity_and_hasher(self.block_size, FnvBuildHasher::default());
            let mut kmer_pos = 0;
            loop {
                match kmergen.next() {
                    Some(kmer) => {
                        let hashval = fhash(&kmer);
                        *wa.entry(hashval).or_insert(0.) += 1.;
                        // have we reached end of block
                        if kmer_pos == self.block_size - 1 {
                            break;
                        }
                        else {
                            kmer_pos += 1;   
                        }
                    },
                    None => { 
                        // end of sequence we have to finish sketching this (possibly not full block)
                        break;
                    },
                }
            }  // end loop
            pminhasha.hash_weigthed_idxmap(&wa);
            let siga = pminhasha.get_signature();
            let current_block = BlockSketched{ numseq:numseq as u32, numblock: numblock as u32, sketch:siga.clone()};
            sketch.push(current_block);
        } // end of for numblock
        //
        return (numseq,sketch);
}  // end of sketch_sequence_in_blocks


    /// sketch in blocks a pack of sequences. Use Rayon
    fn blocksketch_sequences<F>(& self, pack_seq : &[(u32,&Sequence)], fhash : &F)  -> Vec<BlockSketchedSeq> 
        where  F : Fn(&Kmer32bit) -> u32 + Sync + Send {
        //        
        // we sketch in sequence in parallel
        let block_sketched : Vec<BlockSketchedSeq> = pack_seq.into_par_iter().map(|(i, seq)| self.blocksketch_sequence(*i as usize, seq, fhash)).collect();
        return block_sketched;
    }

    /// dump a whole pack of sketches relative to a set of sequences split in blocks and sketched
    fn dump_blocks(&self, out : &mut dyn Write, seqblocks : &Vec<BlockSketchedSeq>) {
        for i in 0..seqblocks.len() {
            // dump numseq
            out.write(&seqblocks[i].0.to_le_bytes()).unwrap();
            let nbblock = seqblocks[i].1.len();
            // dump blocks
            for j in 0..nbblock {
                (seqblocks[i].1)[j].dump(out);
            }
        }
    }  // end of dump_blocks


} // end implementation block for BlockSeqSketcher





/// a block will cover kmer beginning in [i*blockSize : (i+1)*blocksize] so se have some kmers in common between adjacent blocks
pub struct BlockSketched {
    numseq: u32,
    numblock : u32,
    sketch: Vec<u32>
}


impl BlockSketched {
    /// allocator
    pub fn new(numseq: u32, numblock:u32, sketch_size : u32) -> BlockSketched {
        let sketch = Vec::<u32>::with_capacity(sketch_size as usize);
        BlockSketched{numseq, numblock, sketch}
    }

    pub fn get_skech_slice(&self) -> &[u32] {
        &self.sketch
    }
    // dump 
    fn dump(& self, out : &mut dyn Write) {
        out.write(& self.numseq.to_le_bytes()).unwrap();
        out.write(& self.numblock.to_le_bytes()).unwrap();
        for i in 0..self.sketch.len() {
            out.write(& self.sketch[i].to_le_bytes()).unwrap();
        }
    } // end of dump

}  // end of impl block for BlockSketched

// The whole list of BlockSketched must be dumped and need to be reloaded from Julia


#[allow(dead_code)]
/// structure to reload a file consisting of sketch blocks
/// each read should load a whole BlockSketchedSeq
struct SigBlockSketchFileReader {
    fname:String,
    /// signature size in bytes. 4 for u32, 8 for u64
    sig_size : u8,
    /// the number of sketch by object hashed
    sketch_size:usize,
    /// size of kmers used in sketching.
    kmer_size : u8,
    /// block size inside a sequence
    block_size : u32,
    /// read buffer 
    signature_buf:io::BufReader<fs::File>
}  // end of struct SigBlockSketchFileReader




impl SigBlockSketchFileReader {
    /// initialize the fields fname, sketch_size, kmer_size and allocates signature_buf but signatures will be read by next.
    #[allow(dead_code)]
    pub fn new(fname:&String) -> Result<SigBlockSketchFileReader, String> {
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
        //
        // check magic
        //
        io_res = signature_buf.read_exact(&mut buf_u32);
        if io_res.is_err() {
            println!("SigBlockSketchFileReader could no read magic");
            return Err(String::from("SigBlockSketchFileReader could no read magic"));
        }
        let magic = u32::from_le_bytes(buf_u32);
        if magic != MAGIC_BLOCKSIG_DUMP {
            println!("file {} is not a dump of signature", fname);
            return Err(String::from("file is not a dump of signature"));
        }
        //
        // read sig_size
        //
        io_res = signature_buf.read_exact(&mut buf_u32);
        if io_res.is_err() {
            println!("SigBlockSketchFileReader could no read sketch_size");
            return Err(String::from("SigBlockSketchFileReader could no read sketch_size"));
        }
        let sig_size = u32::from_le_bytes(buf_u32);
        if sig_size != 4 {
            println!("SigBlockSketchFileReader could no read sketch_size");
            return Err(String::from("SigBlockSketchFileReader , sig_size != 4 not yet implemented"));            
        }
        //
        // read sketch_size
        //
        io_res = signature_buf.read_exact(&mut buf_u32);
        if io_res.is_err() {
            println!("SigBlockSketchFileReader could no read sketch_size");
            return Err(String::from("SigBlockSketchFileReader could no read sketch_size"));
        }
        let sketch_size = u32::from_le_bytes(buf_u32);
        trace!("read sketch size {}", sketch_size);
        //
        // check kmer_size
        //
        io_res = signature_buf.read_exact(&mut buf_u32);
        if io_res.is_err() {
            println!("SigBlockSketchFileReader could no read kmer_size");
            return Err(String::from("SigBlockSketchFileReader could no read kmer_size"));
        }
        let kmer_size = u32::from_le_bytes(buf_u32);
        trace!("read kmer_size {}", kmer_size);
        //
        // read blcok size
        //
        io_res = signature_buf.read_exact(&mut buf_u32);
        if io_res.is_err() {
            println!("SigBlockSketchFileReader could no read block_size");
            return Err(String::from("SigBlockSketchFileReader could no read block_size"));
        }
        let block_size = u32::from_le_bytes(buf_u32);
        trace!("read block_size {}", block_size);
        //
        Ok(SigBlockSketchFileReader{fname: fname.clone() , sig_size: sig_size as u8 , sketch_size: sketch_size as usize, 
                kmer_size: kmer_size as u8, block_size, signature_buf})
    } // end of new

} // end impl block for SigBlockSketchFileReader





// ==================================================================================
//  Distance between Blocks
// ==================================================================================


pub struct DistBlockSketched<'a> {
    block1 : &'a BlockSketched,
    block2 : &'a BlockSketched
}

impl <'a> DistBlockSketched<'a> {
    pub fn new(block1: &'a BlockSketched, block2: &'a BlockSketched) -> DistBlockSketched<'a> {
        DistBlockSketched{block1, block2}
    }
}   // end of impl DistBlockSketched

// define a distance between BlockSketched
// it is 1. if the blocks comes from the same sequence and a jaccard distance computed by probminhash in the other case.



impl <'a> Distance<u32> for  DistBlockSketched<'a> {

    fn eval(&self, va: &[u32], vb:&[u32]) -> f32 {
        // set maximal distance inside a sequence as we want to pair reads
        if self.block1.numseq == self.block2.numseq {
            return 1.;
        }
        //
        assert_eq!(va.len(), vb.len());
        //
        let mut nb_diff = 0u32;
        for i in 0..va.len() {
            if va[i] != va[i] {
                nb_diff += 1;
            }
        }
        //
        return nb_diff as f32/ va.len() as f32;
    }  // end of eval

} // end implementation Distance for DistBlockSketched


//==============================================================================

mod test {

use super::*;

#[allow(dead_code)]
    fn log_init_test() {
        let mut builder = env_logger::Builder::from_default_env();
        //    builder.filter_level(LevelFilter::Trace);
        let _ = builder.is_test(true).try_init();
    }


#[test]
    fn test_block_sketch() {
        // define a closure for our hash function
        let kmer_revcomp_hash_fn = | kmer : &Kmer32bit | -> u32 {
            let canonical =  kmer.reverse_complement().min(*kmer);
            let hashval = probminhash::invhash::int32_hash(canonical.0);
            hashval
        }; 
        //
        let seqstra = String::from("TCAAAGGGAAACATTCAAAATCAGTATGCGCCCGTTCAGTTACGTATTGCTCTCGCCGTAGGCCTAATGAGATGGGCTGGGTACAGAG");
        let seqa = Sequence::new(seqstra.as_bytes(),2);
        //
        let seqstrb = String::from("TCAAAGCGTCGTATAGCCGGAAACATTCAAAATCAGTATGCGCCCGTTCAGTTACGTATTGCTCTCGCTAATGAGATGGGCTGGGTACAGAG");
        let seqb = Sequence::new(seqstrb.as_bytes(),2);

        let block_size = 10;
        let kmer_size = 4;
        let sketch_size = 6;

        let sketcher = BlockSeqSketcher::new(block_size, kmer_size, sketch_size);
        let sketcha = sketcher.blocksketch_sequence(1, &seqa, &kmer_revcomp_hash_fn);
        let sketchb = sketcher.blocksketch_sequence(2, &seqb, &kmer_revcomp_hash_fn);
    } // end of test_block_sketch

}  // end of module test