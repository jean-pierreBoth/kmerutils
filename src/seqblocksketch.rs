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

use std::io::{Write};


use crate::nohasher::*;

use crate::kmergenerator::*;

use rayon::prelude::*;

type FnvIndexMap<K, V> = IndexMap<K, V, FnvBuildHasher>;

use probminhash::probminhasher::*;
use hnsw_rs::prelude::*;

pub struct BlockSeqSketcher {
    block_size : usize,
    kmer_size : usize,
    sketch_size : usize
}


impl BlockSeqSketcher {


    /// skecth a sequence and returns a Vector of BlockSketched
    fn blocksketch_sequence<F>(& self, numseq: usize, seq : &Sequence, fhash : &F)  -> Vec<BlockSketched> 
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
        return sketch;
}  // end of sketch_sequence_in_blocks


    /// sketch in blocks a pack of sequences. Use Rayon
    fn blocksketch_sequences<F>(& self, pack_seq : &[(u32,&Sequence)], fhash : &F)  -> Vec<Vec<BlockSketched>> 
        where  F : Fn(&Kmer32bit) -> u32 + Sync + Send {
        //        
        // we sketch in sequence in parallel
        let block_sketched : Vec::<Vec<BlockSketched>> = pack_seq.into_par_iter().map(|(i, seq)| self.blocksketch_sequence(*i as usize, seq, fhash)).collect();
        return block_sketched;
    }

    /// dump a whole pack of sketches relative to a set of sequences split in blocks and sketched
    fn dump_blocks(&self, out : &mut dyn Write, seqblocks : &Vec::<Vec<BlockSketched>>) {
        for i in 0..seqblocks.len() {
            let nbblock = seqblocks[i].len();
            for j in 0..nbblock {
                seqblocks[i][j].dump(out);
            }
        } //

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

    // dump 
    fn dump(& self, out : &mut dyn Write) {
        out.write(& self.numseq.to_le_bytes()).unwrap();
        out.write(& self.numblock.to_le_bytes()).unwrap();
        for i in 0..self.sketch.len() {
            out.write(& self.sketch[i].to_le_bytes()).unwrap();
        }
    } // end of dump

}

// The whole list of BlockSketched must be dumped and need to be reloaded from Julia



// define a distance between BlockSketched
// it is 1. if the blocks comes from the same sequence and a jaccard distance computed by probminhash in the other case.

pub struct DistBlockSketched<'a> {
    block1 : &'a BlockSketched,
    block2 : &'a BlockSketched
}

impl <'a> DistBlockSketched<'a> {
    pub fn new(block1: &'a BlockSketched, block2: &'a BlockSketched) -> DistBlockSketched<'a> {
        DistBlockSketched{block1, block2}
    }
}   // end of impl DistBlockSketched




impl <'a> Distance<u32> for  DistBlockSketched<'a> {

    fn eval(&self, va: &[u32], vb:&[u32]) -> f32 {
        if self.block1.numseq != self.block2.numseq {
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



#[allow(dead_code)]
    fn log_init_test() {
        let mut builder = env_logger::Builder::from_default_env();
    //    builder.filter_level(LevelFilter::Trace);
        let _ = builder.is_test(true).try_init();
    }

}