//! This file split a sequence in blocks and sketch each block.I t is then possible to compute 
//! Jaccard Probability Distance between blocks.
//! This module is adaptef to long reads with variable length.
//! 

#![allow(dead_code)]



use log::*;



#[allow(unused_imports)]
use std::hash::{BuildHasher, BuildHasherDefault, Hasher, Hash};
use indexmap::{IndexMap};
#[allow(unused_imports)]
use std::collections::HashMap;

use fnv::{FnvBuildHasher};

use std::io;
use std::io::{Write,Read, ErrorKind};

use std::fs;
use std::fs::OpenOptions;

use crate::nohasher::*;

use typename::TypeName;

use crate::kmergenerator::*;

use rayon::prelude::*;

type FnvIndexMap<K, V> = IndexMap<K, V, FnvBuildHasher>;

use probminhash::probminhasher::*;
use hnsw_rs::prelude::*;


const MAGIC_BLOCKSIG_DUMP : u32 = 0xceabbadd;



/// a block will cover kmer beginning in [i*blockSize : (i+1)*blocksize] so se have some kmers in common between adjacent blocks
/// The type do not implement Copy. So we must take care of using references to avoid cloning
#[derive(Clone, TypeName)]
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



/// For each i, skectch[i] is a vector of size 1!! A bit cumbersome but we need this for Hnsw
/// As BlockSketched do not satisfay Copy, there is a move out when we want to allocate a Vec<BlockSketched>
/// to send it do Distance
pub struct BlockSketchedSeq {
    numseq : usize, 
    pub sketch : Vec<Vec<BlockSketched>>,
}



/// This structure gathers parameters used to sketch sequences in blocks
pub struct BlockSeqSketcher {
    sig_size : u8,
    block_size : usize,
    kmer_size : usize,
    sketch_size : usize
}



impl BlockSeqSketcher {
    /// defines blcok_size, kmer_size, and sketch_size (i.e number of u32 or u64 used in sketching each block)
    pub fn new(block_size : usize, kmer_size: usize, sketch_size:usize) -> BlockSeqSketcher {
        BlockSeqSketcher{sig_size: 4 as u8, block_size, kmer_size, sketch_size}
    }
    /// skecth a sequence and returns a Vector of BlockSketched
    pub fn blocksketch_sequence<F>(& self, numseq: usize, seq : &Sequence, fhash : &F)  -> BlockSketchedSeq 
        where F : Fn(&Kmer32bit) -> u32 + Sync + Send {
        //
        let nb_blocks = if seq.size() % self.block_size == 0 { seq.size() / self.block_size} else  { 1 + seq.size() / self.block_size};
        // estimate number of block to preallocated result
        let mut sketch = Vec::<Vec<BlockSketched>>::with_capacity(nb_blocks);
        //
        // get a kmer generator generate all kmers include in range arg. dependance upon kmer_size 
        // 
        let mut kmergen = KmerSeqIterator::<Kmer32bit>::new(self.kmer_size as u8, &seq);
        kmergen.set_range(0, seq.size()).unwrap();
        let mut wa : FnvIndexMap::<u32,f64> = FnvIndexMap::with_capacity_and_hasher(self.block_size, FnvBuildHasher::default());
        // loop on blocks and compute signature of each block
        for numblock in 0..nb_blocks {
            let mut pminhasha = ProbMinHash3a::<u32,NoHashHasher>::new(self.sketch_size, 0);
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
            let current_block = vec![BlockSketched{ numseq:numseq as u32, numblock: numblock as u32, sketch:siga.clone()}];
            sketch.push(current_block);
            wa.clear();
        } // end of for numblock
        //
        return BlockSketchedSeq{numseq,sketch};
}  // end of sketch_sequence_in_blocks


    /// sketch in blocks a pack of sequences. Use Rayon
    pub fn blocksketch_sequences<F>(& self, pack_seq : &[(u32,&Sequence)], fhash : &F)  -> Vec<BlockSketchedSeq> 
        where  F : Fn(&Kmer32bit) -> u32 + Sync + Send {
        //        
        // we sketch in sequence in parallel
        let block_sketched : Vec<BlockSketchedSeq> = pack_seq.into_par_iter().map(|(i, seq)| self.blocksketch_sequence(*i as usize, seq, fhash)).collect();
        return block_sketched;
    }

    /// dump a whole pack of sketches relative to a set of sequences split in blocks and sketched
    pub fn dump_blocks(&self, out : &mut dyn Write, seqblocks : &Vec<BlockSketchedSeq>) {
        for i in 0..seqblocks.len() {
            // dump numseq
            out.write(&seqblocks[i].numseq.to_le_bytes()).unwrap();
            let nbblock = seqblocks[i].sketch.len();
            // dump blocks
            for j in 0..nbblock {
                let block = &(seqblocks[i].sketch)[j][0];
                block.dump(out);
            }
        }
    }  // end of dump_blocks


    /// initialize dump file. Nota we intialize with size of key signature : 4 bytes.  
    /// 
    /// Format of file is :
    /// -  MAGIC_SIG_DUMP as u32
    /// -  sig_size 4 or 8 dumped as u32 according to type of signature Vec<u32> or Vec<u64>
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
        let sketch_size_u32 = self.sketch_size as u32;
        let kmer_size_u32 = self.kmer_size as u32;
        let mut sigbuf : io::BufWriter<fs::File> = io::BufWriter::with_capacity(1_000_000_000, dumpfile);
        sigbuf.write(& MAGIC_BLOCKSIG_DUMP.to_le_bytes()).unwrap();
        sigbuf.write(& self.sig_size.to_le_bytes()).unwrap();
        sigbuf.write(& sketch_size_u32.to_le_bytes()).unwrap();
        sigbuf.write(& kmer_size_u32.to_le_bytes()).unwrap();
        sigbuf.write(& self.block_size.to_le_bytes()).unwrap();
        //
        return sigbuf;
    }  // end of create_signature_dump

} // end implementation block for BlockSeqSketcher




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
        // read block size
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
    /// emulates iterator API.
    /// Return next object's signature (a Vec<u32> ) if any, None otherwise.
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


} // end impl block for SigBlockSketchFileReader





// ==================================================================================
//  Distance between Blocks
// ==================================================================================

#[derive(TypeName, Default)]
pub struct DistBlockSketched {
}


// define a distance between BlockSketched
// it is 1. if the blocks comes from the same sequence and a jaccard distance computed by probminhash in the other case.
// Point in hnsw_rs will have as data a slice [BlockSketched] of length 1.
// This way we emulate an eval function which could have been  eval(&Obj1, &Obj2) ...


impl  Distance<BlockSketched> for  DistBlockSketched {

    fn eval(&self, va: &[BlockSketched], vb:&[BlockSketched]) -> f32 {
        //
        assert!(va.len() == 1 && vb.len() == 1);
        // set maximal distance inside a sequence as we want to pair reads
        if va[0].numseq == vb[0].numseq {
            return 1.;
        }
        //
        assert_eq!(va[0].sketch.len(), vb[0].sketch.len());
        //
        let sklen = va[0].sketch.len();
        let mut nb_diff = 0u32;
        for i in 0..sklen {
            if va[0].sketch[i] != vb[0].sketch[i] {
                nb_diff += 1;
            }
        }
        //
        return nb_diff as f32/sklen as f32;
    }  // end of eval

} // end implementation Distance for DistBlockSketched


//==============================================================================

mod test {

#[allow(unused_imports)]
use super::*;

#[allow(dead_code)]
    fn log_init_test() {
        let mut builder = env_logger::Builder::from_default_env();
        //    builder.filter_level(LevelFilter::Trace);
        let _ = builder.is_test(true).try_init();
    }


#[test]
    fn test_block_sketch() {
        log_init_test();
        // define a closure for our hash function
        let kmer_revcomp_hash_fn = | kmer : &Kmer32bit | -> u32 {
            let canonical =  kmer.reverse_complement().min(*kmer);
            let hashval = probminhash::invhash::int32_hash(canonical.0);
            hashval
        }; 
        //
        let seqstra = String::from("TCAAAGGGAAACATTCAAAATCAGTATGCGCCCGTTCAGTTACGTATTGCTCTCGCCGTAGGCCTAATGAGATGGGCTGGGTACAGAG");
        let seqa = Sequence::new(seqstra.as_bytes(),2);
        // seqstrb is seqstr with 2 small block of T's inside that were modified
        let seqstrb = String::from("TCAAAGGGAAATTTTTTTCATTCAAAATCAGTATGCGCCCGTTCAGTTACGTATTGCTCTCGCCGTAGGCCTAATGATTTTTTTGATGGGCTGGGTACAGAG");
        let seqb = Sequence::new(seqstrb.as_bytes(),2);

        let block_size = 10;
        let kmer_size = 3;
        let sketch_size = 6;

        let sketcher = BlockSeqSketcher::new(block_size, kmer_size, sketch_size);
        let sketcha = sketcher.blocksketch_sequence(1, &seqa, &kmer_revcomp_hash_fn);
        let sketchb = sketcher.blocksketch_sequence(2, &seqb, &kmer_revcomp_hash_fn);
        // check number of blocks sketch obtained
        println!("sketcha has number o blocks = {:?}",  sketcha.sketch.len());
        println!("sketchb has number o blocks = {:?}",  sketchb.sketch.len());
        // check of distance computations
        let mydist = DistBlockSketched{};
        assert_eq!(mydist.eval(&sketcha.sketch[0], &sketcha.sketch[0]), 1.);
        //
        let dist_1 = mydist.eval(&sketcha.sketch[0], &sketchb.sketch[0]);
        println!("dist_1 = {:?}", dist_1);
        log::info!("dist_1 = {:?}", dist_1);
        //
        let dist_2 = mydist.eval(&sketcha.sketch[1], &sketchb.sketch[1]);
        println!("dist_2 = {:?}", dist_2);
        log::info!("dist_2 = {:?}", dist_2);
    } // end of test_block_sketch

}  // end of module test