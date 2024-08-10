//! This module splits a sequence in blocks and sketch each block. It is then possible to compute
//! Jaccard Probability Distance between blocks among various sequences.    
//! This module is adapted to long reads with variable length.
//!

#![allow(clippy::unnecessary_unwrap)]

use log::*;

use serde::{Deserialize, Serialize};

#[allow(unused_imports)]
use std::hash::{BuildHasher, BuildHasherDefault, Hash, Hasher};

use fnv::{FnvBuildHasher, FnvHashMap};

use std::io;
use std::io::{ErrorKind, Read, Write};

use std::fs;
use std::fs::OpenOptions;

use crate::nohasher::*;

use crate::base::kmergenerator::*;

use rayon::prelude::*;

use hnsw_rs::prelude::*;
use probminhash::probminhasher::*;

/// signature for dump of probminhash with sequences split in blocks
const MAGIC_BLOCKSIG_DUMP: u32 = 0xceabbadd;

/// a block will kmers beginning in [i*blockSize : (i+1)*blocksize] so se have some kmers in common between adjacent blocks.  
/// The type do not implement Copy. So we must take care of using references to avoid cloning
#[derive(Clone, Serialize, Deserialize)]
pub struct BlockSketched {
    numseq: u32,
    numblock: u32,
    sketch: Vec<u32>,
}

impl BlockSketched {
    /// allocator
    pub fn new(numseq: u32, numblock: u32, sketch_size: u32) -> BlockSketched {
        let sketch = Vec::<u32>::with_capacity(sketch_size as usize);
        BlockSketched {
            numseq,
            numblock,
            sketch,
        }
    }

    pub fn get_skech_slice(&self) -> &[u32] {
        &self.sketch
    }
    // dump
    fn dump(&self, out: &mut dyn Write) {
        out.write_all(&self.numseq.to_le_bytes()).unwrap();
        out.write_all(&self.numblock.to_le_bytes()).unwrap();
        for i in 0..self.sketch.len() {
            out.write_all(&self.sketch[i].to_le_bytes()).unwrap();
        }
    } // end of dump
} // end of impl block for BlockSketched

// The whole list of BlockSketched must be dumped and need to be reloaded from Julia

/// For each i, skectch\[i\] is a vector of size 1!! A bit cumbersome but we need this for Hnsw
/// as BlockSketched do not satisfay Copy, there is a move out when we want to allocate a Vec\<BlockSketched\>
/// to send it to Distance
pub struct BlockSketchedSeq {
    numseq: usize,
    pub sketch: Vec<Vec<BlockSketched>>,
}

/// This structure gathers parameters used to sketch sequences in blocks
pub struct BlockSeqSketcher {
    sig_size: u8,
    block_size: usize,
    kmer_size: usize,
    sketch_size: usize,
}

impl BlockSeqSketcher {
    /// defines blcok_size, kmer_size, and sketch_size (i.e number of u32 or u64 used in sketching each block)
    pub fn new(block_size: usize, kmer_size: usize, sketch_size: usize) -> BlockSeqSketcher {
        BlockSeqSketcher {
            sig_size: 4_u8,
            block_size,
            kmer_size,
            sketch_size,
        }
    }
    /// skecth a sequence and returns a Vector of BlockSketched
    pub fn blocksketch_sequence<F>(
        &self,
        numseq: usize,
        seq: &Sequence,
        fhash: &F,
    ) -> BlockSketchedSeq
    where
        F: Fn(&Kmer32bit) -> u32 + Sync + Send,
    {
        //
        assert!(seq.size() > 0);
        let nb_blocks = if seq.size() % self.block_size == 0 {
            seq.size() / self.block_size
        } else {
            1 + seq.size() / self.block_size
        };
        assert!(nb_blocks > 0);
        // estimate number of block to preallocated result
        let mut sketch = Vec::<Vec<BlockSketched>>::with_capacity(nb_blocks);
        //
        // get a kmer generator generate all kmers include in range arg. dependance upon kmer_size
        //
        let mut kmergen = KmerSeqIterator::<Kmer32bit>::new(self.kmer_size as u8, seq);
        kmergen.set_range(0, seq.size()).unwrap();
        let mut wa: FnvHashMap<u32, f64> =
            FnvHashMap::with_capacity_and_hasher(self.block_size, FnvBuildHasher::default());
        // loop on blocks and compute signature of each block
        for numblock in 0..nb_blocks {
            let mut pminhasha = ProbMinHash3a::<u32, NoHashHasher>::new(self.sketch_size, 0);
            let mut kmer_pos = 0;
            while let Some(kmer) = kmergen.next() {
                let hashval = fhash(&kmer);
                *wa.entry(hashval).or_insert(0.) += 1.;
                // have we reached end of block
                if kmer_pos == self.block_size - 1 {
                    break;
                } else {
                    kmer_pos += 1;
                }
            } // end loop
            pminhasha.hash_weigthed_hashmap(&wa);
            let siga = pminhasha.get_signature();
            let current_block = vec![BlockSketched {
                numseq: numseq as u32,
                numblock: numblock as u32,
                sketch: siga.clone(),
            }];
            sketch.push(current_block);
            wa.clear();
        } // end of for numblock
          //
        BlockSketchedSeq { numseq, sketch }
    } // end of sketch_sequence_in_blocks

    /// sketch in blocks a pack of sequences. Use Rayon
    pub fn blocksketch_sequences<F>(
        &self,
        pack_seq: &[(u32, &Sequence)],
        fhash: &F,
    ) -> Vec<BlockSketchedSeq>
    where
        F: Fn(&Kmer32bit) -> u32 + Sync + Send,
    {
        //
        // we sketch in sequence in parallel
        let block_sketched: Vec<BlockSketchedSeq> = pack_seq
            .into_par_iter()
            .map(|(i, seq)| self.blocksketch_sequence(*i as usize, seq, fhash))
            .collect();
        block_sketched
    }

    /// dump a whole pack of sketches relative to a set of sequences split in blocks and sketched
    /// we dump a magic as a u32, num of sequence as u32, nbblocks of a sequnece as u32
    /// and the dump of each block by BlockSketched::dump()
    pub fn dump_blocks(&self, out: &mut dyn Write, seqblocks: &[BlockSketchedSeq]) {
        for seqblock in seqblocks {
            // dump numseq
            let seqnum = seqblock.numseq as u32;
            out.write_all(&seqnum.to_le_bytes()).unwrap();
            let nbblock_u32 = seqblock.sketch.len() as u32;
            assert!(nbblock_u32 > 0);
            // dump number of blocks for sequence of current BlockSketchedSeq
            out.write_all(&nbblock_u32.to_le_bytes()).unwrap();
            // dump blocks
            for j in 0..seqblock.sketch.len() {
                let block = &(seqblock.sketch)[j][0];
                block.dump(out);
            }
        }
    } // end of dump_blocks

    /// initialize dump file. Nota we intialize with size of key signature : 4 bytes.  
    ///
    /// Format of file is :
    /// -  MAGIC_BLOCKSIG_DUMP as u32
    /// -  sig_size 4 or 8 dumped as u32 according to type of signature Vec\<u32\> or Vec\<u64\>
    /// -  sketch_size  : length of vecteur dumped as u32
    /// -  kmer_size    : as u32
    /// -  block_size   : as u32
    ///
    pub fn create_signature_dump(&self, dumpfname: &String) -> io::BufWriter<fs::File> {
        let dumpfile_res = OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true)
            .open(dumpfname);
        //
        let dumpfile = if dumpfile_res.is_ok() {
            dumpfile_res.unwrap()
        } else {
            println!("cannot open {}", dumpfname);
            std::process::exit(1);
        };
        let sketch_size_u32 = self.sketch_size as u32;
        let kmer_size_u32 = self.kmer_size as u32;
        let blocksize_u32 = self.block_size as u32;
        let mut sigbuf: io::BufWriter<fs::File> =
            io::BufWriter::with_capacity(1_000_000_000, dumpfile);
        // dump 17 bytes
        sigbuf
            .write_all(&MAGIC_BLOCKSIG_DUMP.to_le_bytes())
            .unwrap();
        sigbuf.write_all(&self.sig_size.to_le_bytes()).unwrap();
        sigbuf.write_all(&sketch_size_u32.to_le_bytes()).unwrap();
        sigbuf.write_all(&kmer_size_u32.to_le_bytes()).unwrap();
        sigbuf.write_all(&blocksize_u32.to_le_bytes()).unwrap();
        //
        sigbuf
    } // end of create_signature_dump
} // end implementation block for BlockSeqSketcher

/// structure to reload a file consisting of sketch blocks
/// each read should load a whole BlockSketchedSeq
pub struct SigBlockSketchFileReader {
    _fname: String,
    /// signature size in bytes. 4 for u32, 8 for u64
    sig_size: u8,
    /// the number of sketch by object hashed
    sketch_size: usize,
    /// size of kmers used in sketching.
    kmer_size: u8,
    /// block size inside a sequence
    block_size: u32,
    /// read buffer
    signature_buf: io::BufReader<fs::File>,
} // end of struct SigBlockSketchFileReader

impl SigBlockSketchFileReader {
    /// initialize the fields fname, sketch_size, kmer_size and allocates signature_buf but signatures will be read by next.
    pub fn new(fname: &String) -> Result<SigBlockSketchFileReader, String> {
        let dumpfile_res = OpenOptions::new().read(true).open(fname);
        let dumpfile = if dumpfile_res.is_ok() {
            dumpfile_res.unwrap()
        } else {
            println!("cannot open {}", fname);
            std::process::exit(1);
        };
        let mut signature_buf: io::BufReader<fs::File> =
            io::BufReader::with_capacity(1_000_000_000, dumpfile);
        let mut buf_u32 = [0u8; 4];
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
            return Err(String::from(
                "SigBlockSketchFileReader could no read sketch_size",
            ));
        }
        let sig_size = u32::from_le_bytes(buf_u32);
        if sig_size != 4 {
            println!("SigBlockSketchFileReader could no read sketch_size");
            return Err(String::from(
                "SigBlockSketchFileReader , sig_size != 4 not yet implemented",
            ));
        }
        //
        // read sketch_size
        //
        io_res = signature_buf.read_exact(&mut buf_u32);
        if io_res.is_err() {
            println!("SigBlockSketchFileReader could no read sketch_size");
            return Err(String::from(
                "SigBlockSketchFileReader could no read sketch_size",
            ));
        }
        let sketch_size = u32::from_le_bytes(buf_u32);
        trace!("read sketch size {}", sketch_size);
        //
        // check kmer_size
        //
        io_res = signature_buf.read_exact(&mut buf_u32);
        if io_res.is_err() {
            println!("SigBlockSketchFileReader could no read kmer_size");
            return Err(String::from(
                "SigBlockSketchFileReader could no read kmer_size",
            ));
        }
        let kmer_size = u32::from_le_bytes(buf_u32);
        trace!("read kmer_size {}", kmer_size);
        //
        // read block size
        //
        io_res = signature_buf.read_exact(&mut buf_u32);
        if io_res.is_err() {
            println!("SigBlockSketchFileReader could no read block_size");
            return Err(String::from(
                "SigBlockSketchFileReader could no read block_size",
            ));
        }
        let block_size = u32::from_le_bytes(buf_u32);
        trace!("read block_size {}", block_size);
        //
        Ok(SigBlockSketchFileReader {
            _fname: fname.clone(),
            sig_size: sig_size as u8,
            sketch_size: sketch_size as usize,
            kmer_size: kmer_size as u8,
            block_size,
            signature_buf,
        })
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

    /// returns the size in splitting sequences before sketching
    pub fn get_block_size(&self) -> usize {
        self.block_size as usize
    }
    /// emulates iterator API.
    /// Return next sequence blocks signature (a Vec\<u32\> for each block) if any, None otherwise.
    pub fn next(&mut self) -> Option<Vec<Vec<u32>>> {
        // must read sequence num and number of blocks which were dumped as u32 in little endian
        let numseq: u32;
        let mut buf_u32 = [0u8; 4];
        let io_res = self.signature_buf.read_exact(&mut buf_u32);
        if io_res.is_err() {
            println!("cannot read sequence num");
            match io_res.err().unwrap().kind() {
                ErrorKind::UnexpectedEof => return None,
                _ => {
                    println!("an unexpected error occurred reading signature buffer");
                    std::process::exit(1);
                }
            }
        } else {
            numseq = u32::from_le_bytes(buf_u32);
        }
        //
        let io_res = self.signature_buf.read_exact(&mut buf_u32);
        let nbblock: u32 = if io_res.is_err() {
            println!("cannot read number of blocks for sequence {} ", numseq);
            std::process::exit(1);
        } else {
            u32::from_le_bytes(buf_u32)
        };
        // now we have everything
        let nb_bytes = self.sketch_size * std::mem::size_of::<u32>();
        let mut buf: Vec<u8> = (0..nb_bytes).map(|_| 0u8).collect();
        let mut sig = Vec::<Vec<u32>>::with_capacity(nbblock as usize);
        for _ in 0..nbblock as usize {
            let io_res = self.signature_buf.read_exact(buf.as_mut_slice());
            //
            if io_res.is_err() {
                println!("an unexpected error occurred reading signature buffer");
                std::process::exit(1);
            } else {
                let mut sigblock = Vec::<u32>::with_capacity(self.sketch_size);
                for j in 0..self.sketch_size {
                    // split buf in chunks of 4 bytes
                    // CAVEAT the nightly as_chunks experimental API would do this properly
                    buf_u32.copy_from_slice(&buf[j * 4..4 * (j + 1)]);
                    sigblock.push(u32::from_le_bytes(buf_u32));
                }
                sig.push(sigblock);
            }
        }
        Some(sig)
    } // end of next
} // end impl block for SigBlockSketchFileReader

// ==================================================================================
//  Distance between Blocks
// ==================================================================================

/// Defines a distance between BlockSketched with value in \[0-1\]
///
/// Its value is 1. if the blocks comes from the same sequence and a jaccard distance computed by probminhash in the other case.  
/// Point in hnsw_rs will have as data a slice [BlockSketched] of length 1 as required by Hnsw interface.
// This way we emulate an eval function which could have been  eval(&Obj1, &Obj2) ...
#[derive(Default)]
pub struct DistBlockSketched {}

impl Distance<BlockSketched> for DistBlockSketched {
    fn eval(&self, va: &[BlockSketched], vb: &[BlockSketched]) -> f32 {
        //
        assert!(va.len() == 1 && vb.len() == 1);
        // set maximal distance inside a sequence as we want to pair reads
        if va[0].numseq == vb[0].numseq {
            return 1.;
        }
        //
        let nb_diff = distance_jaccard_serial(&va[0].sketch, &vb[0].sketch);
        //        let nb_diff = distance_jaccard_u32_16_simd(&va[0].sketch, &vb[0].sketch);
        //
        nb_diff as f32 / va[0].sketch.len() as f32
    } // end of eval
} // end implementation Distance for DistBlockSketched

#[inline]
fn distance_jaccard_serial(va: &[u32], vb: &[u32]) -> u32 {
    assert_eq!(va.len(), vb.len());
    let dist = va.iter().zip(vb.iter()).filter(|t| t.0 != t.1).count();
    dist as u32
} // end of distance_jaccard_serial

//==============================================================================

#[cfg(test)]
mod tests {

    #[allow(unused_imports)]
    use super::*;
    #[allow(unused_imports)]
    use rand::distributions::{Distribution, Uniform};

    fn log_init_test() {
        let mut builder = env_logger::Builder::from_default_env();
        //    builder.filter_level(LevelFilter::Trace);
        let _ = builder.is_test(true).try_init();
    }

    #[test]
    fn test_block_32bit_sketch() {
        log_init_test();
        // define a closure for our hash function
        let kmer_revcomp_hash_fn = |kmer: &Kmer32bit| -> u32 {
            let canonical = kmer.reverse_complement().min(*kmer);

            probminhash::invhash::int32_hash(canonical.0)
        };
        //
        let seqstra = String::from("TCAAAGGGAAACATTCAAAATCAGTATGCGCCCGTTCAGTTACGTATTGCTCTCGCCGTAGGCCTAATGAGATGGGCTGGGTACAGAG");
        let seqa = Sequence::new(seqstra.as_bytes(), 2);
        // seqstrb is seqstr with 2 small block of T's inside that were modified
        let seqstrb = String::from("TCAAAGGGAAATTTTTTTCATTCAAAATCAGTATGCGCCCGTTCAGTTACGTATTGCTCTCGCCGTAGGCCTAATGATTTTTTTGATGGGCTGGGTACAGAG");
        let seqb = Sequence::new(seqstrb.as_bytes(), 2);

        let block_size = 10;
        let kmer_size = 3;
        let sketch_size = 6;

        let sketcher = BlockSeqSketcher::new(block_size, kmer_size, sketch_size);
        let sketcha = sketcher.blocksketch_sequence(1, &seqa, &kmer_revcomp_hash_fn);
        let sketchb = sketcher.blocksketch_sequence(2, &seqb, &kmer_revcomp_hash_fn);
        // check number of blocks sketch obtained
        println!("sketcha has number of blocks = {:?}", sketcha.sketch.len());
        println!("sketchb has number of blocks = {:?}", sketchb.sketch.len());
        // check of distance computations
        let mydist = DistBlockSketched {};
        // distance between blocks of the same seq is set to 1 in dist function! (See comment)
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
} // end of module test
