//! This module is dedicated to dump reload of signature obtained by
//! the probminhash family algorithms.
//! 


use std::io;
use std::io::{Write,Read};
use std::fs;
use std::fs::OpenOptions;

#[allow(unused_imports)]
use log::*;

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


// initialize dump file. Nota we intialize with size of key signature : 4 bytes
pub fn create_signature_dump(dumpfname:String, kmer_size : u8, sketch_size: usize) -> io::BufWriter<fs::File> {
    let dumpfile_res = OpenOptions::new().write(true).create(true).truncate(true).open(&dumpfname);
    let dumpfile;
    if dumpfile_res.is_ok() {
        dumpfile = dumpfile_res.unwrap();
    } else {
        println!("cannot open {}", dumpfname);
        std::process::exit(1);
    }
    let mut sigbuf : io::BufWriter<fs::File> = io::BufWriter::with_capacity(1_000_000_000, dumpfile);
    sigbuf.write(& MAGIC_SIG_DUMP.to_le_bytes()).unwrap();
    sigbuf.write(& sketch_size.to_le_bytes()).unwrap();
    sigbuf.write(& kmer_size.to_le_bytes()).unwrap();
    sigbuf.write(& 4u32.to_le_bytes()).unwrap();     // We dump signature encoded on 4 bytes!
    //
    return sigbuf;
}

/// structure to reload a file consisting of sketch
struct SigSketchFileReader {
    fname:String,
    /// the number of sketch by object hashed
    sketch_size:usize,
    /// number of kmers used in sketching.
    kmer_size : u8,
    /// read buffer 
    signature_buf:io::BufReader<fs::File>
}

impl SigSketchFileReader {
    /// initialize the fields fname, sketch_size, kmer_size and allocates signature_buf but signatures will be read by next.
    pub fn new(fname:String) -> Result<SigSketchFileReader, String> {
        let dumpfile_res = OpenOptions::new().read(true).create(true).open(&fname);
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
        // read sketch_size
        io_res = signature_buf.read_exact(&mut buf_u32);
        if io_res.is_err() {
            println!("SigSketchFileReader could no read sketch_size");
            return Err(String::from("SigSketchFileReader could no read sketch_size"));
        }
        let sketch_size = u32::from_le_bytes(buf_u32);
        trace!("read sketch size {}", sketch_size);
        // check kmer_size
        io_res = signature_buf.read_exact(&mut buf_u32);
        if io_res.is_err() {
            println!("SigSketchFileReader could no read kmer_size");
            return Err(String::from("SigSketchFileReader could no read kmer_size"));
        }
        let kmer_size = u32::from_le_bytes(buf_u32);
        trace!("read kmer_size {}", kmer_size);
        // read size of signature item (4 or 8 depending on type)
        io_res = signature_buf.read_exact(&mut buf_u32);
        if io_res.is_err() {
            println!("SigSketchFileReader could no read sig_size");
            return Err(String::from("SigSketchFileReader could no read sig_size"));
        }
        // signature base size
        io_res = signature_buf.read_exact(&mut buf_u32);
        if io_res.is_err() {
            println!("SigSketchFileReader could no read sig_size");
            return Err(String::from("SigSketchFileReader could no read sig_size"));
        }
        let sig_size = u32::from_le_bytes(buf_u32);
        trace!("read sig_size {}", sig_size);
        if sig_size != 4 {
            warn!("u64 signature type not yet implemented");
            std::process::exit(1);
        }
        //
        Ok(SigSketchFileReader{fname,sketch_size: sketch_size as usize, kmer_size: kmer_size as u8, signature_buf})
    } // end of new

    /// emulates iterator API. Return next object's signature if any, None otherwise.
    pub fn next(&mut self) -> Option<Vec<u32> > {

        None
    }
     
} // end of impl SigSketchFile



