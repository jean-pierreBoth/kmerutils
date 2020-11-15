//! This module is dedicated to dump and reload of signature obtained by
//! the probminhash family algorithms.
//! 

#![allow(dead_code)]

use std::io;
use std::io::{Write,Read};
use std::io::{ErrorKind};

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


/// initialize dump file. Nota we intialize with size of key signature : 4 bytes.  
/// 
/// Format of file is :
/// -  MAGIC_SIG_DUMP as u32
/// -  sig_size 4 or 8 dumped as u32 according to type of signature Vec<u32> or Vec<u64>
/// -  sketch_size  : length of vecteur dumped as u32
/// -  kmer_size    : as u32
pub fn create_signature_dump(dumpfname:String, kmer_size : u8, sketch_size: usize) -> io::BufWriter<fs::File> {
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
    /// emulates iterator API. Return next object's signature (a Vec<u32> ) if any, None otherwise.
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



#[cfg(test)]
mod tests {
    use super::*;


#[allow(dead_code)]
fn log_init_test() {
    let mut builder = env_logger::Builder::from_default_env();
    //    builder.filter_level(LevelFilter::Trace);
    let _ = builder.is_test(true).try_init();
}



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





}  // end of mod tests