//! This module is dedicated to dump reload of signature obtained by
//! the probminhash family algorithms.
//! 


use std::io;
use ::std::io::{Write};
use ::std::fs;
use ::std::fs::OpenOptions;

const MAGIC_SIG_DUMP : u32 = 0xceabeadd;

// CAVEAT should go to serde/bson

// dumps in an open write buffer a vector os signatures
pub fn dump_signatures_block<D>(signatures : &Vec<Vec<D>>, out : &mut dyn Write) -> io::Result<()> {
    for i in 0..signatures.len() {
        let ptr_usize = signatures[i].as_ptr();
        let vec_u8 = unsafe {
            let ptr_u8 = std::mem::transmute::<*const D, *mut u8>(ptr_usize);
            Vec::from_raw_parts(ptr_u8, signatures[i].len() * std::mem::size_of::<D>(), 
                                        signatures[i].capacity() * std::mem::size_of::<D>())
            };
        out.write(vec_u8.as_slice())?;
    }  // end of for i
    //
    return Ok(());
} // end of dump_signatures


pub fn init_signature_dump(dumpfname:String, kmer_size : u8, sketch_size: usize) -> io::BufWriter<fs::File> {
    let dumpfile_res = OpenOptions::new().write(true).create(true).truncate(true).open(&dumpfname);
    let dumpfile;
    if dumpfile_res.is_ok() {
        dumpfile = dumpfile_res.unwrap();
    } else {
        println!("cannot open {}", dumpfname);
        std::process::exit(1);
    }
    let mut sigbuf : io::BufWriter<fs::File> = io::BufWriter::with_capacity(1_000_000_000, dumpfile);
    sigbuf.write(& MAGIC_SIG_DUMP.to_be_bytes()).unwrap();
    sigbuf.write(& sketch_size.to_be_bytes()).unwrap();
    sigbuf.write(& kmer_size.to_be_bytes()).unwrap();
    sigbuf.write(& 4u32.to_be_bytes()).unwrap();     // We dump signature encoded on 4 bytes!
    //
    return sigbuf;
}