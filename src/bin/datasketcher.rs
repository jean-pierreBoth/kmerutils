//! a tool to sketch sequences of a data file
//! The algorithm used for sketching are probminhash3a and probminhash
//! usage datasketcher -m "method" -o outputfilename -s size
//! - -s size gives the size of signature to use for each sequence.
//!     it depends upon the size of sequence to sketch and the precision needed for further jaccard distance estimation


#[allow(unused_imports)]
use log::*;
#[allow(unused_imports)]
use log::Level::{Debug,Trace};
use env_logger::{Builder};

use clap::{App, Arg};
use time::*;
use ::std::process;
use std::path::Path;




// install a logger facility
fn init_log() -> u64 {
    env_logger::Builder::from_default_env().init();
    println!("\n ************** initializing logger *****************\n");    
    return 1;
}

/// format of file dump
///  - magic 
///  - size of signature as u64
///  - size of component of the signature. As it is a hashed kmer represented by a u32 or a u64,  it is 4 or 8 bytes. 
///  then for each sequence  a vector of size the signature and of type u32 or u64


fn main() {

    init_log();

    let matches = App::new("datasketcher")       
            .arg(Arg::with_name("file")
            .long("file")
            .short("f")
            .takes_value(true)
            .help("expecting a fastq file"))
            .get_matches();

    let fname;

    if matches.is_present("file") {
        fname = matches.value_of("file").ok_or("bad value").unwrap().parse::<String>().unwrap();
        println!("got filename , {}", fname);
    }
    else {
        println!("-f filename is mandatory");
        println!(" usage qualityloade -f name --wawelet (or -w) --port (-p) num");
        process::exit(1);
    }

    let path = Path::new(&fname);
    let f_info_res = path.metadata();
    match f_info_res {
        Ok(meta) => {
            let filesize = meta.len();
            info!("sketching file {}   size : {}", fname, filesize);         
        },
        Err(_e) => {
            error!("file does not exist: {:?}", fname);
            process::exit(1);
        },
    }
//    let start_t = Instant::now();
    let mut reader = needletail::parse_fastx_file(&path).expect("expecting valid filename");

}