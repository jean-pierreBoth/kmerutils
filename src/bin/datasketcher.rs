//! a tool to sketch sequences of a data file
//! The algorithm used for sketching are probminhash3a and probminhash2
//! usage datasketcher -m "method" -o outputfilename -s sketch_size  -k kmer_size
//! - -s  sketch_size gives the size of signature to use for each sequence.
//!     it depends upon the size of sequence to sketch and the precision needed for further jaccard distance estimation
//! - -k kmer_size gives the size of kmer to use 
//!

#[allow(unused_imports)]
use log::*;
#[allow(unused_imports)]
use log::Level::{Debug,Trace, Info};
#[allow(unused_imports)]
use env_logger::{Builder};

use clap::{App, Arg};
use time::*;
use ::std::process;
use std::path::Path;

use kmerutils::sequence::*;
use kmerutils::kmergenerator::*;
use kmerutils::jaccardweight::*;

use needletail::FastxReader;

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

    let matches = Box::new(App::new("datasketcher")       
                    .arg(Arg::with_name("file")
                        .long("file")
                        .short("f")
                        .takes_value(true)
                        .help("expecting a fastq file"))
                    .arg(Arg::with_name("sketch_size")
                        .long("sketch")
                        .short("s")
                        .takes_value(true)
                        .help("expecting sketch size as usize"))
                    .arg(Arg::with_name("kmer_size")
                        .long("kmer")
                        .short("k")
                        .takes_value(true)
                        .help("expecting a kmer size"))
                ) 
                .get_matches();

    let fname;
    let kmer_size;
    let sketch_size;
    // check for all necessary args
    if matches.is_present("file") {
        fname = matches.value_of("file").ok_or("bad value").unwrap().parse::<String>().unwrap();
        println!("got filename , {}", fname);
    }
    else {
        println!("-f filename is mandatory");
        println!(" usage : seqsketcher -f name --sketch (or -s)  s_size --kmer (-k) k_size");
        process::exit(1);
    }
    if matches.is_present("sketch_size") {
        sketch_size = matches.value_of("sketch_size").ok_or("bad value").unwrap().parse::<usize>().unwrap();
        println!("got sketch_size , {}", sketch_size);
    }
    else {
        println!("--sketch is mandatory");
        println!(" usage : seqsketcher -f name --sketch (or -s)  s_size --kmer (-k) k_size");
        process::exit(1);
    }    
    if matches.is_present("kmer_size") {
        kmer_size = matches.value_of("kmer_size").ok_or("bad value").unwrap().parse::<u8>().unwrap();
        println!("got kmer_size , {}", kmer_size);
    }
    else {
        println!("--kmer is mandatory");
        println!(" usage : seqsketcher -f name --sketch (or -s)  s_size --kmer (-k) k_size");
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
    let start_t = Instant::now();
    let mut reader = needletail::parse_fastx_file(&path).expect("expecting valid filename");
    let blocksize = 1000;
    let mut nbseq = 0;
    //

    let kmer_revcomp_hash_fn = | kmer : &Kmer32bit | -> u32 {
        let canonical =  kmer.reverse_complement().min(*kmer);
        let hashval = probminhash::invhash::int32_hash(canonical.0);
        hashval
    };  

    loop {
        let sequenceblock = readblockseq(& mut reader, blocksize);
        if sequenceblock.len() == 0 {
            break;
        }
        nbseq += sequenceblock.len();
        // do the computation
        let signatures = sketch_probminhash3a_kmer32bit(&sequenceblock, sketch_size, kmer_size, kmer_revcomp_hash_fn);
        // dump the signature
    }
    // possibly sketch the last incomplete block

    //
    let elapsed_t = start_t.elapsed().whole_seconds();
    println!(" number of sequences loaded {} ", nbseq);
    println!(" elapsed time (s) in sketching data file {} ", elapsed_t);
}


// read a blcok of nbseq sequences
fn readblockseq(reader : &mut Box<dyn FastxReader>, nbseq : usize) -> Vec<Sequence> {
    //
    let mut veqseq = Vec::<Sequence>::with_capacity(nbseq);
    let mut nb_bad_sequence = 0;
    //
    while let Some(record) = reader.next() {
        // read by block of blocksize sequence to benefit from theading of sketching
        let seqrec = record.expect("invalid record");
        let nb_bad = count_non_acgt(&seqrec.seq());
        if nb_bad == 0 {
            nb_bad_sequence += 1;
            continue;
        }
        // sequence consists in acgt we store  it , compressing in 2 bits/base
        let newseq = Sequence::new(&seqrec.seq(), 2);
        if veqseq.len() < veqseq.capacity() {
            veqseq.push(newseq.clone());
        }
        if veqseq.len() == nbseq {
            break;
        }
    } // end while
    // get here if we have enough sequences or reached end of file
    if nb_bad_sequence > 0 && log_enabled!(Info) {
        info!(" number of non acgt sequences {} ", nb_bad_sequence);
    }
    return veqseq;
}  // end of readblockseq