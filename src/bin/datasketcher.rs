//! a tool to sketch sequences of a data file.  
//! The algorithm used for sketching are probminhash3a and probminhash2
//! 
//! usage :  
//!    datasketcher -f fastqfile -s sketch_size  -k kmer_size ... 
//! -   -s  sketch_size gives the size of signature to use for each sequence.
//!     it depends upon the size of sequence to sketch and the precision needed for further jaccard distance estimation
//! -   -k kmer_size gives the size of kmer to use 
//! -   -b block_size (--block) to do a sketch of sequence by blocks of size block_size, else the whole sequence is sketched in one pass
//! -   -d  dumpname
//! -   ann to get hnsw embedding
//!         - nb (-n) for number of neighbours desired for future use
//!    
//! The format of signature dump file is documented in modules seqblocksketch and seqsketchjaccard 
//!

use log::*;
#[allow(unused_imports)]
use log::Level::{Debug,Trace, Info};

use clap::{App, Arg, SubCommand};
use std::time::*;
use ::std::process;
use std::path::Path;
use std::io::prelude::*;


use kmerutils::base::{sequence::*, kmergenerator::*};
use kmerutils::sketching::*;
use kmerutils::sketching::seqblocksketch::{BlockSeqSketcher, BlockSketched, DistBlockSketched};
use needletail::FastxReader;


use hnsw_rs::prelude::*;
use hnsw_rs::api::AnnT;

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

    let _ = init_log();
    debug!("entering data sketcher, checking log trace");

    let matches = Box::new(App::new("datasketcher")       
                    .arg(Arg::with_name("file")
                        .long("file")
                        .short('f')
                        .takes_value(true)
                        .help("expecting a fastq file"))
                    .arg(Arg::with_name("sketch_size")
                        .long("sketch")
                        .short('s')
                        .takes_value(true)
                        .help("expecting sketch size as usize"))
                    .arg(Arg::with_name("kmer_size")
                        .long("kmer")
                        .short('k')
                        .takes_value(true)
                        .help("expecting a kmer size"))
                    .arg(Arg::with_name("dumpfile")
                        .long("dumpfile")
                        .short('d')
                        .takes_value(true)
                        .help("expecting name of dumpfile for signature"))
                    .arg(Arg::with_name("block_size")
                        .long("block_size")
                        .short('b')
                        .takes_value(true)
                        .help("-b for blocksize if sketching by block"))
                    .subcommand(SubCommand::with_name("ann")
                        .about("ann parameters")
                        .arg(Arg::with_name("nbng")
                                .long("nb")
                                .short('n')
                                .takes_value(true)
                                .help("expecting number of neighbours"))
                    )
                ) 
                .get_matches();

    let fname;       // fasta filename
    let kmer_size;
    let sketch_size;
    let dumpfname;   // for dump of sketches
    //
    let mut do_ann = false;
    let mut block_size : usize = 0;
    let mut sketch_block = false;
    let mut nbng = 0;        // for number of neighbours
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
    //
    if matches.is_present("sketch_size") {
        sketch_size = matches.value_of("sketch_size").ok_or("bad value").unwrap().parse::<usize>().unwrap();
        println!("got sketch_size , {}", sketch_size);
    }
    else {
        println!("--sketch is mandatory");
        println!(" usage : seqsketcher -f name --sketch (or -s)  s_size --kmer (-k) k_size  --block (-b) b_size");
        process::exit(1);
    } 
    // block size
    if matches.is_present("block_size") {
        sketch_block = true;
        block_size = matches.value_of("block_size").ok_or("bad value").unwrap().parse::<usize>().unwrap();
        println!("got block_size , {}", block_size);
    }
    // kmer options
    if matches.is_present("kmer_size") {
        kmer_size = matches.value_of("kmer_size").ok_or("bad value").unwrap().parse::<u8>().unwrap();
        println!("got kmer_size , {}", kmer_size);
    }
    else {
        println!("--kmer is mandatory");
        println!(" usage : seqsketcher -f name --sketch (or -s)  s_size --kmer (-k) k_size");
        process::exit(1
        );
    }  
    // dumpfile optins
    if matches.is_present("dumpfile") {
        dumpfname = matches.value_of("dumpfile").ok_or("bad value").unwrap().parse::<String>().unwrap();
        println!("got dumpfile  , {}", dumpfname);
    }
    else {
        println!("--dumpfile is mandatory");
        println!(" usage : seqsketcher -f name --sketch (or -s)  s_size --kmer (-k) k_size --dumfile (-d) fname");
        process::exit(1);
    }
    // ann asked for
    if let Some(ann_match) = matches.subcommand_matches("ann") {
        println!("got ann command");
        do_ann = true;
        if ann_match.is_present("nbng") {
            println!("got nbng arg");
            let nbng_decoded = ann_match.value_of("nbng").unwrap().parse::<usize>().unwrap();
            nbng = nbng_decoded as u8;
            println!("got nbng {}", nbng);
        }      
    }
    if !matches.args_present() {
        println!(" got no subcommand!");
        log::error!(" got no subcommand!");
    }

    //
    //
    let mut hnsw_opt_seq : Option<Hnsw<u32, DistHamming> > = None;
    let mut hnsw_opt_seqblock : Option<Hnsw<BlockSketched, DistBlockSketched> > = None;
    //
    if do_ann {
        // The fact is that  1. - probminhasher::compute_probminhash_jaccard(va, vb) as f32 is Hamming!
        // except for a multiplicative factor i.e the length of slices!!
        // it is also possible to use a closure defining it 
        /* let mydist_closure = | va : &[u32] , vb: &[u32] |  -> f32  {
            let mut nbdiff = 0;
            for i in 0..va.len() {
                if va[i] != vb[i] {
                    nbdiff += 1;
                }
            }
            1. - (nbdiff as f32)/va.len() as f32
        };
        let my_dist = DistFn::<u32>::new(Box::new(mydist_closure)); */
        println!("initializing hnsw");
        let max_nb_conn = 48.min(3 * nbng as usize);
        let ef_search = 200;
        log::info!("setting max nb conn to : {:?}", max_nb_conn);
        log::info!("setting ef_search to : {:?}", ef_search);
        if !sketch_block {
            hnsw_opt_seq = Some(Hnsw::<u32, DistHamming>::new(max_nb_conn , 700000, 16, ef_search, DistHamming{}));
        } 
        else {
            hnsw_opt_seqblock = Some(Hnsw::<BlockSketched, DistBlockSketched>::new(max_nb_conn , 700000, 16, ef_search, DistBlockSketched{}));
        }
    }  // end if we must do ann

    //
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
    let sequence_pack = if sketch_block { 5000 } else { 10000};
    // dumping info
    log::info!("sketching sequences by pack size {:?}", sequence_pack);
    if sketch_block {
        log::info!("sketching sequences by blocks of size {:?}", block_size);
    }
    else {
        log::info!("sketching sequences in one block each");
    }

    //
    let kmer_revcomp_hash_fn = | kmer : &Kmer32bit | -> u32 {
        let canonical =  kmer.reverse_complement().min(*kmer);
        let hashval = probminhash::invhash::int32_hash(canonical.0);
        hashval
    };
    //
    // create file to dump signature
    //
    let mut sigbuf;
    if !sketch_block {
        log::info!("allocating whole sequence SeqSketcher");
        let sketcher = seqsketchjaccard::SeqSketcher::new(kmer_size as usize, sketch_size);
        sigbuf = sketcher.create_signature_dump(&dumpfname);
    }
    else {
        let sketcher = BlockSeqSketcher::new(block_size, kmer_size as usize, sketch_size);
        sigbuf = sketcher.create_signature_dump(&dumpfname);
    }
    //
    // now we work : read, sketch by block or not, dump, and possibly embed in hnsw 
    //
    let mut nbseq = 0;
    loop {
        let sequencegroup = readblockseq(& mut reader, sequence_pack);
        if sequencegroup.len() == 0 {
            break;
        }
        // do the computation
        if !sketch_block {
            log::info!("sketching entire sequences with probminhash3a algorithm");
            let sketcher = seqsketchjaccard::SeqSketcher::new(kmer_size as usize, sketch_size);
            let sequencegroup_ref : Vec<&Sequence> = sequencegroup.iter().map(|s| s).collect();
            let signatures = sketcher.sketch_probminhash3a_compressedkmer(&sequencegroup_ref, kmer_revcomp_hash_fn);
            trace!("got nb signatures vector {} ", signatures.len());
            // dump the signature
            let resd = seqsketchjaccard::dump_signatures_block_u32(&signatures, &mut sigbuf);
            if !resd.is_ok() {
                println!("\n error occurred dumping signatures");
            }
            if do_ann {
                // insert in hnsw. Must take references and associate an id.
                let mut data_for_hnsw = Vec::<(&Vec<u32>, usize)>::with_capacity(signatures.len());
                for i in 0..signatures.len() {
                    data_for_hnsw.push((&signatures[i], nbseq+i));
                }
                hnsw_opt_seq.as_mut().unwrap().parallel_insert(&data_for_hnsw);
            }
        }
        else {
            // sketching by blocks 
            // we have sequences from [nbseq..nbseq+sequencegroup.len()] (end excluded recall it is different from Julia)
            let blocksketcher = BlockSeqSketcher::new(block_size, kmer_size as usize, sketch_size);
            // transform type from Vec<Sequence> to [(u32, &Sequence)]
            let mut tosketch = Vec::<(u32, &Sequence)>::with_capacity(sequencegroup.len());
            for i in 0..sequencegroup.len() {
                tosketch.push(( (nbseq+i) as u32, &sequencegroup[i]));
            }
            let signatures = blocksketcher.blocksketch_sequences(&tosketch, &kmer_revcomp_hash_fn);
            log::trace!("got nb signatures blocks {} ", signatures.len());
            // dump the signature
            blocksketcher.dump_blocks(&mut sigbuf, &signatures);
            if do_ann {
                // must do parallel insertion
                let nb_blocks_by_seq = 10; // ??
                // how many blocks have we ?
                let mut block_rank : usize = 0;
                // recall : the internal vector is of size 1.
                let mut data_for_hnsw = Vec::<(&Vec<BlockSketched>, usize)>::with_capacity(nb_blocks_by_seq * signatures.len());
                for i in 0..signatures.len() {
                    for j in 0..signatures[i].sketch.len() {
                        data_for_hnsw.push((&signatures[i].sketch[j], block_rank));
                        block_rank += 1;
                    }
                }
                log::debug!("sending (nb seq , nb blocks) in hnsw {:?} , {:?}", signatures.len() , block_rank);
                hnsw_opt_seqblock.as_mut().unwrap().parallel_insert(&data_for_hnsw);
            }  // end of if do_ann          
        }
        // if ann is asked for, do it now by block
        nbseq += sequencegroup.len();
        if nbseq % 1000 == 0 {
            println!(" nbseq loaded : {} ", nbseq);
        }
    }  // end whole loop
    //
    sigbuf.flush().unwrap();
    let elapsed_t = start_t.elapsed().as_secs();
    println!(" number of sequences loaded {} ", nbseq);
    println!(" elapsed time (s) in sketching [inserting in hnsw] data file {} ", elapsed_t);

    if do_ann {
        let res_dump;
        if !sketch_block {
            // dump layer information for user on stdout
            hnsw_opt_seq.as_ref().unwrap().dump_layer_info();
            // dumping hnsw
            let mut hnswname = dumpfname.clone();
            hnswname.push_str("-ann");
            println!(" dumping sketch hnsw in {:?} files", hnswname);
            res_dump = hnsw_opt_seq.as_mut().unwrap().file_dump(&hnswname);
       } else {
            // ann and sketch by block
            hnsw_opt_seqblock.as_ref().unwrap().dump_layer_info();
            let mut hnswname = dumpfname.clone();
            hnswname.push_str("-ann");
            println!(" dumping block sketch hnsw in {:?} files", hnswname);
            res_dump = hnsw_opt_seqblock.as_mut().unwrap().file_dump(&hnswname);
        }
        if res_dump.is_ok() {
            println!(" hnsw dump suceeded");
        }
        else {
            println!(" hnsw dump failed");
        }
    }  // end if do_ann
} // end main


// read a block of nbseq sequences from a fastx record
fn readblockseq(reader : &mut Box<dyn FastxReader>, nbseq : usize) -> Vec<Sequence> {
    //
    trace!("entering in readblockseq");
    let mut veqseq = Vec::<Sequence>::with_capacity(nbseq);
    let mut nb_bad_sequence = 0;
    //
    while let Some(record) = reader.next() {
        // read by block of blocksize sequence to benefit from theading of sketching
        let seqrec = record.expect("invalid record");
        let nb_bad = count_non_acgt(&seqrec.seq());
        if nb_bad > 0 {
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
    trace!("returning from readblockseq , nb seq : {} ", veqseq.len());
    //
    return veqseq;
}  // end of readblockseq

