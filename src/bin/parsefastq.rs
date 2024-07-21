//! fastq parser contains subcommands to compute stats/ or counting kmers.   
//! usage parsefastq -f datafile --bits m (or -b m)  [kmer --count n | -b n] [--csize l] [-t j]
//! -   --bits (or -b)   m : m gives the number of bits per base 2,4 or 8 (8 means no compression)  
//! -   --csize (or -c)  l : size of bloom filter   
//! -   --thread (or -t) j : j gives the number of threads to use  
//! -   --unique (or -u)   : option to get a dump of unique kmers (default is nodump)  




use clap::{Arg, Command, ArgAction};





// general use
#[doc(no_inline)]
use std::process;


// usage parsefastq -f datafile --bits m (or -b m)  [kmer --count n | -b n] [--csize l] [-t j]
// subcommand are given as module call without any -, options are given with -- for long syntax or - for short syntax 



// our modules
use kmerutils::io::*;
use kmerutils::base::{kmercount::*, alphabet::get_ac_from_tg};
use kmerutils::parsearg::*;
use kmerutils::statutils::*;


// for logging (debug mostly, switched at compile time in cargo.toml)
use env_logger::{Builder};

// install a logger facility
fn init_log() -> u64 {
    Builder::from_default_env().init();
    println!("\n ************** initializing logger *****************\n");    
    1
}


fn main() {
    
    let _ = init_log();


    let mut parse_args: ParseFastqArgs = Default::default();
    let mut ret_times_args:ReturnTimesArgs = Default::default();
    let mut kmer_count_args:KmerArgs = Default::default();
    
    let matches = Command::new("parsefastq")
        .arg(Arg::new("bits")
            .long("bits")
            .short('b')
            .value_parser(clap::value_parser!(u8))
            .action(ArgAction::Set)
            .default_missing_value("8")
            .help("nb bits by base : can be 2 , 4 or 8 (compression 4 , 2 , 1)"))
        .arg(Arg::new("file")
            .long("file")
            .short('f')
            .required(true)
            .value_parser(clap::value_parser!(String))
            .action(ArgAction::Set)
            .help("expecting a fastq file"))
        .subcommand(Command::new("ret")
                .about("return times option")
                .arg(Arg::new("base")
                    .long("base")
                    .short('b')
                    .value_parser(clap::value_parser!(char))
                    .action(ArgAction::Set)
                    .help("ask for return times for a given base AT/CG"))
        )
        .subcommand(Command::new("kmer")
            .about("kmer options")
            .arg(Arg::new("count")
                    .long("count")
                    .short('c')
                    .help("kmer subcommand : option to count kmer"))
            .arg(Arg::new("thread")
                    .short('t')
                    .value_parser(clap::value_parser!(usize))
                    .action(ArgAction::Set)
                    .help(" to tell number of thread to be used, -t n"))
            .arg(Arg::new("ksize")
                    .long("ksize")
                    .short('s')
                    .value_parser(clap::value_parser!(usize))
                    .action(ArgAction::Set)
                    .help(" to tell size of kmer to generate, -s n"))
            .arg(Arg::new("csize")
                    .long("csize")
                    .short('c')
                    .value_parser(clap::value_parser!(usize))
                    .action(ArgAction::Set)
                    .help(" to tell size in bits (8 or 16) of counter to use in bloom filter, -s n"))                    
            .arg(Arg::new("unique")
                    .long("unique")
                    .short('u')
                    .help("kmer subcommand : option to determine unique kmer"))
    )
    .get_matches();
    //

    parse_args.nb_bits_by_base = 8;

    if matches.contains_id("bits") {
        println!("got bits option");
        let nb_bits = *matches.get_one::<u8>("bits").unwrap_or(&8u8);
        match nb_bits {
            2 | 4 | 8 => println!("setting nb_bits to {}", nb_bits),
            _  => { println!("setting nb_bits to  possible values : 2, 4 or 8");
                    std::panic!(" bad value for nb_bitssetting nb_bits to  possible values : 2, 4 or 8");
                  },
        }
        println!("got comp option , comp = {}", nb_bits);
        parse_args.nb_bits_by_base = nb_bits;
    }
    if matches.contains_id("file") {
        println!("file option");
        let fname = matches.get_one::<String>("file").unwrap();
        println!("got filename , {}", fname);
        parse_args.filename = fname.clone();
    }
    else {
        println!("-f filename is mandatory");
        process::exit(1);
    }

    if let Some(kmer_match) = matches.subcommand_matches("kmer") {
        println!("got kmer subcommand");
        if kmer_match.contains_id("count") {
            println!(" got subcommand count option");
            kmer_count_args.kmer_task = KmerProcessing::Counting;
        }
        if kmer_match.contains_id("thread") {
            let nb_threads = *kmer_match.get_one::<usize>("thread").unwrap();
            kmer_count_args.nb_threads = nb_threads;
        }
        if kmer_match.contains_id("ksize") {
            let k_size = *kmer_match.get_one::<usize>("ksize").expect("expecting kmer size");
            kmer_count_args.kmer_size = k_size;
        }                        
        if kmer_match.contains_id("csize") {
            let c_size = *kmer_match.get_one::<usize>("csize").expect("expect nb bits for bloom counter");
            kmer_count_args.counter_size = c_size;
        }
        if kmer_match.contains_id("unique") {
            println!(" got subcommand unicity option");
            kmer_count_args.kmer_task = KmerProcessing::Unicity;
        }
        parse_args.kmer_args = kmer_count_args;
    };
        // for return times command
    if let Some(ret_match) =  matches.subcommand_matches("ret") { 
        println!("got ret subcommand");
        ret_times_args.to_do = true;
        if ret_match.contains_id("base") {
            println!(" got subcommand ret base affectation");
            let base_char = *ret_match.get_one::<char>("base").expect("base A or C");
            let base = get_ac_from_tg(base_char as u8);
            if base != b'A' && base != b'C' {
                println!(" bad base for return times analyze : {} ", base as char);
                println!(" base must be b'A' or b'C' ");
                std::process::exit(1);
            }
            else {
                ret_times_args.searched_base = base;       
            }
        }
        parse_args.ret_times_args_opt = Some(ret_times_args);
    };
     
    if !matches.args_present() {
        println!(" got no subcommand!");
        log::error!(" got no subcommand!");
    }

    
    
    let seqvec = match parse_with_needletail(parse_args.clone()) {
        Ok(sthing) => sthing,
        Err(e) => {
            println!("coud not parse file : {}", e);
            process::exit(1);
        }
    };
    //
    println!("got nb read : {}  ", seqvec.len());
    // for long read data a maxreadlen = 1000000 should be OK.
    let maxreadlen = 1000000;
    log::info!("setting max read length for size histogram to : {}", maxreadlen);
    let prec = 3;
    let base_distribution_res= get_base_count_par(&seqvec, maxreadlen, prec);
    match base_distribution_res {
        Some(base_distribution) => { let _res= base_distribution.ascii_dump_readlen_distribution("readlen.histo");
                                    },
        _                       => std::process::exit(1),
    }
    //
    match parse_args.kmer_args.kmer_task {
        KmerProcessing::Counting => {
            let mut multi_kmer_file_name = parse_args.filename.clone();
            multi_kmer_file_name.push_str(".multi_kmer.bin");
            if let Some(pos) = multi_kmer_file_name.rfind('/') {
                multi_kmer_file_name = multi_kmer_file_name.split_off(pos+1);
            }
            println!("dumping multiple kmers in file : {} ", multi_kmer_file_name);
            //
            if parse_args.kmer_args.kmer_size <= 32 && parse_args.kmer_args.kmer_size > 16 {
                let kmer_counter : KmerCounterPool<Kmer64bit> = count_kmer_threaded_one_to_many(&seqvec,
                                                                                                   parse_args.kmer_args.nb_threads,
                                                                                                   parse_args.kmer_args.counter_size,
                                                                                                   parse_args.kmer_args.kmer_size);
                let _res = threaded_dump_kmer_counter(&kmer_counter, &multi_kmer_file_name, &seqvec,
                                                      parse_args.kmer_args.nb_threads, parse_args.kmer_args.kmer_size);
            }
            else if parse_args.kmer_args.kmer_size == 16 {
                let kmer_counter : KmerCounterPool<Kmer16b32bit> = count_kmer_threaded_one_to_many(&seqvec,
                                                                                                   parse_args.kmer_args.nb_threads,
                                                                                                   parse_args.kmer_args.counter_size,
                                                                                                   parse_args.kmer_args.kmer_size);
                let _res = threaded_dump_kmer_counter(&kmer_counter, &multi_kmer_file_name, &seqvec,
                                                      parse_args.kmer_args.nb_threads, parse_args.kmer_args.kmer_size);
            }
            else if parse_args.kmer_args.kmer_size <= 14 {
                let kmer_counter : KmerCounterPool<Kmer32bit> = count_kmer_threaded_one_to_many(&seqvec,
                                                                                                parse_args.kmer_args.nb_threads,
                                                                                                parse_args.kmer_args.counter_size,
                                                                                                parse_args.kmer_args.kmer_size);
                let _res = threaded_dump_kmer_counter(&kmer_counter, &multi_kmer_file_name, &seqvec,
                                                      parse_args.kmer_args.nb_threads, parse_args.kmer_args.kmer_size);
            };
        },
        KmerProcessing::Unicity => {
            let kmer_counter = filter1_kmer_16b32bit(&seqvec);
            let mut kmer1_file_name = parse_args.filename.clone();
            kmer1_file_name.push_str(".once_kmer.bin");
            if let Some(pos) = kmer1_file_name.rfind('/') {
                kmer1_file_name = kmer1_file_name.split_off(pos+1);
            }
            println!("dumping unique kmers in file : {} ", kmer1_file_name);
            let _res = kmer_counter.dump_in_file_once_kmer16b32bit(&kmer1_file_name, &seqvec);
        }
        KmerProcessing::None => (),
    }
    

    
//    generate_all_kmer32(&seqvec);
    
}
    
