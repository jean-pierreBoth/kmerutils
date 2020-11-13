//! fastq parser , a kind of driver for subcommands.
//! contains subcommands to compute stats/ or counting kmers. 
//! usage parsefastq -f datafile --bits m (or -b m)  [kmer --count n | -b n] [--csize l] [-t j]
//! - -b  m : m gives the number of bits per base 2,4 or 8 (8 means no compression)
//! - -t  j : j gives the number of threads to use




use clap::{App, Arg, SubCommand};





// general use
#[doc(no_inline)]
use std::process;


// usage parsefastq -f datafile --bits m (or -b m)  [kmer --count n | -b n] [--csize l] [-t j]
// subcommand are given as module call without any -, options are given with -- for long syntax or - for short syntax 



// our modules
use kmerutils::io::*;
use kmerutils::kmercount::*;
use kmerutils::parsearg::*;
use kmerutils::statutils::*;
use kmerutils::alphabet::get_ac_from_tg;    

fn main() {
    
    if cfg!(verbose_1 = "1") {
        println!(" appel main : ");
    }

    let mut parse_args: ParseFastqArgs = Default::default();
    let mut ret_times_args:ReturnTimesArgs = Default::default();
    let mut kmer_count_args:KmerArgs = Default::default();
    
    let matches = App::new("parsefastq")
        .arg(Arg::with_name("bits")
             .long("bits")
             .short("b")
             .takes_value(true)
             .default_value("8")
             .help("nb bits by base : can be 2 , 4 or 8 (compression 4 , 2 , 1)"))
        .arg(Arg::with_name("file")
             .long("file")
             .short("f")
             .takes_value(true)
             .help("expecting a fastq file"))
        .subcommand(SubCommand::with_name("ret")
                    .about("return times option")
                    .arg(Arg::with_name("base")
                         .takes_value(true)
                         .long("base")
                         .short("b")
                         .help("ask for return times for a given base AT/CG"))
        )
        .subcommand(SubCommand::with_name("kmer")
                    .about("kmer options")
                    .arg(Arg::with_name("count")
                         .long("count")
                         .short("c")
                         .takes_value(false)
                         .help("kmer subcommand : option to count kmer"))
                    .arg(Arg::with_name("thread")
                         .short("t")
                         .takes_value(true)
                         .help(" to tell number of thread to be used, -t n"))
                    .arg(Arg::with_name("ksize")
                         .long("ksize")
                         .short("s")
                         .takes_value(true)
                         .help(" to tell size of kmer to generate, -s n"))
                    .arg(Arg::with_name("csize")
                         .long("csize")
                         .short("c")
                         .takes_value(true)
                         .help(" to tell size of counter to use in bloom filter, -s n"))                    
                    .arg(Arg::with_name("unique")
                         .long("unique")
                         .short("u")
                         .takes_value(false)
                         .help("kmer subcommand : option to determine unique kmer"))
        )
        .get_matches();
    //

    let mut nb_bits: u8 = 8;
    parse_args.nb_bits_by_base = nb_bits;
    let fname;

    if matches.is_present("bits") {
        println!("got bits option");
        nb_bits = matches.value_of("bits").ok_or("bad value").unwrap().parse::<u8>().unwrap_or(8 as u8);
        match nb_bits {
            2 | 4 | 8 => println!("setting nb_bits to {}", nb_bits),
            _  => println!("setting nb_bits to : possible values : 2, 4 or 8"),
        }
        println!("got comp option , comp = {}", nb_bits);
        parse_args.nb_bits_by_base = nb_bits;
    }
    if matches.is_present("file") {
        println!("file option");
        fname = matches.value_of("file").ok_or("bad value").unwrap().parse::<String>().unwrap();
        println!("got filename , {}", fname);
        parse_args.filename = fname;
    }
    else {
        println!("-f filename is mandatory");
        process::exit(1);
    }

    match matches.subcommand() {
        ("kmer", Some(kmer_match)) => {
            println!("got kmer subcommand");
            if kmer_match.is_present("count") {
                println!(" got subcommand count option");
                kmer_count_args.kmer_task = KmerProcessing::Counting;
            }
            if kmer_match.is_present("thread") {
                let nb_threads = kmer_match.value_of("thread").unwrap().parse::<usize>().unwrap();
                kmer_count_args.nb_threads = nb_threads;
            }
            if kmer_match.is_present("ksize") {
                let k_size = kmer_match.value_of("ksize").unwrap().parse::<usize>().unwrap();
                kmer_count_args.kmer_size = k_size;
            }                        
            if kmer_match.is_present("csize") {
                let c_size = kmer_match.value_of("csize").unwrap().parse::<usize>().unwrap();
                kmer_count_args.counter_size = c_size;
            }
            if kmer_match.is_present("unique") {
                println!(" got subcommand unicity option");
                kmer_count_args.kmer_task = KmerProcessing::Unicity;
            }
            parse_args.kmer_args = kmer_count_args;
        },
        // for return times command
        ("ret", Some(ret_match)) => {
            println!("got ret subcommand");
            ret_times_args.to_do = true;
            if ret_match.is_present("base") {
                println!(" got subcommand ret base affectation");
                let base_char = ret_match.value_of("base").unwrap().parse::<char>().unwrap();
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
        },
     
        ("", None)               => println!("no subcommand at all"),
        _                        => unreachable!(),
    }

    let seqvec;
    match parse_with_needletail(parse_args.clone()) {
        Ok(sthing) => seqvec = sthing,
        Err(e) => {
            println!("coud not parse file : {}", e);
            process::exit(1);
        }
    }
    //
    println!("got nb read : {}  ", seqvec.len());
    //
    let maxreadlen = 1000000;
    let prec = (maxreadlen as f64).log10() as usize;
    let base_distribution_res= get_base_count_par(&seqvec, maxreadlen, prec);
    match base_distribution_res {
        Some(base_distribution) => { let _res= base_distribution.ascii_dump_readlen_distribution(&"readlen.histo");
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
                let kmer_counter : KmerCounterPool<Kmer16b32bit> = count_kmer_threaded_one_to_many(&seqvec,
                                                                                                   parse_args.kmer_args.nb_threads,
                                                                                                   parse_args.kmer_args.counter_size,
                                                                                                   parse_args.kmer_args.kmer_size);
                let _res = threaded_dump_kmer_counter(&kmer_counter, &multi_kmer_file_name, &seqvec,
                                                      parse_args.kmer_args.nb_threads as usize, parse_args.kmer_args.kmer_size);
            }
            else if parse_args.kmer_args.kmer_size == 16 {
                let kmer_counter : KmerCounterPool<Kmer16b32bit> = count_kmer_threaded_one_to_many(&seqvec,
                                                                                                   parse_args.kmer_args.nb_threads,
                                                                                                   parse_args.kmer_args.counter_size,
                                                                                                   parse_args.kmer_args.kmer_size);
                let _res = threaded_dump_kmer_counter(&kmer_counter, &multi_kmer_file_name, &seqvec,
                                                      parse_args.kmer_args.nb_threads as usize, parse_args.kmer_args.kmer_size);
            }
            else if parse_args.kmer_args.kmer_size <= 14 {
                let kmer_counter : KmerCounterPool<Kmer32bit> = count_kmer_threaded_one_to_many(&seqvec,
                                                                                                parse_args.kmer_args.nb_threads,
                                                                                                parse_args.kmer_args.counter_size,
                                                                                                parse_args.kmer_args.kmer_size);
                let _res = threaded_dump_kmer_counter(&kmer_counter, &multi_kmer_file_name, &seqvec,
                                                      parse_args.kmer_args.nb_threads as usize, parse_args.kmer_args.kmer_size);
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
    
