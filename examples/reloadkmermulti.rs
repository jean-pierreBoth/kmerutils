use std::process;

#[macro_use]
extern crate lazy_static;

// for logging (debug mostly, switched at compile time in cargo.toml)

use clap::{App, Arg};

use std::path::Path;

// our modules
extern crate kmerutils;

use kmerutils::kmercount::*;

lazy_static! {
    #[allow(dead_code)]
    static ref LOG: u64 = {
        let res = init_log();
        res
    };
}

// install a logger facility
fn init_log() -> u64 {
    simple_logger::init().unwrap();
    println!("\n ************** initializing logger *****************\n");    
    return 1;
}



fn main() {
    //
    // the reference to LOG will force the call to lazy_static! call to init_log to get LOG initialized.
    //
    if *LOG != 1 {
        println!(" LOG = {:?}", *LOG);
    }
    //
    let filename;
    //
    let matches = App::new("reloadkmermulti")
        .arg(Arg::with_name("file")
             .long("file")
             .short("f")
             .takes_value(true)
             .help("expecting dumped file .bin"))
        .get_matches();


    if matches.is_present("file") {
        filename = matches.value_of("file").ok_or("bad value").unwrap().parse::<String>().unwrap();
        println!("got filename , {}", filename);
    }
    else {
        println!("-f filename is mandatory");
        println!(" usage qualityloade -f name");
        process::exit(1);
    }
    {
        let path = Path::new(&filename);
        let f_info_res = path.metadata();
        let _filesize:u64;
        match f_info_res {
            Ok(meta) => {
                _filesize = meta.len();           
            },
            Err(_e) => {
                println!("file does not exist: {:?}", filename);
                process::exit(1);
            },
        }
    }
    let _res = KmerCountReload::load_multiple_kmers_from_file(&String::from(filename));
}
