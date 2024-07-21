use std::process;

#[macro_use]
extern crate lazy_static;

use env_logger::{Builder};

// for logging (debug mostly, switched at compile time in cargo.toml)

use clap::{Arg, Command, ArgAction};

use std::path::Path;

// our modules
extern crate kmerutils;

use kmerutils::base::kmercount::*;

lazy_static! {
    #[allow(dead_code)]
    static ref LOG: u64 = {
        
        init_log()
    };
}

// install a logger facility
fn init_log() -> u64 {
    Builder::from_default_env().init();
    println!("\n ************** initializing logger *****************\n");    
    1
}



fn main() {
    //
    // the reference to LOG will force the call to lazy_static! call to init_log to get LOG initialized.
    //
    if *LOG != 1 {
        println!(" LOG = {:?}", *LOG);
    }
    //
    let matches = Command::new("reloadkmermulti")
        .arg(Arg::new("file")
            .long("file")
            .short('f')
            .required(true)
            .action(ArgAction::Set)
            .value_parser(clap::value_parser!(String))
            .help("expecting dumped file .bin"))
        .get_matches();


    let filename = matches.get_one::<String>("file").unwrap();
    println!("got filename , {}", filename);

    {
        let path = Path::new(filename);
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
    let _res = KmerCountReload::load_multiple_kmers_from_file(filename);
}
