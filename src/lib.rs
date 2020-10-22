extern crate rand;

// for logging (debug mostly, switched at compile time in cargo.toml)
#[macro_use]
extern crate lazy_static;



// basic stuff
pub mod alphabet;
pub mod sequence;
pub mod statutils;
pub mod io;
pub mod quality;
pub mod parsearg;


// quality client-sever stuff

pub mod qualclient;
pub mod qserverclient;

// hashing stuff

pub mod minhash;
pub mod nohasher;
pub mod mininvhashkmer;
#[macro_use]
pub mod nthash;


pub mod seqsketcher;

// kmer stuff

pub mod kmer;
pub mod kmergenerator;
pub mod kmercount;

// contig generation

pub mod anchor;
pub mod redisbase;

pub mod prelude;


lazy_static! {
    #[allow(dead_code)]
    pub static ref LOG: u64 = {
        let res = init_log();
        res
    };
}
// install a logger facility
// set RUST_LOG to trace, warn debug off ....
fn init_log() -> u64 {
    env_logger::Builder::from_default_env().init();
    println!("\n ************** initializing logger from env *****************\n");    
    return 1;
}
