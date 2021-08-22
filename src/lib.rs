extern crate rand;

// for logging (debug mostly, switched at compile time in cargo.toml)
#[macro_use]
extern crate lazy_static;



// basic stuff bases definitions, kmer , kmer generation sequences,
pub mod base;



pub mod statutils;
pub mod io;
pub mod parsearg;


// quality stuff
pub mod quality;

// hashing stuff

pub mod hashed;
pub mod nohasher;

// sketching methods
pub mod sketching;


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
