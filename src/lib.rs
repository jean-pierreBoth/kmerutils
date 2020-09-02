extern crate rand;

// for logging (debug mostly, switched at compile time in cargo.toml)
extern crate log;
extern crate simple_logger;

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
pub mod superminhash; 
pub mod nohasher;
pub mod mininvhashkmer;
#[macro_use]
pub mod nthash;

pub mod invhash;

pub mod seqsketcher;

// kmer stuff

pub mod kmer;
pub mod kmergenerator;
pub mod kmercount;

// contig generation

pub mod anchor;
pub mod redisbase;

pub mod prelude;

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    // initialize once log system for tests.
    fn init_log() {
        let _res = simple_logger::init();
    }
}  // end of tests
