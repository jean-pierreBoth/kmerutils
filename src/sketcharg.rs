// describe sketching paramaeters


use serde::{Deserialize, Serialize};

/// Specify of we skecth using the Probminhash or SuperMinHash algorithms
#[derive(Copy,Clone,Serialize,Deserialize)]
pub enum SketchAlgo {
    PROB3A,
    SUPER,
}
// This is redundant with struct Sketcher for DNA case and RNA case, but it makes
// possible the factorization of all parameters

#[derive(Copy,Clone,Serialize,Deserialize)]
pub struct SeqSketcherParams {
    kmer_size : usize,
    sketch_size : usize,
    algo : SketchAlgo,  
}


impl SeqSketcherParams {
    /// 
    pub fn new(kmer_size: usize, sketch_size : usize, algo : SketchAlgo) -> Self {
        SeqSketcherParams{kmer_size, sketch_size, algo}
    }

    /// returns kmer size
    pub fn get_kmer_size(&self) -> usize {
        self.kmer_size
    }

    /// return sketch size
    pub fn get_sketch_size(&self) -> usize {
        self.sketch_size
    }  

    /// get get sketching algorithm PROB or SUPER
    pub fn get_algo(&self) -> SketchAlgo {
        self.algo
    }
}  // end of SeqSketcherParams

//==========================================================================================