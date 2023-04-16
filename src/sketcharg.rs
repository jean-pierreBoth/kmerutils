//! describe sketching paramaters, kmer size, sketching size and sketching algorithms.


use std::io::{BufReader, BufWriter };

use std::fs::OpenOptions;
use std::path::{Path, PathBuf};


use serde::{Deserialize, Serialize};
use serde_json::{to_writer};


/// specify if we process DNA sequence or AA sequences
#[derive(Copy,Clone,Serialize,Deserialize,Debug)]
pub enum DataType {
    DNA,
    AA,
}

/// Specify which algo we use for sketching :  Probminhash or SuperMinHash or Hyperloglog (SetSketch) algorithms.  
/// - PROB3A is the value for asking ProbMinHashh3a, 
/// - SUPER for first version SuperMinHash.(f64 signature)
/// - SUPER2 for second version of superminhash (u64 signature)
/// - HLL for SetSketch based on hyperloglog (u16, u32 or more signature)
#[derive(Copy,Clone,Serialize,Deserialize,Debug)]
pub enum SketchAlgo {
    PROB3A,
    SUPER,
    SUPER2,
    HLL,
}
// This is redundant with struct Sketcher for DNA case and RNA case, but it makes
// possible the factorization of all parameters

/// describe sketching paramaters, kmer size, sketching size and sketching algorithms.
/// 
#[derive(Copy,Clone,Serialize,Deserialize,Debug)]
pub struct SeqSketcherParams {
    kmer_size : usize,
    sketch_size : usize,
    algo : SketchAlgo, 
    data_t : DataType
}


impl SeqSketcherParams {
    /// 
    pub fn new(kmer_size: usize, sketch_size : usize, algo : SketchAlgo, data_t: DataType) -> Self {
        SeqSketcherParams{kmer_size, sketch_size, algo, data_t}
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

    /// returns AA or DNA type for data
    pub fn get_data_t(&self) -> DataType {
        self.data_t
    }

    /// serialized dump
    pub fn dump_json(&self, filename : &String) -> Result<(), String> {
        //
        let filepath = PathBuf::from(filename.clone());
        //
        log::info!("dumping sketching parameters in json file : {}", filename);
        //
        let fileres = OpenOptions::new().write(true).create(true).truncate(true).open(&filepath);
        if fileres.is_err() {
            log::error!("SeqSketcher dump : dump could not open file {:?}", filepath.as_os_str());
            println!("SeqSketcher dump: could not open file {:?}", filepath.as_os_str());
            return Err("SeqSketcher dump failed".to_string());
        }
        // 
        let mut writer = BufWriter::new(fileres.unwrap());
        let _ = to_writer(&mut writer, &self).unwrap();
        //
        Ok(())
    } // end of dump


    /// reload from a json dump
    pub fn reload_json(dirpath : &Path) -> Result<SeqSketcherParams, String> {
        log::info!("in reload_json");
        //
        let filepath = dirpath.join("sketchparams_dump.json");
        let fileres = OpenOptions::new().read(true).open(&filepath);
        if fileres.is_err() {
            log::error!("Sketcher reload_json : reload could not open file {:?}", filepath.as_os_str());
            println!("Sketcher reload_json: could not open file {:?}", filepath.as_os_str());
            return Err("Sketcher reload_json could not open file".to_string());            
        }
        //
        let loadfile = fileres.unwrap();
        let reader = BufReader::new(loadfile);
        let sketch_params:SeqSketcherParams = serde_json::from_reader(reader).unwrap();
        //
        log::info!("SeqSketcher reload, kmer_size : {}, sketch_size : {}", 
            sketch_params.get_kmer_size(), sketch_params.get_sketch_size());     
        //
        Ok(sketch_params)
    } // end of reload_json

}  // end of SeqSketcherParams


//==========================================================================================