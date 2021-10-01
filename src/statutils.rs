//! This module computes basic statistics on base distribution, and read length distribution.

// It gathers some heuristics on how to pairs blocks of read before going to sketching


use ::histogram::*;


use ::ndarray::*;


// rayon usage
use ::rayon::prelude::*;


// general use
use ::std::result;
use ::std::fs::OpenOptions;
use ::std::io::prelude::*;
//
use std::io;

// our includes

use crate::base::sequence::*;


// should be allocated via ndarray::Dim([sz as usize, 4]) so a sz by 4 matrix
// rust use row order as opposed to fortran which is row column

/// We keep in ReadBaseDistribution main statistics :
///  - readsizehisto  
///    an histogram of read length distribution
///   
///  - acgt_distribution  
///  For each  percentage and each base we store the fraction of reads in acgt_distribution.  
///  So  acgt_distribution[percent, cbase] = f means there is a fraction f of reads where cbase occurs at level percent%.   
///  base columns are in order a,c,g,t. For a given column we have the fraction of reads in which the base occurs at each
///  percentage. 
/// 

pub struct ReadBaseDistribution {
    pub readsizehisto : histogram::Histogram,
    /// upper value of histo
    upper_histo : usize,
    /// count number of values greater than upper , histo trackable value
    pub histo_out : usize,
    pub non_acgt : usize,
    /// array of dimensions 101,4.  101 lines for percentage can be in 0..100 included!  
    /// columns are A C G T in thid order
    pub acgt_distribution : Array2<f64>,
}



impl ReadBaseDistribution {
    /// readmaxsize is the maximum length represented in the histogram  
    /// prec : is the precision for histogram representation1,2 or 3
    pub fn new(readmaxsize : usize, prec : usize, d: ndarray::Ix2) -> ReadBaseDistribution {
        // TODO check that prec < log10(readmaxsize)
        // record histogram length from 1 to readmaxsize with slot of size readmaxsize/10**prec
        // lowest value arg in init must be >= 1
        assert!(prec >= 1, "precision for histogram construction should range >= 1");
        let histo = Histogram::configure().max_value(readmaxsize as u64).precision(prec as u32).build().unwrap();
        let array : Array2<f64> = Array2::zeros(d);
        ReadBaseDistribution{readsizehisto : histo, upper_histo : readmaxsize, histo_out : 0,
                             non_acgt : 0 as usize , acgt_distribution : array}
    }

    
    /// dump result in file of name name.
    /// first line is number of non acgt
    /// Then for each slot from 0 to 100 include , we get number of reads
    /// where a base (a,c,g,t) occurs at a given percentage.
    /// format is 4 column for bases in order a,c,g,t
    /// This file can be reloaded by Julia module BaseDistribution
    
    pub fn ascii_dump_acgt_distribution(&self, name : &String) -> result::Result<(), io::Error> {
        // open filename
        let fileres = OpenOptions::new().write(true).create(true).truncate(true).open(name);
        match fileres {
            Ok(mut file) => {
                //
                let (nbrow,_) = self.acgt_distribution.dim();
                for i in 0..nbrow {
                    write!(file, "{} {}  {}  {} \n", self.acgt_distribution[[i,0]] , self.acgt_distribution[[i,1]] ,
                           self.acgt_distribution[[i,2]] , self.acgt_distribution[[i,3]])?;  // ? get directly to io::Error !
                }
                file.flush()?;
                return Ok(());
            },
            Err(e) => {
                println!("could not open file {}", name);
                return Err(e);
            }
        } // end match
    } // end of ascii_dump_acgt_distribution


    /// dump a list of 1000 couples (readsize, nbreads)
    /// read length are are obtained at regular quantiles and 
    /// abscisses (read length) for j = 0..1000 are computed as first read length greater than maxlen * j / nbpoints
    pub fn ascii_dump_readlen_distribution(&self, name : &str) -> result::Result<(), io::Error> {
        //
        let nb_entries = self.readsizehisto.entries() as usize;
        let maxlen : usize;
        match self.readsizehisto.maximum() {
           Ok(maxlen_ok ) =>  { 
                                maxlen = maxlen_ok as usize;
                                println!("ascii_dump_readlen_distribution nb_entries {}  maxlen {}", nb_entries, maxlen);},
           Err(str) => {    println!("error : {}", String::from(str)); 
                            println!("Error : ascii_dump_readlen_distribution nb_entries {}", nb_entries);
                            return Err(std::io::Error::new(std::io::ErrorKind::Other, "histogram error!"));
                        }
        };
        //
        if nb_entries < 100 {
            println!("Error : ascii_dump_readlen_distribution nb_entries too small : {}", nb_entries);
            return Err(std::io::Error::new(std::io::ErrorKind::Other, "histogram error!"));
        }
        let nbslot = nb_entries / 100;
        //
        let mut readsize : Vec<u64> = (0..(nbslot+1)).map(|_| 0u64).collect();
        for i in 0..(nbslot+1) {
            readsize[i] =  self.readsizehisto.percentile(100. * (i as f64/ nbslot as f64) as f64).unwrap();
        }
        let nb_points = 1000usize;
        let mut nb_read_vec = Vec::<usize>::with_capacity(nb_points);
        let mut abscisse = Vec::<usize>::with_capacity(nb_points);
        //
        let mut first_i = 0;
        let mut current_i = 0;
        for j in 0..nb_points {
            // find first i such that readsize[i] >= maxlen * (j/nbslot) 
            let threshold = (maxlen * j)/nb_points;
            while readsize[current_i]  < (threshold as u64)  && current_i < nbslot {
                current_i +=1;
            }
            if current_i < nbslot && current_i > first_i {
                let nb_in_slot = ((current_i - first_i) * nb_entries as usize) / nbslot; 
                nb_read_vec.push(nb_in_slot);
                abscisse.push(readsize[current_i] as usize);
            }
            first_i = current_i;
        }
        //
        let fileres = OpenOptions::new().write(true).create(true).truncate(true).open(name);
        match fileres {
            Ok(mut file) => {
                //
                for i in 0..nb_read_vec.len() {
                    write!(file, "{}  {} \n",abscisse[i] , nb_read_vec[i])?;  // ? get directly to io::Error !
                }
                file.flush()?;
                return Ok(());
            },
            Err(e) => {
                println!("could not open file {}", name);
                return Err(e);
            }
        } // end match
    } // end of ascii_readlen_distribution

    /// record read len in histogram keeping track of number of values outside histogram
    fn record_read_len(&mut self, sz : usize) -> () {
        log::trace!("record_read_len sz : {}", sz);
        let res = self.readsizehisto.increment(sz as u64);
        if !res.is_ok() {
            self.histo_out += 1;
        }
    } // end of record_read_len

} // end of implementation for BaseDistribution



/// computes base distribution statistics(in one thread) and dump result in file of name "bases.histo-1thread"
/// maxreadlen is maximum read length in histogram of read length distribution.
#[allow(dead_code)]
fn get_base_count(seq_array : &Vec<Sequence > , maxreadlen:usize) {
    //
    println!(" in get_base_count");
    let start_t = std::time::Instant::now();
    //
    let prec = (maxreadlen as f64).log10() as usize;
    let v_ref = &seq_array;
    let low = 0;
    let up = v_ref.len();
    let mut base_distribution = get_base_count_by_slice(v_ref.get(low..up).unwrap(), maxreadlen, prec).unwrap();       
    (*base_distribution).acgt_distribution *= 1./(seq_array.len() as f64);
    //
    let elapsed_t = start_t.elapsed().as_secs();
    println!(" elapsed time (s) in get_base_count {} ", elapsed_t);
    //
    base_distribution.ascii_dump_acgt_distribution(&String::from("bases.histo-1thread")).unwrap();
} // end of get_base_count



/// takes a slice on sequences and sendback base distribution on a channel
/// maxreadlen is the maximum size we record in histogram. It must be adapted to read length distribution
/// to get a significant historgram


fn get_base_count_by_slice(seq_array : &[Sequence] , maxreadlen : usize, prec : usize) -> result::Result< Box<ReadBaseDistribution> , ()> {
    //
    log::trace!(" in get_base_count_by_slice len : {} ", seq_array.len());
    //
    let start_t = std::time::Instant::now();
    // allocate a structure to store percentages, we want to cover an interval [0..100] so sz = 101
    let sz : usize = 101;

    let mut base_distribution : Box<ReadBaseDistribution>  =
        Box::new(ReadBaseDistribution::new(maxreadlen , prec, ndarray::Dim([sz as usize, 4])) );
    //
    let mut histo : Vec<u64> = vec![0,0,0,0];
    let mut nb_bad = 0;
    for i in 0..seq_array.len() {
        let decompressed = seq_array[i].decompress();
        base_distribution.record_read_len(decompressed.len());
        nb_bad = nb_bad + seq_array[i].base_count(&mut histo);
        for j in 0..4 {
            // compute percentage of each base in sequence
            let percent : usize = ((100 * histo[j]) as f64 / decompressed.len() as f64).round() as usize;
            (*base_distribution).acgt_distribution[ [percent ,j ] ] += 1.;
        }     
    }
    (*base_distribution).non_acgt = nb_bad as usize;
    //
    let elapsed_t = start_t.elapsed().as_secs();
    println!(" elapsed time (s) in get_base_count {} ", elapsed_t);
    //
    if nb_bad > 0 {
        log::trace!(" out  get_base_count_par_slice nb_bad : {} ", nb_bad);
    }
    //
    Ok(base_distribution)
}





/// takes a slice on sequences and sendback base distribution on a channel
/// maxreadlen is the maximum size we record in histogram. It must be adapted to read length distribution
/// to get a significant historgram.    
/// Result is dumped in a file named : "bases.histo"
pub fn get_base_count_par (seq_array : &Vec<Sequence> ,  maxreadlen :usize, prec : usize) -> Option<Box<ReadBaseDistribution>> {
    //
    log::info!(" in get_base_count_par");
    let nbthreads : usize = 2;
    let nb_cpus = num_cpus::get_physical();
    log::info!(" in get_base_count_par,  number of cpus found : {} " , nb_cpus);
    let start_t = std::time::Instant::now();
    //
    // allocate nbthreads BaseDistribution.
    let v_ref = &seq_array;   // Vec n implemente pas Copy mais Ref oui!

    let distrib_collector : Vec<result::Result< Box<ReadBaseDistribution> , () >  > = (0..nbthreads).into_par_iter().map(|i|  {
        let low = (v_ref.len() / nbthreads)  * i;
        let up = if i < nbthreads-1 {
                (v_ref.len() / nbthreads)  * (i as usize + 1)}
            else { 
                v_ref.len()
            };
        log::trace!(" in get_base_count_par thread {} , low up {} {} ", i, low , up);
        let res: result::Result< Box<ReadBaseDistribution> , () > = get_base_count_by_slice(v_ref.get(low..up).unwrap(), maxreadlen , prec);
        res
    }).collect();

   
    // now we can reduce result summing all result in first 
    let sz = 101;    
    let mut base_distribution : Box<ReadBaseDistribution>  = Box::new(ReadBaseDistribution::new(maxreadlen , prec, ndarray::Dim([sz as usize, 4])) );
    
    for distrib in &distrib_collector {
        let ref_distrib = distrib.as_ref().unwrap();
        base_distribution.readsizehisto.merge(&ref_distrib.readsizehisto);
        // CAVEAT should check that upper_histo are coherent (the same) they are all initialized to maxreadlen
        base_distribution.upper_histo = maxreadlen;
        base_distribution.histo_out += ref_distrib.histo_out;
        base_distribution.non_acgt += ref_distrib.non_acgt;
        //   documentation suggest something like that should be possible
        base_distribution.acgt_distribution += &(ref_distrib.acgt_distribution);
    } // end of for on distrib_collector
    base_distribution.acgt_distribution *= 1./(seq_array.len() as f64);
    //
    let elapsed_t = start_t.elapsed().as_secs();
    println!(" elapsed time (s) in get_base_count {} ", elapsed_t);
    // dump
    log::info!("nb read outside max size for histogram {}", base_distribution.histo_out);
    base_distribution.ascii_dump_acgt_distribution(&String::from("bases.histo")).unwrap();
    //
    return Some(base_distribution);
} // end of get_base_count_par

