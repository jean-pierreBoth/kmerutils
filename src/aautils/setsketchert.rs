//! provide minimal tool to sketch AA sequences by probminhash3a or superminhash and hyperloglog




use std::marker::PhantomData;

use std::io::{BufReader, BufWriter };


use std::fs::OpenOptions;
use std::fmt::Debug;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};
use serde_json::to_writer;

use fnv::{FnvHashMap, FnvBuildHasher};
use std::hash::BuildHasherDefault;

use num;
use num::{Integer,ToPrimitive, FromPrimitive, Bounded};

use rand_distr::uniform::SampleUniform;

use crate::nohasher::*;

use crate::base::kmertraits::*;
use crate::aautils::kmeraa::*;

use rayon::prelude::*;

use probminhash::{probminhasher::*, superminhasher::SuperMinHash, setsketcher::SetSketcher, setsketcher::SetSketchParams};

use crate::sketcharg::{SeqSketcherParams, SketchAlgo};

/// This trait gathers interface to all sketcher : SuperMinhash, Probminhash3a, Probminhash3, ...  
/// 
/// It is useful when we need to send various sketchers in external functions as a impl Trait.
pub trait SeqSketcherAAT<Kmer> 
    where   Kmer : CompressedKmerT + KmerBuilder<Kmer>,
            KmerGenerator<Kmer> :  KmerGenerationPattern<Kmer> {
    /// Signature type of the sketch algo, f64 or f32 for SuperMinHash, Kmer::Val for ProbMinhashs
    type Sig : Serialize + Clone + Send + Sync;
    //
    fn get_kmer_size(&self) -> usize;
    /// returns the length of the sketch vector we want.
    fn get_sketch_size(&self) -> usize;
    //
    fn get_algo(&self) -> SketchAlgo;
    /// This function receive a vector of concatenated sequences and returns for each sequence a sketch.
    /// So the function returns a vector of Sketches.
    /// F is a hashing function (possibly just extracting Kmer::Val) to apply to kmer before sending to sketcher.
    fn sketch_compressedkmeraa<F>(&self, vseq : &Vec<&SequenceAA>, fhash : F) -> Vec<Vec<Self::Sig> > 
                    where F : Fn(&Kmer) -> Kmer::Val + Send + Sync; 
    /// This function implements the sketching a File of Sequences, 
    /// (The sequence are not concatenated, so we have many sequences) and make one sketch Vector for the sequence collection.
    /// It returns the same signature as sketch_compressedkmer for interface homogeneity (msg system for //) but
    /// but the returned vec has size 1!
    fn sketch_compressedkmeraa_seqs<F>(&self, vseq : &Vec<&SequenceAA>, fhash : F) -> Vec<Vec<Self::Sig> > 
                    where F : Fn(&Kmer) -> Kmer::Val + Send + Sync;  
}


//============================================================================================


// A structure providing ProbMinHash3a sketching for SequenceAA by implementing the generic trait SeqSketcherAAT<Kmer>
#[derive(Serialize,Deserialize,Copy,Clone)]
pub struct ProbHash3aSketch<Kmer> {
    //
    _kmer_marker: PhantomData<Kmer>,
    //
    params : SeqSketcherParams,
}


impl <Kmer> ProbHash3aSketch<Kmer> {

    pub fn new(params : &SeqSketcherParams) -> Self {
        ProbHash3aSketch{_kmer_marker : PhantomData,  params : params.clone()}
    }

} // end of impl ProbHash3aSketch



impl <Kmer> SeqSketcherAAT<Kmer> for ProbHash3aSketch<Kmer> 
        where   Kmer : CompressedKmerT + KmerBuilder<Kmer> + Send + Sync,
                Kmer::Val : num::PrimInt + Send + Sync + Debug + Clone + Serialize,
                KmerGenerator<Kmer> :  KmerGenerationPattern<Kmer> {

    type Sig = Kmer::Val;


    fn get_kmer_size(&self) -> usize {
        self.params.get_kmer_size()
    }

    fn get_sketch_size(&self) -> usize {
        self.params.get_sketch_size()
    }

    fn get_algo(&self) -> SketchAlgo {
        SketchAlgo::PROB3A
    }

    fn sketch_compressedkmeraa<F> (&self, vseq : &Vec<&SequenceAA>, fhash : F) -> Vec<Vec<Self::Sig> > 
            where  F : Fn(&Kmer) -> Kmer::Val + Send + Sync   {
        //
        log::debug!("entering sketch_compressedkmeraa for probminhash");
        //
        let comput_closure = | seqb : &SequenceAA, i:usize | -> (usize,Vec<Kmer::Val>) {
            // if we get very large sequence (many Gb length) we must be cautious on size of hashmap; i.e about number of different kmers!!! 
            let nb_kmer = get_nbkmer_guess(seqb);
            let mut wb : FnvHashMap::<Kmer::Val,u64> = FnvHashMap::with_capacity_and_hasher(nb_kmer, FnvBuildHasher::default());
            let mut kmergen = KmerSeqIterator::<Kmer>::new(self.get_kmer_size(), &seqb);
            kmergen.set_range(0, seqb.size()).unwrap();
            loop {
                match kmergen.next() {
                    Some(kmer) => {
                        let hashval = fhash(&kmer);
                        *wb.entry(hashval).or_insert(0) += 1;
                    },
                    None => break,
                }
            }  // end loop 
            let mut pminhashb = ProbMinHash3a::<Kmer::Val,NoHashHasher>::new(self.get_sketch_size(), 
                <Kmer::Val>::default());
            pminhashb.hash_weigthed_hashmap(&wb);
            let sigb = pminhashb.get_signature();
            // get back from usize to Kmer32bit ?. If fhash is inversible possible, else NO.
            return (i,sigb.clone());
        };
        //
        let sig_with_rank : Vec::<(usize,Vec<Kmer::Val>)> = (0..vseq.len()).into_par_iter().map(|i| comput_closure(vseq[i],i)).collect();
        // re-order from jac_with_rank to jaccard_vec as the order of return can be random!!
        let mut jaccard_vec = Vec::<Vec<Kmer::Val>>::with_capacity(vseq.len());
        for _ in 0..vseq.len() {
            jaccard_vec.push(Vec::new());
        }
        // CAVEAT , boxing would avoid the clone?
        for i in 0..sig_with_rank.len() {
            let slot = sig_with_rank[i].0;
            jaccard_vec[slot] = sig_with_rank[i].1.clone();
        }
        jaccard_vec
    }


    // recall we revceive a vecor of sequences originating from one file, we return a vector of size 1
    // containing a vecor of size sketch size and type Self::Sig
    fn sketch_compressedkmeraa_seqs<F>(&self, vseq : &Vec<&SequenceAA>, fhash : F) -> Vec<Vec<Self::Sig> >
            where   Kmer : CompressedKmerT + KmerBuilder<Kmer> + Send + Sync,
                    F : Fn(&Kmer) -> Kmer::Val + Send + Sync,
                    Kmer::Val : num::PrimInt + Send + Sync + Debug,
                    KmerGenerator<Kmer> :  KmerGenerationPattern<Kmer> {
        //
        log::debug!("entering  sketch_compressedkmeraa_seqs for Probminhash");
        // we must estimate nb kmer to avoid reallocation in FnvHashMap
        let nb_kmer = get_nbkmer_guess_seqs(vseq);
        //
        let mut wb : FnvHashMap::<Kmer::Val,u64> = FnvHashMap::with_capacity_and_hasher(nb_kmer, FnvBuildHasher::default());
        //
        let mut nb_kmer_generated : u64 = 0;
        // we loop on sequences and generate kmer. TODO // on sequences
        for seq in vseq {
            let mut kmergen = KmerSeqIterator::<Kmer>::new(self.get_kmer_size(), seq);
            kmergen.set_range(0, seq.size()).unwrap();
            loop {
                match kmergen.next() {
                    Some(kmer) => {
                        nb_kmer_generated += 1;
                        let hashval = fhash(&kmer);
                        *wb.entry(hashval).or_insert(0) += 1;
                    },
                    None => break,
                }
                if log::log_enabled!(log::Level::Debug) && nb_kmer_generated % 500_000_000 == 0 {
                    log::debug!("nb kmer generated : {:#}", nb_kmer_generated);
                }
            }  // end loop 
        }
        let mut pminhashb : ProbMinHash3a<Kmer::Val, NoHashHasher> = ProbMinHash3a::<Kmer::Val,NoHashHasher>::new(self.get_sketch_size(),
                    <Kmer::Val>::default());
        //
        pminhashb.hash_weigthed_hashmap(&wb);
        let sigb = pminhashb.get_signature();
        //
        let mut v = Vec::<Vec<Self::Sig>>::with_capacity(1);
        v.push(sigb.clone());
        //
        return v;
    }

}  // end of impl SeqSketcherAAT for ProHash3aSketch



//==================================================================================================================



/// A structure providing SuperMinHash sketching for SequenceAA by implementing the generic trait SeqSketcherAAT\<Kmer\>.  
///  The type argument S encodes for f32 or f64 as the SuperMinHash can sketch to f32 or f64
#[derive(Serialize,Deserialize,Copy,Clone)]
pub struct SuperHashSketch<Kmer, S : num::Float> {
    //
    _kmer_marker: PhantomData<Kmer>,
    //
    _sig_marker : PhantomData<S>,
    //
    params : SeqSketcherParams,
}


impl <Kmer, S : num::Float> SuperHashSketch<Kmer, S> {


    pub fn new(params : &SeqSketcherParams) -> Self {
        SuperHashSketch{_kmer_marker : PhantomData, _sig_marker : PhantomData,  params : params.clone()}
    }

} // end of impl ProbHash3aSketch

impl <Kmer, S> SeqSketcherAAT<Kmer> for SuperHashSketch<Kmer, S> 
        where   Kmer : CompressedKmerT + KmerBuilder<Kmer> + Send + Sync,
                Kmer::Val : num::PrimInt + Send + Sync + Debug,
                KmerGenerator<Kmer> :  KmerGenerationPattern<Kmer>,
                S : num::Float + SampleUniform + Send + Sync + Debug + Serialize {

    type Sig = S;

    fn get_kmer_size(&self) -> usize {
        self.params.get_kmer_size()
    }

    fn get_sketch_size(&self) -> usize {
        self.params.get_sketch_size()
    }

    fn get_algo(&self) -> SketchAlgo {
        SketchAlgo::SUPER
    }

    /// a generic implementation of superminhash  against our standard compressed Kmer types.  
    /// Kmer::Val is the base type u32, u64 on which compressed kmer representations relies.
    /// F is a hash function returning morally a u32, usize or u64.  
    /// The argument type of the hashing function F specify the type of Kmer to generate along the sequence.  
    fn sketch_compressedkmeraa<F>(&self, vseq : &Vec<&SequenceAA>, fhash : F) -> Vec<Vec<Self::Sig> >
        where F : Fn(&Kmer) -> Kmer::Val + Send + Sync {
        //
        log::debug!("entering sketch_compressedkmeraa for superminhash");
        //
        let comput_closure = | seqb : &SequenceAA, i:usize | -> (usize,Vec<Self::Sig>) {
            //
            log::trace!(" in sketch_compressedkmeraa (SuperMinHash), closure");
            let mut nb_kmer_generated : u64 = 0;
            //
            let bh = BuildHasherDefault::<NoHashHasher>::default();
            let mut sminhash : SuperMinHash<Self::Sig, Kmer::Val, NoHashHasher>= SuperMinHash::new(self.get_sketch_size(), bh);

            let mut kmergen = KmerSeqIterator::<Kmer>::new(self.get_kmer_size(), &seqb);
            kmergen.set_range(0, seqb.size()).unwrap();
            loop {
                match kmergen.next() {
                    Some(kmer) => {
                        nb_kmer_generated += 1;
                        let hashval = fhash(&kmer);
                        if sminhash.sketch(&hashval).is_err() {
                            log::error!("could not hash kmer : {:?}", kmer.get_uncompressed_kmer());
                            std::panic!("could not hash kmer : {:?}", kmer.get_uncompressed_kmer());
                        }
                    },
                    None => break,
                }
                if log::log_enabled!(log::Level::Debug) && nb_kmer_generated % 500_000_000 == 0 {
                    log::debug!("nb kmer generated : {:#}", nb_kmer_generated);
                }
            }  // end loop 
            let sigb = sminhash.get_hsketch();
            // get back from usize to Kmer32bit ?. If fhash is inversible possible, else NO.
            return (i,sigb.clone());
        };
        //
        let sig_with_rank : Vec::<(usize,Vec<Self::Sig>)> = (0..vseq.len()).into_par_iter().map(|i| comput_closure(vseq[i],i)).collect();
        // re-order from jac_with_rank to jaccard_vec as the order of return can be random!!
        let mut jaccard_vec = Vec::<Vec<Self::Sig>>::with_capacity(vseq.len());
        for _ in 0..vseq.len() {
            jaccard_vec.push(Vec::new());
        }
        // CAVEAT , boxing would avoid the clone?
        for i in 0..sig_with_rank.len() {
            let slot = sig_with_rank[i].0;
            jaccard_vec[slot] = sig_with_rank[i].1.clone();
        }
        jaccard_vec
    } // end of sketch_compressedkmeraa



    fn sketch_compressedkmeraa_seqs<F>(&self, vseq : &Vec<&SequenceAA>, fhash : F) -> Vec<Vec<Self::Sig> >
            where   Kmer : CompressedKmerT + KmerBuilder<Kmer> + Send + Sync,
                    F : Fn(&Kmer) -> Kmer::Val + Send + Sync,
                    Kmer::Val : num::PrimInt + Send + Sync + Debug,
                    KmerGenerator<Kmer> :  KmerGenerationPattern<Kmer> {
        //
        log::debug!("entering  sketch_compressedkmeraa_seqs for SuperMinHashSketch");
        //
        let bh = BuildHasherDefault::<NoHashHasher>::default();
        let mut setsketch : SuperMinHash<Self::Sig, Kmer::Val, NoHashHasher> = SuperMinHash::new(self.get_sketch_size(), bh);
        //
        let mut nb_kmer_generated : u64 = 0;
        // we loop on sequences and generate kmer. TODO // on sequences
        for seq in vseq {
            let mut kmergen = KmerSeqIterator::<Kmer>::new(self.get_kmer_size(), &seq);
            kmergen.set_range(0, seq.size()).unwrap();
            loop {
                match kmergen.next() {
                    Some(kmer) => {
                        nb_kmer_generated += 1;
                        let hashval = fhash(&kmer);
                        if setsketch.sketch(&hashval).is_err() {
                            log::error!("could not hash kmer : {:?}", kmer.get_uncompressed_kmer());
                            std::panic!("could not hash kmer : {:?}", kmer.get_uncompressed_kmer());
                        }
                    },
                    None => break,
                }
                if log::log_enabled!(log::Level::Debug) && nb_kmer_generated % 500_000_000 == 0 {
                    log::debug!("nb kmer generated : {:#}", nb_kmer_generated);
                }
            }  // end loop 
        }
        //
        let mut v: Vec<Vec<S>> = Vec::<Vec<Self::Sig>>::with_capacity(1);
        let sig = setsketch.get_hsketch();
        v.push(sig.clone());
        //
        return v;
    }

} // end of SuperHashSketch


//=========================================================================================================================




/// Defines the maximum number of threads to use in // iteratos in [HyperLogLogSketch]
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct HllSeqsThreading {
    /// max number of threads in iterator
    nb_iter_thread : usize,
    /// number of bases above which we use threading in // iterators
    thread_threshold : usize,
}


impl  HllSeqsThreading {
    pub fn new(nb_iter_thread : usize, thread_threshold : usize) -> Self {
        HllSeqsThreading{nb_iter_thread, thread_threshold}
    }

    // returns the number of intern iterator threads
    pub fn get_nb_iter_threads(&self) -> usize {
        self.nb_iter_thread
    }

    // return the number of base (in a list of sequences) above wich the is threading)
    pub fn get_thread_threshold(&self) -> usize {
        self.thread_threshold
    }
} // end impl HllSeqsThreading


impl Default for HllSeqsThreading {
    fn default() -> Self {
        HllSeqsThreading{ nb_iter_thread : 4, thread_threshold : 10_000_000}
    }
}

//===========================================

/// A structure providing HyperLogLog sketching for SequenceAA by implementing the generic trait SeqSketcherAAT\<Kmer\>.  
///  The type argument S encodes for u16 or u32 but u16 should be sufficient for hyperloglog. 
///  Currently HyperLogLogSketch uses rayon parallel iterator in sketch_compressedkmer_seqs to dispatch
///  the sketching of sequences into threads and the use a merge strategy.  
/// 
///  **The number of internal thtreads can be bounded by the structure HllSeqsThreading which defines the maximum number of threads blocks in // iterator**.   
///  **The number of threads is defined by (nb_bases/thread_threshold).ilog(3).min(HllSeqsThreading::nb_iter_thread).max(1).**  
///  This is useful if the caller use also multithreading.  

#[derive(Serialize,Deserialize,Copy,Clone)]
pub struct HyperLogLogSketch<Kmer, S : num::Integer> {
    //
    params : SeqSketcherParams,
    // this sketcher needs its particular parameters
    hll_params : SetSketchParams,
    //
    hll_threads : HllSeqsThreading,
    //
    _kmer_marker: PhantomData<Kmer>,
    //
    _sig_marker : PhantomData<S>,
}  // end of HyperLogLogSketch



impl <Kmer, S : Integer> HyperLogLogSketch<Kmer, S> {
    pub fn new(seq_params : &SeqSketcherParams, hll_params : SetSketchParams, hll_threads : HllSeqsThreading) -> Self {
        HyperLogLogSketch{params : seq_params.clone(), hll_params, hll_threads, _kmer_marker : PhantomData, _sig_marker: PhantomData}
    }

        // building block for sketch_compressedkmer_seqs. sketch a list of sequence and return a sketch to merge!
    pub fn sketch_compressedkmer_seqs_block<F>(&self, vseq : &[&SequenceAA], fhash : F) -> SetSketcher<S, Kmer::Val, NoHashHasher>
        where   Kmer : CompressedKmerT + KmerBuilder<Kmer> + Send + Sync,
                F : Fn(&Kmer) -> Kmer::Val + Send + Sync,
                Kmer::Val : num::PrimInt + Send + Sync + Debug,
                KmerGenerator<Kmer> :  KmerGenerationPattern<Kmer>,
                S : Integer + Bounded + Copy + Clone + FromPrimitive + ToPrimitive + Send + Sync + Debug + Serialize {
    //
    log::debug!("entering  sketch_compressedkmer_seqs_block for AA HyperLogLogSketch");
    //
    let bh = BuildHasherDefault::<NoHashHasher>::default();
    let mut setsketch : SetSketcher<S, Kmer::Val, NoHashHasher>= SetSketcher::new(self.hll_params, bh);
    //
    let mut nb_kmer_generated : u64 = 0;
    // we loop on sequences and generate kmer. TODO // on sequences
    for seq in vseq {
        let mut kmergen = KmerSeqIterator::<Kmer>::new(self.get_kmer_size(), &seq);
        kmergen.set_range(0, seq.size()).unwrap();
        loop {
            match kmergen.next() {
                Some(kmer) => {
                    nb_kmer_generated += 1;
                    let hashval = fhash(&kmer);
                    if setsketch.sketch(&hashval).is_err() {
                        log::error!("could not hash kmer : {:?}", kmer.get_uncompressed_kmer());
                        std::panic!("could not hash kmer : {:?}", kmer.get_uncompressed_kmer());
                    }
                },
                None => break,
            }
            if log::log_enabled!(log::Level::Debug) && nb_kmer_generated % 500_000_000 == 0 {
                log::debug!("nb kmer generated : {:#}", nb_kmer_generated);
            }
        }  // end loop 
    }
    //
    return setsketch;
}


} // en of impl HyperLogLogSketch


impl <Kmer,S> SeqSketcherAAT<Kmer> for HyperLogLogSketch<Kmer, S> 
        where   Kmer : CompressedKmerT + KmerBuilder<Kmer> + Send + Sync,
                Kmer::Val : num::PrimInt + Send + Sync + Debug,
                KmerGenerator<Kmer> :  KmerGenerationPattern<Kmer>,
                S : Integer + Bounded + Copy + Clone + FromPrimitive + ToPrimitive + Send + Sync + Debug + Serialize  {

    type Sig = S;

    fn get_kmer_size(&self) -> usize {
        self.params.get_kmer_size()
    }

    /// returns sketch size
    fn get_sketch_size(&self) -> usize {
        self.params.get_sketch_size()
    }

    fn get_algo(&self) -> SketchAlgo {
        SketchAlgo::HLL
    }


    /// a generic implementation of superminhash  against our standard compressed Kmer types.  
    /// Kmer::Val is the base type u32, u64 on which compressed kmer representations relies.
    /// F is a hash function returning morally a u32, usize or u64.  
    /// The argument type of the hashing function F specify the type of Kmer to generate along the sequence.  
    fn sketch_compressedkmeraa<F>(&self, vseq : &Vec<&SequenceAA>, fhash : F) -> Vec<Vec<Self::Sig> >
        where F : Fn(&Kmer) -> Kmer::Val + Send + Sync {
        //
        log::debug!("entering sketch_compressedkmeraa for setsketch");

        let comput_closure = | seqb : &SequenceAA, i:usize | -> (usize,Vec<Self::Sig>) {
            //
            log::debug!(" in sketch_compressedkmeraa, closure");
            let mut nb_kmer_generated : u64 = 0;
            //
            let bh = BuildHasherDefault::<NoHashHasher>::default();
            let mut setsketch : SetSketcher<Self::Sig, Kmer::Val, NoHashHasher>= SetSketcher::new(self.hll_params, bh);

            let mut kmergen = KmerSeqIterator::<Kmer>::new(self.get_kmer_size(), &seqb);
            kmergen.set_range(0, seqb.size()).unwrap();
            loop {
                match kmergen.next() {
                    Some(kmer) => {
                        nb_kmer_generated += 1;
                        let hashval = fhash(&kmer);
                        if setsketch.sketch(&hashval).is_err() {
                            log::error!("could not hash kmer : {:?}", kmer.get_uncompressed_kmer());
                            std::panic!("could not hash kmer : {:?}", kmer.get_uncompressed_kmer());
                        }
                    },
                    None => break,
                }
                if log::log_enabled!(log::Level::Debug) && nb_kmer_generated % 500_000_000 == 0 {
                    log::debug!("nb kmer generated : {:#}", nb_kmer_generated);
                }
            }  // end loop 
            let sigb = setsketch.get_signature();
            // get back from usize to Kmer32bit ?. If fhash is inversible possible, else NO.
            return (i,sigb.clone());
        };
        //
        let sig_with_rank : Vec::<(usize,Vec<Self::Sig>)> = (0..vseq.len()).into_par_iter().map(|i| comput_closure(vseq[i],i)).collect();
        // re-order from jac_with_rank to jaccard_vec as the order of return can be random!!
        let mut jaccard_vec = Vec::<Vec<Self::Sig>>::with_capacity(vseq.len());
        for _ in 0..vseq.len() {
            jaccard_vec.push(Vec::new());
        }
        // CAVEAT , boxing would avoid the clone?
        for i in 0..sig_with_rank.len() {
            let slot = sig_with_rank[i].0;
            jaccard_vec[slot] = sig_with_rank[i].1.clone();
        }
        jaccard_vec
    } // end of sketch_compressedkmeraa


    fn sketch_compressedkmeraa_seqs<F>(&self, vseq : &Vec<&SequenceAA>, fhash : F) -> Vec<Vec<Self::Sig> >
        where   Kmer : CompressedKmerT + KmerBuilder<Kmer> + Send + Sync,
                F : Fn(&Kmer) -> Kmer::Val + Send + Sync,
                Kmer::Val : num::PrimInt + Send + Sync + Debug,
                KmerGenerator<Kmer> :  KmerGenerationPattern<Kmer> {
        //
        log::debug!("entering  sketch_compressedkmeraa_seqs for AA setskecth");
        // Now we will try to dispatch work. We use hll for large sequence (>= 10^7) otherwise we propably should use probminhash
        let thread_threshold = self.hll_threads.get_thread_threshold();
        const BASE_LOG : usize = 3;
        let total_size = vseq.iter().fold(0, |acc, s| acc + s.size() );
        if total_size <= BASE_LOG * thread_threshold {
            let sketch =  self.sketch_compressedkmer_seqs_block(vseq, fhash);
            let mut v_sketch = Vec::<Vec<Self::Sig>>::new();
            v_sketch.push(sketch.get_signature().clone());
            return v_sketch;
        };
        // we must split work in equal parts.  At most 4 , we do not have so many threads. The threading must be treated at a higher level 
        // to correctly dispatch tasks.
        let nb_sequences = vseq.len();
        // we are sure that nb_blocks will be greater than 1.
        let nb_thread_max = self.hll_threads.get_nb_iter_threads();
        let nb_blocks : usize = nb_thread_max.min((total_size / thread_threshold).ilog(BASE_LOG) as usize).max(1);
        let block_size = nb_sequences / nb_blocks;
        log::debug!("total_size , block_size : {}", block_size);
        // TODO: compute frontiers when sequences do not have equal length
        let mut frontiers = Vec::<usize>::with_capacity(nb_blocks+1);
        for i in 0..nb_blocks {
            if i == 0 {
                frontiers.push(0);
            }
            else {
                frontiers.push(nb_sequences.min(i * block_size));
            }
        }
        frontiers.push(nb_sequences);
        //
        let v_sketch : Vec<SetSketcher<S, Kmer::Val, NoHashHasher> > = (0..nb_blocks).into_par_iter().map(|i| self.sketch_compressedkmer_seqs_block(&vseq[frontiers[i]..frontiers[i+1]], &fhash)).collect();
        // we allocate a sketcher that will contain the union. Signature is initialized to 0.
        let bh = BuildHasherDefault::<NoHashHasher>::default();
        let mut setsketch : SetSketcher<S, Kmer::Val, NoHashHasher>= SetSketcher::new(self.hll_params, bh);
        // now we can merge signatures
        for sketch in v_sketch {
            let res = setsketch.merge(&sketch);
            if res.is_err() {
                log::error!("an error occurred in merging signatures");
                std::panic!("an error occurred in merging signatures");
            }
        }
        //
        let sig = setsketch.get_signature();
        let mut v = Vec::<Vec<Self::Sig>>::with_capacity(1);
        v.push(sig.clone());
        //
        return v;
    } // end of sketch_compressedkmeraa_seqs

} // end of impl SeqSketcherT


//============================================================================================

// TODO this should be factorized with DNA case.

// This structure (deprecated, prefer ProbHash3aSketch and SuperHashSketch) describes the kmer size used in computing sketches and the number of sketch we want.
/// It gathers methods for sketch_superminhash, sketch_probminhash3a 
#[derive(Serialize,Deserialize,Copy,Clone)]
pub struct SeqSketcher {
    kmer_size : usize,
    sketch_size : usize
}  // end of SeqSketcher


impl SeqSketcher {
    /// 
    pub fn new(kmer_size: usize, sketch_size : usize) -> Self {
        SeqSketcher{kmer_size, sketch_size}
    }

    /// returns kmer size
    pub fn get_kmer_size(&self) -> usize {
        self.kmer_size
    }

    /// return sketch size
    pub fn get_sketch_size(&self) -> usize {
        self.sketch_size
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
    pub fn reload_json(dirpath : &Path) -> Result<SeqSketcher, String> {
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
        let sketch_params:SeqSketcher = serde_json::from_reader(reader).unwrap();
        //
        log::info!("SeqSketcher reload, kmer_size : {}, sketch_size : {}", 
            sketch_params.get_kmer_size(), sketch_params.get_sketch_size());     
        //
        Ok(sketch_params)
    } // end of reload_json





    /// A generic version of sketching with probminhash3a on compressed kmer for amino acids
    pub fn sketch_probminhash3a<'b, Kmer : CompressedKmerT + KmerBuilder<Kmer>, F>(&self, vseq : &'b Vec<&SequenceAA>, fhash : F) -> Vec<Vec<Kmer::Val> >
        where F : Fn(&Kmer) -> Kmer::Val + Send + Sync,
              Kmer::Val : num::PrimInt + Send + Sync + Debug,
              KmerGenerator<Kmer> :  KmerGenerationPattern<Kmer> {

        let comput_closure = | seqb : &SequenceAA, i:usize | -> (usize,Vec<Kmer::Val>) {
            // if we get very large sequence (many Gb length) we must be cautious on size of hashmap; i.e about number of different kmers!!! 
            let nb_kmer = get_nbkmer_guess(seqb);
            let mut wb : FnvHashMap::<Kmer::Val,u64> = FnvHashMap::with_capacity_and_hasher(nb_kmer, FnvBuildHasher::default());
            let mut kmergen = KmerSeqIterator::<Kmer>::new(self.kmer_size, &seqb);
            kmergen.set_range(0, seqb.size()).unwrap();
            loop {
                match kmergen.next() {
                    Some(kmer) => {
                        let hashval = fhash(&kmer);
                        *wb.entry(hashval).or_insert(0) += 1;
                    },
                    None => break,
                }
            }  // end loop 
            let mut pminhashb = ProbMinHash3a::<Kmer::Val,NoHashHasher>::new(self.sketch_size, 
                <Kmer::Val>::default());
            pminhashb.hash_weigthed_hashmap(&wb);
            let sigb = pminhashb.get_signature();
            // get back from usize to Kmer32bit ?. If fhash is inversible possible, else NO.
            return (i,sigb.clone());
        };
        //
        let sig_with_rank : Vec::<(usize,Vec<Kmer::Val>)> = (0..vseq.len()).into_par_iter().map(|i| comput_closure(vseq[i],i)).collect();
        // re-order from jac_with_rank to jaccard_vec as the order of return can be random!!
        let mut jaccard_vec = Vec::<Vec<Kmer::Val>>::with_capacity(vseq.len());
        for _ in 0..vseq.len() {
            jaccard_vec.push(Vec::new());
        }
        // CAVEAT , boxing would avoid the clone?
        for i in 0..sig_with_rank.len() {
            let slot = sig_with_rank[i].0;
            jaccard_vec[slot] = sig_with_rank[i].1.clone();
        }
        jaccard_vec
    }  // end of sketch_probminhash3a


   //  Superminhash

    /// a generic implementation of superminhash  against our standard compressed Kmer types.  
    /// Kmer::Val is the base type u32, u64 on which compressed kmer representations relies.
    /// F is a hash function returning morally a u32, usize or u64.  
    /// The argument type of the hashing function F specify the type of Kmer to generate along the sequence.  
    pub fn sketch_superminhash<'b, Kmer : CompressedKmerT + KmerBuilder<Kmer>, F>(&self, vseq : &'b Vec<&SequenceAA>, fhash : F) -> Vec<Vec<f64> >
        where F : Fn(&Kmer) -> Kmer::Val + Send + Sync,
              Kmer::Val : num::PrimInt + Send + Sync + Debug,
              KmerGenerator<Kmer> :  KmerGenerationPattern<Kmer> {
        //
        log::debug!("entering sketch_superminhash_compressedkmer");
        //
        let comput_closure = | seqb : &SequenceAA, i:usize | -> (usize,Vec<f64>) {
            //
            log::debug!(" in sketch_superminhash_compressedkmer, closure");
            //
            let bh = BuildHasherDefault::<fnv::FnvHasher>::default();
            // generic arg is here type sent to sketching
            let mut sminhash : SuperMinHash<f64, Kmer::Val, fnv::FnvHasher>= SuperMinHash::new(self.sketch_size, bh);

            let mut kmergen = KmerSeqIterator::<Kmer>::new(self.kmer_size, &seqb);
            kmergen.set_range(0, seqb.size()).unwrap();
            loop {
                match kmergen.next() {
                    Some(kmer) => {
                        let hashval = fhash(&kmer);
                        if sminhash.sketch(&hashval).is_err() {
                            log::error!("could not hash kmer : {:?}", kmer.get_uncompressed_kmer());
                            std::panic!("could not hash kmer : {:?}", kmer.get_uncompressed_kmer());
                        }
                    },
                    None => break,
                }
            }  // end loop 
            let sigb = sminhash.get_hsketch();
            // get back from usize to Kmer32bit ?. If fhash is inversible possible, else NO.
            return (i,sigb.clone());
        };
        //
        let sig_with_rank : Vec::<(usize,Vec<f64>)> = (0..vseq.len()).into_par_iter().map(|i| comput_closure(vseq[i],i)).collect();
        // re-order from jac_with_rank to jaccard_vec as the order of return can be random!!
        let mut jaccard_vec = Vec::<Vec<f64>>::with_capacity(vseq.len());
        for _ in 0..vseq.len() {
            jaccard_vec.push(Vec::new());
        }
        // CAVEAT , boxing would avoid the clone?
        for i in 0..sig_with_rank.len() {
            let slot = sig_with_rank[i].0;
            jaccard_vec[slot] = sig_with_rank[i].1.clone();
        }
        jaccard_vec
    } // end of sketch_superminhash

} // end of SeqSketcher (AA case)



//=========================================================


#[cfg(test)]
mod tests {

use super::*;

use crate::sketcharg::{SeqSketcherParams, SketchAlgo, DataType};

use std::str::FromStr;

    fn log_init_test() {
        let mut builder = env_logger::Builder::from_default_env();
        //    builder.filter_level(LevelFilter::Trace);
        let _ = builder.is_test(true).try_init();
    }

    #[test]
    fn test_seqaa_probminhash_64bit() {
        log_init_test();
        //
        log::debug!("test_seqaa_probminhash");
        //
        let str1 = "MTEQIELIKLYSTRILALAAQMPHVGSLDNPDASAMKRSPLCGSKVTVDVIMQNGKITFDGFEVLAPASEYKNRHASILLSLDATAEACASIAAQNSA";
        // The second string is the first half of the first repeated
        let str2 = "MTEQIELIKLYSTRILALAAQMPHVGSLDNPDASAMKRSPLCGSKVMTEQIELIKLYSTRILALAAQMPHVGSLDNPDASAMKRSPLCGSKV";

        let seq1 = SequenceAA::from_str(str1).unwrap();
        let seq2 = SequenceAA::from_str(str2).unwrap();
        let vseq = vec![&seq1, &seq2];
        let kmer_size = 5;
        let sketch_size = 400;
        let sketcher = SeqSketcher::new(kmer_size, sketch_size);
        let nb_alphabet_bits = Alphabet::new().get_nb_bits();
        // we need a hash function from u128 to f64
        let kmer_hash_fn = | kmer : &KmerAA64bit | -> <KmerAA64bit as CompressedKmerT>::Val {
            let mask : <KmerAA64bit as CompressedKmerT>::Val = num::NumCast::from::<u64>((0b1 << nb_alphabet_bits*kmer.get_nb_base()) - 1).unwrap();
            let hashval = kmer.get_compressed_value() & mask;
            hashval
        };
        let mask : u64 = num::NumCast::from::<u64>((0b1 << nb_alphabet_bits*kmer_size as u8) - 1).unwrap();
        log::debug!("mask = {:b}", mask);
        //
        log::info!("calling sketch_probminhash3a_compressedKmerAA64bit");
        let signatures = sketcher.sketch_probminhash3a(&vseq, kmer_hash_fn); 
        // get distance between the 2 strings  
        let sig1 = &signatures[0];
        let sig2 = &signatures[1];
        //
        let inter : u64 = sig1.iter().zip(sig2.iter()).map(|(a,b)| if a==b {1} else {0}).sum();
        let dist = inter as f64/sig1.len() as f64;
        log::info!("inter : {:?} length {:?} jaccard distance {:?}", inter, sig1.len(), dist );
        assert!( (dist-0.5).abs() < 1./10.);
    } // end of test_seqaa_probminhash_64bit


    #[test]
    fn test_seqaa_probminhash_trait_64bit() {
        log_init_test();
        //
        log::debug!("test_seqaa_probminhash");
        //
        let str1 = "MTEQIELIKLYSTRILALAAQMPHVGSLDNPDASAMKRSPLCGSKVTVDVIMQNGKITFDGFEVLAPASEYKNRHASILLSLDATAEACASIAAQNSA";
        // The second string is the first half of the first repeated
        let str2 = "MTEQIELIKLYSTRILALAAQMPHVGSLDNPDASAMKRSPLCGSKVMTEQIELIKLYSTRILALAAQMPHVGSLDNPDASAMKRSPLCGSKV";

        let seq1 = SequenceAA::from_str(str1).unwrap();
        let seq2 = SequenceAA::from_str(str2).unwrap();
        let vseq = vec![&seq1, &seq2];
        let kmer_size = 5;
        let sketch_size = 800;
        let sketch_args = SeqSketcherParams::new(kmer_size, sketch_size, SketchAlgo::PROB3A, DataType::AA);
        let sketcher = ProbHash3aSketch::<KmerAA64bit>::new(&sketch_args);
        let nb_alphabet_bits = Alphabet::new().get_nb_bits();
        // we need a hash function from u128 to f64
        let kmer_hash_fn = | kmer : &KmerAA64bit | -> <KmerAA64bit as CompressedKmerT>::Val {
            let mask : <KmerAA64bit as CompressedKmerT>::Val = num::NumCast::from::<u64>((0b1 << nb_alphabet_bits*kmer.get_nb_base()) - 1).unwrap();
            let hashval = kmer.get_compressed_value() & mask;
            hashval
        };
        let mask : u64 = num::NumCast::from::<u64>((0b1 << nb_alphabet_bits*kmer_size as u8) - 1).unwrap();
        log::debug!("mask = {:b}", mask);
        //
        log::info!("calling sketch_compressedkmeraa for ProbHash3aSketch::<KmerAA64bit>");
        let signatures = sketcher.sketch_compressedkmeraa(&vseq, kmer_hash_fn); 
        // get distance between the 2 strings  
        let sig1 = &signatures[0];
        let sig2 = &signatures[1];
        //
        let inter : u64 = sig1.iter().zip(sig2.iter()).map(|(a,b)| if a==b {1} else {0}).sum();
        let dist = inter as f64/sig1.len() as f64;
        log::info!("inter : {:?} length {:?} jaccard distance {:?}", inter, sig1.len(), dist );
        assert!( (dist-0.5).abs() < 1./10.);
    } // end of test_seqaa_probminhash_64bit


    #[test]
    fn test_seqaa_superminhash_trait_64bit() {
        log_init_test();
        //
        log::debug!("test_seqaa_superminhash_trait_64bit");
        //
        let str1 = "MTEQIELIKLYSTRILALAAQMPHVGSLDNPDASAMKRSPLCGSKVTVDVIMQNGKITFDGFEVLAPASEYKNRHASILLSLDATAEACASIAAQNSA";
        // The second string is the first half of the first repeated
        let str2 = "MTEQIELIKLYSTRILALAAQMPHVGSLDNPDASAMKRSPLCGSKVMTEQIELIKLYSTRILALAAQMPHVGSLDNPDASAMKRSPLCGSKV";

        let seq1 = SequenceAA::from_str(str1).unwrap();
        let seq2 = SequenceAA::from_str(str2).unwrap();
        let vseq = vec![&seq1, &seq2];
        let kmer_size = 5;
        let sketch_size = 800;
        let sketch_args = SeqSketcherParams::new(kmer_size, sketch_size, SketchAlgo::PROB3A, DataType::AA);
        let nb_alphabet_bits = Alphabet::new().get_nb_bits();
        let mask : u64 = num::NumCast::from::<u64>((0b1 << nb_alphabet_bits*kmer_size as u8) - 1).unwrap();
        log::debug!("mask = {:b}", mask);
        //
        // we need a hash function from u128 to f64
        let kmer_hash_fn = | kmer : &KmerAA64bit | -> <KmerAA64bit as CompressedKmerT>::Val {
            let mask : <KmerAA64bit as CompressedKmerT>::Val = num::NumCast::from::<u64>((0b1 << nb_alphabet_bits*kmer.get_nb_base()) - 1).unwrap();
            let hashval = kmer.get_compressed_value() & mask;
            hashval
        };
        // first we sketch with SuperHashSketch<f64>
        log::info!("calling sketch_compressedkmeraa for SuperHashSketch::<KmerAA64bit, f64>");
        let sketcher_f64 = SuperHashSketch::<KmerAA64bit, f64>::new(&sketch_args);
        let signatures = sketcher_f64.sketch_compressedkmeraa(&vseq, kmer_hash_fn); 
        // get distance between the 2 strings  
        let sig1 = &signatures[0];
        let sig2 = &signatures[1];
        //
        let inter : u64 = sig1.iter().zip(sig2.iter()).map(|(a,b)| if a==b {1} else {0}).sum();
        let dist = inter as f64/sig1.len() as f64;
        log::info!("SuperHashSketch::<KmerAA64bit, f64> inter : {:?} length {:?} jaccard distance {:?}", inter, sig1.len(), dist );
        assert!( (dist-0.5).abs() < 1./10.);
        //
        // now we sketch with SuperHashSketch<f32>
        let sketcher_f32 = SuperHashSketch::<KmerAA64bit, f32>::new(&sketch_args);
        let signatures = sketcher_f32.sketch_compressedkmeraa(&vseq, kmer_hash_fn); 
        // get distance between the 2 strings  
        let sig1 = &signatures[0];
        let sig2 = &signatures[1];
        //
        let inter : u64 = sig1.iter().zip(sig2.iter()).map(|(a,b)| if a==b {1} else {0}).sum();
        let dist = inter as f64/sig1.len() as f64;
        log::info!("SuperHashSketch::<KmerAA64bit, f32> inter : {:?} length {:?} jaccard distance {:?}", inter, sig1.len(), dist );
        assert!( (dist-0.5).abs() < 1./10.);
    } // end of test_seqaa_superminhash_trait_64bit



    #[test]
    fn test_seqaa_probminhash_32bit() {
        log_init_test();
        //
        log::debug!("test_seqaa_probminhash");
        //
        let str1 = "MTEQIELIKLYSTRILALAAQMPHVGSLDNPDASAMKRSPLCGSKVTVDVIMQNGKITFDGFEVLAPASEYKNRHASILLSLDATAEACASIAAQNSA";
        // The second string is the first half of the first repeated
        let str2 = "MTEQIELIKLYSTRILALAAQMPHVGSLDNPDASAMKRSPLCGSKVMTEQIELIKLYSTRILALAAQMPHVGSLDNPDASAMKRSPLCGSKV";

        let seq1 = SequenceAA::from_str(str1).unwrap();
        let seq2 = SequenceAA::from_str(str2).unwrap();
        let vseq = vec![&seq1, &seq2];
        let kmer_size = 5;
        let sketch_size = 400;
        let sketcher = SeqSketcher::new(kmer_size, sketch_size);
        let nb_alphabet_bits = Alphabet::new().get_nb_bits();
        // we need a hash function from u128 to f64
        let kmer_hash_fn = | kmer : &KmerAA32bit | -> <KmerAA32bit as CompressedKmerT>::Val {
            let mask : <KmerAA32bit as CompressedKmerT>::Val = num::NumCast::from::<u32>((0b1 << nb_alphabet_bits*kmer.get_nb_base()) - 1).unwrap();
            let hashval = kmer.get_compressed_value() & mask;
            hashval
        };
        let mask : u64 = num::NumCast::from::<u32>((0b1 << nb_alphabet_bits*kmer_size as u8) - 1).unwrap();
        log::debug!("mask = {:b}", mask);
        //
        log::info!("calling sketch_probminhash3a_compressedKmerAA32bit");
        let signatures = sketcher.sketch_probminhash3a(&vseq, kmer_hash_fn); 
        // get distance between the 2 strings  
        let sig1 = &signatures[0];
        let sig2 = &signatures[1];
        //
        let inter : u64 = sig1.iter().zip(sig2.iter()).map(|(a,b)| if a==b {1} else {0}).sum();
        let dist = inter as f64/sig1.len() as f64;
        log::info!("inter : {:?} length {:?} jaccard distance {:?}", inter, sig1.len(), dist );
        assert!( (dist-0.5).abs() < 1./10.);
    } // end of test_seqaa_probminhash_32bit



    #[test]
    fn test_seqaa_probminhash_gen() {
        log_init_test();
        //
        log::debug!("test_seqaa_probminhash");
        //
        let str1 = "MTEQIELIKLYSTRILALAAQMPHVGSLDNPDASAMKRSPLCGSKVTVDVIMQNGKITFDGFEVLAPASEYKNRHASILLSLDATAEACASIAAQNSA";
        // The second string is the first half of the first repeated
        let str2 = "MTEQIELIKLYSTRILALAAQMPHVGSLDNPDASAMKRSPLCGSKVMTEQIELIKLYSTRILALAAQMPHVGSLDNPDASAMKRSPLCGSKV";

        let seq1 = SequenceAA::from_str(str1).unwrap();
        let seq2 = SequenceAA::from_str(str2).unwrap();
        let vseq = vec![&seq1, &seq2];
        let kmer_size = 5;
        let sketch_size = 400;
        let sketcher = SeqSketcher::new(kmer_size, sketch_size);
        let nb_alphabet_bits = Alphabet::new().get_nb_bits();
        // we need a hash function from u128 to f64
        let kmer_hash_fn = | kmer : &KmerAA32bit | -> <KmerAA32bit as CompressedKmerT>::Val {
            let mask : <KmerAA32bit as CompressedKmerT>::Val = num::NumCast::from::<u32>((0b1 << nb_alphabet_bits*kmer.get_nb_base()) - 1).unwrap();
            let hashval = kmer.get_compressed_value() & mask;
            hashval
        };
        let mask : u64 = num::NumCast::from::<u32>((0b1 << nb_alphabet_bits*kmer_size as u8) - 1).unwrap();
        log::debug!("mask = {:b}", mask);
        //
        log::info!("calling sketch_probminhash3a_compressed_kmeraa for KmerAA32bit");
        let signatures = sketcher.sketch_probminhash3a(&vseq, kmer_hash_fn); 
        // get distance between the 2 strings  
        let sig1 = &signatures[0];
        let sig2 = &signatures[1];
        //
        let inter : u64 = sig1.iter().zip(sig2.iter()).map(|(a,b)| if a==b {1} else {0}).sum();
        let dist = inter as f64/sig1.len() as f64;
        log::info!("inter : {:?} length {:?} jaccard distance {:?}", inter, sig1.len(), dist );
        assert!( (dist-0.5).abs() < 1./10.);
    } // end of test_seqaa_probminhash_32bit


}  // end of mod tests in aautils::seqsketchjaccard