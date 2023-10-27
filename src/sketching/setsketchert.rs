//! module defining trait SeqSketcherT and implementers for probminhash3a optdensminhash and hyperloglog skecther

/* 
- TODO:   in sketch_compressedkmer_seqs compute frontiers when sequences do not have equal length
- TODO:   In SuperMinHas::sketch_compressedkmer_seqs we could // with a mutex on setsketcher.

*/

use std::marker::PhantomData;

use std::fmt::Debug;


use std::hash::{BuildHasherDefault, Hasher};

use serde::{Deserialize, Serialize};

use fnv::{FnvHashMap, FnvBuildHasher};

use num;
use num::{Integer, ToPrimitive, FromPrimitive, Unsigned, Bounded};

use rand_distr::uniform::SampleUniform;

use crate::nohasher::*;

use crate::base::{kmer::*, kmergenerator::*, kmergenerator::KmerSeqIteratorT};

use super::nbkmerguess::*;


use rayon::prelude::*;

use crate::sketcharg::{SeqSketcherParams, SketchAlgo};

use probminhash::{probminhasher::*, superminhasher::SuperMinHash, densminhash::*, setsketcher::SetSketcher, setsketcher::SetSketchParams};


#[cfg(feature="sminhash2")]
use probminhash::superminhasher2::SuperMinHash2;



/// This trait gathers interface to all sketcher : SuperMinhash, Probminhash3a, Probminhash3, ...  
/// 
/// It is useful when we need to send various sketchers in external functionsas as a impl Trait.
pub trait SeqSketcherT<Kmer> 
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
    /// This function receive a vector of (possibly concatenated) sequences and returns for each sequence a sketch.  
    /// The function returns a vector of Sketches (one for each sequence).
    /// F is a hashing function (possibly just extracting Kmer::Val) to apply to kmer before sending to sketcher.
    fn sketch_compressedkmer<F>(&self, vseq : &Vec<&Sequence>, fhash : F) -> Vec<Vec<Self::Sig> > 
                    where F : Fn(&Kmer) -> Kmer::Val + Send + Sync;
    /// This function implements the sketching a file of Sequences, 
    /// (The sequence are not concatenated, so we have many sequences) and make one sketch Vector for the sequence collection (for the file).  
    /// **It returns the same signature as sketch_compressedkmer for interface homogeneity (same msg system for //)
    /// but the returned intern vec has size 1!**
    fn sketch_compressedkmer_seqs<F>(&self, vseq : &Vec<&Sequence>, fhash : F) -> Vec<Vec<Self::Sig> > 
                    where F : Fn(&Kmer) -> Kmer::Val + Send + Sync;                
} // end of SeqSketcherT<Kmer>




/// A structure providing ProbMinHash3a sketching implementing the generic trait SeqSketcherT\<Kmer\>.  
/// 
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



impl <Kmer> SeqSketcherT<Kmer> for ProbHash3aSketch<Kmer> 
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

    fn sketch_compressedkmer<F> (&self, vseq : &Vec<&Sequence>, fhash : F) -> Vec<Vec<Self::Sig> > 
            where  F : Fn(&Kmer) -> Kmer::Val + Send + Sync   {
        //
        log::debug!("entering sketch_probminhash3a_compressedkmer");
        //
        let comput_closure = | seqb : &Sequence, i:usize | -> (usize,Vec<Kmer::Val>) {
            // if we get very large sequence (many Gb length) we must be cautious on size of hashmap; i.e about number of different kmers!!! 
            let nb_kmer = get_nbkmer_guess(seqb);
            let mut wb : FnvHashMap::<Kmer::Val,u64> = FnvHashMap::with_capacity_and_hasher(nb_kmer, FnvBuildHasher::default());
            let mut kmergen = KmerSeqIterator::<Kmer>::new(self.get_kmer_size() as u8, &seqb);
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
        log::debug!("exiting sketch_probminhash3a_compressedkmer");
        jaccard_vec
    }



    // This functin implement the sketching a File of Sequences, (The sequence are not concatenated, so we have many sequences) and make one sketch Vector 
    fn sketch_compressedkmer_seqs<F>(&self, vseq : &Vec<&Sequence>, fhash : F) -> Vec<Vec<Self::Sig> > 
        where   Kmer : CompressedKmerT + KmerBuilder<Kmer> + Send + Sync,
                F : Fn(&Kmer) -> Kmer::Val + Send + Sync,
                Kmer::Val : num::PrimInt + Send + Sync + Debug,
                KmerGenerator<Kmer> :  KmerGenerationPattern<Kmer> {
        //
        log::debug!("entering sketch_compressedkmer_seqs for ProHash3aSketch");
        //
        // we must estimate nb kmer to avoid reallocation in FnvHashMap
        let nb_kmer = get_nbkmer_guess_seqs(vseq);
        //
        let mut wb : FnvHashMap::<Kmer::Val,u64> = FnvHashMap::with_capacity_and_hasher(nb_kmer, FnvBuildHasher::default());
        //
        let mut nb_kmer_generated : u64 = 0;
        // we loop on sequences and generate kmer. TODO // on sequences
        for seq in vseq {
            let mut kmergen = KmerSeqIterator::<Kmer>::new(self.get_kmer_size() as u8, &seq);
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
    } // end of sketch_compressedkmer_seqs

}  // end of impl SeqSketcherT for ProHash3aSketch


//=========================================================================================================

///
///  A structure providing SuperMinHash sketching implementing the generic trait SeqSketcherT\<Kmer\>.  
///  The type argument S encodes for f32 or f64 as the SuperMinHash can sketch to f32 or f64
#[derive(Serialize,Deserialize,Copy,Clone)]
pub struct SuperHashSketch<Kmer, S: num::Float> {
    //
    _kmer_marker: PhantomData<Kmer>,
    //
    _sig_marker: PhantomData<S>,
    //
    params : SeqSketcherParams,
}


impl <Kmer, S : num::Float> SuperHashSketch<Kmer,S> {


    pub fn new(params : &SeqSketcherParams) -> Self {
        SuperHashSketch{_kmer_marker : PhantomData, _sig_marker: PhantomData,  params : params.clone()}
    }

} // end of impl SuperHashSketch




impl <Kmer,S> SeqSketcherT<Kmer> for SuperHashSketch<Kmer, S> 
        where   Kmer : CompressedKmerT + KmerBuilder<Kmer> + Send + Sync,
                Kmer::Val : num::PrimInt + Send + Sync + Debug,
                KmerGenerator<Kmer> :  KmerGenerationPattern<Kmer>,
                S : num::Float + SampleUniform + Send + Sync + Debug + Serialize  {

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
    fn sketch_compressedkmer<F>(&self, vseq : &Vec<&Sequence>, fhash : F) -> Vec<Vec<Self::Sig> >
        where F : Fn(&Kmer) -> Kmer::Val + Send + Sync {
        //
        log::debug!("entering sketch_superminhash_compressedkmer");
        //
        let comput_closure = | seqb : &Sequence, i:usize | -> (usize,Vec<Self::Sig>) {
            //
            log::debug!(" in sketch_compressedkmer, closure");
            let mut nb_kmer_generated : u64 = 0;
            //
            let bh = BuildHasherDefault::<NoHashHasher>::default();
            let mut sminhash : SuperMinHash<Self::Sig, Kmer::Val, NoHashHasher>= SuperMinHash::new(self.get_sketch_size(), bh);

            let mut kmergen = KmerSeqIterator::<Kmer>::new(self.get_kmer_size() as u8, &seqb);
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
    } // end of sketch_compressedkmer


    fn sketch_compressedkmer_seqs<F>(&self, vseq : &Vec<&Sequence>, fhash : F) -> Vec<Vec<Self::Sig> >
            where   Kmer : CompressedKmerT + KmerBuilder<Kmer> + Send + Sync,
                    F : Fn(&Kmer) -> Kmer::Val + Send + Sync,
                    Kmer::Val : num::PrimInt + Send + Sync + Debug,
                    KmerGenerator<Kmer> :  KmerGenerationPattern<Kmer> {
        //
        log::debug!("entering  sketch_compressedkmer_seqs for SuperMinHashSketch");
        //
        let bh = BuildHasherDefault::<NoHashHasher>::default();
        let mut setsketch : SuperMinHash<Self::Sig, Kmer::Val, NoHashHasher> = SuperMinHash::new(self.get_sketch_size(), bh);
        //
        let mut nb_kmer_generated : u64 = 0;
        // we loop on sequences and generate kmer. TODO // on sequences
        for seq in vseq {
            let mut kmergen = KmerSeqIterator::<Kmer>::new(self.get_kmer_size() as u8, &seq);
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
        let mut v = Vec::<Vec<Self::Sig>>::with_capacity(1);
        let sig = setsketch.get_hsketch();
        v.push(sig.clone());
        //
        return v;
    }



} // end of SuperHashSketch

//====================================================================================


///  A structure providing Optimal Densification MinHash (OptDensMinHash in probminhash crate) sketching implementing the generic trait SeqSketcherT\<Kmer\>.  
///  The type argument S encodes for f32 or f64 as for SuperMinHash
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct OptDensHashSketch<Kmer, S: num::Float> {
    //
    _kmer_marker: PhantomData<Kmer>,
    //
    _sig_marker: PhantomData<S>,
    //
    params : SeqSketcherParams,
}


impl <Kmer, S : num::Float> OptDensHashSketch<Kmer,S> {
    pub fn new(params : &SeqSketcherParams) -> Self {
        OptDensHashSketch{_kmer_marker : PhantomData, _sig_marker: PhantomData,  params : params.clone()}
    }
}  // end of OptDensMinHashSketch



impl <Kmer,S> SeqSketcherT<Kmer> for OptDensHashSketch<Kmer, S> 
        where   Kmer : CompressedKmerT + KmerBuilder<Kmer> + Send + Sync,
                Kmer::Val : num::PrimInt + Send + Sync + Debug,
                KmerGenerator<Kmer> :  KmerGenerationPattern<Kmer>,
                S : num::Float + SampleUniform + Send + Sync + Debug + Serialize  {

    type Sig = S;

    fn get_kmer_size(&self) -> usize {
        self.params.get_kmer_size()
    }

    fn get_sketch_size(&self) -> usize {
        self.params.get_sketch_size()
    }

    fn get_algo(&self) -> SketchAlgo {
        SketchAlgo::OPTDENS
    }

    fn sketch_compressedkmer<F>(&self, vseq : &Vec<&Sequence>, fhash : F) -> Vec<Vec<Self::Sig> >
    where F : Fn(&Kmer) -> Kmer::Val + Send + Sync {
        //
       log::debug!("entering OptDensHashSketch::sketch_compressedkmer");
              //
            let comput_closure = | seqb : &Sequence, i:usize | -> (usize,Vec<Self::Sig>) {
            //
            log::debug!(" in sketch_compressedkmer, closure");
            let mut nb_kmer_generated : u64 = 0;
            //
            let bh = BuildHasherDefault::<NoHashHasher>::default();
            let mut sminhash : OptDensMinHash<Self::Sig, Kmer::Val, NoHashHasher>= OptDensMinHash::new(self.get_sketch_size(), bh);

            let mut kmergen = KmerSeqIterator::<Kmer>::new(self.get_kmer_size() as u8, &seqb);
            kmergen.set_range(0, seqb.size()).unwrap();
            loop {
                match kmergen.next() {
                    Some(kmer) => {
                        nb_kmer_generated += 1;
                        let hashval = fhash(&kmer);
                        sminhash.sketch(&hashval);
                    },
                    None => break,
                }
                if log::log_enabled!(log::Level::Debug) && nb_kmer_generated % 500_000_000 == 0 {
                    log::debug!("nb kmer generated : {:#}", nb_kmer_generated);
                }
            }  // end loop 
            // do not forget to close sketching (it calls densification!)
            sminhash.end_sketch();
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
        //
        jaccard_vec
    } // end of sketch_compressedkmer


    fn sketch_compressedkmer_seqs<F>(&self, vseq : &Vec<&Sequence>, fhash : F) -> Vec<Vec<Self::Sig> >
            where   Kmer : CompressedKmerT + KmerBuilder<Kmer> + Send + Sync,
                    F : Fn(&Kmer) -> Kmer::Val + Send + Sync,
                    Kmer::Val : num::PrimInt + Send + Sync + Debug,
                    KmerGenerator<Kmer> :  KmerGenerationPattern<Kmer> {
        //
        log::debug!("entering  OptDensHashSketch::sketch_compressedkmer_seqs");
        //
        let bh = BuildHasherDefault::<NoHashHasher>::default();
        let mut setsketch : OptDensMinHash<Self::Sig, Kmer::Val, NoHashHasher> = OptDensMinHash::new(self.get_sketch_size(), bh);
        //
        let mut nb_kmer_generated : u64 = 0;
        // we loop on sequences and generate kmer. TODO // on sequences
        for seq in vseq {
            let mut kmergen = KmerSeqIterator::<Kmer>::new(self.get_kmer_size() as u8, &seq);
            kmergen.set_range(0, seq.size()).unwrap();
            loop {
                match kmergen.next() {
                    Some(kmer) => {
                        nb_kmer_generated += 1;
                        let hashval = fhash(&kmer);
                        setsketch.sketch(&hashval);
                    },
                    None => break,
                }
                if log::log_enabled!(log::Level::Debug) && nb_kmer_generated % 500_000_000 == 0 {
                    log::debug!("nb kmer generated : {:#}", nb_kmer_generated);
                }
            }  // end loop 
        }
        //
        let mut v = Vec::<Vec<Self::Sig>>::with_capacity(1);
        // do not forget to close sketching (it calls densification!)
        setsketch.end_sketch();
        let sig = setsketch.get_hsketch();
        v.push(sig.clone());
        //
        return v;
    } // end of sketch_compressedkmer_seqs

} // end of impl SeqSketcherT<Kmer> for OptDensHashSketch


//=====================================================================================

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

//
///  A structure providing SetSketcher (HyperLogLog) sketching implementing the generic trait SeqSketcherT\<Kmer\>.  
///  The type argument S encodes for u16 , u32 or u64 as the SetSketcher can sketch to u16,  u32 or u64
/// 
///  Currently HyperLogLogSketch uses rayon parallel iterator in sketch_compressedkmer_seqs to dispatch
///  the sketching of sequences into threads and the use a merge strategy.  
///  **The number of internal thtreads can be bounded by the structure HllSeqsThreading which defines the maximum number of threads blocks in // iterator**.   
///  **The number of threads is defined by (nb_bases/thread_threshold).ilog(3).min(HllSeqsThreading::nb_iter_thread).max(1).**  
///  This is useful if the caller use also multithreading.  

#[derive(Serialize,Deserialize,Copy,Clone)]
pub struct HyperLogLogSketch<Kmer, S: num::Integer> {
    //
    params : SeqSketcherParams,
    // this sketcher needs its particular parameters
    hll_params : SetSketchParams,
    //
    hll_threads : HllSeqsThreading,
    //
    _kmer_marker: PhantomData<Kmer>,
    //
    _sig_marker: PhantomData<S>,

} // end of HyperLogLogSketch


impl <Kmer, S : Integer> HyperLogLogSketch<Kmer, S> {
    pub fn new(seq_params : &SeqSketcherParams, hll_params : SetSketchParams, hll_threads : HllSeqsThreading) -> Self {
        HyperLogLogSketch{params : seq_params.clone(), hll_params, hll_threads, 
                _kmer_marker :  PhantomData, _sig_marker: PhantomData}
    }

    // building block for sketch_compressedkmer_seqs. sketch a list of sequence and return a sketch to merge!
    pub fn sketch_compressedkmer_seqs_block<F>(&self, vseq : &[&Sequence], fhash : F) -> SetSketcher<S, Kmer::Val, NoHashHasher>
            where   Kmer : CompressedKmerT + KmerBuilder<Kmer> + Send + Sync,
                    F : Fn(&Kmer) -> Kmer::Val + Send + Sync,
                    Kmer::Val : num::PrimInt + Send + Sync + Debug,
                    KmerGenerator<Kmer> :  KmerGenerationPattern<Kmer>,
                    S : Integer + Bounded + Copy + Clone + FromPrimitive + ToPrimitive + Send + Sync + Debug + Serialize {
        //
        log::trace!("entering  sketch_compressedkmer_seqs_block for HyperLogLogSketch");
        //
        let bh = BuildHasherDefault::<NoHashHasher>::default();
        let mut setsketch : SetSketcher<S, Kmer::Val, NoHashHasher>= SetSketcher::new(self.hll_params, bh);
        //
        let mut nb_kmer_generated : u64 = 0;
        // we loop on sequences and generate kmer. TODO // on sequences
        for seq in vseq {
            let mut kmergen = KmerSeqIterator::<Kmer>::new(self.get_kmer_size() as u8, &seq);
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


impl <Kmer,S> SeqSketcherT<Kmer> for HyperLogLogSketch<Kmer, S> 
        where   Kmer : CompressedKmerT + KmerBuilder<Kmer> + Send + Sync,
                Kmer::Val : num::PrimInt + Send + Sync + Debug,
                KmerGenerator<Kmer> :  KmerGenerationPattern<Kmer>,
                S : Integer + Bounded + Copy + Clone + FromPrimitive + ToPrimitive + Send + Sync + Debug + Serialize  {

    type Sig = S;

    fn get_kmer_size(&self) -> usize {
        self.params.get_kmer_size()
    }

    fn get_sketch_size(&self) -> usize {
        self.params.get_sketch_size()
    }

    fn get_algo(&self) -> SketchAlgo {
        SketchAlgo::HLL
    }

    // This funtions sketch a list of of Sequence and returns a Sketch vector for each one.
    // In fact each sequence is a file that was concatenated in a sequence.
    fn sketch_compressedkmer<F>(&self, vseq : &Vec<&Sequence>, fhash : F) -> Vec<Vec<Self::Sig> >
        where F : Fn(&Kmer) -> Kmer::Val + Send + Sync {
        //
        log::debug!("entering sketch_compressedkmer for HyperLogLogSketch");
        //
        let comput_closure = | seqb : &Sequence, i:usize | -> (usize,Vec<Self::Sig>) {
            //
            log::debug!(" in sketch_compressedkmer, closure");
            let mut nb_kmer_generated : u64 = 0;
            //
            let bh = BuildHasherDefault::<NoHashHasher>::default();
            let mut setsketch : SetSketcher<Self::Sig, Kmer::Val, NoHashHasher>= SetSketcher::new(self.hll_params, bh);

            let mut kmergen = KmerSeqIterator::<Kmer>::new(self.get_kmer_size() as u8, &seqb);
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
            let sigb = setsketch.get_signature().clone();
            // closure in function in // iter, we drop explicitly
            drop(setsketch);
            // get back from usize to Kmer32bit ?. If fhash is inversible possible, else NO.
            return (i,sigb);
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
    } // end of sketch_compressedkmer


    // This function implement the sketching a File of Sequences.
    // By analogy to sketch_compressedkmer it returns a Vec<Vec<Self::Sig>>, but the resulting vector contains only one Vec<Self::Sig>
    // Algo:
    // The sequence are not concatenated, so we have many sequences. We dispatch sequences to sketch_compressedkmer_seqs_block
    // by parallelizing and merge sketch Vector.  
    // 
    fn sketch_compressedkmer_seqs<F>(&self, vseq : &Vec<&Sequence>, fhash : F) -> Vec<Vec<Self::Sig>>
        where F : Fn(&Kmer) -> Kmer::Val + Send + Sync  {
        //
        if log::log_enabled!(log::Level::Debug) {
            log::debug!("entering sketch_compressedkmer_seqs for HyperLogLogSketch");
            log::debug!("memory  : {:?}", memory_stats::memory_stats().unwrap());
        }
        // Now we will try to dispatch work. We use hll for large sequence (>= 10^7) otherwise we propably should use probminhash
        let thread_threshold = self.hll_threads.get_thread_threshold();
        const BASE_LOG : usize = 3;
        let total_size = vseq.iter().fold(0, |acc, s| acc + s.size() );
        if total_size <= BASE_LOG * thread_threshold {
            log::debug!("  calling directly sketch_compressedkmer_seqs_block, total size : {}", total_size);
            let sketch =  self.sketch_compressedkmer_seqs_block(vseq, fhash);
            let mut v_sketch = Vec::<Vec<Self::Sig>>::with_capacity(1);
            v_sketch.push(sketch.get_signature().clone());
            drop(sketch);
            if log::log_enabled!(log::Level::Debug) {
                log::debug!("exiting sketch_compressedkmer_seqs for HyperLogLogSketch");
                log::debug!("memory  : {:?}", memory_stats::memory_stats().unwrap());
            }
            return v_sketch;
        }
        // we must split work in equal parts.  A few threads , we do not have so many threads. The threading must be treated at a higher level 
        // to correctly dispatch tasks.
        let nb_sequences = vseq.len();
        // now we are sure that nb_blocks is >= 1 !!! 
        let nb_thread_max = self.hll_threads.get_nb_iter_threads();
        let nb_blocks : usize = nb_thread_max.min((total_size / thread_threshold).ilog(BASE_LOG) as usize);
        // we want blocks to be the same size as far as possible
        let block_size = (nb_sequences as f64 / nb_blocks as f64).round() as usize;
        log::debug!("nb_base : {}, nb_seq {}, block_size (nb_seq/thread): {}, nb_blocks : {}", total_size, nb_sequences, block_size, nb_blocks);
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
        for sketch in &v_sketch {
            let res = setsketch.merge(&sketch);
            if res.is_err() {
                log::error!("an error occurred in merging signatures");
                std::panic!("an error occurred in merging signatures");
            }
        }
        //
        let sig = setsketch.get_signature().clone();
        let mut v = Vec::<Vec<Self::Sig>>::with_capacity(1);
        v.push(sig);
        // explicit drop to monitor memory
        drop(v_sketch);
        drop(setsketch);
        if log::log_enabled!(log::Level::Debug) {
            log::debug!("exiting sketch_compressedkmer_seqs for HyperLogLogSketch");
            log::debug!("memory  : {:?}", memory_stats::memory_stats().unwrap());
        }
        //
        return v;
    } // end of sketch_compressedkmer_seqs


} // end of impl for HyperLogLogSketch






//=====================================================================================
///
///  A structure providing SuperMinHash2 sketching implementing the generic trait SeqSketcherT\<Kmer\>.  
///  The type argument S encodes for u32 or u64 as the SuperMinHash2 can sketch to u32 or u64
#[cfg(feature="sminhash2")]
#[derive(Clone)]
pub struct SuperHash2Sketch<Kmer, S: Integer  + Unsigned, H : Hasher + Default> {
    //
    _kmer_marker: PhantomData<Kmer>,
    //
    _sig_marker: PhantomData<S>,
    //
    build_hasher : BuildHasherDefault<H>,
    //
    params : SeqSketcherParams,
}

#[cfg(feature="sminhash2")]
impl <Kmer, S :  Integer  + Unsigned,  H : Hasher + Default> SuperHash2Sketch<Kmer,S, H> {


    pub fn new(params : &SeqSketcherParams, build_hasher: BuildHasherDefault<H>) -> Self {
        SuperHash2Sketch{_kmer_marker : PhantomData, _sig_marker: PhantomData,  build_hasher, params : params.clone()}
    }

} // end of impl ProbHash3aSketch
#[cfg(feature="sminhash2")]
impl <Kmer,S, H> SeqSketcherT<Kmer> for SuperHash2Sketch<Kmer, S, H> 
        where   Kmer : CompressedKmerT + KmerBuilder<Kmer> + Send + Sync,
                Kmer::Val : num::PrimInt + Send + Sync + Debug,
                KmerGenerator<Kmer> :  KmerGenerationPattern<Kmer>,
                H : Hasher + Default,
                S :Integer  + Unsigned + ToPrimitive + FromPrimitive + Bounded + Copy + Clone + Send + Sync + Serialize + std::fmt::Debug {

    type Sig = S;

    fn get_kmer_size(&self) -> usize {
        self.params.get_kmer_size()
    }

    fn get_sketch_size(&self) -> usize {
        self.params.get_sketch_size()
    }

    fn get_algo(&self) -> SketchAlgo {
        SketchAlgo::SUPER2
    }

    /// a generic implementation of superminhash  against our standard compressed Kmer types.  
    /// Kmer::Val is the base type u32, u64 on which compressed kmer representations relies.
    /// F is a hash function returning morally a u32, usize or u64.  
    /// The argument type of the hashing function F specify the type of Kmer to generate along the sequence.  
    fn sketch_compressedkmer<F>(&self, vseq : &Vec<&Sequence>, fhash : F) -> Vec<Vec<Self::Sig> >
        where F : Fn(&Kmer) -> Kmer::Val + Send + Sync {
        //
        log::debug!("entering sketch_superminhash2_compressedkmer");
        //
        let comput_closure = | seqb : &Sequence, i:usize | -> (usize,Vec<Self::Sig>) {
            //
            log::debug!(" in sketch_compressedkmer, closure");
            let mut nb_kmer_generated : u64 = 0;
            //
            let mut sminhash : SuperMinHash2<Self::Sig, Kmer::Val, H>= SuperMinHash2::new(self.get_sketch_size(), self.build_hasher.clone());

            let mut kmergen = KmerSeqIterator::<Kmer>::new(self.get_kmer_size() as u8, &seqb);
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
    } // end of sketch_compressedkmer




    #[cfg(feature="sminhash2")]
    fn sketch_compressedkmer_seqs<F>(&self, vseq : &Vec<&Sequence>, fhash : F) -> Vec<Vec<Self::Sig> >
            where   Kmer : CompressedKmerT + KmerBuilder<Kmer> + Send + Sync,
                    F : Fn(&Kmer) -> Kmer::Val + Send + Sync,
                    Kmer::Val : num::PrimInt + Send + Sync + Debug,
                    KmerGenerator<Kmer> :  KmerGenerationPattern<Kmer> {
        //
        log::debug!("entering  sketch_compressedkmer_seqs for SuperHash2Sketch");
        //
        let bh = BuildHasherDefault::<NoHashHasher>::default();
        let mut setsketch : SuperMinHash2<Self::Sig, Kmer::Val, NoHashHasher> = SuperMinHash2::new(self.get_sketch_size(), bh);
        //
        let mut nb_kmer_generated : u64 = 0;
        // we loop on sequences and generate kmer. TODO // on sequences
        for seq in vseq {
            let mut kmergen = KmerSeqIterator::<Kmer>::new(self.get_kmer_size() as u8, &seq);
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
        let mut v = Vec::<Vec<Self::Sig>>::with_capacity(1);
        let sig = setsketch.get_hsketch();
        v.push(sig.clone());
        //
        return v;
    }


} // end of SuperHash2Sketch

