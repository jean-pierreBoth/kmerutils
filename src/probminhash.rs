//! Implementation of ProbMinHash as described in O. Ertl

use log::{trace};

use std::fmt::{Debug};

use rand::distributions::{Distribution,Uniform};
use rand::prelude::*;
use rand_xoshiro::Xoshiro256PlusPlus;

use wyhash::*;

/// Structure for defining exponential sampling of parameter lambda restricted with support restricted
/// to unit interval [0,1).
/// Specially adapted for ProbminHash3 and 4.
// All comments follow notations in Ertl article
#[derive(Clone, Copy, Debug)]
pub struct ExpRestricted01 {
    /// parameter of exponential
    lambda : f64,
    c1 : f64,
    // abciss of point for which A3 is under exponential
    c2 : f64,
    c3 : f64,
    /// we build upon a uniform [0,1) sampling
    unit_range : Uniform<f64>,
} // end of struct ExpRestricted01


impl ExpRestricted01  {
    pub fn new(lambda : f64) -> Self {
        let c1 = (lambda.exp() - 1.) / lambda;
        let c2 = (2./(1. + (-lambda).exp())).ln()/ lambda;
        let c3 = (1. - (-lambda).exp()) / lambda;
        ExpRestricted01{lambda, c1, c2, c3, unit_range:Uniform::<f64>::new(0.,1.)}
    }

    /// return lambda parameter of exponential
    pub fn get_lambda(&self) -> f64 {
        self.lambda
    }
}


impl Distribution<f64> for  ExpRestricted01  {
    fn sample<R : Rng + ?Sized>(&self, rng: &mut R) -> f64 {
        let mut x = self.c1 * rng.sample(&self.unit_range);
        if x < 1.  { return x }
        loop {
            // check if we can sample in A3
            x = rng.sample(&self.unit_range);
            if x < self.c2 { return x}
            // 
            let mut y = 0.5 * rng.sample(&self.unit_range);
            if y > 1. - x {
                // transform a point in A5 to a point in A6
                x = 1. - x;
                y = 1. - y;
            }
            if x <= self.c3 * (1. - y) { return x }
            if self.c1 * y <= (1. - x) { return x }
            if y * self.c1 * self.lambda <= (self.lambda * (1.- x)).exp_m1() { return x }
        }        
    } // end sample
} 

// structure to keep track of max values in hash set
// adapted from class MaxValueTracker
struct MaxValueTracker {
    m : usize,
    // last_index = 2*m-2. max of array is at slot last_index
    last_index : usize,
    // dimensioned to m hash functions
    values : Vec<f64>
}


impl MaxValueTracker {
    pub fn new(m:usize) -> Self {
        let last_index  = ((m << 1) - 2) as usize;  // 0-indexation for the difference with he paper, lastIndex = 2*m-2
        let vlen = last_index+1;
        let values : Vec::<f64> = (0..vlen).map( |_| f64::INFINITY).collect();
        MaxValueTracker{m, last_index, values}
    }

 
    
    // update slot k with value value
    // 0 indexation imposes some changes with respect to the the algo 4 of the paper
    // parent of k is m + (k/2)
    // and accordingly
    // sibling ok k is k+1 if k even, k-1 else so it is given by bitxor(k,1)
    fn update(&mut self, k:usize, value:f64) {
        assert!(k < self.m);
        trace!("\n max value tracker update k, value , value at k {} {} {} ", k, value, self.values[k]);
        let mut current_value = value;
        let mut current_k = k;
        let mut more = false;
        if current_value < self.values[current_k] {
            more = true;
        }
        
        while more {
            trace!("mxvt update k value {} {}", current_k, current_value);
            self.values[current_k] = current_value;
            let pidx = self.m + (current_k/2) as usize;   // m + upper integer value of k/2 beccause of 0 based indexation
            if pidx > self.last_index {
                break;
            }
            let siblidx = current_k^1;      // get sibling index of k with numeration beginning at 0
            assert!(self.values[siblidx] <= self.values[pidx]);
            assert!(self.values[current_k] <= self.values[pidx]);
            //
            if self.values[siblidx] >= self.values[pidx] && self.values[current_k] >= self.values[pidx] {
                break;     // means parent current and sibling are equals no more propagation needed
            }
            // now either self.values[siblidx] <self.values[pidx] or current_value < self.values[pidx]
            trace!("propagating current_value {} sibling  {} ? ", current_value, self.values[siblidx]);
            //
            if current_value < self.values[siblidx] {
                trace!("     propagating sibling value {} to parent {}", self.values[siblidx], pidx);
                current_value = self.values[siblidx];
            }
            else {
                trace!("     propagating current_value {} to parent {}", current_value, pidx);   
            }
            current_k = pidx;
            if current_value >= self.values[current_k]  {
                more = false;
            }
        }
    } // end of update function 


    /// return the maximum value maintained in the data structure
    pub fn get_max_value(&self) -> f64 {
        return self.values[self.last_index]
    }

    #[allow(dead_code)]
    pub fn get_parent_slot(&self, slot : usize) -> usize {
        assert!(slot <= self.m);
        return self.m + (slot/2) as usize   // m + upper integer value of k/2 beccause of 0 based indexation
    }

    /// get value MaxValueTracker at slot
    #[allow(dead_code)]
    pub fn get_value(&self, slot: usize) -> f64 {
        self.values[slot]
    } // end of get_value

    #[allow(dead_code)]
    pub fn dump(&self) {
        println!("\n\nMaxValueTracker dump : ");
        for i in 0..self.values.len() {
            println!(" i  value   {}   {} ", i , self.values[i]);
        }
    } // end of dump
} // end of impl MaxValueTracker

/// A Trait to define association of a weight to an object.
/// Typically we could implement Trait WeightedSet for an IndexMap<Object, f64> giving a weight to each object
/// or 
pub trait WeightedSet {
    type Object;
    fn get_weight(&self, obj:&Self::Object) -> f64;
}


/// implementation of the algorithm ProbMinHash3 as described in Etrl Paper
/// D must be convertible injectively into a use for random generator initialization hence the requirement Into<usize>
/// If all data are directly referred to by an index D is the index (usize)
pub struct ProbMinHash3<D> 
            where D:Copy+Eq+Into<usize>+Debug   {
    m : usize,
    /// field to keep track of max hashed values
    maxvaluetracker : MaxValueTracker,
    /// a exponential law restricted to interval [0., 1)
    exp01 : ExpRestricted01,
    ///  final signature of distribution. allocated to size m
    signature : Vec<D>,
} // end of struct ProbMinHash3


impl<D> ProbMinHash3<D> 
            where D:Copy+Eq+Debug+Into<usize> {

    pub fn new(nbhash:usize, initobj : D) -> Self {
        assert!(nbhash >= 2);
        let lambda = ((nbhash as f64)/((nbhash - 1) as f64)).ln();
        let h_signature = (0..nbhash).map( |_| initobj).collect();
        ProbMinHash3{m:nbhash, maxvaluetracker: MaxValueTracker::new(nbhash as usize), 
                    exp01:ExpRestricted01::new(lambda), signature:h_signature}
    } // end of new
    

    /// incrementally adds an item in hash signature.
    /// It is the building block of the computation, but this method 
    /// does not check for unicity of id added in hash computation.  
    /// It is user responsability to enforce that. See method hashWSet
    pub fn hash_item(&mut self, id:D, weight:f64) {
        assert!(weight > 0.);
        trace!("hash_item : id {:?}  weight {} ", id, weight);
        let winv = 1./weight;
        let unif0m = Uniform::<usize>::new(0, self.m);
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(id.into() as u64);
        let mut h = winv * self.exp01.sample(&mut rng);
        let mut i = 1;
        let mut qmax = self.maxvaluetracker.get_max_value();
        while h < qmax {
            let k = unif0m.sample(&mut rng);
            assert!(k < self.m);
            if h < self.maxvaluetracker.values[k] {
                self.signature[k] = id;
                // 
                self.maxvaluetracker.update(k, h);
                qmax = self.maxvaluetracker.get_max_value();
            }
            h = winv * i as f64;
            i = i + 1;
            if h >= qmax {
                break;
            }
            h = h + winv * self.exp01.sample(&mut rng);
            trace!("hash_item :  i h qmax =  {}   {}   {} ", i, h, qmax);
        }
    } // end of hash_item

    /// return final signature.
    pub fn get_signature(&self) -> &Vec<D> {
        return &self.signature
    }

    /// hash data when given by an iterable WeightedSet
    pub fn hash_wset<T>(&mut self, data: &mut T)
        where T: WeightedSet<Object=D> + Iterator<Item=D> {
            while let Some(obj) = &data.next() {
                let weight = data.get_weight(&obj);
                self.hash_item(*obj, weight);
            }
    } // end of hash method


}  // end of impl ProbMinHash3



/// fisher Yates random permutation generation, with lazy generation
/// of an array of size n
pub struct FYshuffle {
    m: usize,
    rng : Xoshiro256PlusPlus,
    /// uniform distribution on [0,1)
    unif_01 : Uniform<f64>,
    ///
    v : Vec<usize>,
    ///
    lastidx : usize,
}

impl FYshuffle {
    pub fn new(m: usize) -> FYshuffle {
        let v : Vec<usize> = (0..m).map(|x| x).collect();
        // seed Xoshiro256PlusPlus with WyRng...
        let mut rng_init = WyRng::default();
        FYshuffle{m:m, rng:Xoshiro256PlusPlus::seed_from_u64(rng_init.next_u64()), unif_01: Uniform::<f64>::new(0., 1.), v : v, lastidx:m}
    }

    // See https://www.geeksforgeeks.org/generate-a-random-permutation-of-1-to-n/
    // and The algorithm design manual S.Skiena P.458
    pub fn next(&mut self) -> usize {
        if self.lastidx >= self.m {
            self.lastidx = 0;
        }
        let xsi = self.unif_01.sample(&mut self.rng);
        // sample between self.lastidx (included) and self.m (excluded)
        let idx = self.lastidx + (xsi * (self.m - self.lastidx) as f64) as usize;
        let val = self.v[idx];
        self.v[idx] = self.v[self.lastidx];
        self.v[self.lastidx] = val;
        self.lastidx += 1;
        val
    }

    pub fn reset(&mut self) {
        self.lastidx = 0;
    }

    pub fn get_values(&self) -> &Vec<usize> {
        &self.v
    }

}  // end of impl FYshuffle


//=================================================================


#[cfg(test)]
mod tests {

use log::*;

#[allow(dead_code)]
fn log_init() {
    let _ = env_logger::builder().is_test(true).try_init();
}

use super::*;


    #[test]
    // This tests exponential random sampling in [0,1) 
    // by comparing theoretical mean and estimated mean and checking for deviation 
    // with nb_sampled = 1_000_000_000 we get 
    // mu_th 0.4585059174632017 mean 0.45850733816056904  sigma  0.000009072128699429336 
    // test = (mu_th - mean)/sigma = -0.15660022189165437
    // But as it needs some time we set nb_sampled to 10_000_000. 
    // test is often negative, so mu_th is approximated by above ?. to check
    fn test_exp01() {
        log_init();
        //
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(234567 as u64);
        let mut xsi;
        let lambda = 0.5f64;
        let mut mu_th = (-lambda - 1.) * (-lambda).exp() + 1.;
        mu_th =  mu_th / (lambda * (1. - (-lambda).exp()));
        //
        let nb_sampled = 10_000_000;
        let mut sampled = Vec::<f64>::with_capacity(nb_sampled);
        let exp01 = ExpRestricted01::new(lambda);
        //  
        for _ in 0..nb_sampled {
            xsi = exp01.sample(&mut rng);
            sampled.push(xsi);
        }
        let sum = sampled.iter().fold(0., |acc, x| acc +x);
        let mean = sum / nb_sampled as f64;
        //
        let mut s2 = sampled.iter().fold(0., |acc, x| acc +(x-mean)*(x-mean));
        s2 = s2 / (nb_sampled - 1)  as f64;
        // 
        println!("mu_th {} mean {}  sigma  {} ", mu_th, mean, (s2/nb_sampled as f64).sqrt());
        let test = (mu_th - mean) / (s2/nb_sampled as f64).sqrt();
        println!("test {}", test);
    }
    #[test]    
    // This test stores random values in a MaxValueTracker and check for max at higher end of array
    fn test_max_value_tracker() {
        log_init();
        //
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(45678 as u64);

        let nbhash = 10;
        let unif_01 = Uniform::<f64>::new(0., 1.);
        let unif_m = Uniform::<usize>::new(0, nbhash);

        let mut tracker = MaxValueTracker::new(nbhash);
        //
        let mut vmax = 0f64;
        let loop_size = 500;
        //
        for _ in 0..loop_size {
            let k = unif_m.sample(&mut rng);
            assert!(k < nbhash);
            let xsi = unif_01.sample(&mut rng);
            vmax = vmax.max(xsi);
            tracker.update(k,xsi);
            // check equality of max
            assert!( !( vmax > tracker.get_max_value() && vmax < tracker.get_max_value()) );
            // check for sibling and their parent coherence
        }
        // check for sibling and their parent coherence
        for i in 0..nbhash {
            let sibling = i^1;
            let sibling_value = tracker.get_value(sibling);
            let i_value = tracker.get_value(i);
            let pidx = tracker.get_parent_slot(i);
            let pidx_value = tracker.get_value(pidx);
            assert!(sibling_value <=  pidx_value && i_value <= pidx_value);
            assert!( !( sibling_value > pidx_value  &&   i_value >  pidx_value) );
        }
        assert!(!( vmax > tracker.get_max_value()  && vmax < tracker.get_max_value() ));
        tracker.dump();
    } // end of test_probminhash_count_range_intersection

    #[test] 
    // This test checks that with equal weights we fall back to Jaccard estimate
    fn test_probminhash3_count_intersection_equal_weights() {
        //
        log_init();
        //
        debug!("test_probminhash_count_intersection_equal_weights");
        println!("test_probminhash_count_intersection_equal_weights");
        // we construct 2 ranges [a..b] [c..d], with a<b, b < d, c<d sketch them and compute jaccard.
        // we should get something like max(b,c) - min(b,c)/ (b-a+d-c)
        //
        let set_size = 100;
        let nbhash = 10;
        //
        // choose weights for va and vb elements
        let mut wa = Vec::<f64>::with_capacity(set_size);
        let mut wb = Vec::<f64>::with_capacity(set_size);
        // initialize wa, weight 20 up to 130
        for i in 0..set_size {
            if i < 70 {
                wa.push(20.);
            }
            else {
                wa.push(0.);
            }
        }
        // initialize wb weight 10 above 70
        for i in 0..set_size {
            if i < 50 {
                wb.push(0.);
            }
            else {
                wb.push(10.);
            }
        }        
        // compute Jp as in 
        let mut jp = 0.;
        for i in 0..set_size {
            if wa[i] > 0. && wb[i] > 0. {
                let mut den = 0.;
                for j in 0..set_size {
                    den += (wa[j]/wa[i]).max(wb[j]/wb[i]);
                }
                jp += 1./den;
            }
        }
        trace!("Jp = {} ",jp);
        // probminhash 
        trace!("\n\n hashing wa");
        let mut waprobhash = ProbMinHash3::new(nbhash, 0);
        for i in 0..set_size {
            if wa[i] > 0. {
                waprobhash.hash_item(i, wa[i]);
            }
        }
        waprobhash.maxvaluetracker.dump();
        //
        trace!("\n\n hashing wb");
        let mut wbprobhash = ProbMinHash3::new(nbhash,0 );
        for i in 0..set_size {
            if wb[i] > 0. {
                wbprobhash.hash_item(i, wb[i]);
            }
        }        
        let siga = waprobhash.get_signature();
        let sigb = wbprobhash.get_signature();
        let mut inter = 0;
        for i in 0..siga.len() {
            if siga[i] == sigb[i] {
                inter += 1;
            }
        }
        //
        waprobhash.maxvaluetracker.dump();
        wbprobhash.maxvaluetracker.dump();
        //
        info!("jp = {} , inter / card = {} ", jp, inter as f64/siga.len() as f64);
    } // end of test_prob_count_intersection



#[test] 
// This test checks JaccardProbability with unequal weights inside sets
fn test_probminhash3_count_intersection_unequal_weights() {
    //
    log_init();
    //
    println!("test_probminhash_count_intersection_unequal_weights");
    debug!("test_probminhash_count_intersection_unequal_weights");
    // we construct 2 ranges [a..b] [c..d], with a<b, b < d, c<d sketch them and compute jaccard.
    // we should get something like max(b,c) - min(b,c)/ (b-a+d-c)
    //
    let set_size = 100;
    let nbhash = 50;
    //
    // choose weights for va and vb elements
    let mut wa = Vec::<f64>::with_capacity(set_size);
    let mut wb = Vec::<f64>::with_capacity(set_size);
    // initialize wa, weight 20 up to 130
    for i in 0..set_size {
        if i < 70 {
            wa.push(2. * i as f64);
        }
        else {
            wa.push(0.);
        }
    }
    // initialize wb weight 10 above 70
    for i in 0..set_size {
        if i < 50 {
            wb.push(0.);
        }
        else {
            wb.push( (i as f64).powi(4));
        }
    }        
    // compute Jp as in 
    let mut jp = 0.;
    for i in 0..set_size {
        if wa[i] > 0. && wb[i] > 0. {
            let mut den = 0.;
            for j in 0..set_size {
                den += (wa[j]/wa[i]).max(wb[j]/wb[i]);
            }
            jp += 1./den;
        }
    }
    trace!("Jp = {} ",jp);
    // probminhash 
    trace!("\n\n hashing wa");
    let mut waprobhash = ProbMinHash3::new(nbhash, 0);
    for i in 0..set_size {
        if wa[i] > 0. {
            waprobhash.hash_item(i, wa[i]);
        }
    }
    waprobhash.maxvaluetracker.dump();
    //
    trace!("\n\n hashing wb");
    let mut wbprobhash = ProbMinHash3::new(nbhash, 0);
    for i in 0..set_size {
        if wb[i] > 0. {
            wbprobhash.hash_item(i, wb[i]);
        }
    }        
    let siga = waprobhash.get_signature();
    let sigb = wbprobhash.get_signature();
    let mut inter = 0;
    for i in 0..siga.len() {
        if siga[i] == sigb[i] {
            inter += 1;
        }
    }
    //
    waprobhash.maxvaluetracker.dump();
    wbprobhash.maxvaluetracker.dump();
    //
    info!("jp = {} , inter / card = {} ", jp, inter as f64/siga.len() as f64);
} // end of test_prob_count_intersection


#[test]
// We check we have a unifom distribution of values at each rank of v
// variance is 5
fn test_fyshuffle() {

    log_init();

    let m = 4;
    let mut fypermut = FYshuffle::new(m);
    let nb_permut = 500000;
    let mut freq : Vec<usize> = (0..m).map(|_| 0).collect();

    for _ in 0..nb_permut {
        fypermut.next();
        let v = fypermut.get_values();
        for k in 0..v.len() {
            freq[k] += v[k];
        }
    }

    let th_freq = 1.5;
    let th_var = 5./4.;
    let sigma = (th_var/ (nb_permut as f64)).sqrt();
    for i in 0..freq.len() {
        let rel_error = ((freq[i] as f64)/ (nb_permut as f64) - th_freq)/ sigma;
        trace!(" slot i {} , rel error = {}", i, rel_error);
        assert!( rel_error.abs() < 3.)
    }
    trace!("  freq = {:?}", freq);
} // end of test_fyshuffle




}  // end of module tests