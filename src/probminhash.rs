//! Implementation of ProbMinHash as described in O. Ertl

use rand::distributions::{Distribution,Uniform};
use rand::prelude::*;
use rand_xoshiro::Xoshiro256PlusPlus;

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
        let mut current_value = value;
        let mut current_k = k;
        while current_value < self.values[k] {
            self.values[k] = current_value;
            let pidx = self.m + (current_k/2) as usize;   // m + upper integer value of k/2 beccause of 0 based indexation
            if pidx > self.last_index {
                break;
            }
            let siblidx = current_k^1;      // get sibling index of k with numerationbeginning at 0
            if self.values[siblidx] >= self.values[pidx] {
                break;                      // means parent stores the value of sibling, no more propagation needed
            }
            // now we now parent stores the value of index k
            if value < self.values[siblidx] {
                // if value is less than its sibling , we must set the parent to sibling value and propagate
                current_value = self.values[siblidx];
            }
            current_k = pidx;
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
    }
} // end of impl MaxValueTracker

/// A Trait to define association of a weight to an object.
/// Typically we could implement Trait WeightedSet for an IndexMap<Object, f64> giving a weight to each object
/// or 
pub trait WeightedSet {
    type Object;
    fn get_weight(&self, obj:&Self::Object) -> f64;
}


pub struct ProbMinHash3 {
    m : usize,
    /// field to keep track of max hashed values
    maxvaluetracker : MaxValueTracker,
    /// a expoential law restricted to interval [0., 1)
    exp01 : ExpRestricted01,
    ///  final signature of distribution. allocated to size m
    signature : Vec<usize>,
} // end of struct ProbMinHash3


impl ProbMinHash3 {
    pub fn new(nbhash:usize) -> Self {
        assert!(nbhash >= 2);
        let lambda = ((nbhash as f64)/((nbhash - 1) as f64)).ln();
        let h_signature = Vec::<usize>::with_capacity(nbhash as usize);
        ProbMinHash3{m:nbhash, maxvaluetracker: MaxValueTracker::new(nbhash as usize), 
                    exp01:ExpRestricted01::new(lambda), signature:h_signature}
    } // end of new
    

    /// incrementally adds an item in hash signature.
    /// It is the building block of the computation, but this method 
    /// does not check for unicity of id added in hash computation.  
    /// It is user responsability to enforce that. See method hashWSet
    pub fn hash_item(&mut self, id:usize, weight:f64) {
        assert!(weight > 0.);
        let winv = 1./weight;
        let unif0m = Uniform::<usize>::new(0, self.m);
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(id as u64);
        let mut h = winv * self.exp01.sample(&mut rng);
        let mut i = 1;
        let mut qmax = self.maxvaluetracker.get_max_value();
        while h < qmax {
            let k = unif0m.sample(&mut rng);
            if h < self.maxvaluetracker.values[k] {
                self.maxvaluetracker.values[k] = h;
                self.signature[k] = id;
                // 
                self.maxvaluetracker.update(k, h);
                qmax = self.maxvaluetracker.get_max_value();
            }
            i = i + 1;
            h = winv * (i-1) as f64;
            if h >= qmax {
                break;
            }
            h = h + winv * self.exp01.sample(&mut rng);
        }
    } // end of hash_item

    /// return final signature
    pub fn get_signature(&self) -> &Vec<usize> {
        return &self.signature
    }

    /// hash data when given by an iterable WeightedSet
    pub fn hash_wset<T>(&mut self, data: &mut T)
        where T: WeightedSet<Object=usize> + Iterator<Item=usize> {
            while let Some(obj) = &data.next() {
                let weight = data.get_weight(&obj);
                self.hash_item(*obj, weight);
            }
    } // end of hash method


}  // end of impl ProbMinHash3


//=================================================================


#[cfg(test)]
mod tests {

use super::*;

    #[test]    
    // This test stores random values in a MaxValueTracker and check for max at higher end of array
    fn test_max_value_tracker() {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(45678 as u64);

        let nbhash = 100;
        let unif_01 = Uniform::<f64>::new(0., 1.);
        let unif_m = Uniform::<usize>::new(0, nbhash);

        let mut tracker = MaxValueTracker::new(nbhash);
        //
        let mut vmax = 0f64;
        let loop_size = 5000;
        //
        for _ in 0..loop_size {
            let k = unif_m.sample(&mut rng);
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
    }

}  // end of module tests