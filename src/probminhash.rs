//! Implementation of ProbMinHash as described in O. Ertl

use wyhash::WyHash;
use rand::distributions;

/// Structure for defining exponential sampling of parameter lambda restricted with support restricted
/// to unit interval [0,1).
/// Specially adapted for ProbminHash3 and 4.
// All comments follow notations in Ertl article
#[derive(Clone, Copy, Debug)]
pub struct Exp_restricted_01 {
    /// parameter of exponential
    lambda : f64,
    c1 : f64,
    // abciss of point for which A3 is under exponential
    c2 : f64,
    c3 : f64,
    unit_range : Uniform
} // end of struct Exp_restricted_01


impl Exp_restricted_01  {
    pub fn new(lambda : f64) {
        let c1 = (lambda.exp() - 1.) / lambda;
        let c2 = (2./(1. + (-lambda).exp())).ln()/ lambda;
        let c3 = (1. - (-lambda).exp()) / lambda;
        Exp_restricted_01(lambda, c1, c2, c3, WyRng::default())
    }
}


impl Distribution<f64> for  Exp_restricted_01  {
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
            if x <= c3 * (1. - y) { return x }
            if c1 * y <= (1. - x) { return x }
            if y * c1 * self.lambda <= (self.lambda * (1.- x)).exp() - 1 { return x }
        }

        
    } // end sample
} 