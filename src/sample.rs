//! Draw samples

use crate::pcr::{simulate_tree, Reaction, SparseTree};
use crate::sequence::{sequence, Sequencer};
use ndarray::prelude::*;
use ndarray::concatenate;
use ndarray_npy::write_npy;
use rand::{seq::IndexedRandom, Rng};
use rand_distr::Zipf;

/// Parametrises sampling
pub struct Sampler {
    /// Total number of observations
    pub observed: u32,

    /// Distribution observations per species are sampled from
    pub distribution: Zipf<f32>,
}

/// Function to create Sampler
impl Sampler {
    /// Create sampler holding read depth target and parametrised Zipfian distribution
    pub fn new(observed: u32, max_obs: f32, exponent: f32) -> Sampler {
        Sampler {
            observed,
            distribution: Zipf::new(max_obs, exponent).unwrap(),
        }
    }
}

/// Draw from sampler and convert to integer
fn draw<R: Rng + ?Sized>(sampler: &Sampler, rng: &mut R) -> u32 {
    rng.sample(sampler.distribution).round() as u32
}

/// Join the node id, pcr mutation, and sequencing mutation arrays
fn join(arr1: Array2<u32>, arr2: Array2<u32>, arr3: Array1<u32>) -> Array2<u32> {
    concatenate![Axis(1), arr1, arr2, arr3.insert_axis(Axis(1))]
}

/// Run library sampler and write out each UMI's observations directly to file
pub fn sample<R: Rng + ?Sized>(
    sampler: &Sampler,
    reaction: &Reaction,
    sequencer: &Sequencer,
    rng: &mut R,
    out: &str,
) {
    // Loop while observed < N
    let mut umi_number = 0;
    let mut observation_count = 0;
    while observation_count < sampler.observed {
        println!("Simulating UMI number {}...", &umi_number);
        println!("Observations: {}/{}...", &observation_count, &sampler.observed);
        let umi_obs = draw(sampler, rng);
        umi_number += 1;
        observation_count += umi_obs;

        // Simulate reaction and sequencing
        let tree = simulate_tree(reaction, umi_obs, rng);
        let sequencer_errors = sequence(
            &tree
                .observations
                .index_axis(ndarray::Axis(1), tree.observations.shape()[1] - 1),
            &reaction,
            sequencer,
            rng,
        );
        
        // Join arrays and write out
        let umi_array = join(tree.observations, tree.mutations, sequencer_errors);
        let filename = format!("{out}{umi_number}.npy");
        println!("Writing UMI number {} to {}...", &umi_number, &filename);
        let _ = write_npy(filename, &umi_array);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand_chacha::rand_core::SeedableRng;
    use rand_chacha::ChaCha8Rng;

    #[test]
    fn test_sampler() {
        let mut rng = ChaCha8Rng::seed_from_u64(927); // Draws 11.0 from Zipf defined below
        let observed = 1_000_000;
        let max_obs = 1e3;
        let exponent = 1.5;
        let sampler = Sampler::new(observed, max_obs, exponent);
        assert_eq!(sampler.observed, observed);
        let draw = rng.sample(sampler.distribution);
        assert_eq!(draw, 11.0);
    }

    #[test]
    fn test_draw() {
        let mut rng = ChaCha8Rng::seed_from_u64(927);
        let observed = 4;
        let max_obs = 1e3;
        let exponent = 1.5;
        let sampler = Sampler::new(observed, max_obs, exponent);
        let draw = draw(&sampler, &mut rng);
        assert_eq!(draw, 11);
    }

    #[test]
    fn test_join() {
        let arr1 = array![[1, 2], [3, 4]];
        let arr2 = array![[1, 2], [3, 4]];
        let arr3 = array![5, 6];
        let joined = join(arr1, arr2, arr3);
        assert_eq!(joined, array![[1, 2, 1, 2, 5], [3, 4, 3, 4, 6]])
    }
}
