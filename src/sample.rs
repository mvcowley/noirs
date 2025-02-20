//! Draw samples

use rand::Rng;
use rand_distr::Zipf;

/// Parametrises sampling
pub struct Sampler {
    /// Total number of observations
    pub observed: u32,

    /// Distribution observations per species are sampled from
    pub distribution: Zipf<f32>
}

impl Sampler {
    fn new(observed: u32, max_obs: f32, exponent: f32) -> Sampler {
        Sampler {
            observed,
            distribution: Zipf::new(max_obs, exponent).unwrap()
        }
    }
}

// fn draw<R: Rng + ?Sized>(sampler: Sampler, rng: &mut R) -> Vec<u32>

#[cfg(test)]
mod tests {
    use super::*;
    use rand_chacha::rand_core::SeedableRng;
    use rand_chacha::ChaCha8Rng;

    #[test]
    fn test_sampler() {
        let mut rng = ChaCha8Rng::seed_from_u64(927); // Draws 11.0 from Zipf defined below
        let observed = 1000000;
        let max_obs = 1e3;
        let exponent = 1.5;
        let sampler = Sampler::new(observed, max_obs, exponent);
        assert_eq!(sampler.observed, observed);
        let draw = rng.sample(sampler.distribution);
        assert_eq!(draw, 11.0);
    }
}

