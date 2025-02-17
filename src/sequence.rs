//! Error profiles for sequencing

use ndarray::prelude::*;

use crate::pcr::PcrParameters;
use rand::Rng;
use rand_distr::{Binomial, Distribution};

/// Parameters of sequencing technology chosen
pub struct Sequencer {
    // Error rate
    error: f64,
}

pub fn sequence<R: Rng + ?Sized>(
    library: &Array1<u32>,
    reaction: &PcrParameters,
    sequencer: Sequencer,
    rng: &mut R,
) -> Array1<u32> {
    let bin = Binomial::new(reaction.sites as u64, sequencer.error).unwrap();
    library.mapv(|_| bin.sample(rng).try_into().unwrap())
}


#[cfg(test)]
mod tests {
    use super::*;
    use rand_chacha::rand_core::SeedableRng;
    use rand_chacha::ChaCha8Rng;

    #[test]
    fn test_sequence() {
        let library = array![4, 5, 3];
        let reaction = PcrParameters {
            sites: 3,
            efficiencies: vec![0.95; 2],
            errors: vec![0.5; 2],
        };
        let sequencer = Sequencer {
            error: 0.5
        };
        let mut rng = ChaCha8Rng::seed_from_u64(927); // Draws [2, 1, 3, 2, 1, 0] from binomial with n=3 and p=0.5
        let sequence_errors = sequence(&library, &reaction, sequencer, &mut rng);
        assert_eq!(sequence_errors, array![2, 1, 3]);
    }
}
