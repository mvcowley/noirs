//! Generate noisy adaptive immune receptor data
mod amplicon;
mod fastx;
mod noise;
mod parse;
mod pcr;
mod sample;
mod sequence;
mod zipf;

use rand_chacha::rand_core::SeedableRng;
use rand_chacha::ChaCha8Rng;

fn main() {
    // let mut rng = ChaCha8Rng::seed_from_u64(987); // run 1
    // let mut rng = ChaCha8Rng::seed_from_u64(392); // run 2
    let mut rng = ChaCha8Rng::seed_from_u64(654); // run 3

    // Setup sampler
    let total_observations = 1_000_000;
    let max_observations = 1e3;
    let exponent = 1.5;
    let sampler = sample::Sampler::new(total_observations, max_observations, exponent);

    // Setup PCR reaction
    let reaction = pcr::Reaction {
        sites: 12,
        efficiencies: vec![0.8; 30],
        errors: vec![0.001; 30],
    };

    // Setup sequencer
    let sequencer = sequence::Sequencer { error: 0.005 };

    // Sample library
    sample::sample(&sampler, &reaction, &sequencer, &mut rng, "../out/noirs_out_3_high_error/");
}
