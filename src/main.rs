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
    let mut rng = ChaCha8Rng::seed_from_u64(987); // Draws [0.24346048, true, false, true]
    let efficiencies = vec![0.8; 30];
    let tree = pcr::simulate_tree(efficiencies, 100000, &mut rng);
    println!("{:?}", tree.matrix)
}
