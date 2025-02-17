//! Generate noisy adaptive immune receptor data
mod amplicon;
mod fastx;
mod noise;
mod parse;
mod pcr;
mod sample;
mod sequence;
mod zipf;

use ndarray_npy::write_npy;
use rand_chacha::rand_core::SeedableRng;
use rand_chacha::ChaCha8Rng;

fn main() {
    let mut rng = ChaCha8Rng::seed_from_u64(987); // Draws [0.24346048, true, false, true]
    let reaction = pcr::PcrParameters {
        sites: 300,
        efficiencies: vec![0.8; 30],
        errors: vec![0.0001; 30],
    };
    let tree = pcr::simulate_tree(reaction, 100, &mut rng);
    println!("{:?}", tree.observations);
    let _ = write_npy("../out/feat-mut_observations.npy", &tree.observations);
    println!("{:?}", tree.mutations);
    let _ = write_npy("../out/feat-mut_mutations.npy", &tree.mutations);
}
