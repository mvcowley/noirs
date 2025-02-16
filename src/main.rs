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
use rand_distr::{Binomial, Distribution};
use rand::Rng;
use ndarray_npy::write_npy;

fn main() {
    let mut rng = ChaCha8Rng::seed_from_u64(987); // Draws [0.24346048, true, false, true]
    println!("{:?}", rng.random::<f32>());
    println!("{:?}", rng.random_bool(0.5));
    println!("{:?}", rng.random_bool(0.5));
    println!("{:?}", rng.random_bool(0.5));
    let bin = Binomial::new(3, 0.5).unwrap();
    println!("{:?}", bin.sample(&mut rng));
    println!("{:?}", bin.sample(&mut rng));
    // println!("{:?}", bin.sample(&mut rng));
    // let efficiencies = vec![0.8; 30];
    // let tree = pcr::simulate_tree(efficiencies, 100, &mut rng);
    // println!("{:?}", tree.observations);
    // let _ = write_npy("../out/molecule.npy", &tree.observations);
}
