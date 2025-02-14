//! Generate noisy adaptive immune receptor data
mod amplicon;
mod fastx;
mod noise;
mod parse;
mod pcr;
mod sample;
mod sequence;
mod zipf;

use rand::Rng;
use rand_chacha::rand_core::SeedableRng;
use rand_chacha::ChaCha8Rng;

fn main() {
    let mut rng = ChaCha8Rng::seed_from_u64(987); // Draws [0.24346048, true, false, true]
    println!("{:?}", rng.random::<f32>());
    println!("{:?}", rng.random_bool(0.5));
    println!("{:?}", rng.random_bool(0.5));
    println!("{:?}", rng.random_bool(0.5));
    println!("{:?}", rng.random::<f32>());
    println!("{:?}", rng.random::<f32>());
    println!("{:?}", rng.random_bool(0.5));
    println!("{:?}", rng.random_bool(0.5));
    println!("{:?}", rng.random_bool(0.5));
}
