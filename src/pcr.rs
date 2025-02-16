//! Coalescent PCR simulator

use indexmap::IndexMap;
use ndarray::prelude::*;
use rand::Rng;

/// SparseTree with matrix field that will be populated by node IDs
pub struct SparseTree {
    ///  An array representing the path through the sparse binary tree of
    ///  N observations of the sequencer. u32 as there can be (2^31)-1 leaf nodes
    pub observations: Array2<u32>,

    /// An array representing the mutations accumulated by each observation during evolution.
    /// As even 1 error per round is unlikely for a given amplicon, u8 is used.
    pub mutations: Array2<u8>
}

/// Functions to create and update the SparseTree matrix
impl SparseTree {

    /// Create a matrix of ones so that binary tree path calculation is easy e.g. 1*2 = 2
    fn new(reads: &u32, rounds: &Vec<f32>) -> SparseTree {
        let axis1 = rounds.len() + 1;
        SparseTree {
            observations: Array2::<u32>::ones((usize::try_from(*reads).unwrap(), axis1).f()),
            mutations: Array2::<u8>::zeros((usize::try_from(*reads).unwrap(), axis1).f()),
        }
    }
}

/// Drop duplicates from node vector
fn get_uniques(current_nodes: &Vec<u32>) -> Vec<u32> {
    let mut nodes = current_nodes.clone();
    nodes.sort();
    nodes.dedup();
    nodes
}

/// Evolve unique nodes to next round
fn evolve_nodes<R: Rng + ?Sized>(
    unique_nodes: &Vec<u32>,
    efficiency: f32,
    rng: &mut R,
) -> Array1<u32> {
    let rand_vec: Vec<f32> = (0..unique_nodes.len()).map(|_| rng.random()).collect();
    let rand_arr = Array1::from_vec(rand_vec);
    let evo_arr = rand_arr.map(|x| ((*x < efficiency) as u32) + 1);
    let new_nodes = Array::from_vec(unique_nodes.to_vec()) * evo_arr;
    new_nodes
}

/// Returns a 1 for any non-zero number, or zero for zero
fn is_non_zero(n: u32) -> u32 {
    (n != 0) as u32
}

/// Updates SparseTree object with the results of the next PCR round
fn evolve_tree<R: Rng + ?Sized>(tree: &mut SparseTree, round: usize, efficiency: f32, error: f64, rng: &mut R) {
    let current_nodes = tree.observations.index_axis(Axis(1), round);
    let unique_nodes = get_uniques(&current_nodes.to_vec());
    let evolved_nodes = evolve_nodes(&unique_nodes, efficiency, rng);
    let evolve_map: IndexMap<u32, u32> = unique_nodes
        .iter()
        .zip(evolved_nodes.iter())
        .map(|(&orig, &evolved)| (orig, evolved))
        .collect();
    let updated_nodes = current_nodes.mapv(|node| {
        let new = *evolve_map.get(&node).unwrap();
        let evolved = is_non_zero(new - &node);
        new + (evolved * (rng.random_bool(0.5) as u32))
    });
    let mutations = current_nodes.mapv(|node| {
        let new = *evolve_map.get(&node).unwrap();
        let evolved = is_non_zero(new - &node);
        // TODO: Allow for multiple mutations per round: draw from distribution
        (evolved * (rng.random_bool(error) as u32)) as u8
    });
    tree.observations
        .slice_mut(s![.., (round + 1) as usize])
        .assign(&updated_nodes);
    tree.mutations
        .slice_mut(s![.., (round + 1) as usize])
        .assign(&mutations);
}

/// Create and evolve a SparseTree with `reads` through `rounds`
pub fn simulate_tree<R: Rng + ?Sized>(
    efficiencies: Vec<f32>,
    errors: Vec<f64>,
    reads: u32,
    rng: &mut R,
) -> SparseTree {
    let mut tree = SparseTree::new(&reads, &efficiencies);
    for round in 0..efficiencies.len() {
        evolve_tree(&mut tree, round, efficiencies[round], errors[round], rng);
    }
    tree
}

#[cfg(test)]
mod tests {

    use super::*;
    use rand_chacha::rand_core::SeedableRng;
    use rand_chacha::ChaCha8Rng;

    #[test]
    fn test_new_tree() {
        let reads: u32 = 2;
        let efficiencies: Vec<f32> = vec![0.95, 0.95];
        let tree = SparseTree::new(&reads, &efficiencies);
        assert_eq!(tree.observations, array![[1, 1, 1], [1, 1, 1]]);
        assert_eq!(tree.mutations, array![[0, 0, 0], [0, 0, 0]]);
    }

    #[test]
    fn test_get_uniques() {
        let nodes: Vec<u32> = vec![3, 2, 2, 1];
        let unique_nodes = get_uniques(&nodes);
        assert_eq!(unique_nodes, vec![1, 2, 3])
    }

    #[test]
    fn test_evolve_nodes() {
        let unique_nodes: Vec<u32> = vec![1, 2, 3];
        let efficiency: f32 = 0.5;
        let mut rng = ChaCha8Rng::seed_from_u64(927); // Draws [0.4124102, 0.7353993, 0.8147212]
        let evolved_nodes = evolve_nodes(&unique_nodes, efficiency, &mut rng);
        assert_eq!(evolved_nodes, array![2, 2, 3])
    }

    #[test]
    fn test_is_non_zero() {
        assert_eq!(is_non_zero(0), 0);
        assert_eq!(is_non_zero(1), 1);
        assert_eq!(is_non_zero(2), 1);
    }

    // #[test]
    // fn test_evolve_tree_rep_success() {
    //     let efficiencies = vec![0.5, 0.5];
    //     let mut tree = SparseTree::new(&3, &efficiencies);
    //     let mut rng = ChaCha8Rng::seed_from_u64(987); // Draws [0.24346048, true, false, true]
    //     let round = 0;
    //     evolve_tree(&mut tree, round, efficiencies[round], &mut rng);
    //     assert_eq!(tree.observations, array![[1, 3, 1], [1, 2, 1], [1, 3, 1]]);
    // }

    // #[test]
    // fn test_evolve_tree_rep_fail() {
    //     let efficiencies = vec![0.2, 0.5];
    //     let mut tree = SparseTree::new(&3, &efficiencies);
    //     let mut rng = ChaCha8Rng::seed_from_u64(987); // Draws [0.24346048, true, false, true]
    //     let round = 0;
    //     evolve_tree(&mut tree, round, efficiencies[round], &mut rng);
    //     assert_eq!(tree.observations, array![[1, 1, 1], [1, 1, 1], [1, 1, 1]]);
    // }

//     #[test]
//     fn test_simulate_tree() {
//         let efficiencies = vec![0.8, 0.8];
//         let mut rng = ChaCha8Rng::seed_from_u64(987);
//         // Draws [0.24346048, true, false, true, 0.84027797, 0.7916585, false, true, false]
//         let tree = simulate_tree(efficiencies, 3, &mut rng);
//         assert_eq!(tree.observations, array![[1, 3, 6], [1, 2, 2], [1, 3, 6]]);
//     }
}
