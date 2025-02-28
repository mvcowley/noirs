//! Coalescent PCR simulator

use indexmap::IndexMap;
use ndarray::prelude::*;
use rand::Rng;
use rand_distr::{Binomial, Distribution};
use std::collections::HashSet;

/// Holds parameters of polymerase chain reaction
pub struct Reaction {
    /// Length of the amplified molecule
    pub sites: u16,
    /// Vector of the reaction success probabilities of each cycle of PCR
    pub efficiencies: Vec<f32>,
    /// Vector of the error probabilities of each cycle of PCR
    pub errors: Vec<f64>,
}

/// SparseTree with matrix field that will be populated by node IDs
pub struct SparseTree {
    ///  An array representing the path through the sparse binary tree of
    ///  N observations of the sequencer. u32 as there can be (2^31)-1 leaf nodes
    pub observations: Array2<u32>,

    /// An array representing the mutations accumulated by each observation during evolution.
    /// As even 1 error per cycle is unlikely for a given amplicon, u8 would be ok, but u32 helps
    /// for later merge with other data for a given UMI
    pub mutations: Array2<u32>,
}

/// Functions to create and update the SparseTree matrix
impl SparseTree {
    /// Create a matrix of ones so that binary tree path calculation is easy e.g. 1*2 = 2
    fn new(reads: &u32, reaction: &Reaction) -> SparseTree {
        let cycles = reaction.efficiencies.len() + 1;
        SparseTree {
            observations: Array2::<u32>::ones((usize::try_from(*reads).unwrap(), cycles).f()),
            mutations: Array2::<u32>::zeros((usize::try_from(*reads).unwrap(), cycles).f()),
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

/// Check per unique node if reaction cycle is a success
fn check_reaction_success<R: Rng + ?Sized>(
    unique_nodes: &Vec<u32>,
    efficiency: f32,
    rng: &mut R,
) -> Array1<u32> {
    let rand_vec: Vec<f32> = (0..unique_nodes.len()).map(|_| rng.random()).collect();
    let rand_arr = Array1::from_vec(rand_vec);
    rand_arr.map(|x| (*x < efficiency) as u32)
}

/// Evolve nodes according to reaction successes
fn evolve_nodes(unique_nodes: &Vec<u32>, reaction_successes: &Array1<u32>) -> Array1<u32> {
    Array::from_vec(unique_nodes.to_vec()) * (reaction_successes + 1)
}

/// Returns a 1 for any non-zero number, or zero for zero
fn is_non_zero(n: u32) -> u32 {
    (n != 0) as u32
}

/// Simulates the number of mutations which occur per evolved node
fn simulate_mutations<R: Rng + ?Sized>(
    updated_nodes: &Array1<u32>,
    unique_nodes: &Vec<u32>,
    cycle: usize,
    reaction: &Reaction,
    rng: &mut R,
) -> Array1<u32> {
    let bin = Binomial::new(reaction.sites as u64, reaction.errors[cycle]).unwrap();
    let updated_uniques = get_uniques(&updated_nodes.to_vec());
    let set_current: HashSet<&u32> = unique_nodes.iter().collect();
    let evolved_uniques: Vec<&u32> = updated_uniques
        .iter()
        .filter(|x| !set_current.contains(x))
        .collect();
    let unique_mutations: Vec<u32> = evolved_uniques
        .iter()
        .map(|_| bin.sample(rng).try_into().unwrap())
        .collect();
    let mutation_map: IndexMap<&u32, u32> = evolved_uniques
        .iter()
        .zip(unique_mutations.iter())
        .map(|(&node, &mutation)| (node, mutation))
        .collect();
    updated_nodes.mapv(|node| *mutation_map.get(&node).unwrap_or_else(|| &0))
}

/// Updates SparseTree object with the results of the next PCR cycle
fn evolve_tree<R: Rng + ?Sized>(
    tree: &mut SparseTree,
    cycle: usize,
    reaction: &Reaction,
    rng: &mut R,
) {
    let current_nodes = tree.observations.index_axis(Axis(1), cycle);
    let unique_nodes = get_uniques(&current_nodes.to_vec());
    let reaction_successes =
        check_reaction_success(&unique_nodes, reaction.efficiencies[cycle], rng);
    let evolved_nodes = evolve_nodes(&unique_nodes, &reaction_successes);
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
    tree.observations
        .slice_mut(s![.., (cycle + 1) as usize])
        .assign(&updated_nodes);

    let mutations = simulate_mutations(&updated_nodes, &unique_nodes, cycle, reaction, rng);
    tree.mutations
        .slice_mut(s![.., (cycle + 1) as usize])
        .assign(&mutations);
}

/// Create and evolve a SparseTree with `reads` through `cycles`
pub fn simulate_tree<R: Rng + ?Sized>(
    reaction: &Reaction,
    reads: u32,
    rng: &mut R,
) -> SparseTree {
    let mut tree = SparseTree::new(&reads, &reaction);
    for cycle in 0..reaction.efficiencies.len() {
        evolve_tree(&mut tree, cycle, reaction, rng);
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
        let reaction = Reaction {
            sites: 12,
            efficiencies: vec![0.95; 2],
            errors: vec![0.00001; 2],
        };
        let tree = SparseTree::new(&reads, &reaction);
        assert_eq!(tree.observations, array![[1, 1, 1], [1, 1, 1]]);
        assert_eq!(tree.mutations, array![[0, 0, 0], [0, 0, 0]]);
    }

    #[test]
    fn test_get_uniques() {
        let nodes: Vec<u32> = vec![3, 2, 2, 1];
        let unique_nodes = get_uniques(&nodes);
        assert_eq!(unique_nodes, vec![1, 2, 3]);
    }

    #[test]
    fn test_check_reaction_success() {
        let unique_nodes: Vec<u32> = vec![1, 2, 3];
        let efficiency: f32 = 0.5;
        let mut rng = ChaCha8Rng::seed_from_u64(927); // Draws [0.4124102, 0.7353993, 0.8147212]
        let reaction_sucesses = check_reaction_success(&unique_nodes, efficiency, &mut rng);
        assert_eq!(reaction_sucesses, array![1, 0, 0]);
    }

    #[test]
    fn test_evolve_nodes() {
        let unique_nodes: Vec<u32> = vec![1, 2, 3];
        let reaction_successes = array![1, 0, 0];
        let evolved_nodes = evolve_nodes(&unique_nodes, &reaction_successes);
        assert_eq!(evolved_nodes, array![2, 2, 3]);
    }

    #[test]
    fn test_simulate_mutations() {
        let updated_nodes = array![4, 5, 3];
        let unique_nodes = vec![2, 3];
        let reaction = Reaction {
            sites: 3,
            efficiencies: vec![0.95; 2],
            errors: vec![0.5; 2],
        };
        let mut rng = ChaCha8Rng::seed_from_u64(927); // Draws [2, 1, 3, 2, 1, 0] from binomial with n=3 and p=0.5
        let cycle = 0;
        let mutated_nodes =
            simulate_mutations(&updated_nodes, &unique_nodes, cycle, &reaction, &mut rng);
        assert_eq!(mutated_nodes, array![2, 1, 0]);
    }

    #[test]
    fn test_is_non_zero() {
        assert_eq!(is_non_zero(0), 0);
        assert_eq!(is_non_zero(1), 1);
        assert_eq!(is_non_zero(2), 1);
    }

    #[test]
    fn test_evolve_tree_rep_success() {
        let reaction = Reaction {
            sites: 3,
            efficiencies: vec![0.5; 2],
            errors: vec![0.5; 2],
        };
        let mut tree = SparseTree::new(&3, &reaction);
        let mut rng = ChaCha8Rng::seed_from_u64(987); // Draws [0.24346048, true, false, true,
                                                      // 2, 2]
        let cycle = 0;
        evolve_tree(&mut tree, cycle, &reaction, &mut rng);
        assert_eq!(tree.observations, array![[1, 3, 1], [1, 2, 1], [1, 3, 1]]);
        assert_eq!(tree.mutations, array![[0, 2, 0], [0, 2, 0], [0, 2, 0]]);
    }

    #[test]
    fn test_evolve_tree_rep_fail() {
        let reaction = Reaction {
            sites: 3,
            efficiencies: vec![0.2; 2],
            errors: vec![0.5; 2],
        };
        let mut tree = SparseTree::new(&3, &reaction);
        let mut rng = ChaCha8Rng::seed_from_u64(987); // Draws [0.24346048, true, false, true,
                                                      // 2, 2]
        let cycle = 0;
        evolve_tree(&mut tree, cycle, &reaction, &mut rng);
        assert_eq!(tree.observations, array![[1, 1, 1], [1, 1, 1], [1, 1, 1]]);
        assert_eq!(tree.mutations, array![[0, 0, 0], [0, 0, 0], [0, 0, 0]]);
    }

    #[test]
    fn test_simulate_tree() {
        let mut rng = ChaCha8Rng::seed_from_u64(987);
        let reaction = Reaction {
            sites: 3,
            efficiencies: vec![0.3; 2],
            errors: vec![0.5; 2],
        };
        // Draws [0.24346048, true, false, true, 2, 2, 0.06048143, 0.4701972, false, true,
        // false, 1]
        let tree = simulate_tree(&reaction, 3, &mut rng);
        assert_eq!(tree.observations, array![[1, 3, 3], [1, 2, 5], [1, 3, 3]]);
        assert_eq!(tree.mutations, array![[0, 2, 0], [0, 2, 1], [0, 2, 0]]);
    }
}
