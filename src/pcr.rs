//! Error profiles for pcr

use indexmap::IndexMap;

use ndarray::prelude::*;
use ndarray_rand::{rand_distr, RandomExt};
use rand::Rng;

/// SparseTree with matrix field that will be populated by node IDs
struct SparseTree {
    /// matrix is an array of u32 as there can be (2^31)-1 leaf nodes
    matrix: Array2<u32>,
}

/// Functions to create and update the SparseTree matrix
impl SparseTree {
    /// Create a matrix of ones so that binary tree path calculation is easy e.g. 1*2 = 2
    fn new(reads: &u32, rounds: &Vec<f32>) -> SparseTree {
        let axis1 = rounds.len() + 1;
        SparseTree {
            matrix: Array2::<u32>::ones((usize::try_from(*reads).unwrap(), axis1).f()),
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
fn evolve_nodes(unique_nodes: &Vec<u32>, efficiency: f32) -> Array1<u32> {
    let mut rng = rand::thread_rng();
    let rand_arr = Array::random(unique_nodes.len(), rand_distr::Uniform::new(0., 1.));
    let evo_arr =
        rand_arr.map(|x| ((*x < efficiency) as u32) * 2 + (1 * (rng.gen_bool(0.5) as u32)));
    let new_nodes = Array::from_vec(unique_nodes.to_vec()) * evo_arr;
    new_nodes
}

/// Updates SparseTree object with the results of the next PCR round
fn evolve_tree(tree: &mut SparseTree, round: usize, efficiency: f32) {
    let current_nodes = tree.matrix.index_axis(Axis(0), round);
    let unique_nodes = get_uniques(&current_nodes.to_vec());
    let evolved_nodes = evolve_nodes(&unique_nodes, efficiency);
    let evolve_map: IndexMap<u32, u32> = unique_nodes
        .iter()
        .zip(evolved_nodes.iter())
        .map(|(&orig, &evolved)| (orig, evolved))
        .collect();
    let updated_nodes = current_nodes.mapv(|node| *evolve_map.get(&node).unwrap());
    tree.matrix
        .slice_mut(s![(round + 1) as usize, ..])
        .assign(&updated_nodes);
}

/// Create and evolve a SparseTree with `reads` through `rounds`
pub fn simulate_tree(efficiencies: Vec<f32>, reads: u32) -> SparseTree {
    let mut tree = SparseTree::new(&reads, &efficiencies);
    for round in 0..efficiencies.len() {
        evolve_tree(&mut tree, round, efficiencies[round]);
    }
    tree
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn initialise_tree() {
        let rounds: Vec<f32> = vec![0.95, 0.95];
        let reads: u32 = 2;
        let tree = SparseTree::new(&reads, &rounds);
        assert_eq!(tree.matrix, array![[1, 1], [1, 1]])
    }

    // #[test]
    // fn create_tree() {
    //     let rounds = vec![0.95, 0.95, 0.95];
    //     let reads = 5;
    //     let tree = construct_tree(rounds, reads);
    // }
    //
    //
}
