//! Error profiles for pcr

use indexmap::IndexMap;
use ndarray::prelude::*;
use rand::Rng;

/// SparseTree with matrix field that will be populated by node IDs
pub struct SparseTree {
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

/// Updates SparseTree object with the results of the next PCR round
fn evolve_tree<R: Rng + ?Sized>(tree: &mut SparseTree, round: usize, efficiency: f32, rng: &mut R) {
    let current_nodes = tree.matrix.index_axis(Axis(0), round);
    let unique_nodes = get_uniques(&current_nodes.to_vec());
    let evolved_nodes = evolve_nodes(&unique_nodes, efficiency, rng);
    let evolve_map: IndexMap<u32, u32> = unique_nodes
        .iter()
        .zip(evolved_nodes.iter())
        .map(|(&orig, &evolved)| (orig, evolved))
        .collect();
    let updated_nodes = current_nodes
        .mapv(|node| *evolve_map.get(&node).unwrap() + ((rng.random::<f32>() < 0.5) as u32) * 1);
    tree.matrix
        .slice_mut(s![(round + 1) as usize, ..])
        .assign(&updated_nodes);
}

/// Create and evolve a SparseTree with `reads` through `rounds`
pub fn simulate_tree<R: Rng + ?Sized>(
    efficiencies: Vec<f32>,
    reads: u32,
    rng: &mut R,
) -> SparseTree {
    let mut tree = SparseTree::new(&reads, &efficiencies);
    for round in 0..efficiencies.len() {
        evolve_tree(&mut tree, round, efficiencies[round], rng);
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
        assert_eq!(tree.matrix, array![[1, 1, 1], [1, 1, 1]])
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
        let mut rng = ChaCha8Rng::seed_from_u64(9274); // Draws [0.5790517, 0.3029861, 0.11362618]
        let evolved_nodes = evolve_nodes(&unique_nodes, efficiency, &mut rng);
        assert_eq!(evolved_nodes, array![1, 4, 6])
    }

    #[test]
    fn test_evolve_tree() {}
}
