//! Error profiles for pcr

use ndarray::prelude::*;
use ndarray_rand::{rand_distr, RandomExt};

/// SparseTree with matrix field that will be populated by node IDs
struct SparseTree {
    /// matrix is an array of u32 as there can be (2^31)-1 leaf nodes
    matrix: Array2<u32>,
}

/// Functions to create and update the SparseTree matrix
impl SparseTree {
    /// Create a matrix of ones so that binary tree path calculation is easy e.g. 1*2 = 2
    fn new(rounds: &Vec<f32>, reads: &u32) -> SparseTree {
        let axis1 = rounds.len() + 1;
        SparseTree {
            matrix: Array2::<u32>::ones((usize::try_from(*reads).unwrap(), axis1).f()),
        }
    }
    /// Update a slice of the matrix with the next step of the simulation
    fn update(&mut self, read: u32, path: Array1<u32>) {
        self.matrix
            .slice_mut(s![usize::try_from(read).unwrap(), ..])
            .assign(&path);
    }
}

// If next round for other reads has the same node id, transcript did not participate in this round
// of PCR.
// fn in_next_round(next_round: ArrayView1<u64>, current_node: &u64) -> bool {
//     next_round.iter().any(|&x| x == *current_node)
// }
//
// fn calc_node(current_node: u64, child: f32) -> u64 {
//     let next_node: u64 = current_node * 2;
//     if child < 0.5 {
//         next_node
//     } else {
//         next_node + 1
//     }
// }

/// Drop duplicates from node vector
fn get_uniques(current_nodes: &Vec<u32>) -> Vec<u32> {
    let mut nodes = current_nodes.clone();
    nodes.sort();
    nodes.dedup();
    nodes
}

/// Evolve unique nodes to next round
fn evolve_nodes(unique_nodes: Vec<u32>, efficiency: f32) -> Array1<u32> {
    let rand_arr = Array::random(unique_nodes.len(), rand_distr::Uniform::new(0., 1.));
    let evo_arr = rand_arr.map(|x| ((*x < efficiency) as u32) * 2);
    let new_nodes = Array::from_vec(unique_nodes) * evo_arr;
    new_nodes
}

pub fn evolve_tree(tree: &mut SparseTree, round: u8, efficiency: f32) {
    let current_nodes = tree.matrix.index_axis(Axis(0), round as usize).to_vec();
    let unique_nodes = get_uniques(&current_nodes);
    let evolved_nodes = evolve_nodes(unique_nodes, efficiency); 
}

// pub fn trace_path(
//     tree: &SparseTree,
//     efficiencies: &Vec<f32>,
// ) -> ArrayBase<OwnedRepr<u64>, Dim<[usize; 1]>> {
//     let mut rng = thread_rng();
//     let mut path: Array1<u64> = Array1::<u64>::ones(efficiencies.len());
//     for (i, efficiency) in efficiencies.iter().enumerate() {
//         let next_round = tree.matrix.slice(s![.., i + 1]);
//         let current_node = path[[i]];
//         if in_next_round(next_round, &current_node) {
//             path[[i + 1]] = current_node;
//         } else {
//             let replicate: f32 = rng.gen();
//             if replicate > *efficiency {
//                 let child: f32 = rng.gen();
//                 let next_node: u64 = calc_node(current_node, child);
//                 path[[i + 1]] = next_node;
//             } else {
//                 path[[i + 1]] = current_node;
//             }
//         }
//     }
//     path
// }
//
pub fn simulate_tree(rounds: Vec<f32>, reads: u32) -> SparseTree {
    let mut tree = SparseTree::new(&rounds, &reads);
    for read in 0..=reads {
        let path = trace_path(&tree, &rounds);
        tree.update(read, path);
    }
    tree
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn initialise_tree() {
        let rounds = vec![0.95, 0.95];
        let reads = 2;
        let tree = SparseTree::new(&rounds, &reads);
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
