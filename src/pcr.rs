//! Error profiles for pcr

use ndarray::{array, s, Array1, Array2, ArrayView1, ArrayView2};
use rand::{thread_rng, Rng};

struct SparseTree {
    matrix: Array2<u64>,
}

// SparseTree::matrix is an array of u64 as there can be 2^30 leaf nodes
// Matrix is ones so that binary tree path calculation is easy e.g. 1*2 = 2
impl SparseTree {
    fn new(rounds: &Vec<f32>, reads: &u32) -> Self {
        let axis1 = rounds.len() + 1;
        Self {
            matrix: Array2::<u64>::ones((reads as int, axis1 as int)),
        }
    }
    fn update(&mut self, read: u32, path: Array1<u64>) {
        self.matrix.slice(s![read as int]).assign_to(path);
    }
}

// If next round for other reads has the same node id, transcript did not participate in this round
// of PCR.
fn skip_round(next_round: ArrayView2<'_,u64>, current_node: &u64) -> bool {
    next_round.iter().any(|&x| x == *current_node)
}

fn calc_node(current_node: u64, child: f32) -> u64 {
    let next_node: u64 = current_node * 2;
    if child < 0.5 {
        next_node
    } else {
        next_node + 1
    }
}

pub fn trace_path(tree: &SparseTree, rounds: &Vec<f32>) -> Array1<u64> {
    let mut rng = thread_rng();
    let mut path: Array1<u64> = Array1::<u64>::ones(rounds.len());
    for (i, efficiency) in rounds.iter().enumerate() {
        let next_round = tree.matrix.slice(s![.., i + 1]);
        let current_node = path[[i]];
        let replicate: f32 = rng.gen();
        if skip_round(next_round, &current_node) {
            path[[i + 1]] = current_node;
        } else {
            if replicate > *efficiency {
                let child: f32 = rng.gen();
                let next_node: u64 = calc_node(current_node, child);
                path[[i + 1]] = next_node;
            } else {
                path[[i + 1]] = current_node;
            }
        }
    }
    path
}

pub fn construct_tree(rounds: Vec<f32>, reads: u32) -> SparseTree {
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
