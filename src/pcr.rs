//! Error profiles for pcr
//! SparseTree::matrix is an array of u64 as there can be 2^30 leaf nodes

use ndarray::{s, Array1, Array2};

struct SparseTree {
    matrix: Array2<u64>,
}

impl SparseTree {
    fn new(rounds: Vec<f32>, reads: u32) -> Self {
        Self {
            matrix: Array2::<u64>::default((reads, rounds.len())),
        }
    }
    fn update(&mut self, read: u32, path: ndarray::Array1<u64>) {
        self.matrix.slice(s![read, ...]).assign_to(path);
    }
}

pub fn trace_path(tree: &SparseTree, &rounds: Vec<f32>) -> Array1<u64> {}

pub fn construct_tree(rounds: Vec<f32>, reads: u32) -> SparseTree {
    let mut tree = SparseTree::new(rounds, reads);
    for read in 0..=reads {
        let path = trace_path(&tree, &rounds);
        tree.update(read, path);
    }
    tree
}
