//! Error profiles for pcr

use ndarray::{s, Array1, Array2};
use rand::{Rng, thread_rng};


struct SparseTree {
    matrix: Array2<u64>,
}

// SparseTree::matrix is an array of u64 as there can be 2^30 leaf nodes
// Matrix is ones so that binary tree path calculation is easy e.g. 1*2 = 2
impl SparseTree {
    fn new(rounds: Vec<f32>, reads: u32) -> Self {
        Self {
            matrix: Array2::<u64>::ones((reads, rounds.len() + 1)),
        }
    }
    fn update(&mut self, read: u32, path: Array1<u64>) {
        self.matrix.slice(s![read, ...]).assign_to(path);
    }
}

fn skip_round(next_state: &ArrayView2<u64> &current_node: &u64) -> bool {
    
}

fn calc_node(current_node: u64, child: f32) -> u64 {
    let next_node: u64 = current_node * 2;
    if child < 0.5 {
        next_node
    }
    else {
        next_node + 1
    }
}

pub fn trace_path(tree: &SparseTree, &rounds: Vec<f32>) -> Array1<u64> {
    let mut rng = thread_rng();
    let mut path: Array1<u64> = Array1::<u64>::ones(rounds.len())
    for (i, efficiency) in rounds.iter().enumerate() {
        let next_state = tree.matrix.slice(s![.., i + 1]);
        let current_node: u64 = path.slice(s![i]);
        let replicate: f32 = rng.gen();
        if skip_round(&next_state, &current_node) {
            path.slice_mut(s![i +1]).assign(&current_node)
        }
        else {
            if replicate > efficiency {
                let child: f32 = rng.gen();
                let next_node: u64 = calc_node(current_node, child);
                path.slice_mut(s![i + 1]).assign(&next_node)
            }
            else {
                path.slice_mut(s![i +1]).assign(&current_node)
            }
        }
    }
    path
}

pub fn construct_tree(rounds: Vec<f32>, reads: u32) -> SparseTree {
    let mut tree = SparseTree::new(rounds, reads);
    for read in 0..=reads {
        let path = trace_path(&tree, &rounds);
        tree.update(read, path);
    }
    tree
}
