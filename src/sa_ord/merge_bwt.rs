use std::collections::HashMap;
use std::hash::Hash;
use std::clone::Clone;
use std::fmt::Debug;
use bitvec::vec::BitVec;
use itertools::Itertools;


pub fn generate_f_vec<T: Ord + Hash + Clone + Debug>(
    bwts: &[&[T]],
) -> HashMap<T, usize>{
    let mut num_occurrences: HashMap<T, usize> = HashMap::new();
    for &bwt in bwts.iter() {
        for c in bwt {
            if let Some(i) = num_occurrences.get_mut(c) {
                *i += 1;
            } else {
                num_occurrences.insert(c.clone(), 1);
            }
        }
    }
    let ordered_chars = num_occurrences.keys().sorted().cloned().collect::<Vec<T>>();
    let mut total = 1usize;
    let mut f_vec = HashMap::with_capacity(num_occurrences.len());
    for c in ordered_chars {
        f_vec.insert(c, total);
        total += num_occurrences[&c];
    }
    f_vec
}

pub trait BwtMerger<T: Ord + Hash + Clone + Debug> {
    type BinTrackingType where Self::BinTrackingType: S
    /// Smallest character in the alphabet being considered
    const TERMINATOR_CHAR: T;

    fn generate_z_vec(
        bwts: &[&[T]], 
        previous_z_vec: &[Self::BinTrackingType],
        chunk_size: usize,
    ) -> Vec<Self::BinTrackingType>
    where usize: std::convert::From<Self::BinTrackingType> + std::convert::TryInto<Self::BinTrackingType> + Clone
    {
        let mut f_vec = generate_f_vec(bwts);
        let mut k_i = vec![0; bwts.len()];
        
        let total_len = bwts.iter().map(|&x| x.len()).sum();
        let new_z_vec: Vec<Self::BinTrackingType> = vec![0; total_len];
        for k in 0..total_len {
            let b = usize::from(previous_z_vec[k]);
            let bwt = bwts[b];
            k_i[b] += 1;
            let c = &bwt[k_i[b]];
            let j = if *c == Self::TERMINATOR_CHAR {
                b
            } else {
                f_vec[&c] += 1;
                f_vec[&c]
            };
            new_z_vec[j] = b;
        }
        new_z_vec
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let result = 2 + 2;
        assert_eq!(result, 4);
    }
}


