//! FROZEN pre-optimization implementation of `sa_ord::sais`, snapshotted from git commit 1df94d4
//! (HEAD at the time of writing). It is used purely as an oracle by the equivalence tests.
//! DO NOT "improve" or refactor this file.
//!
//! The only deviations from the original source are:
//!  * imports (`sa_ord::errors` instead of the local `errors` module),
//!  * lint allowances so the verbatim code is warning free in the test crate,
//!  * the public entry point `reference_sais` wraps the original recursive function (renamed
//!    `sais_impl`) which additionally carries a `depth` argument, purely so tests can observe
//!    recursion depth via `reference_max_depth`. No algorithmic change.
#![allow(dead_code, clippy::all, clippy::pedantic, clippy::nursery)]

use bitvec::vec::BitVec;
use core::hash::Hash;
use itertools::Itertools;
use sa_ord::errors::{Error, Result};
use std::cell::Cell;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::fmt::Debug;

thread_local! {
    static MAX_DEPTH: Cell<usize> = Cell::new(0);
}

/// Frozen reference implementation (oracle).
pub fn reference_sais<S: Ord + Clone + Hash + Debug>(text: &[S]) -> Result<Vec<usize>> {
    sais_impl(text, 1)
}

/// Runs the reference and returns the maximum recursion depth reached (1 = no recursion,
/// 2 = one recursive call on the reduced string, ...).
pub fn reference_max_depth<S: Ord + Clone + Hash + Debug>(text: &[S]) -> usize {
    MAX_DEPTH.with(|d| d.set(0));
    let _ = reference_sais(text);
    MAX_DEPTH.with(|d| d.get())
}

/// Position types (L or S).
#[derive(Debug)]
struct SuffixTypes {
    /// `BitVec`, with 1-bits indicating S type and 0 bits L-type
    types: BitVec,
}

impl SuffixTypes {
    fn new<S: Ord + Clone + Hash + Debug>(text: &[S]) -> Self {
        let text_len = text.len();
        let mut types: BitVec = BitVec::with_capacity(text_len);
        for _ in 0..text_len - 1 {
            types.push(false);
        }

        types.push(true);

        for text_index in (0..text_len - 1).rev() {
            if text[text_index] == text[text_index + 1] {
                // if the characters are equal, the next position determines
                // the lexicographical order
                let v = types[text_index + 1];
                types.set(text_index, v);
            } else {
                types.set(text_index, text[text_index] < text[text_index + 1]);
            }
        }
        SuffixTypes { types }
    }

    /// Check if p is S-position.
    fn is_stype(&self, p: usize) -> bool {
        self.types[p]
    }

    /// Check if p is L-position.
    fn is_ltype(&self, p: usize) -> bool {
        !(self.types[p])
    }

    /// Check if p is LMS-position.
    fn is_leftmost_stype(&self, p: usize) -> bool {
        p != 0 && self.is_stype(p) && self.is_ltype(p - 1)
    }
}

#[derive(Debug, Clone, Hash, PartialEq, Eq)]
enum SortableSuffixType<S: Ord + Clone + Hash + Debug> {
    /// 'Smaller' type, indicating that the suffix is smaller than the suffix to its right
    SType(S),
    /// 'Larger' type, indicating that the suffix is larger than the suffix to its right
    LType(S),
}

impl<S: Ord + Clone + Hash + Debug> PartialOrd for SortableSuffixType<S> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<S: Ord + Clone + Hash + Debug> Ord for SortableSuffixType<S> {
    fn cmp(&self, other: &Self) -> Ordering {
        match (self, other) {
            (SortableSuffixType::SType(x), SortableSuffixType::LType(y)) => match x.cmp(y) {
                Ordering::Equal => Ordering::Greater,
                cmp => cmp,
            },
            (SortableSuffixType::LType(x), SortableSuffixType::SType(y)) => match x.cmp(y) {
                Ordering::Equal => Ordering::Less,
                cmp => cmp,
            },
            (
                SortableSuffixType::SType(x) | SortableSuffixType::LType(x),
                SortableSuffixType::SType(y) | SortableSuffixType::LType(y),
            ) => x.cmp(y),
        }
    }
}

/// Finds the start and end indices of the LMS substrings in the provided arra of suffix classifications
fn find_lms_substring_indices(suffix_classifications: &SuffixTypes) -> Vec<(usize, usize)> {
    let lms_num = (0..suffix_classifications.types.len())
        .filter(|idx| suffix_classifications.is_leftmost_stype(*idx))
        .count();
    let mut wstrings: Vec<(usize, usize)> = Vec::with_capacity(lms_num);
    let mut start = None;
    for current_index in 0..suffix_classifications.types.len() {
        if suffix_classifications.is_leftmost_stype(current_index) {
            if let Some(current_start) = start {
                wstrings.push((current_start, current_index + 1));
            }
            start = Some(current_index);
        }
    }
    wstrings.push((suffix_classifications.types.len() - 1, suffix_classifications.types.len()));
    wstrings
}

/// Counts the number of characters of each value and type (L or S).
fn count_text_bin_sizes<S: Ord + Clone + Hash + Debug>(
    text: &[S],
    suffix_classifications: &SuffixTypes,
) -> HashMap<SortableSuffixType<S>, usize> {
    let mut bin_sizes: HashMap<SortableSuffixType<S>, usize> = HashMap::new();
    for (index, suffix_type) in suffix_classifications.types.iter().enumerate() {
        let char_var = if *suffix_type {
            SortableSuffixType::SType(text[index].clone())
        } else {
            SortableSuffixType::LType(text[index].clone())
        };
        if let Some(x) = bin_sizes.get_mut(&char_var) {
            *x += 1;
        } else {
            bin_sizes.insert(char_var, 1);
        }
    }
    bin_sizes
}

/// Derives the bin offsets for the suffix array given the list of bin sizes.
fn get_bin_offsets<S: Ord + Clone + Hash + Debug>(
    bin_sizes: &HashMap<SortableSuffixType<S>, usize>,
    forward_edge: bool,
) -> HashMap<SortableSuffixType<S>, usize> {
    let ordered_chars = bin_sizes.keys().cloned().sorted().collect::<Vec<SortableSuffixType<S>>>();
    let mut running_counter = 0;
    let mut bin_offsets: HashMap<SortableSuffixType<S>, usize> =
        HashMap::with_capacity(ordered_chars.len());

    for x_type in &ordered_chars {
        if let Some(x) = bin_sizes.get(x_type) {
            if forward_edge {
                bin_offsets.insert(x_type.clone(), running_counter);
                running_counter += *x;
            } else {
                running_counter += *x;
                bin_offsets.insert(x_type.clone(), running_counter - 1);
            }
        }
    }
    bin_offsets
}

/// Reduces the provided text into a reduced substring in which LMS substrings are represented by
/// their indexing in the ordering of the LMS substrings.
fn get_reduced_substring<S: Ord + Clone + Hash + Debug>(
    text: &[S],
    w_slice_indices: &[(usize, usize)],
    stypes: &SuffixTypes,
    bin_sizes: &HashMap<SortableSuffixType<S>, usize>,
) -> Result<(bool, Vec<usize>)> {
    let mut sort_sa = Vec::with_capacity(text.len());
    for _ in 0..text.len() {
        sort_sa.push(-1);
    }
    let mut bin_edges = get_bin_offsets(bin_sizes, false);
    for (start_idx, _) in w_slice_indices {
        let x = bin_edges.get_mut(&SortableSuffixType::SType(text[*start_idx].clone())).ok_or(
            Error::ValueNotInBinEdges {
                val: format!("{:?}", SortableSuffixType::SType(text[*start_idx].clone())),
            },
        )?;

        sort_sa[*x] = isize::try_from(*start_idx).map_err(|e| Error::TryFromIntError { e })?;
        if *x > 0 {
            *x -= 1;
        }
    }

    suffix_array_modification_propagation(&mut sort_sa, text, stypes, bin_sizes, true)?;
    suffix_array_modification_propagation(&mut sort_sa, text, stypes, bin_sizes, false)?;

    let validated_sa = sort_sa
        .iter()
        .map(|i| usize::try_from(*i).map_err(|e| Error::TryFromIntError { e }))
        .collect::<Result<Vec<usize>>>()?;

    let substring_start_to_ends: HashMap<usize, usize> = w_slice_indices.iter().copied().collect();
    let mut unique_substring_counter = 0;
    let mut reduced_substring_mapping: HashMap<usize, usize> =
        HashMap::with_capacity(w_slice_indices.len());

    let mut previous_substring: Option<&[S]> = None;
    for text_index in &validated_sa {
        if stypes.is_leftmost_stype(*text_index) {
            let substring_end = substring_start_to_ends
                .get(text_index)
                .ok_or(Error::LmsIndexError { beginning: *text_index })?;
            let current_substring = &text[*text_index..*substring_end];
            if previous_substring.is_some() && previous_substring != Some(current_substring) {
                unique_substring_counter += 1;
            }
            previous_substring = Some(current_substring);
            reduced_substring_mapping.insert(*text_index, unique_substring_counter);
        }
    }
    let reduced_substring = w_slice_indices
        .iter()
        .map(|(lms_substring_start, _)| {
            Ok(*reduced_substring_mapping
                .get(lms_substring_start)
                .ok_or(Error::LmsReductionError { text_index: *lms_substring_start })?)
        })
        .collect::<Result<Vec<usize>>>()?;

    Ok((reduced_substring.len() > unique_substring_counter + 1, reduced_substring))
}

/// Modifies the provided suffix array by doing either a forward or reverse pass.
fn suffix_array_modification_propagation<S: Ord + Clone + Hash + Debug>(
    sa: &mut Vec<isize>,
    text: &[S],
    stypes: &SuffixTypes,
    bin_sizes: &HashMap<SortableSuffixType<S>, usize>,
    forward_pass: bool,
) -> Result<()> {
    let mut bin_edges = get_bin_offsets(bin_sizes, forward_pass);

    let iter_order: Box<dyn Iterator<Item = usize>> =
        if forward_pass { Box::new(0..sa.len()) } else { Box::new((0..sa.len()).rev()) };
    for index in iter_order {
        let mut val = sa[index];
        if val == -1 {
            continue;
        };
        if val == 0 {
            val = isize::try_from(sa.len()).map_err(|e| Error::TryFromIntError { e })? - 1;
        } else {
            val -= 1;
        };
        let sortable_suffix_type = if forward_pass {
            SortableSuffixType::LType(
                text[usize::try_from(val).map_err(|e| Error::TryFromIntError { e })?].clone(),
            )
        } else {
            SortableSuffixType::SType(
                text[usize::try_from(val).map_err(|e| Error::TryFromIntError { e })?].clone(),
            )
        };
        let check_function_result = if forward_pass {
            stypes.is_ltype(usize::try_from(val).map_err(|e| Error::TryFromIntError { e })?)
        } else {
            stypes.is_stype(usize::try_from(val).map_err(|e| Error::TryFromIntError { e })?)
        };
        if check_function_result {
            let edge = bin_edges.get_mut(&sortable_suffix_type).ok_or_else(|| {
                Error::ValueNotInBinEdges { val: format!("{:?}", sortable_suffix_type) }
            })?;
            sa[*edge] = val;
            if forward_pass {
                *edge += 1;
            } else if *edge > 0 {
                *edge -= 1;
            }
        }
    }
    Ok(())
}

/// Performs SA-IS (original `sais`, with an added `depth` argument for instrumentation).
fn sais_impl<S: Ord + Clone + Hash + Debug>(text: &[S], depth: usize) -> Result<Vec<usize>> {
    MAX_DEPTH.with(|d| d.set(d.get().max(depth)));
    let mut sa: Vec<isize> = vec![-1; text.len()];
    // Step 1 - Classify all suffixes in the text with SType + LeftmostSType / LType label
    let stypes = SuffixTypes::new(text);
    // Step 2 - Divide the buckets into S- and L-Buckets on. TheL-Buckets in front of the S-Buckets.
    let bin_sizes = count_text_bin_sizes(text, &stypes);

    // Sort the LeftmostSType strings lexicographically and write them in their respective S-Buckets.
    let w_slice_indices = find_lms_substring_indices(&stypes);
    let (duplicates, reduced_substring) =
        get_reduced_substring(text, &w_slice_indices, &stypes, &bin_sizes)?;

    let mut reduced_suffix_array: Vec<usize>;
    if duplicates {
        reduced_suffix_array = sais_impl(&reduced_substring, depth + 1)?;
    } else {
        reduced_suffix_array = vec![0; reduced_substring.len()];
        for (original_index, &sort_index) in reduced_substring.iter().enumerate() {
            reduced_suffix_array[sort_index] = original_index;
        }
    }
    let mut bin_edges = get_bin_offsets(&bin_sizes, false);
    for reduced_string_index in reduced_suffix_array.iter().rev() {
        let w_slice_start = w_slice_indices[*reduced_string_index].0;
        let s_char = text[w_slice_start].clone();
        let x = bin_edges.get_mut(&SortableSuffixType::SType(s_char)).ok_or(
            Error::ValueNotInBinEdges {
                val: format!("{:?}", SortableSuffixType::SType(text[w_slice_start].clone())),
            },
        )?;

        sa[*x] = isize::try_from(w_slice_start).map_err(|e| Error::TryFromIntError { e })?;
        if *x > 0 {
            *x -= 1;
        }
    }

    suffix_array_modification_propagation(&mut sa, text, &stypes, &bin_sizes, true)?;
    suffix_array_modification_propagation(&mut sa, text, &stypes, &bin_sizes, false)?;

    return sa
        .iter()
        .map(|i| usize::try_from(*i).map_err(|e| Error::TryFromIntError { e }))
        .collect::<Result<Vec<usize>>>();
}
