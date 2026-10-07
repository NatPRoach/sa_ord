/// This module is a pure (attempted idiomatic) Rust implementation of the SA-IS suffix tree
/// algorithm for slices of Ord + Clone + Hash + Debug types
pub mod errors;

use bitvec::vec::BitVec;
use core::hash::Hash;
use errors::{Error, Result};
#[cfg(test)]
use itertools::Itertools;
#[cfg(test)]
use std::cmp::Ordering;
use std::collections::HashMap;
use std::fmt::Debug;

/// Position types (L or S). Implementation of Suffix Types `BitVec` partially borrowed from
/// [Rust Bio](https://github.com/rust-bio/rust-bio/blob/43d198e87c69867fad39f07374d3fdcd9d879729/src/data_structures/suffix_array.rs)
#[derive(Debug)]
struct SuffixTypes {
    /// `BitVec`, with 1-bits indicating S type and 0 bits L-type
    types: BitVec,
}

impl SuffixTypes {
    /// Instantiates a new `SuffixTypes` object. Characterizes the provided text into L or S type
    /// suffixes and stores these types in a `BitVec`
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

/// An object that allows for sorting of suffix types based on L or S substring properties for the
/// purposes of ordering 'bins' of the suffix array as described in the SA-IS implementation paper.
/// `Bins` encodes the same ordering as plain indices; this type is kept as its test reference.
#[cfg(test)]
#[derive(Debug, Clone, Hash, PartialEq, Eq)]
enum SortableSuffixType<S: Ord + Clone + Hash + Debug> {
    /// 'Smaller' type, indicating that the suffix is smaller than the suffix to its right
    SType(S),
    /// 'Larger' type, indicating that the suffix is larger than the suffix to its right
    LType(S),
}

#[cfg(test)]
impl<S: Ord + Clone + Hash + Debug> PartialOrd for SortableSuffixType<S> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[cfg(test)]
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
#[cfg(test)]
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
#[cfg(test)]
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

/// The suffix array bins ("buckets") of a text, laid out as plain vectors so the induction passes
/// never hash or clone characters.
///
/// Bins are numbered in suffix array order: the bin of a position is `2 * rank + is_stype`, where
/// `rank` is the rank of its character among the distinct characters of the text. This puts the
/// L-type bin of each character directly before its S-type bin, matching the ordering of
/// `SortableSuffixType`.
#[derive(Debug)]
struct Bins {
    /// Bin of each text position
    position_bins: Vec<usize>,
    /// Exclusive end offset of each bin in the suffix array
    ends: Vec<usize>,
}

impl Bins {
    /// Assigns every position of the text to its bin and counts the bin sizes.
    fn new<S: Ord + Hash>(text: &[S], suffix_classifications: &SuffixTypes) -> Self {
        // Number the distinct characters in order of first appearance, then rank them.
        let mut first_seen: HashMap<&S, usize> = HashMap::new();
        let mut position_bins: Vec<usize> = text
            .iter()
            .map(|character| {
                let next_id = first_seen.len();
                *first_seen.entry(character).or_insert(next_id)
            })
            .collect();
        let mut distinct: Vec<(&S, usize)> = first_seen.into_iter().collect();
        distinct.sort_unstable_by_key(|(character, _)| *character);
        let mut rank_of_id = vec![0; distinct.len()];
        for (rank, (_, id)) in distinct.iter().enumerate() {
            rank_of_id[*id] = rank;
        }
        for id in &mut position_bins {
            *id = rank_of_id[*id];
        }
        Self::from_ranks(position_bins, distinct.len(), suffix_classifications)
    }

    /// Builds the bins of a text that is already given as character ranks in
    /// `0..alphabet_size`, such as a reduced string, without hashing.
    fn from_ranks(
        ranks: Vec<usize>,
        alphabet_size: usize,
        suffix_classifications: &SuffixTypes,
    ) -> Self {
        let mut position_bins = ranks;
        let mut ends = vec![0; 2 * alphabet_size];
        for (index, bin) in position_bins.iter_mut().enumerate() {
            *bin = 2 * *bin + usize::from(suffix_classifications.is_stype(index));
            ends[*bin] += 1;
        }
        let mut running_counter = 0;
        for end in &mut ends {
            running_counter += *end;
            *end = running_counter;
        }
        Bins { position_bins, ends }
    }

    /// Offset of the first slot of every bin.
    fn heads(&self) -> Vec<usize> {
        let mut heads = Vec::with_capacity(self.ends.len());
        heads.push(0);
        heads.extend_from_slice(&self.ends[..self.ends.len() - 1]);
        heads
    }

    /// Offset of the last slot of every bin (empty bins are never written to).
    fn tails(&self) -> Vec<usize> {
        self.ends.iter().map(|end| end.saturating_sub(1)).collect()
    }
}

/// Reduces the provided text into a reduced substring in which LMS substrings are represented by
/// their indexing in the ordering of the LMS substrings.
#[cfg(test)]
fn get_reduced_substring<S: Ord + Clone + Hash + Debug>(
    text: &[S],
    w_slice_indices: &[(usize, usize)],
    stypes: &SuffixTypes,
    _bin_sizes: &HashMap<SortableSuffixType<S>, usize>,
) -> Result<(bool, Vec<usize>)> {
    reduce_lms_substrings(text, w_slice_indices, stypes, &Bins::new(text, stypes))
}

/// Reduces the provided text into a reduced substring in which LMS substrings are represented by
/// their indexing in the ordering of the LMS substrings.
fn reduce_lms_substrings<S: Ord + Clone + Hash + Debug>(
    text: &[S],
    w_slice_indices: &[(usize, usize)],
    stypes: &SuffixTypes,
    bins: &Bins,
) -> Result<(bool, Vec<usize>)> {
    let mut sort_sa = vec![-1; text.len()];
    let mut bin_edges = bins.tails();
    for (start_idx, _) in w_slice_indices {
        let x = &mut bin_edges[bins.position_bins[*start_idx]];

        sort_sa[*x] = isize::try_from(*start_idx).map_err(|e| Error::TryFromIntError { e })?;
        if *x > 0 {
            *x -= 1;
        }
    }

    induce_pass(&mut sort_sa, stypes, bins, true)?;
    induce_pass(&mut sort_sa, stypes, bins, false)?;

    let validated_sa = sort_sa
        .iter()
        .map(|i| usize::try_from(*i).map_err(|e| Error::TryFromIntError { e }))
        .collect::<Result<Vec<usize>>>()?;

    // Both indexed by text position: the end of the LMS substring starting there, and its name.
    let mut substring_ends: Vec<Option<usize>> = vec![None; text.len()];
    for &(start, end) in w_slice_indices {
        substring_ends[start] = Some(end);
    }
    let mut unique_substring_counter = 0;
    let mut reduced_substring_mapping: Vec<Option<usize>> = vec![None; text.len()];

    let mut previous_substring: Option<&[S]> = None;
    for text_index in &validated_sa {
        if stypes.is_leftmost_stype(*text_index) {
            let substring_end = substring_ends[*text_index]
                .ok_or(Error::LmsIndexError { beginning: *text_index })?;
            let current_substring = &text[*text_index..substring_end];
            if previous_substring.is_some() && previous_substring != Some(current_substring) {
                unique_substring_counter += 1;
            }
            previous_substring = Some(current_substring);
            reduced_substring_mapping[*text_index] = Some(unique_substring_counter);
        }
    }
    let reduced_substring = w_slice_indices
        .iter()
        .map(|(lms_substring_start, _)| {
            reduced_substring_mapping[*lms_substring_start]
                .ok_or(Error::LmsReductionError { text_index: *lms_substring_start })
        })
        .collect::<Result<Vec<usize>>>()?;

    Ok((reduced_substring.len() > unique_substring_counter + 1, reduced_substring))
}

/// Modifies the provided suffix array by doing either a forward or reverse pass, corresponding to
/// steps 2 or 3 in the induction of SA from SA of the reduced substring described in the SA-IS
/// paper.
///
/// # Errors
/// Can theoretically return an Err if something goes wrong with the internals of this crate, but
/// never should
#[cfg(test)]
fn suffix_array_modification_propagation<S: Ord + Clone + Hash + Debug>(
    sa: &mut [isize],
    text: &[S],
    stypes: &SuffixTypes,
    _bin_sizes: &HashMap<SortableSuffixType<S>, usize>,
    forward_pass: bool,
) -> Result<()> {
    induce_pass(sa, stypes, &Bins::new(text, stypes), forward_pass)
}

/// Modifies the provided suffix array by doing either a forward (L-type) or reverse (S-type)
/// induction pass, using the precomputed bins of the text.
///
/// # Errors
/// Can theoretically return an Err if something goes wrong with the internals of this crate, but
/// never should
fn induce_pass(
    sa: &mut [isize],
    stypes: &SuffixTypes,
    bins: &Bins,
    forward_pass: bool,
) -> Result<()> {
    let mut bin_edges = if forward_pass { bins.heads() } else { bins.tails() };

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
        let position = usize::try_from(val).map_err(|e| Error::TryFromIntError { e })?;
        let check_function_result =
            if forward_pass { stypes.is_ltype(position) } else { stypes.is_stype(position) };
        if check_function_result {
            let edge = &mut bin_edges[bins.position_bins[position]];
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

/// Performs SA-IS (Suffix Array - Induced Sorting) of a provided text.
///
/// # Errors
///
/// Will return `Error::TryFromIntError` if there are more LMS substrings than there are positive
/// numbers in isize (if so you likely have bigger problems).
/// May also Error if the last character of the text is not a sentinel character (i.e. the
/// unambiguous smallest character in the text).
pub fn sais<S: Ord + Clone + Hash + Debug>(text: &[S]) -> Result<Vec<usize>> {
    // Step 1 - Classify all suffixes in the text with SType + LeftmostSType / LType label
    let stypes = SuffixTypes::new(text);
    // Step 2 - Divide the buckets into S- and L-Buckets on. TheL-Buckets in front of the S-Buckets.
    let bins = Bins::new(text, &stypes);
    induced_sort(text, &stypes, &bins)
}

/// Performs the remaining SA-IS steps for a text whose suffix types and bins are already known.
///
/// # Errors
///
/// Same as `sais`.
fn induced_sort<S: Ord + Clone + Hash + Debug>(
    text: &[S],
    stypes: &SuffixTypes,
    bins: &Bins,
) -> Result<Vec<usize>> {
    let mut sa: Vec<isize> = vec![-1; text.len()];

    // Sort the LeftmostSType strings lexicographically and write them in their respective S-Buckets.
    let w_slice_indices = find_lms_substring_indices(stypes);
    let (duplicates, reduced_substring) =
        reduce_lms_substrings(text, &w_slice_indices, stypes, bins)?;

    let mut reduced_suffix_array: Vec<usize>;
    if duplicates {
        // The reduced string already consists of ranks, so its bins need no hashing.
        let reduced_stypes = SuffixTypes::new(&reduced_substring);
        let alphabet_size = reduced_substring.iter().max().map_or(0, |name| name + 1);
        let reduced_bins =
            Bins::from_ranks(reduced_substring.clone(), alphabet_size, &reduced_stypes);
        reduced_suffix_array = induced_sort(&reduced_substring, &reduced_stypes, &reduced_bins)?;
    } else {
        reduced_suffix_array = vec![0; reduced_substring.len()];
        for (original_index, &sort_index) in reduced_substring.iter().enumerate() {
            reduced_suffix_array[sort_index] = original_index;
        }
    }
    let mut bin_edges = bins.tails();
    for reduced_string_index in reduced_suffix_array.iter().rev() {
        let w_slice_start = w_slice_indices[*reduced_string_index].0;
        let x = &mut bin_edges[bins.position_bins[w_slice_start]];

        sa[*x] = isize::try_from(w_slice_start).map_err(|e| Error::TryFromIntError { e })?;
        if *x > 0 {
            *x -= 1;
        }
    }

    induce_pass(&mut sa, stypes, bins, true)?;
    induce_pass(&mut sa, stypes, bins, false)?;

    sa.iter()
        .map(|i| usize::try_from(*i).map_err(|e| Error::TryFromIntError { e }))
        .collect::<Result<Vec<usize>>>()
}

#[cfg(test)]
mod tests {
    use super::*;
    use bitvec::bitvec;
    use bitvec::prelude::*;

    #[test]
    fn test_suffix_types() {
        let suffix_types = SuffixTypes::new("mmiissiissiippii$".as_bytes());

        let expected_largest_types = bitvec![1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0];
        let expected_smallest_types = bitvec![0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1];
        let expected_lms_type = bitvec![0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1];
        assert_eq!(
            suffix_types.types.len(),
            expected_largest_types.len(),
            "Suffix Types length ({}) didn't match expectation ({})",
            suffix_types.types.len(),
            expected_largest_types.len()
        );
        for idx in 0..suffix_types.types.len() {
            assert_eq!(suffix_types.is_ltype(idx), expected_largest_types[idx], "Suffix Type at position {} is L type was expected to resolve to {}, resolved to {}", idx, expected_largest_types[idx], suffix_types.is_ltype(idx));
            assert_eq!(suffix_types.is_stype(idx), expected_smallest_types[idx], "Suffix Type at position {} is S type was expected to resolve to {}, resolved to {}", idx, expected_smallest_types[idx], suffix_types.is_ltype(idx));
            assert_eq!(suffix_types.is_leftmost_stype(idx), expected_lms_type[idx], "Suffix Type at position {} is LMS type was expected to resolve to {}, resolved to {}", idx, expected_lms_type[idx], suffix_types.is_leftmost_stype(idx));
        }
    }

    #[test]
    fn test_w_slice_indices() {
        let text = "mmiissiissiippii$".as_bytes();
        let suffix_types = SuffixTypes::new(text);
        let w_slice_indices = find_lms_substring_indices(&suffix_types);

        assert_eq!(w_slice_indices, vec![(2, 7), (6, 11), (10, 17), (16, 17)]);
        let w_slices: Vec<&[u8]> =
            w_slice_indices.iter().map(|&(start, end)| &text[start..end]).collect();

        assert_eq!(
            w_slices,
            vec!["iissi".as_bytes(), "iissi".as_bytes(), "iippii$".as_bytes(), "$".as_bytes()]
        );
    }

    #[test]
    fn test_count_text_bin_sizes() {
        let m = "m".as_bytes()[0];
        let i = "i".as_bytes()[0];
        let s = "s".as_bytes()[0];
        let p = "p".as_bytes()[0];
        let terminator = "$".as_bytes()[0];

        let text = "mmiissiissiippii$".as_bytes();
        let suffix_types = SuffixTypes::new(text);
        let text_bin_sizes = count_text_bin_sizes(text, &suffix_types);

        assert_eq!(
            text_bin_sizes,
            HashMap::from([
                (SortableSuffixType::LType(m), 2),
                (SortableSuffixType::LType(i), 2),
                (SortableSuffixType::SType(i), 6),
                (SortableSuffixType::LType(s), 4),
                (SortableSuffixType::LType(p), 2),
                (SortableSuffixType::SType(terminator), 1),
            ])
        );
    }

    #[test]
    fn test_get_bin_offsets() {
        let m = "m".as_bytes()[0];
        let i = "i".as_bytes()[0];
        let s = "s".as_bytes()[0];
        let p = "p".as_bytes()[0];
        let terminator = "$".as_bytes()[0];

        let text = "mmiissiissiippii$".as_bytes();
        let suffix_types = SuffixTypes::new(text);
        let text_bin_sizes = count_text_bin_sizes(text, &suffix_types);
        let forward_bin_offsets = get_bin_offsets(&text_bin_sizes, true);
        let reverse_bin_offsets = get_bin_offsets(&text_bin_sizes, false);

        assert_eq!(
            forward_bin_offsets,
            HashMap::from([
                (SortableSuffixType::SType(terminator), 0),
                (SortableSuffixType::LType(i), 1),
                (SortableSuffixType::SType(i), 3),
                (SortableSuffixType::LType(m), 9),
                (SortableSuffixType::LType(p), 11),
                (SortableSuffixType::LType(s), 13),
            ])
        );
        assert_eq!(
            reverse_bin_offsets,
            HashMap::from([
                (SortableSuffixType::SType(terminator), 0),
                (SortableSuffixType::LType(i), 2),
                (SortableSuffixType::SType(i), 8),
                (SortableSuffixType::LType(m), 10),
                (SortableSuffixType::LType(p), 12),
                (SortableSuffixType::LType(s), 16),
            ])
        );
    }

    #[test]
    fn test_get_reduced_substring() {
        let text = "mmiissiissiippii$".as_bytes();
        let suffix_types = SuffixTypes::new(text);
        let text_bin_sizes = count_text_bin_sizes(text, &suffix_types);
        let w_slice_indices = find_lms_substring_indices(&suffix_types);
        let (duplicates, reduced_substring) =
            get_reduced_substring(text, &w_slice_indices, &suffix_types, &text_bin_sizes).unwrap();
        assert!(duplicates, "Reduced substring of 'mmiissiissiippii$' has duplicates not identified by our approach");
        assert_eq!(reduced_substring, vec![2, 2, 1, 0]);
    }

    #[test]
    fn test_suffix_array_modification_propagation() {
        let text = "mmiissiissiippii$".as_bytes();
        let suffix_types = SuffixTypes::new(text);
        let text_bin_sizes = count_text_bin_sizes(text, &suffix_types);
        let mut initial_suffix_array =
            vec![16, -1, -1, -1, -1, -1, 10, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1];
        suffix_array_modification_propagation(
            &mut initial_suffix_array,
            text,
            &suffix_types,
            &text_bin_sizes,
            true,
        )
        .unwrap();
        assert_eq!(
            initial_suffix_array,
            vec![16, 15, 14, -1, -1, -1, 10, 6, 2, 1, 0, 13, 12, 9, 5, 8, 4]
        );
        suffix_array_modification_propagation(
            &mut initial_suffix_array,
            text,
            &suffix_types,
            &text_bin_sizes,
            false,
        )
        .unwrap();
        assert_eq!(
            initial_suffix_array,
            vec![16, 15, 14, 10, 6, 2, 11, 7, 3, 1, 0, 13, 12, 9, 5, 8, 4]
        );
    }

    #[test]
    fn full_test() {
        let paper_example = "mmiissiissiippii$".as_bytes();
        assert_eq!(
            sais(paper_example).unwrap(),
            [16, 15, 14, 10, 6, 2, 11, 7, 3, 1, 0, 13, 12, 9, 5, 8, 4].to_vec()
        );
    }

    #[test]
    fn wholly_repetitive_completes() {
        let homopolymer = [1_u32, 1, 1, 1, 1, 1, 1, 0];
        assert_eq!(sais(&homopolymer).unwrap(), [7_usize, 6, 5, 4, 3, 2, 1, 0].to_vec());
    }
}
