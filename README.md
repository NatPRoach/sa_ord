# SA-Ord
Pure Rust implementation of the [SA-IS suffix array algorithm](https://ieeexplore.ieee.org/document/4976463)
for slices of types that impl `Ord` + `Hash` + `Clone` + `Debug`

The last character of the slice must be the smallest character in the text in order for the algorithm
to work properly. This can be achieved by adding a terminator character.

Usage:
```rust
use sa_ord::sais

let paper_example = "mmiissiissiippii$".as_bytes();
assert_eq!(
    sais(paper_example).unwrap(),
    [16, 15, 14, 10, 6, 2, 11, 7, 3, 1, 0, 13, 12, 9, 5, 8, 4].to_vec()
);
```