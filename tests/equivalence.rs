//! Equivalence tests for `sa_ord::sais`.
//!
//! The current implementation is compared against
//!  * `reference::reference_sais` - a frozen snapshot of the pre-optimization implementation,
//!  * a naive oracle (sort of `0..n` by suffix) for small inputs,
//!  * the `suffix` crate (DNA strings, after stripping the terminator suffix).
//!
//! Every text handed to `sais` ends with a unique, strictly smallest terminator. Texts are
//! generated as `u32` "bodies" with values in `1..=k` and then mapped into each element type
//! under test (the terminator is always the minimum of the type).

mod reference;

use reference::{reference_max_depth, reference_sais};
use sa_ord::sais;
use std::fmt::Debug;
use std::hash::Hash;
use std::panic::{catch_unwind, AssertUnwindSafe};

// ---------------------------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------------------------

/// Deterministic xorshift64* RNG.
struct Rng(u64);

impl Rng {
    fn new(seed: u64) -> Self {
        // Mix the seed so that small seeds don't give correlated streams; never zero.
        let mut r = Rng(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) ^ 0xD1B5_4A32_D192_ED03);
        if r.0 == 0 {
            r.0 = 1;
        }
        for _ in 0..4 {
            r.next();
        }
        r
    }

    fn next(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    /// Uniform-ish value in `0..n` (n > 0).
    fn below(&mut self, n: usize) -> usize {
        ((self.next() >> 11) % n as u64) as usize
    }
}

/// Naive suffix array: stable sort of `0..n` by suffix.
fn naive<S: Ord>(text: &[S]) -> Vec<usize> {
    let mut idx: Vec<usize> = (0..text.len()).collect();
    idx.sort_by(|&a, &b| text[a..].cmp(&text[b..]));
    idx
}

const NAIVE_LIMIT: usize = 2000;

/// Core assertion: `sais == reference` and (for small inputs) `sais == naive`.
fn check<S: Ord + Clone + Hash + Debug>(text: &[S], label: &str) {
    if text.len() == 1 {
        // Known quirk of the frozen implementation: a terminator-only text yields Err (the lone
        // position is S-type but not *leftmost* S-type, so it is missing from the LMS mapping).
        // The new implementation must keep returning Err (see `terminator_only`).
        assert!(sais(text).is_err(), "[{label}] n=1: expected Err like the reference");
        assert!(reference_sais(text).is_err(), "[{label}] n=1: reference drifted");
        return;
    }
    let got = sais(text).unwrap_or_else(|e| panic!("[{label}] sais returned Err: {e}"));
    let want =
        reference_sais(text).unwrap_or_else(|e| panic!("[{label}] reference returned Err: {e}"));
    assert_eq!(got, want, "[{label}] sais != reference (len {})", text.len());
    if text.len() <= NAIVE_LIMIT {
        assert_eq!(got, naive(text), "[{label}] sais != naive (len {})", text.len());
    }
}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
enum Tok {
    Sentinel,
    A,
    B,
    C,
    D,
    E,
}

fn tok(v: u32) -> Tok {
    match v {
        0 => Tok::Sentinel,
        1 => Tok::A,
        2 => Tok::B,
        3 => Tok::C,
        4 => Tok::D,
        5 => Tok::E,
        _ => unreachable!("Tok only supports values 0..=5"),
    }
}

/// Runs `check` for `body ++ [terminator]` in every element type that can represent it.
/// `body` values must lie in `1..=65535`.
fn check_all_types(body: &[u32], label: &str) {
    let max = body.iter().copied().max().unwrap_or(0);
    assert!(body.iter().all(|&v| (1..=65535).contains(&v)), "bad body value");
    let term = |mut v: Vec<u32>| {
        v.push(0);
        v
    };
    let t = term(body.to_vec());

    // u32
    check(&t, &format!("{label}/u32"));

    // u8
    if max <= 255 {
        let text: Vec<u8> = t.iter().map(|&v| v as u8).collect();
        check(&text, &format!("{label}/u8"));
    }

    // i32 with negatives, terminator i32::MIN
    let text: Vec<i32> = t.iter().map(|&v| if v == 0 { i32::MIN } else { v as i32 - 3 }).collect();
    check(&text, &format!("{label}/i32"));

    // char ('\0' terminator, body mapped monotonically into the supplementary planes)
    let text: Vec<char> = t
        .iter()
        .map(|&v| if v == 0 { '\0' } else { char::from_u32(0x10000 + v).expect("valid char") })
        .collect();
    check(&text, &format!("{label}/char"));

    // (u8, u8)
    let text: Vec<(u8, u8)> = t.iter().map(|&v| ((v >> 8) as u8, v as u8)).collect();
    check(&text, &format!("{label}/tuple"));

    // enum with sentinel first variant
    if max <= 5 {
        let text: Vec<Tok> = t.iter().map(|&v| tok(v)).collect();
        check(&text, &format!("{label}/enum"));
    }

    // String tokens, "" terminator (zero padded so lexicographic == numeric)
    let text: Vec<String> =
        t.iter().map(|&v| if v == 0 { String::new() } else { format!("w{v:05}") }).collect();
    check(&text, &format!("{label}/string"));
}

/// Random body of length `len` over alphabet `1..=k`.
fn random_body(rng: &mut Rng, len: usize, k: usize) -> Vec<u32> {
    (0..len).map(|_| rng.below(k) as u32 + 1).collect()
}

/// Random body made of runs of repeated symbols.
fn runs_body(rng: &mut Rng, len: usize, k: usize, max_run: usize) -> Vec<u32> {
    let mut v = Vec::with_capacity(len);
    while v.len() < len {
        let c = rng.below(k) as u32 + 1;
        let r = rng.below(max_run) + 1;
        for _ in 0..r.min(len - v.len()) {
            v.push(c);
        }
    }
    v
}

fn fib_body(len: usize) -> Vec<u32> {
    let (mut a, mut b): (Vec<u32>, Vec<u32>) = (vec![1], vec![1, 2]);
    while b.len() < len {
        let mut n = b.clone();
        n.extend_from_slice(&a);
        a = b;
        b = n;
    }
    b.truncate(len);
    b
}

fn periodic_body(pattern: &[u32], len: usize) -> Vec<u32> {
    (0..len).map(|i| pattern[i % pattern.len()]).collect()
}

const ALPHABETS: [usize; 5] = [1, 2, 4, 26, 256];

// ---------------------------------------------------------------------------------------------
// Tiny sizes
// ---------------------------------------------------------------------------------------------

#[test]
fn terminator_only() {
    check_all_types(&[], "n=1");
    // Behaviour record: the frozen implementation returns Err for a text of length 1.
    assert!(sais(&[0u8]).is_err());
    assert!(reference_sais(&[0u8]).is_err());
}

#[test]
fn sizes_two_and_three_exhaustive() {
    // All bodies of length 1 and 2 over alphabet 1..=3 (texts of length 2 and 3).
    for a in 1..=3u32 {
        check_all_types(&[a], "n=2");
        for b in 1..=3u32 {
            check_all_types(&[a, b], "n=3");
        }
    }
}

#[test]
fn small_exhaustive_binary_and_ternary() {
    // Every body up to length 9 over {1,2} and up to length 6 over {1,2,3}.
    for len in 0..=9usize {
        for code in 0..(1u32 << len) {
            let body: Vec<u32> = (0..len).map(|i| ((code >> i) & 1) + 1).collect();
            // u32 + u8 only, keeps this fast
            let mut t = body.clone();
            t.push(0);
            check(&t, "exhaustive2/u32");
            let t8: Vec<u8> = t.iter().map(|&v| v as u8).collect();
            check(&t8, "exhaustive2/u8");
        }
    }
    for len in 0..=6usize {
        for code in 0..3u32.pow(len as u32) {
            let mut c = code;
            let mut t: Vec<u32> = (0..len)
                .map(|_| {
                    let d = c % 3 + 1;
                    c /= 3;
                    d
                })
                .collect();
            t.push(0);
            check(&t, "exhaustive3/u32");
        }
    }
}

#[test]
fn classic_examples() {
    check(b"mmiissiissiippii$", "mississippi-ish");
    check(b"banana$", "banana");
    check(b"abracadabra$", "abracadabra");
    check(b"$", "just terminator");
    check(b"a$", "one char");
}

// ---------------------------------------------------------------------------------------------
// Random texts, all alphabets, all element types
// ---------------------------------------------------------------------------------------------

#[test]
fn random_small_all_alphabets_all_types() {
    for &k in &ALPHABETS {
        for seed in 0..60u64 {
            let mut rng = Rng::new(seed * 1000 + k as u64);
            let len = rng.below(120);
            let body = random_body(&mut rng, len, k);
            check_all_types(&body, &format!("rand k={k} seed={seed} len={len}"));
        }
    }
}

#[test]
fn random_medium_up_to_5000() {
    // u32 and u8 only (the heavy-weight types are covered at small sizes). Sizes straddle the
    // naive-oracle limit of 2000.
    let sizes = [200usize, 777, 1999, 2000, 2001, 3500, 5000];
    for &k in &ALPHABETS {
        for (i, &len) in sizes.iter().enumerate() {
            for seed in 0..3u64 {
                let mut rng = Rng::new(7_000 + seed * 100 + i as u64 * 10 + k as u64);
                let body = random_body(&mut rng, len - 1, k);
                let mut t = body.clone();
                t.push(0);
                check(&t, &format!("rand-med/u32 k={k} len={len} seed={seed}"));
                if k <= 255 {
                    let t8: Vec<u8> = t.iter().map(|&v| v as u8).collect();
                    check(&t8, &format!("rand-med/u8 k={k} len={len} seed={seed}"));
                }
            }
        }
    }
}

#[test]
fn random_with_long_runs() {
    for &k in &[2usize, 4, 26] {
        for seed in 0..25u64 {
            let mut rng = Rng::new(31_000 + seed * 7 + k as u64);
            let len = rng.below(400) + 1;
            let body = runs_body(&mut rng, len, k, 40);
            check_all_types(&body, &format!("runs k={k} seed={seed} len={len}"));
        }
    }
    // A few bigger ones with very long runs.
    for seed in 0..5u64 {
        let mut rng = Rng::new(99_000 + seed);
        let body = runs_body(&mut rng, 1800, 3, 300);
        check_all_types(&body, &format!("long-runs seed={seed}"));
    }
}

#[test]
fn random_wide_alphabet_sparse_values() {
    // Values scattered over 1..=65535 (exercises u32/i32/char/tuple/string with big gaps).
    for seed in 0..30u64 {
        let mut rng = Rng::new(555 + seed);
        let len = rng.below(200);
        let body: Vec<u32> = (0..len).map(|_| rng.below(65535) as u32 + 1).collect();
        check_all_types(&body, &format!("sparse seed={seed}"));
    }
}

// ---------------------------------------------------------------------------------------------
// Degenerate / repetitive
// ---------------------------------------------------------------------------------------------

#[test]
fn homopolymer() {
    for len in [0usize, 1, 2, 3, 10, 100, 1000, 1999] {
        check_all_types(&vec![1; len], &format!("homopolymer len={len}"));
    }
}

#[test]
fn periodic_ab_abc() {
    for k in [1usize, 2, 3, 5, 10, 50, 333] {
        check_all_types(&periodic_body(&[1, 2], 2 * k), &format!("(ab)^{k}"));
        check_all_types(&periodic_body(&[1, 2], 2 * k + 1), &format!("(ab)^{k}a"));
        check_all_types(&periodic_body(&[1, 2, 3], 3 * k), &format!("(abc)^{k}"));
        check_all_types(&periodic_body(&[1, 2, 3], 3 * k + 2), &format!("(abc)^{k}ab"));
        check_all_types(&periodic_body(&[2, 1], 2 * k), &format!("(ba)^{k}"));
        check_all_types(&periodic_body(&[3, 2, 1], 3 * k), &format!("(cba)^{k}"));
    }
}

#[test]
fn fibonacci_words() {
    for len in [0usize, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987, 1597] {
        check_all_types(&fib_body(len), &format!("fib len={len}"));
    }
    // Also fib-1 / fib+1 to break exact-boundary luck.
    for len in [6usize, 7, 12, 14, 232, 234, 986, 988] {
        check_all_types(&fib_body(len), &format!("fib len={len}"));
    }
}

#[test]
fn strictly_descending_and_ascending() {
    for m in [1usize, 2, 3, 5, 26, 100, 255, 256, 1000, 1999] {
        let desc: Vec<u32> = (1..=m as u32).rev().collect();
        let asc: Vec<u32> = (1..=m as u32).collect();
        check_all_types(&desc, &format!("descending m={m}"));
        check_all_types(&asc, &format!("ascending m={m}"));
    }
}

#[test]
fn all_distinct_permutations() {
    // Every suffix of a permutation starts with a distinct symbol, so the reduced string has no
    // duplicates (the non-recursive branch).
    for m in [1usize, 2, 3, 4, 5, 7, 26, 100, 255, 1000, 1500] {
        for seed in 0..4u64 {
            let mut rng = Rng::new(4242 + seed * 17 + m as u64);
            let mut p: Vec<u32> = (1..=m as u32).collect();
            for i in (1..p.len()).rev() {
                let j = rng.below(i + 1);
                p.swap(i, j);
            }
            check_all_types(&p, &format!("perm m={m} seed={seed}"));
        }
    }
}

// ---------------------------------------------------------------------------------------------
// Large inputs: compared against the reference only (naive is skipped above NAIVE_LIMIT)
// ---------------------------------------------------------------------------------------------

fn check_large(body: &[u32], label: &str) {
    let mut t = body.to_vec();
    t.push(0);
    check(&t, &format!("{label}/u32"));
    if body.iter().all(|&v| v <= 255) {
        let t8: Vec<u8> = t.iter().map(|&v| v as u8).collect();
        check(&t8, &format!("{label}/u8"));
    }
}

/// Large (10^5) inputs. If this takes more than ~2s in a debug build it is `#[ignore]`d; run with
/// `cargo +stable test -- --ignored`.
#[test]
#[ignore = "large inputs (1e5), slow in debug builds"]
fn large_random() {
    for &k in &[2usize, 4, 26, 256] {
        let mut rng = Rng::new(123_456 + k as u64);
        let body = random_body(&mut rng, 100_000 - 1, k);
        check_large(&body, &format!("large rand k={k}"));
    }
}

#[test]
#[ignore = "large inputs (1e5), slow in debug builds"]
fn large_degenerate() {
    check_large(&vec![1; 99_999], "large homopolymer");
    check_large(&periodic_body(&[1, 2], 99_999), "large (ab)^k");
    check_large(&periodic_body(&[1, 2, 3], 99_999), "large (abc)^k");
    check_large(&fib_body(99_999), "large fib");
    let mut rng = Rng::new(8);
    check_large(&runs_body(&mut rng, 99_999, 4, 500), "large runs");
    let desc: Vec<u32> = (1..=50_000u32).rev().collect();
    check_large(&desc, "large descending");
    let asc: Vec<u32> = (1..=50_000u32).collect();
    check_large(&asc, "large ascending");
}

// ---------------------------------------------------------------------------------------------
// Cross-check against the `suffix` crate
// ---------------------------------------------------------------------------------------------

/// For `s` without terminator, `sais(s + "$")[1..] == SuffixTable::new(s).table()`: the
/// terminator suffix (position n) is always first, and since '$' < 'A' < 'C' < 'G' < 'T' the
/// remaining order is that of the plain suffix array of `s`.
fn suffix_crate_check(s: &str, label: &str) {
    let mut with_term = s.as_bytes().to_vec();
    with_term.push(b'$');
    let ours = sais(&with_term).expect("sais ok");
    assert_eq!(ours[0], s.len(), "[{label}] terminator suffix must come first");
    let table = suffix::SuffixTable::new(s);
    let theirs: Vec<usize> = table.table().iter().map(|&x| x as usize).collect();
    assert_eq!(&ours[1..], &theirs[..], "[{label}] differs from suffix crate");
}

#[test]
fn suffix_crate_random_dna() {
    for seed in 0..80u64 {
        let mut rng = Rng::new(2024 + seed);
        let len = rng.below(1500);
        let s: String = (0..len).map(|_| b"ACGT"[rng.below(4)] as char).collect();
        suffix_crate_check(&s, &format!("dna seed={seed} len={len}"));
    }
}

#[test]
fn suffix_crate_dna_edge_cases() {
    // Empty `s` gives the length-1 text "$", for which sais returns Err (see `terminator_only`).
    assert!(sais(b"$").is_err());
    assert!(suffix::SuffixTable::new("").table().is_empty());
    suffix_crate_check("A", "single");
    suffix_crate_check(&"A".repeat(500), "homopolymer");
    suffix_crate_check(&"AC".repeat(300), "(AC)^k");
    suffix_crate_check(&"ACG".repeat(200), "(ACG)^k");
    suffix_crate_check("GATTACA", "gattaca");
    let fib: String = fib_body(1000).iter().map(|&v| if v == 1 { 'A' } else { 'C' }).collect();
    suffix_crate_check(&fib, "fib");
    let mut rng = Rng::new(1);
    let runs: String =
        runs_body(&mut rng, 3000, 4, 60).iter().map(|&v| b"ACGT"[v as usize - 1] as char).collect();
    suffix_crate_check(&runs, "runs");
}

// ---------------------------------------------------------------------------------------------
// Invalid input: behaviour record
// ---------------------------------------------------------------------------------------------

#[derive(Debug, PartialEq, Eq)]
enum Outcome {
    Ok(Vec<usize>),
    Err,
    Panic,
}

fn classify<F: FnOnce() -> sa_ord::errors::Result<Vec<usize>>>(f: F) -> Outcome {
    match catch_unwind(AssertUnwindSafe(f)) {
        Ok(Ok(v)) => Outcome::Ok(v),
        Ok(Err(_)) => Outcome::Err,
        Err(_) => Outcome::Panic,
    }
}

/// CURRENT (frozen) behaviour on invalid input, as observed when this test was written, on a
/// debug build: `b"banana"` (last char not smallest) and `b"a$b$"` (last char smallest but not
/// unique) both return `Ok` with the (accidentally correct) suffix array listed below.
///
/// The new implementation must fall in the same class (Ok with identical output / Err / panic).
fn expected_invalid(name: &str) -> Outcome {
    match name {
        // Observed behaviours (filled in from a run against the frozen reference).
        "banana" => Outcome::Ok(vec![5, 3, 1, 0, 4, 2]),
        "a$b$" => Outcome::Ok(vec![3, 1, 0, 2]),
        _ => unreachable!(),
    }
}

#[test]
fn invalid_input_behaviour_matches_reference() {
    let cases: Vec<(&str, &[u8])> = vec![
        ("banana", b"banana"),
        ("a$b$", b"a$b$"),
        ("aa", b"aa"),
        ("ba", b"ba"),
        ("ab", b"ab"),
        ("abab", b"abab"),
        ("$$", b"$$"),
        ("a$$", b"a$$"),
        ("zzzz", b"zzzz"),
        ("b$a$", b"b$a$"),
        ("aab$b$", b"aab$b$"),
        ("abcabc", b"abcabc"),
        ("cba", b"cba"),
    ];
    for (name, text) in cases {
        let want = classify(|| reference_sais(text));
        let got = classify(|| sais(text));
        eprintln!("invalid input {name:?}: reference = {want:?}, new = {got:?}");
        assert_eq!(got, want, "invalid input {name:?}: behaviour class changed");
    }
}

#[test]
fn invalid_input_recorded_behaviour() {
    for (name, text) in [("banana", &b"banana"[..]), ("a$b$", &b"a$b$"[..])] {
        let want = expected_invalid(name);
        assert_eq!(classify(|| reference_sais(text)), want, "reference drifted for {name:?}");
        assert_eq!(classify(|| sais(text)), want, "sais differs for {name:?}");
    }
}

#[test]
fn empty_input_behaviour_matches_reference() {
    // The reference computes `text.len() - 1`, which panics on overflow in debug builds. In
    // release builds it would wrap, so only assert the class where the behaviour is defined.
    if cfg!(debug_assertions) {
        let empty: &[u8] = &[];
        let want = classify(|| reference_sais(empty));
        let got = classify(|| sais(empty));
        eprintln!("empty input: reference = {want:?}, new = {got:?}");
        assert_eq!(got, want);
    }
}

// ---------------------------------------------------------------------------------------------
// Recursion depth
// ---------------------------------------------------------------------------------------------

/// Confirms that the suite exercises the recursive branch, including >= 2 levels of recursion
/// (depth 3 in the reference's 1-based counting: top level + two recursive calls).
#[test]
fn recursion_depth_is_exercised() {
    let mut t = periodic_body(&[1, 2], 2000);
    t.push(0);
    let d_ab = reference_max_depth(&t);

    let mut t = fib_body(1597);
    t.push(0);
    let d_fib = reference_max_depth(&t);

    let mut t = vec![1u32; 100];
    t.push(0);
    let d_homo = reference_max_depth(&t);

    let mut t = periodic_body(&[1, 2, 3], 3000);
    t.push(0);
    let d_abc = reference_max_depth(&t);

    eprintln!(
        "depths: (ab)^1000={d_ab} fib1597={d_fib} homopolymer100={d_homo} (abc)^1000={d_abc}"
    );
    // Depth is 1-based: 1 = no recursion, 3 = recursion twice (>= 2 recursive levels).
    assert!(d_fib >= 3, "fib word should recurse at least 2 levels, got depth {d_fib}");
    assert!(d_ab >= 2, "(ab)^k should recurse, got depth {d_ab}");
    assert!(d_abc >= 2, "(abc)^k should recurse, got depth {d_abc}");
    // And the recursive inputs agree with all oracles.
    check_all_types(&fib_body(1597), "fib-depth");
}
