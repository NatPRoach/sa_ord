//! Benchmarks for `sa_ord::sais` over a DNA-like `u8` alphabet, a small `Copy` enum alphabet and
//! non-`Copy` `String` tokens, plus a comparison against the `suffix` crate.
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use rand::{rngs::StdRng, Rng, SeedableRng};
use std::hint::black_box;

const SIZES: [usize; 3] = [10_000, 100_000, 1_000_000];
const SEED: u64 = 0x5EED_1234_ABCD_9876;

/// Small fixed-size `Copy` alphabet. `Sentinel` sorts first and is used only as the terminator.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum Nucleotide {
    Sentinel,
    A,
    C,
    G,
    T,
}

/// Random ACGT text of length `len` (including a trailing `$`, which is smaller than any base).
fn dna_text(len: usize) -> Vec<u8> {
    let mut rng = StdRng::seed_from_u64(SEED);
    let mut text: Vec<u8> = (0..len - 1).map(|_| b"ACGT"[rng.random_range(0..4)]).collect();
    text.push(b'$');
    text
}

/// Enum equivalent of the `dna_text` input. Producing the same bases in the same order in a different datatype.
fn enum_text(len: usize) -> Vec<Nucleotide> {
    const BASES: [Nucleotide; 4] = [Nucleotide::A, Nucleotide::C, Nucleotide::G, Nucleotide::T];
    let mut rng = StdRng::seed_from_u64(SEED);
    let mut text: Vec<Nucleotide> = (0..len - 1).map(|_| BASES[rng.random_range(0..4)]).collect();
    text.push(Nucleotide::Sentinel);
    text
}

/// Random tokens from a 64-word vocabulary, terminated by `""` (smaller than every word).
fn string_text(len: usize) -> Vec<String> {
    let vocab: Vec<String> = (0..64).map(|i| format!("word{:02}", i)).collect();
    let mut rng = StdRng::seed_from_u64(SEED);
    let mut text: Vec<String> =
        (0..len - 1).map(|_| vocab[rng.random_range(0..vocab.len())].clone()).collect();
    text.push(String::new());
    text
}

/// Prefix of the Fibonacci word over {a, b} of length `len` (including a trailing `$`).
fn fibonacci_text(len: usize) -> Vec<u8> {
    let (mut prev, mut cur) = (b"a".to_vec(), b"ab".to_vec());
    while cur.len() < len - 1 {
        let next = [cur.as_slice(), prev.as_slice()].concat();
        prev = std::mem::replace(&mut cur, next);
    }
    cur.truncate(len - 1);
    cur.push(b'$');
    cur
}

fn configure(group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>, n: usize) {
    group.throughput(Throughput::Elements(n as u64));
    group.sample_size(if n >= 1_000_000 { 10 } else { 20 });
}

fn bench_sais(c: &mut Criterion) {
    let mut group = c.benchmark_group("dna_u8");
    for &n in &SIZES {
        configure(&mut group, n);
        let text = dna_text(n);
        group.bench_with_input(BenchmarkId::new("sa_ord", n), &text, |b, t| {
            b.iter(|| sa_ord::sais(black_box(t)).unwrap())
        });
        // `suffix` terminates the text itself, so give it the bases without our `$`.
        let text = String::from_utf8(text[..n - 1].to_vec()).unwrap();
        group.bench_with_input(BenchmarkId::new("suffix_crate", n), &text, |b, t| {
            b.iter(|| suffix::SuffixTable::new(black_box(t.as_str())))
        });
    }
    group.finish();

    let mut group = c.benchmark_group("enum_copy");
    for &n in &SIZES {
        configure(&mut group, n);
        let text = enum_text(n);
        group.bench_with_input(BenchmarkId::new("sa_ord", n), &text, |b, t| {
            b.iter(|| sa_ord::sais(black_box(t)).unwrap())
        });
    }
    group.finish();

    let mut group = c.benchmark_group("string_tokens");
    for &n in &SIZES {
        configure(&mut group, n);
        let text = string_text(n);
        group.bench_with_input(BenchmarkId::new("sa_ord", n), &text, |b, t| {
            b.iter(|| sa_ord::sais(black_box(t)).unwrap())
        });
    }
    group.finish();

    // Highly repetitive text recurses deeply, which stresses the `usize` reduced levels.
    let mut group = c.benchmark_group("fibonacci_u8");
    let n = 100_000;
    configure(&mut group, n);
    let text = fibonacci_text(n);
    group.bench_with_input(BenchmarkId::new("sa_ord", n), &text, |b, t| {
        b.iter(|| sa_ord::sais(black_box(t)).unwrap())
    });
    group.finish();
}

criterion_group!(benches, bench_sais);
criterion_main!(benches);
