use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use ie_coloring::{Graph, baselines, ie};
use std::hint::black_box;
use std::time::Duration;

fn configured() -> Criterion {
    Criterion::default()
        .sample_size(10)
        .warm_up_time(Duration::from_millis(500))
        .measurement_time(Duration::from_secs(3))
}

/// Exact chromatic number on G(n, 1/2) as n grows.
fn exact_vs_n(c: &mut Criterion) {
    let mut group = c.benchmark_group("chromatic_number/gnp_0.5");
    for n in [10, 12, 14, 16, 18, 20, 22] {
        let g = Graph::random(n, 0.5, n as u64);
        group.bench_with_input(BenchmarkId::new("ie", n), &g, |b, g| b.iter(|| ie::chromatic_number(black_box(g))));
        group.bench_with_input(BenchmarkId::new("ie_unbracketed", n), &g, |b, g| {
            b.iter(|| ie::chromatic_number_unbracketed(black_box(g)))
        });
        if n <= 20 {
            group.bench_with_input(BenchmarkId::new("ie_polyspace", n), &g, |b, g| {
                b.iter(|| ie::chromatic_number_polyspace(black_box(g)))
            });
            group
                .bench_with_input(BenchmarkId::new("lawler", n), &g, |b, g| b.iter(|| baselines::lawler(black_box(g))));
        }
        group.bench_with_input(BenchmarkId::new("dsatur_bnb", n), &g, |b, g| {
            b.iter(|| baselines::dsatur_exact(black_box(g)))
        });
    }
    group.finish();
}

/// Exact chromatic number on 18-vertex G(n, p) as the density varies.
fn exact_vs_density(c: &mut Criterion) {
    let mut group = c.benchmark_group("chromatic_number/n18_density");
    for p in [0.1, 0.3, 0.5, 0.7, 0.9] {
        let g = Graph::random(18, p, 99);
        group.bench_with_input(BenchmarkId::new("ie", p), &g, |b, g| b.iter(|| ie::chromatic_number(black_box(g))));
        group.bench_with_input(BenchmarkId::new("ie_polyspace", p), &g, |b, g| {
            b.iter(|| ie::chromatic_number_polyspace(black_box(g)))
        });
        group.bench_with_input(BenchmarkId::new("lawler", p), &g, |b, g| b.iter(|| baselines::lawler(black_box(g))));
        group.bench_with_input(BenchmarkId::new("dsatur_bnb", p), &g, |b, g| {
            b.iter(|| baselines::dsatur_exact(black_box(g)))
        });
    }
    group.finish();
}

/// Structured graphs that are hard for clique bounds / greedy reasoning.
fn exact_structured(c: &mut Criterion) {
    let mut group = c.benchmark_group("chromatic_number/structured");
    let graphs = [
        ("mycielski5_n23", Graph::mycielski(5)),
        ("queen5x5_n25", Graph::queen(5, 5)),
        ("queen4x5_n20", Graph::queen(4, 5)),
    ];
    for (name, g) in &graphs {
        group.bench_function(BenchmarkId::new("ie", name), |b| b.iter(|| ie::chromatic_number(black_box(g))));
        group
            .bench_function(BenchmarkId::new("dsatur_bnb", name), |b| b.iter(|| baselines::dsatur_exact(black_box(g))));
        if g.n() <= 20 {
            group.bench_function(BenchmarkId::new("lawler", name), |b| b.iter(|| baselines::lawler(black_box(g))));
        }
    }
    group.finish();
}

/// Heuristic colourings for reference (they are not exact).
fn heuristics(c: &mut Criterion) {
    let mut group = c.benchmark_group("heuristics/gnp_0.5");
    for n in [20, 40, 64] {
        let g = Graph::random(n, 0.5, n as u64);
        let order: Vec<usize> = (0..n).collect();
        group.bench_with_input(BenchmarkId::new("greedy", n), &g, |b, g| {
            b.iter(|| baselines::greedy(black_box(g), &order))
        });
        group.bench_with_input(BenchmarkId::new("largest_first", n), &g, |b, g| {
            b.iter(|| baselines::largest_first(black_box(g)))
        });
        group.bench_with_input(BenchmarkId::new("dsatur", n), &g, |b, g| b.iter(|| baselines::dsatur(black_box(g))));
    }
    group.finish();
}

/// Proposition 3 (number of colourings) and Proposition 4 (chromatic sum).
fn counting(c: &mut Criterion) {
    let mut group = c.benchmark_group("counting/gnp_0.5");
    for n in [10, 12, 14, 16, 18] {
        let g = Graph::random(n, 0.5, n as u64);
        group.bench_with_input(BenchmarkId::new("chromatic_polynomial", n), &g, |b, g| {
            b.iter(|| ie::chromatic_polynomial(black_box(g)))
        });
        if n <= 14 {
            group.bench_with_input(BenchmarkId::new("chromatic_sum", n), &g, |b, g| {
                b.iter(|| ie::chromatic_sum(black_box(g)))
            });
        }
    }
    group.finish();
}

criterion_group! {
    name = benches;
    config = configured();
    targets = exact_vs_n, exact_vs_density, exact_structured, heuristics, counting
}
criterion_main!(benches);
