//! Side-by-side comparison of the inclusion–exclusion algorithms with
//! classical exact and heuristic colourings. Prints a Markdown table.
//!
//! ```text
//! cargo run --release --example compare [max_n]
//! ```

use ie_coloring::baselines::{self, color_count};
use ie_coloring::{Graph, ie};
use std::time::Instant;

/// Runs `f` and returns its result with the wall-clock time in milliseconds.
fn timed<T>(f: impl FnOnce() -> T) -> (T, f64) {
    let start = Instant::now();
    let out = f();
    (out, start.elapsed().as_secs_f64() * 1e3)
}

fn cell(result: Option<(usize, f64)>, chi: usize) -> String {
    match result {
        None => "–".into(),
        Some((value, ms)) => {
            assert_eq!(value, chi, "exact algorithms disagree");
            format!("{ms:.1}")
        }
    }
}

fn main() {
    ie_coloring::modular::primes(); // build the prime table before timing anything
    let max_n: usize = std::env::args().nth(1).map_or(24, |s| s.parse().expect("max_n must be a number"));
    let mut instances: Vec<(String, Graph)> = Vec::new();
    for n in (10..=max_n).step_by(2) {
        for p in [0.2, 0.5, 0.8] {
            instances.push((format!("G({n},{p})"), Graph::random(n, p, (n as u64) << 8 | (p * 10.0) as u64)));
        }
    }
    instances.push(("Petersen".into(), Graph::petersen()));
    instances.push(("Grötzsch (M4)".into(), Graph::mycielski(4)));
    instances.push(("Mycielski M5".into(), Graph::mycielski(5)));
    instances.push(("Queen 4x5".into(), Graph::queen(4, 5)));
    instances.push(("Queen 5x5".into(), Graph::queen(5, 5)));
    instances.retain(|(_, g)| g.n() <= max_n.max(25));

    println!("Times in milliseconds; `–` = skipped (too large for that algorithm).\n");
    println!(
        "| graph | n | m | χ | IE (Prop. 1) | IE polyspace (Prop. 2) | Lawler | DSATUR B&B | greedy | largest-first | DSATUR |"
    );
    println!("|---|---|---|---|---|---|---|---|---|---|---|");
    for (name, g) in &instances {
        let n = g.n();
        let (chi, ie_ms) = timed(|| ie::chromatic_number(g));
        let poly = (n <= 22).then(|| timed(|| ie::chromatic_number_polyspace(g)));
        let lawler = (n <= 22).then(|| timed(|| baselines::lawler(g)));
        let bnb = timed(|| color_count(&baselines::dsatur_exact(g)));
        let order: Vec<usize> = (0..n).collect();
        println!(
            "| {name} | {n} | {} | {chi} | {ie_ms:.1} | {} | {} | {} | {} | {} | {} |",
            g.edge_count(),
            cell(poly, chi),
            cell(lawler, chi),
            cell(Some(bnb), chi),
            color_count(&baselines::greedy(g, &order)),
            color_count(&baselines::largest_first(g)),
            color_count(&baselines::dsatur(g)),
        );
    }
}
