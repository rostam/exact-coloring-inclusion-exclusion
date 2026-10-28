use ie_coloring::baselines::{self, brute_force_chromatic_number, color_count};
use ie_coloring::graph::{Graph, iter_bits};
use ie_coloring::ie::{self, Exec};

/// Graphs with known chromatic numbers.
fn named_graphs() -> Vec<(&'static str, Graph, usize)> {
    vec![
        ("empty0", Graph::new(0), 0),
        ("empty1", Graph::new(1), 1),
        ("empty7", Graph::new(7), 1),
        ("K2", Graph::complete(2), 2),
        ("K7", Graph::complete(7), 7),
        ("K12", Graph::complete(12), 12),
        ("P9", Graph::path(9), 2),
        ("C5", Graph::cycle(5), 3),
        ("C8", Graph::cycle(8), 2),
        ("C13", Graph::cycle(13), 3),
        ("W6", Graph::wheel(6), 4),
        ("W7", Graph::wheel(7), 3),
        ("K3,4", Graph::complete_bipartite(3, 4), 2),
        ("Petersen", Graph::petersen(), 3),
        ("Grotzsch", Graph::mycielski(4), 4),
        ("Mycielski5", Graph::mycielski(5), 5),
        ("Queen3x3", Graph::queen(3, 3), 5),
        ("Queen4x4", Graph::queen(4, 4), 5),
        ("Queen5x5", Graph::queen(5, 5), 5),
    ]
}

/// A spread of random graphs: sizes 1..=13 and densities 0.1..=0.9.
fn random_graphs() -> impl Iterator<Item = Graph> {
    (1..=13).flat_map(|n| (1..=9).map(move |d| Graph::random(n, d as f64 / 10.0, (n * 100 + d) as u64)))
}

#[test]
fn chromatic_number_of_named_graphs() {
    for (name, g, chi) in named_graphs() {
        assert_eq!(ie::chromatic_number(&g), chi, "{name}");
        assert_eq!(ie::chromatic_number_unbracketed(&g), chi, "{name} (unbracketed)");
        assert_eq!(color_count(&baselines::dsatur_exact(&g)), chi, "{name} (DSATUR B&B)");
        if g.n() <= 23 {
            assert_eq!(ie::chromatic_number_polyspace(&g), chi, "{name} (polyspace)");
            assert_eq!(baselines::lawler(&g), chi, "{name} (Lawler)");
        }
    }
}

#[test]
fn all_exact_algorithms_agree_on_random_graphs() {
    for g in random_graphs() {
        let chi = brute_force_chromatic_number(&g);
        assert_eq!(ie::chromatic_number(&g), chi, "{:?}", g.edges());
        assert_eq!(ie::chromatic_number_unbracketed(&g), chi);
        assert_eq!(ie::chromatic_number_polyspace(&g), chi);
        assert_eq!(baselines::lawler(&g), chi);
        let bb = baselines::dsatur_exact(&g);
        assert!(g.is_proper_coloring(&bb));
        assert_eq!(color_count(&bb), chi);
    }
}

#[test]
fn k_colorability_is_a_threshold_at_chi() {
    for g in random_graphs().step_by(7) {
        let chi = brute_force_chromatic_number(&g);
        for k in 0..=g.n() + 1 {
            assert_eq!(ie::is_k_colorable(&g, k), k >= chi, "k = {k}, edges {:?}", g.edges());
        }
    }
}

#[test]
fn heuristics_produce_proper_colorings_no_better_than_chi() {
    for g in random_graphs() {
        let chi = ie::chromatic_number(&g);
        let order: Vec<usize> = (0..g.n()).collect();
        for c in [baselines::greedy(&g, &order), baselines::largest_first(&g), baselines::dsatur(&g)] {
            assert!(g.is_proper_coloring(&c));
            assert!(color_count(&c) >= chi);
        }
        assert!(baselines::greedy_clique(&g).len() <= chi);
    }
}

/// Lemma 2: c_k(G) = Σ_Y (-1)^{n-|Y|} i(Y)^k equals the number of k-tuples of
/// nonempty stable sets whose union is V, computed by enumeration.
#[test]
fn lemma2_counts_covers_exactly() {
    for g in random_graphs().filter(|g| g.n() <= 6) {
        let n = g.n();
        let counts = ie::independent_set_counts(&g);
        let stable: Vec<u64> = (1..1u64 << n).filter(|&s| g.is_independent(s)).collect();
        assert_eq!(counts[(1 << n) - 1] as usize, stable.len());
        for k in 1..=3 {
            let formula: i128 = counts
                .iter()
                .enumerate()
                .map(|(y, &c)| {
                    let sign = if (n as u32 - (y as u64).count_ones()).is_multiple_of(2) { 1 } else { -1 };
                    sign * (c as i128).pow(k as u32)
                })
                .sum();
            let mut brute = 0i128;
            let mut idx = vec![0usize; k];
            'tuples: loop {
                if idx.iter().fold(0, |acc, &i| acc | stable[i]) == g.all() {
                    brute += 1;
                }
                for slot in idx.iter_mut() {
                    *slot += 1;
                    if *slot < stable.len() {
                        continue 'tuples;
                    }
                    *slot = 0;
                }
                break;
            }
            assert_eq!(formula, brute, "k = {k}, edges {:?}", g.edges());
        }
    }
}

#[test]
fn independent_set_counter_matches_table() {
    for g in random_graphs().filter(|g| g.n() >= 8) {
        let counts = ie::independent_set_counts(&g);
        for y in (0..1u64 << g.n()).step_by(37) {
            assert_eq!(ie::count_independent_sets(g.adjacency(), y), counts[y as usize] as u64 + 1);
        }
    }
    // Pure paths and cycles exercise the Fibonacci / Lucas closed forms.
    assert_eq!(ie::count_independent_sets(Graph::path(10).adjacency(), Graph::path(10).all()), 144);
    assert_eq!(ie::count_independent_sets(Graph::cycle(10).adjacency(), Graph::cycle(10).all()), 123);
}

#[test]
fn maximum_independent_set_is_maximum() {
    for g in random_graphs() {
        let mis = ie::maximum_independent_set(g.adjacency(), g.all());
        assert!(g.is_independent(mis));
        let alpha = (0..1u64 << g.n()).filter(|&s| g.is_independent(s)).map(u64::count_ones).max().unwrap();
        assert_eq!(mis.count_ones(), alpha);
    }
}

// ----- number of colourings --------------------------------------------------

/// Chromatic polynomial by deletion–contraction, as coefficients of k^j.
fn deletion_contraction(adj: &[u64]) -> Vec<i128> {
    let n = adj.len();
    let Some(u) = (0..n).find(|&u| adj[u] != 0) else {
        let mut p = vec![0; n + 1];
        p[n] = 1;
        return p;
    };
    let v = adj[u].trailing_zeros() as usize;
    let mut deleted = adj.to_vec();
    deleted[u] &= !(1 << v);
    deleted[v] &= !(1 << u);
    // Contract v into u, then drop v and renumber.
    let mut merged = deleted.clone();
    for w in iter_bits(deleted[v]) {
        merged[u] |= 1 << w;
        merged[w] |= 1 << u;
    }
    let low = (1u64 << v) - 1;
    let contracted: Vec<u64> =
        (0..n).filter(|&w| w != v).map(|w| (merged[w] & low) | ((merged[w] >> (v + 1)) << v)).collect();
    let a = deletion_contraction(&deleted);
    let b = deletion_contraction(&contracted);
    (0..=n).map(|j| a[j] - b.get(j).copied().unwrap_or(0)).collect()
}

fn eval(poly: &[i128], k: i128) -> i128 {
    poly.iter().rev().fold(0, |acc, &c| acc * k + c)
}

#[test]
fn chromatic_polynomial_matches_deletion_contraction() {
    for g in random_graphs().filter(|g| g.n() <= 9) {
        assert_eq!(ie::chromatic_polynomial(&g), deletion_contraction(g.adjacency()), "{:?}", g.edges());
    }
    let pet = Graph::petersen();
    assert_eq!(ie::chromatic_polynomial(&pet), deletion_contraction(pet.adjacency()));
}

#[test]
fn chromatic_polynomial_closed_forms() {
    for n in 1..=12usize {
        // K_n: k (k-1) ... (k-n+1); P_n (a tree): k (k-1)^(n-1).
        let kn = ie::chromatic_polynomial(&Graph::complete(n));
        let tree = ie::chromatic_polynomial(&Graph::path(n));
        for k in 0..=14i128 {
            assert_eq!(eval(&kn, k), (0..n as i128).map(|i| k - i).product::<i128>());
            assert_eq!(eval(&tree, k), k * (k - 1).pow(n as u32 - 1));
            if n >= 3 {
                // C_n: (k-1)^n + (-1)^n (k-1).
                let cyc = ie::chromatic_polynomial(&Graph::cycle(n));
                let sign = if n % 2 == 0 { 1 } else { -1 };
                assert_eq!(eval(&cyc, k), (k - 1).pow(n as u32) + sign * (k - 1));
            }
        }
    }
}

#[test]
fn stable_partitions_of_edgeless_graph_are_stirling_numbers() {
    // π(E_n, r) = S(n, r), Stirling numbers of the second kind.
    let n = 15;
    let mut stirling = vec![vec![0u128; n + 1]; n + 1];
    stirling[0][0] = 1;
    for i in 1..=n {
        for r in 1..=i {
            stirling[i][r] = r as u128 * stirling[i - 1][r] + stirling[i - 1][r - 1];
        }
    }
    assert_eq!(ie::stable_partition_counts(&Graph::new(n)), stirling[n]);
    let mut expect = vec![0u128; 11];
    expect[10] = 1;
    assert_eq!(ie::stable_partition_counts(&Graph::complete(10)), expect);
}

#[test]
fn count_colorings_matches_enumeration() {
    for g in random_graphs().filter(|g| g.n() <= 7) {
        for k in 0..=4u64 {
            let n = g.n();
            let brute = (0..k.pow(n as u32))
                .filter(|&code| {
                    let colour: Vec<usize> = (0..n).map(|v| (code / k.pow(v as u32) % k) as usize).collect();
                    g.is_proper_coloring(&colour)
                })
                .count();
            assert_eq!(ie::count_colorings(&g, k), Some(brute as u128), "k = {k}");
        }
    }
    // Large values still exact: 20 colours of the Petersen graph.
    let pet = Graph::petersen();
    assert_eq!(ie::count_colorings(&pet, 20), Some(eval(&deletion_contraction(pet.adjacency()), 20) as u128));
}

// ----- chromatic sum, optimal colouring, approximation -----------------------

fn brute_force_chromatic_sum(g: &Graph) -> usize {
    let n = g.n();
    let k = g.max_degree() as u64 + 1;
    (0..k.pow(n as u32))
        .filter_map(|code| {
            let colour: Vec<usize> = (0..n).map(|v| (code / k.pow(v as u32) % k) as usize).collect();
            g.is_proper_coloring(&colour).then(|| baselines::coloring_sum(&colour))
        })
        .min()
        .unwrap_or(0)
}

#[test]
fn chromatic_sum_matches_brute_force() {
    for g in random_graphs().filter(|g| g.n() <= 8) {
        assert_eq!(ie::chromatic_sum(&g), brute_force_chromatic_sum(&g), "{:?}", g.edges());
    }
    assert_eq!(ie::chromatic_sum(&Graph::complete(6)), 21);
    assert_eq!(ie::chromatic_sum(&Graph::complete_bipartite(3, 5)), 3 * 2 + 5);
    assert_eq!(ie::chromatic_sum(&Graph::petersen()), brute_force_chromatic_sum(&Graph::petersen()));
}

#[test]
fn optimal_coloring_is_proper_and_optimal() {
    let graphs = named_graphs().into_iter().map(|(_, g, _)| g).filter(|g| g.n() <= 16).chain(random_graphs());
    for g in graphs {
        let colouring = ie::optimal_coloring(&g);
        assert!(g.is_proper_coloring(&colouring), "{:?}", g.edges());
        assert_eq!(color_count(&colouring), ie::chromatic_number(&g));
    }
}

#[test]
fn approximation_is_within_guarantee() {
    for g in random_graphs().filter(|g| g.n() >= 6).chain([Graph::mycielski(5), Graph::queen(4, 5)]) {
        let chi = ie::chromatic_number(&g);
        for eps in [0.0, 0.1, 0.3, 0.7, 1.0, 2.0] {
            let approx = ie::approximate_chromatic_number(&g, eps);
            assert!(
                chi <= approx && approx as f64 <= ((1.0 + eps) * chi as f64).ceil(),
                "eps {eps}: {chi} vs {approx}"
            );
        }
        assert_eq!(ie::approximate_chromatic_number(&g, 0.0), chi);
    }
}

// ----- parallel evaluation ---------------------------------------------------

/// Random graphs large enough that the subset sums span many parallel chunks.
fn chunked_graphs(max_n: usize) -> impl Iterator<Item = Graph> {
    (13..=max_n).flat_map(|n| [0.2, 0.5, 0.8].map(move |p| Graph::random(n, p, n as u64 * 7 + (p * 10.0) as u64)))
}

#[test]
fn parallel_chromatic_number_matches_sequential() {
    for g in chunked_graphs(22) {
        let chi = ie::chromatic_number_with(&g, Exec::Sequential);
        assert_eq!(ie::chromatic_number_with(&g, Exec::Parallel), chi, "{:?}", g.edges());
        assert_eq!(color_count(&baselines::dsatur_exact(&g)), chi);
        if g.n() <= 18 {
            assert_eq!(ie::chromatic_number_polyspace_with(&g, Exec::Parallel), chi);
            assert_eq!(ie::chromatic_number_polyspace_with(&g, Exec::Sequential), chi);
        }
    }
    for g in [Graph::mycielski(5), Graph::queen(5, 5)] {
        assert_eq!(ie::chromatic_number_with(&g, Exec::Parallel), ie::chromatic_number_with(&g, Exec::Sequential));
    }
}

#[test]
fn parallel_counting_matches_sequential() {
    for g in chunked_graphs(17) {
        assert_eq!(
            ie::stable_partition_counts_with(&g, Exec::Parallel),
            ie::stable_partition_counts_with(&g, Exec::Sequential),
            "{:?}",
            g.edges()
        );
    }
    for g in chunked_graphs(14) {
        assert_eq!(ie::chromatic_sum_with(&g, Exec::Parallel), ie::chromatic_sum_with(&g, Exec::Sequential));
    }
}

/// The block-parallel table (recursion on the highest vertex) must equal the
/// original sequential table (recursion on the lowest vertex).
#[test]
fn parallel_table_matches_lowest_vertex_recursion() {
    for g in chunked_graphs(20).chain([Graph::new(16), Graph::complete(15), Graph::mycielski(5)]) {
        let adj = g.adjacency();
        let mut reference = vec![0u32; 1 << g.n()];
        for y in 1..1usize << g.n() {
            let v = y.trailing_zeros() as usize;
            reference[y] = reference[y & (y - 1)] + reference[y & (y - 1) & !(adj[v] as usize)] + 1;
        }
        assert_eq!(ie::independent_set_counts_with(&g, Exec::Parallel), reference, "{:?}", g.edges());
        assert_eq!(ie::independent_set_counts_with(&g, Exec::Sequential), reference);
    }
}

// ----- Petersen graph and generalized Petersen graphs --------------------------

/// Every valid `GP(n, k)` with at most `2 * max_n` vertices.
fn generalized_petersen_graphs(max_n: usize) -> impl Iterator<Item = (usize, usize, Graph)> {
    (3..=max_n).flat_map(|n| (1..n.div_ceil(2)).map(move |k| (n, k, Graph::generalized_petersen(n, k))))
}

/// `GP(n, k)` is bipartite exactly when `n` is even and `k` is odd; otherwise
/// it contains an odd cycle and, being cubic and not `K_4`, has `χ = 3` by
/// Brooks' theorem (Castagna and Prins, 1972).
fn gp_chromatic_number(n: usize, k: usize) -> usize {
    if n.is_multiple_of(2) && !k.is_multiple_of(2) { 2 } else { 3 }
}

/// Named generalized Petersen graphs: (name, n, k, χ, α).
const NAMED_GP: [(&str, usize, usize, usize, u32); 8] = [
    ("cube Q3", 4, 1, 2, 4),
    ("pentagonal prism", 5, 1, 3, 4),
    ("Petersen", 5, 2, 3, 4),
    ("Dürer", 6, 2, 3, 4),
    ("Möbius–Kantor", 8, 3, 2, 8),
    ("dodecahedron", 10, 2, 3, 8),
    ("Desargues", 10, 3, 2, 10),
    ("Nauru", 12, 5, 2, 12),
];

#[test]
fn generalized_petersen_structure() {
    assert_eq!(Graph::petersen(), Graph::generalized_petersen(5, 2));
    for (n, k, g) in generalized_petersen_graphs(15) {
        assert_eq!(g.n(), 2 * n);
        assert_eq!(g.edge_count(), 3 * n, "GP({n},{k})");
        assert!((0..2 * n).all(|v| g.degree(v) == 3), "GP({n},{k}) is cubic");
    }
    // The Petersen graph is triangle-free (girth 5).
    let pet = Graph::petersen();
    for (u, v) in pet.edges() {
        assert_eq!(pet.neighbors(u) & pet.neighbors(v), 0);
    }
    for (name, n, k, _, alpha) in NAMED_GP {
        let g = Graph::generalized_petersen(n, k);
        let brute = (0..1u64 << g.n()).filter(|&s| g.is_independent(s)).map(u64::count_ones).max().unwrap();
        assert_eq!(brute, alpha, "{name}");
        assert_eq!(ie::maximum_independent_set(g.adjacency(), g.all()).count_ones(), alpha, "{name}");
    }
}

#[test]
fn generalized_petersen_chromatic_numbers() {
    for (name, n, k, chi, _) in NAMED_GP {
        assert_eq!(gp_chromatic_number(n, k), chi, "{name}");
    }
    // Up to GP(12, k): 24 vertices, a 64 MiB table.
    for (n, k, g) in generalized_petersen_graphs(12) {
        let chi = gp_chromatic_number(n, k);
        let label = format!("GP({n},{k})");
        assert_eq!(ie::chromatic_number_with(&g, Exec::Parallel), chi, "{label}");
        assert_eq!(ie::chromatic_number_with(&g, Exec::Sequential), chi, "{label}");
        assert_eq!(ie::chromatic_number_unbracketed(&g), chi, "{label}");
        assert!(!ie::is_k_colorable(&g, chi - 1) && ie::is_k_colorable(&g, chi), "{label}");
        assert_eq!(color_count(&baselines::dsatur_exact(&g)), chi, "{label}");
        if 2 * n <= 20 {
            assert_eq!(ie::chromatic_number_polyspace(&g), chi, "{label}");
            assert_eq!(baselines::lawler(&g), chi, "{label}");
        }
        if 2 * n <= 18 {
            assert_eq!(brute_force_chromatic_number(&g), chi, "{label}");
        }
    }
}

#[test]
fn generalized_petersen_colouring_counts() {
    for (n, k, g) in generalized_petersen_graphs(9) {
        // A connected bipartite graph has exactly two 2-colourings.
        let two = if gp_chromatic_number(n, k) == 2 { 2 } else { 0 };
        assert_eq!(ie::count_colorings(&g, 2), Some(two), "GP({n},{k})");
    }
    for (n, k, g) in generalized_petersen_graphs(6) {
        assert_eq!(ie::chromatic_polynomial(&g), deletion_contraction(g.adjacency()), "GP({n},{k})");
        let threes = (0..3u64.pow(2 * n as u32))
            .filter(|&code| {
                let colour: Vec<usize> = (0..2 * n).map(|v| (code / 3u64.pow(v as u32) % 3) as usize).collect();
                g.is_proper_coloring(&colour)
            })
            .count();
        assert_eq!(ie::count_colorings(&g, 3), Some(threes as u128), "GP({n},{k})");
    }
    // The Petersen graph has 120 proper 3-colourings, all equivalent up to
    // permuting the colours: 120 = 3! · 20.
    assert_eq!(ie::count_colorings(&Graph::petersen(), 3), Some(120));
    assert_eq!(ie::stable_partition_counts(&Graph::petersen())[3], 20);
}

#[test]
fn generalized_petersen_sum_optimal_and_approximate() {
    for (n, k, g) in generalized_petersen_graphs(8) {
        let label = format!("GP({n},{k})");
        let chi = gp_chromatic_number(n, k);
        let sum = ie::chromatic_sum(&g);
        if 2 * n <= 10 {
            assert_eq!(sum, brute_force_chromatic_sum(&g), "{label}");
        }
        if chi == 2 {
            // Bipartite and cubic: colour classes of size n, so Σ = n·1 + n·2.
            assert_eq!(sum, 3 * n, "{label}");
        }
        assert_eq!(ie::chromatic_sum_with(&g, Exec::Sequential), sum, "{label}");
    }
    for (n, k, g) in generalized_petersen_graphs(10) {
        let colouring = ie::optimal_coloring(&g);
        assert!(g.is_proper_coloring(&colouring), "GP({n},{k})");
        assert_eq!(color_count(&colouring), gp_chromatic_number(n, k), "GP({n},{k})");
        for eps in [0.0, 0.5, 1.0] {
            let approx = ie::approximate_chromatic_number(&g, eps);
            let chi = gp_chromatic_number(n, k);
            assert!(chi <= approx && approx as f64 <= ((1.0 + eps) * chi as f64).ceil(), "GP({n},{k}) eps {eps}");
        }
    }
}
