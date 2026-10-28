//! Well-known colouring algorithms used as references and for comparison:
//!
//! * [`greedy`] / [`largest_first`] – first-fit heuristics,
//! * [`dsatur`] – Brélaz's DSATUR heuristic (1979),
//! * [`dsatur_exact`] – exact DSATUR-based branch and bound (Brélaz; Brown),
//! * [`lawler`] – Lawler's `O(2.4423^n)` dynamic program over maximal
//!   independent sets (1976), the first entry of Table 1 in the paper,
//! * [`brute_force_chromatic_number`] – plain backtracking, for testing.
//!
//! Colourings are `Vec<usize>` with colours `0, 1, 2, ...`.

use crate::graph::{Graph, full_mask, iter_bits};

/// Number of distinct colours used (colours are assumed to be `0..count`).
pub fn color_count(colouring: &[usize]) -> usize {
    colouring.iter().map(|&c| c + 1).max().unwrap_or(0)
}

/// `Σ_v C(v)` with colours counted from 1, as in the chromatic sum problem.
pub fn coloring_sum(colouring: &[usize]) -> usize {
    colouring.iter().map(|&c| c + 1).sum()
}

/// First-fit greedy colouring in the given vertex order.
pub fn greedy(g: &Graph, order: &[usize]) -> Vec<usize> {
    let mut colour = vec![usize::MAX; g.n()];
    let mut classes: Vec<u64> = Vec::new();
    for &v in order {
        let c = classes.iter().position(|&cls| cls & g.neighbors(v) == 0).unwrap_or(classes.len());
        if c == classes.len() {
            classes.push(0);
        }
        classes[c] |= 1 << v;
        colour[v] = c;
    }
    colour
}

/// Welsh–Powell: greedy in order of decreasing degree.
pub fn largest_first(g: &Graph) -> Vec<usize> {
    let mut order: Vec<usize> = (0..g.n()).collect();
    order.sort_by_key(|&v| std::cmp::Reverse(g.degree(v)));
    greedy(g, &order)
}

/// Brélaz's DSATUR heuristic: repeatedly colour the uncoloured vertex seeing
/// the most distinct colours among its neighbours (ties: most uncoloured
/// neighbours) with the smallest available colour.
pub fn dsatur(g: &Graph) -> Vec<usize> {
    let mut colour = vec![usize::MAX; g.n()];
    let mut classes: Vec<u64> = Vec::new();
    let mut uncoloured = g.all();
    while uncoloured != 0 {
        let v = select_dsatur(g, &classes, uncoloured);
        let c = classes.iter().position(|&cls| cls & g.neighbors(v) == 0).unwrap_or(classes.len());
        if c == classes.len() {
            classes.push(0);
        }
        classes[c] |= 1 << v;
        colour[v] = c;
        uncoloured &= !(1 << v);
    }
    colour
}

fn saturation(g: &Graph, classes: &[u64], v: usize) -> usize {
    classes.iter().filter(|&&cls| cls & g.neighbors(v) != 0).count()
}

fn select_dsatur(g: &Graph, classes: &[u64], uncoloured: u64) -> usize {
    iter_bits(uncoloured)
        .max_by_key(|&v| (saturation(g, classes, v), (g.neighbors(v) & uncoloured).count_ones()))
        .unwrap()
}

/// A clique found greedily (by decreasing degree); its size is a lower bound on χ.
pub fn greedy_clique(g: &Graph) -> Vec<usize> {
    let mut order: Vec<usize> = (0..g.n()).collect();
    order.sort_by_key(|&v| std::cmp::Reverse(g.degree(v)));
    let mut best = Vec::new();
    for &start in &order {
        let mut clique = vec![start];
        let mut candidates = g.neighbors(start);
        for &v in &order {
            if candidates >> v & 1 == 1 {
                clique.push(v);
                candidates &= g.neighbors(v);
            }
        }
        if clique.len() > best.len() {
            best = clique;
        }
    }
    best
}

/// Exact DSATUR branch and bound. Starts from the DSATUR heuristic as the
/// incumbent and a greedy clique (precoloured, to break colour symmetry) as
/// the lower bound, then backtracks over DSATUR vertex choices.
pub fn dsatur_exact(g: &Graph) -> Vec<usize> {
    let n = g.n();
    if n == 0 {
        return Vec::new();
    }
    let mut best = dsatur(g);
    let clique = greedy_clique(g);
    if color_count(&best) == clique.len() {
        return best;
    }
    let mut colour = vec![usize::MAX; n];
    let mut classes = Vec::new();
    let mut uncoloured = g.all();
    for (c, &v) in clique.iter().enumerate() {
        colour[v] = c;
        classes.push(1u64 << v);
        uncoloured &= !(1 << v);
    }
    let mut search = BranchAndBound { g, lower: clique.len(), best_count: color_count(&best), best: &mut best };
    search.recurse(&mut colour, &mut classes, uncoloured);
    best
}

struct BranchAndBound<'a> {
    g: &'a Graph,
    lower: usize,
    best_count: usize,
    best: &'a mut Vec<usize>,
}

impl BranchAndBound<'_> {
    /// Returns `true` once a colouring matching the lower bound is found.
    fn recurse(&mut self, colour: &mut Vec<usize>, classes: &mut Vec<u64>, uncoloured: u64) -> bool {
        if uncoloured == 0 {
            self.best_count = classes.len();
            self.best.clone_from(colour);
            return self.best_count == self.lower;
        }
        if classes.len() >= self.best_count {
            return false;
        }
        let v = select_dsatur(self.g, classes, uncoloured);
        let rest = uncoloured & !(1 << v);
        for c in 0..classes.len() {
            if classes[c] & self.g.neighbors(v) == 0 {
                classes[c] |= 1 << v;
                colour[v] = c;
                let done = self.recurse(colour, classes, rest);
                classes[c] &= !(1 << v);
                if done {
                    return true;
                }
                if classes.len() >= self.best_count {
                    // The incumbent improved below the current palette: prune.
                    break;
                }
            }
        }
        if classes.len() + 1 < self.best_count {
            classes.push(1 << v);
            colour[v] = classes.len() - 1;
            let done = self.recurse(colour, classes, rest);
            classes.pop();
            if done {
                return true;
            }
        }
        colour[v] = usize::MAX;
        false
    }
}

/// Lawler's algorithm: `χ(S) = 1 + min { χ(S - I) : I maximal stable in G[S] }`
/// over all subsets `S` in increasing order. Time `O(2.4423^n)`, space `2^n`.
pub fn lawler(g: &Graph) -> usize {
    let n = g.n();
    assert!(n <= 30, "Lawler's table supports at most 30 vertices");
    let adj = g.adjacency();
    let mut chi = vec![0u8; 1 << n];
    for s in 1..1u64 << n {
        let mut best = u8::MAX;
        maximal_independent_sets(adj, s, 0, s, 0, &mut |i| {
            best = best.min(chi[(s & !i) as usize] + 1);
        });
        chi[s as usize] = best;
    }
    chi[full_mask(n) as usize] as usize
}

/// Bron–Kerbosch with pivoting on the complement graph restricted to `within`:
/// reports every maximal stable set of `G[within]`.
fn maximal_independent_sets(adj: &[u64], within: u64, r: u64, p: u64, x: u64, report: &mut impl FnMut(u64)) {
    if p | x == 0 {
        report(r);
        return;
    }
    let non_nbrs = |v: usize| within & !adj[v] & !(1 << v);
    let pivot = iter_bits(p | x).max_by_key(|&u| (p & non_nbrs(u)).count_ones()).unwrap();
    let (mut p, mut x) = (p, x);
    for v in iter_bits(p & !non_nbrs(pivot)) {
        let nv = non_nbrs(v);
        maximal_independent_sets(adj, within, r | 1 << v, p & nv, x & nv, report);
        p &= !(1 << v);
        x |= 1 << v;
    }
}

/// Smallest `k` admitting a colouring found by plain backtracking. Only for
/// small test graphs.
pub fn brute_force_chromatic_number(g: &Graph) -> usize {
    (0..=g.n()).find(|&k| backtrack_colourable(g, k, 0, &mut vec![0; g.n()])).unwrap()
}

fn backtrack_colourable(g: &Graph, k: usize, v: usize, colour: &mut Vec<usize>) -> bool {
    if v == g.n() {
        return true;
    }
    for c in 0..k {
        if (0..v).all(|u| !g.has_edge(u, v) || colour[u] != c) {
            colour[v] = c;
            if backtrack_colourable(g, k, v + 1, colour) {
                return true;
            }
        }
    }
    false
}
