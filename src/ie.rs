//! The inclusion–exclusion colouring algorithms of Björklund and Husfeldt,
//! "Inclusion–Exclusion Based Algorithms for Graph Colouring" (ECCC TR06-044).
//!
//! # Notation
//!
//! The paper writes `s(X)` for the number of nonempty stable sets of `G` that
//! avoid `X`. We index by the complement instead: `i(Y) = s(V - Y)` is the
//! number of nonempty stable sets *contained in* `Y`. With that, Lemma 2 reads
//!
//! ```text
//! c_k(G) = sum_{Y ⊆ V} (-1)^{n - |Y|} i(Y)^k,
//! ```
//!
//! the number of ways to cover `V` with `k` stable sets, and by Lemma 1
//! `χ(G) = min { k : c_k(G) > 0 }`.
//!
//! All sums are evaluated with the exact multi-modulus scheme from
//! [`crate::modular`], so every answer is deterministic and exact.
//!
//! The `2^n`-term inclusion–exclusion sums are split into chunks of subsets
//! that are summed in parallel with rayon. Modular and wrapping arithmetic is
//! associative and commutative, so the result does not depend on how the
//! chunks are combined: parallel and sequential runs return identical values.
//! The `*_with` functions take an [`Exec`] to choose explicitly.

use std::ops::Range;

use rayon::prelude::*;

use crate::baselines;
use crate::graph::{Graph, full_mask, iter_bits};
use crate::modular::{self, ModEval, Modulus, Montgomery, Wrap64};

/// How the inclusion–exclusion sums over all subsets are evaluated.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Exec {
    Sequential,
    /// Chunks of subsets summed on rayon's global thread pool.
    Parallel,
}

/// Subsets per parallel task: large enough to amortise scheduling, small
/// enough to balance load.
const CHUNK: usize = 1 << 12;

/// Folds `chunk` over `0..len` split into ranges and merges the partial
/// results with `combine`, which must be associative and commutative.
fn fold_subsets<T: Send>(
    len: usize,
    exec: Exec,
    chunk: impl Fn(Range<usize>) -> T + Sync,
    combine: impl Fn(T, T) -> T + Sync,
) -> T {
    if exec == Exec::Sequential || len <= CHUNK {
        return chunk(0..len);
    }
    (0..len.div_ceil(CHUNK))
        .into_par_iter()
        .map(|i| chunk(i * CHUNK..((i + 1) * CHUNK).min(len)))
        .reduce_with(&combine)
        .expect("at least one chunk")
}

/// Largest graph for which the `2^n`-entry table of [`independent_set_counts`]
/// is built (about 4 GiB at `n = 30`).
pub const MAX_TABLE_VERTICES: usize = 30;

/// Largest graph accepted by the polynomial-size-table algorithms
/// ([`stable_partition_counts`], [`chromatic_sum`]). Ordered partition counts
/// are kept in `u128`, which is exact up to `25^25 < 2^128`.
pub const MAX_POLY_VERTICES: usize = 25;

// ----- Proposition 1: chromatic number in time and space 2^n n^O(1) ----------

/// `i(Y)` for every `Y ⊆ V`, indexed by the bitmask of `Y`: the number of
/// nonempty stable sets contained in `Y`.
///
/// This is recursion (5) of the paper, written for `Y = V - X`: pick any
/// `v ∈ Y`; a stable set inside `Y` either avoids `v`, or contains `v` and
/// avoids `N(v)`, or is exactly `{v}`:
///
/// ```text
/// i(Y) = i(Y - v) + i(Y - N[v]) + 1.
/// ```
///
/// We always take `v` to be the highest vertex of `Y`. Then for the block of
/// subsets `Y` with highest vertex `j`, i.e. `Y ∈ [2^j, 2^(j+1))`, both
/// arguments on the right lie below `2^j`. Each block therefore depends only
/// on earlier blocks, and its entries can be filled in parallel. The whole
/// table takes `O(2^n)` time.
pub fn independent_set_counts(g: &Graph) -> Vec<u32> {
    independent_set_counts_with(g, Exec::Parallel)
}

/// [`independent_set_counts`] with an explicit execution strategy.
pub fn independent_set_counts_with(g: &Graph, exec: Exec) -> Vec<u32> {
    let n = g.n();
    assert!(n <= MAX_TABLE_VERTICES, "table algorithms support at most {MAX_TABLE_VERTICES} vertices");
    let mut t = vec![0u32; 1 << n];
    fill_blocks(&mut t, 1, g.adjacency(), exec, |done, low, non_nbrs, out| {
        out[0] = done[low] + done[low & non_nbrs] + 1;
    });
    t
}

/// Fills a table of `2^n` rows of `stride` entries, block by block: the rows
/// `Y = low | 2^j` for `low < 2^j` are computed by
/// `row(done, low, !N(j), out)`, where `done` holds all rows below `2^j`.
/// Row 0 must be initialised by the caller.
fn fill_blocks(
    t: &mut [u32],
    stride: usize,
    adj: &[u64],
    exec: Exec,
    row: impl Fn(&[u32], usize, usize, &mut [u32]) + Sync,
) {
    for (j, &nbrs) in adj.iter().enumerate() {
        let (done, rest) = t.split_at_mut(stride << j);
        let block = &mut rest[..stride << j];
        let non_nbrs = !(nbrs as usize);
        let fill = |first: usize, out: &mut [u32]| {
            for (i, out_row) in out.chunks_exact_mut(stride).enumerate() {
                row(done, first + i, non_nbrs, out_row);
            }
        };
        if exec == Exec::Parallel && 1 << j > CHUNK {
            block.par_chunks_mut(CHUNK * stride).enumerate().for_each(|(c, out)| fill(c * CHUNK, out));
        } else {
            fill(0, block);
        }
    }
}

/// Lemma 2 evaluated modulo `M`: `c_k = Σ_Y (-1)^{n-|Y|} i(Y)^k`.
struct CoverCount<'a> {
    counts: &'a [u32],
    n: usize,
    k: u64,
    exec: Exec,
}

impl ModEval for CoverCount<'_> {
    fn eval<M: Modulus>(&self, m: &M) -> Vec<u64> {
        let acc = fold_subsets(
            self.counts.len(),
            self.exec,
            |range| {
                let mut acc = [0u64; 2]; // [even, odd] parity of n - |Y|
                for y in range {
                    let term = m.pow(m.lift(self.counts[y] as u64), self.k);
                    let parity = (self.n as u32 - y.count_ones()) as usize & 1;
                    acc[parity] = m.add(acc[parity], term);
                }
                acc
            },
            |a, b| [m.add(a[0], b[0]), m.add(a[1], b[1])],
        );
        vec![m.sub(acc[0], acc[1])]
    }
}

/// Number of bits of `x`, so that `x < 2^bit_len(x)`.
fn bit_len(x: u64) -> u32 {
    64 - x.leading_zeros()
}

/// `0 <= c_k <= i(V)^k`, since `c_k` counts some of the `k`-tuples of stable sets.
fn cover_count_bound_bits(total_stable_sets: u64, k: usize) -> u32 {
    (k as u32 * bit_len(total_stable_sets)).max(1)
}

fn covers_exist(counts: &[u32], n: usize, k: usize, exec: Exec) -> bool {
    let e = CoverCount { counts, n, k: k as u64, exec };
    modular::is_nonzero(&e, cover_count_bound_bits(counts[counts.len() - 1] as u64, k))
}

/// Decides whether `G` has a proper `k`-colouring, i.e. whether `c_k(G) > 0`.
pub fn is_k_colorable(g: &Graph, k: usize) -> bool {
    if g.n() == 0 {
        return true;
    }
    if k == 0 {
        return false;
    }
    covers_exist(&independent_set_counts(g), g.n(), k, Exec::Parallel)
}

/// The chromatic number by Proposition 1: build the `i(Y)` table, then binary
/// search for the least `k` with `c_k(G) > 0`. The search is bracketed by
/// `[1, DSATUR upper bound]` instead of `[1, n]`, which only saves iterations.
pub fn chromatic_number(g: &Graph) -> usize {
    chromatic_number_with(g, Exec::Parallel)
}

/// [`chromatic_number`] with an explicit execution strategy.
pub fn chromatic_number_with(g: &Graph, exec: Exec) -> usize {
    if g.n() == 0 {
        return 0;
    }
    let counts = independent_set_counts_with(g, exec);
    let (mut lo, mut hi) = (1, baselines::color_count(&baselines::dsatur(g)));
    while lo < hi {
        let mid = (lo + hi) / 2;
        if covers_exist(&counts, g.n(), mid, exec) {
            hi = mid;
        } else {
            lo = mid + 1;
        }
    }
    lo
}

/// Prop. 1 exactly as stated in the paper: binary search over `[1, n]`
/// without any heuristic upper bound. Exposed so benchmarks can measure the
/// effect of the DSATUR bracket used by [`chromatic_number`].
pub fn chromatic_number_unbracketed(g: &Graph) -> usize {
    if g.n() == 0 {
        return 0;
    }
    let counts = independent_set_counts(g);
    let (mut lo, mut hi) = (1, g.n());
    while lo < hi {
        let mid = (lo + hi) / 2;
        if covers_exist(&counts, g.n(), mid, Exec::Parallel) {
            hi = mid;
        } else {
            lo = mid + 1;
        }
    }
    lo
}

// ----- Proposition 2: polynomial space ---------------------------------------

/// Number of stable sets (including the empty set) of `G[set]`, by branching
/// on a vertex of maximum degree; components of maximum degree at most two
/// (paths and cycles) are counted in closed form. Uses polynomial space.
pub fn count_independent_sets(adj: &[u64], set: u64) -> u64 {
    let mut best = (0, 0u32);
    for v in iter_bits(set) {
        let d = (adj[v] & set).count_ones();
        if d > best.1 {
            best = (v, d);
        }
    }
    if best.1 <= 2 {
        return count_paths_and_cycles(adj, set);
    }
    let v = best.0;
    count_independent_sets(adj, set & !(1 << v)) + count_independent_sets(adj, set & !(adj[v] | 1 << v))
}

/// Stable-set count of a graph whose components are paths and cycles: a path
/// on `m` vertices has `F(m + 2)` stable sets and a cycle has `L(m)`
/// (Fibonacci and Lucas numbers).
fn count_paths_and_cycles(adj: &[u64], set: u64) -> u64 {
    let mut total = 1u64;
    let mut rest = set;
    while rest != 0 {
        let mut comp = rest & rest.wrapping_neg();
        let mut frontier = comp;
        while frontier != 0 {
            let mut next = 0;
            for v in iter_bits(frontier) {
                next |= adj[v];
            }
            frontier = next & set & !comp;
            comp |= frontier;
        }
        let m = comp.count_ones() as usize;
        let degree_sum: u32 = iter_bits(comp).map(|v| (adj[v] & set).count_ones()).sum();
        let is_cycle = degree_sum as usize == 2 * m && m >= 3;
        let (mut a, mut b) = (1u64, 1u64); // F(1), F(2)
        for _ in 0..m {
            (a, b) = (b, a + b);
        }
        // Now a = F(m + 1), b = F(m + 2), so L(m) = F(m - 1) + F(m + 1) = 3a - b.
        total *= if is_cycle { 3 * a - b } else { b };
        rest &= !comp;
    }
    total
}

/// The chromatic number by Proposition 2, in polynomial space: instead of
/// tabulating `i(Y)`, it is recomputed for every `Y` by
/// [`count_independent_sets`].
///
/// Because computing `i(Y)` is now the expensive part, a single pass over all
/// `Y` accumulates `c_k` for every candidate `k` (up to the DSATUR bound) and
/// for enough moduli to decide each one exactly.
pub fn chromatic_number_polyspace(g: &Graph) -> usize {
    chromatic_number_polyspace_with(g, Exec::Parallel)
}

/// [`chromatic_number_polyspace`] with an explicit execution strategy. Each
/// parallel task keeps its own `O(n · moduli)` accumulators, so the space
/// stays polynomial.
pub fn chromatic_number_polyspace_with(g: &Graph, exec: Exec) -> usize {
    let n = g.n();
    assert!(n < 64, "polynomial-space algorithm supports at most 63 vertices");
    if n == 0 {
        return 0;
    }
    let adj = g.adjacency();
    let hi = baselines::color_count(&baselines::dsatur(g));
    let total = count_independent_sets(adj, g.all()) - 1;
    let bound = cover_count_bound_bits(total, hi);
    let mut primes = Vec::new();
    let mut covered = Wrap64.bits();
    for &p in modular::primes() {
        if covered >= bound {
            break;
        }
        primes.push(Montgomery::new(p));
        covered += Montgomery::new(p).bits();
    }

    // sums[j][sign][k] for modulus j (0 is Wrap64, j >= 1 is primes[j - 1]).
    let sums = fold_subsets(
        1 << n,
        exec,
        |range| {
            let mut sums = vec![[vec![0u64; hi + 1], vec![0u64; hi + 1]]; primes.len() + 1];
            for y in range {
                let y = y as u64;
                let x = count_independent_sets(adj, y) - 1;
                let sign = ((n as u32 - y.count_ones()) & 1) as usize;
                add_powers(&Wrap64, x, &mut sums[0][sign]);
                for (m, s) in primes.iter().zip(&mut sums[1..]) {
                    add_powers(m, x, &mut s[sign]);
                }
            }
            sums
        },
        |mut a, b| {
            for (j, (sa, sb)) in a.iter_mut().zip(&b).enumerate() {
                for sign in 0..2 {
                    for (x, &y) in sa[sign].iter_mut().zip(&sb[sign]) {
                        *x = if j == 0 { Wrap64.add(*x, y) } else { primes[j - 1].add(*x, y) };
                    }
                }
            }
            a
        },
    );
    (1..=hi)
        .find(|&k| {
            Wrap64.sub(sums[0][0][k], sums[0][1][k]) != 0
                || primes.iter().zip(&sums[1..]).any(|(m, s)| m.sub(s[0][k], s[1][k]) != 0)
        })
        .expect("c_k > 0 for k = DSATUR bound")
}

/// Adds `x, x^2, ..., x^(len-1)` (mod `m`) to `sums[1..]`.
fn add_powers<M: Modulus>(m: &M, x: u64, sums: &mut [u64]) {
    let x = m.lift(x);
    let mut power = x;
    for s in &mut sums[1..] {
        *s = m.add(*s, power);
        power = m.mul(power, x);
    }
}

// ----- Proposition 3: number of colourings -----------------------------------

/// `s_j(Y)`, the number of stable sets of size `j` (including the empty set at
/// `j = 0`) inside `Y`, for all `Y` and `j <= stride - 1`, laid out row-major
/// with `stride = α(G) + 1`.
fn sized_independent_set_counts(g: &Graph, exec: Exec) -> (Vec<u32>, usize) {
    let n = g.n();
    assert!(n <= MAX_POLY_VERTICES, "at most {MAX_POLY_VERTICES} vertices supported");
    let stride = maximum_independent_set(g.adjacency(), g.all()).count_ones() as usize + 1;
    let mut t = vec![0u32; stride << n];
    t[0] = 1;
    // s_j(Y) = s_j(Y - v) + s_{j-1}(Y - N[v]) for the highest vertex v of Y.
    fill_blocks(&mut t, stride, g.adjacency(), exec, |done, low, non_nbrs, out| {
        let (a, b) = (low * stride, (low & non_nbrs) * stride);
        out[0] = 1;
        for j in 1..stride {
            out[j] = done[a + j] + done[b + j - 1];
        }
    });
    (t, stride)
}

/// `π_G(r)` for `r = 0..=n`: the number of partitions of `V` into exactly `r`
/// nonempty stable sets (unordered), by Lemma 4 and Proposition 3.
///
/// For each `Y`, the dynamic program `A(l, m)` of the paper counts ordered
/// `l`-tuples of nonempty stable sets inside `Y` with total size `m`; then
/// `r! π_G(r) = Σ_Y (-1)^{n-|Y|} A(r, n)`. The sum is taken in wrapping
/// `u128` arithmetic, which is exact because the true value is below `2^128`.
pub fn stable_partition_counts(g: &Graph) -> Vec<u128> {
    stable_partition_counts_with(g, Exec::Parallel)
}

/// [`stable_partition_counts`] with an explicit execution strategy.
pub fn stable_partition_counts_with(g: &Graph, exec: Exec) -> Vec<u128> {
    let n = g.n();
    if n == 0 {
        return vec![1];
    }
    let (t, stride) = sized_independent_set_counts(g, exec);
    let ordered = fold_subsets(
        1 << n,
        exec,
        |range| ordered_partition_sum(&t, stride, n, range),
        |mut a, b| {
            for (x, y) in a.iter_mut().zip(b) {
                *x = x.wrapping_add(y);
            }
            a
        },
    );
    let mut factorial = 1u128;
    let mut out = vec![0u128; n + 1];
    for r in 1..=n {
        factorial *= r as u128;
        debug_assert_eq!(ordered[r] % factorial, 0);
        out[r] = ordered[r] / factorial;
    }
    out
}

/// The part of `Σ_Y (-1)^{n-|Y|} A(r, n, Y)` (for every `r`) contributed by
/// the subsets `Y` in `range`, in wrapping `u128` arithmetic.
fn ordered_partition_sum(t: &[u32], stride: usize, n: usize, range: Range<usize>) -> Vec<u128> {
    let mut ordered = vec![0u128; n + 1];
    let mut cur = vec![0u128; n + 1];
    let mut next = vec![0u128; n + 1];
    for y in range.start.max(1)..range.end {
        let s = &t[y * stride..(y + 1) * stride];
        let alpha = s.iter().rposition(|&c| c != 0).unwrap();
        let negative = (n as u32 - y.count_ones()) & 1 == 1;
        // cur = S(z)^l restricted to degrees [l, min(n, l * alpha)].
        cur.fill(0);
        for j in 1..=alpha.min(n) {
            cur[j] = s[j] as u128;
        }
        for (l, total) in ordered.iter_mut().enumerate().skip(1) {
            let top = (l * alpha).min(n);
            if top < l {
                break;
            }
            if top == n {
                *total = if negative { total.wrapping_sub(cur[n]) } else { total.wrapping_add(cur[n]) };
            }
            // next = cur * S, degrees [l + 1, min(n, (l + 1) * alpha)].
            let new_top = ((l + 1) * alpha).min(n);
            for m in l + 1..=new_top {
                let mut acc = 0u128;
                for j in 1..=alpha.min(m - l) {
                    acc = acc.wrapping_add((s[j] as u128).wrapping_mul(cur[m - j]));
                }
                next[m] = acc;
            }
            std::mem::swap(&mut cur, &mut next);
            next.fill(0);
        }
    }
    ordered
}

/// `p_G(k)`, the number of proper colourings with colours `{1..k}`, as
/// `Σ_r π_G(r) k (k-1) ... (k-r+1)`. Returns `None` on `u128` overflow.
pub fn count_colorings(g: &Graph, k: u64) -> Option<u128> {
    let pi = stable_partition_counts(g);
    let mut total = 0u128;
    let mut falling = 1u128;
    for (r, &p) in pi.iter().enumerate() {
        if r > 0 {
            if k < r as u64 {
                break;
            }
            falling = falling.checked_mul((k - r as u64 + 1) as u128)?;
        }
        total = total.checked_add(p.checked_mul(falling)?)?;
    }
    Some(total)
}

/// Coefficients `a_0, ..., a_n` of the chromatic polynomial
/// `p_G(k) = Σ_j a_j k^j`, obtained by expanding the falling factorials of
/// [`stable_partition_counts`]. Exact for `n <= MAX_POLY_VERTICES`.
pub fn chromatic_polynomial(g: &Graph) -> Vec<i128> {
    let pi = stable_partition_counts(g);
    let n = g.n();
    let mut coeffs = vec![0i128; n + 1];
    let mut falling = vec![0i128; n + 1]; // coefficients of k (k-1) ... (k-r+1)
    falling[0] = 1;
    for (r, &p) in pi.iter().enumerate() {
        if r > 0 {
            // multiply by (k - (r - 1))
            for j in (0..=r).rev() {
                let shifted = if j > 0 { falling[j - 1] } else { 0 };
                falling[j] = shifted.wrapping_sub(falling[j].wrapping_mul(r as i128 - 1));
            }
        }
        for j in 0..=r {
            coeffs[j] = coeffs[j].wrapping_add((p as i128).wrapping_mul(falling[j]));
        }
    }
    coeffs
}

// ----- Proposition 4: chromatic sum -------------------------------------------

/// Weight enumerator of covering colour-class tuples, evaluated modulo `M`.
struct WeightedCovers<'a> {
    table: &'a [u32],
    stride: usize,
    n: usize,
    colours: usize,
    max_weight: usize,
    exec: Exec,
}

impl ModEval for WeightedCovers<'_> {
    fn eval<M: Modulus>(&self, m: &M) -> Vec<u64> {
        let k = self.max_weight;
        let (pos, neg) = fold_subsets(
            1 << self.n,
            self.exec,
            |range| self.partial_sums(m, range),
            |(mut pa, mut na), (pb, nb)| {
                for w in 0..=k {
                    pa[w] = m.add(pa[w], pb[w]);
                    na[w] = m.add(na[w], nb[w]);
                }
                (pa, na)
            },
        );
        pos.iter().zip(&neg).map(|(&a, &b)| m.sub(a, b)).collect()
    }
}

impl WeightedCovers<'_> {
    /// Sums of `Π_i Q_Y(z^i)` over the `Y` in `range`, split by the sign
    /// `(-1)^{n-|Y|}` into `(positive, negative)`.
    fn partial_sums<M: Modulus>(&self, m: &M, range: Range<usize>) -> (Vec<u64>, Vec<u64>) {
        let k = self.max_weight;
        let mut pos = vec![0u64; k + 1];
        let mut neg = vec![0u64; k + 1];
        let mut f = vec![0u64; k + 1];
        let mut g = vec![0u64; k + 1];
        for y in range {
            let s: Vec<u64> =
                self.table[y * self.stride..(y + 1) * self.stride].iter().map(|&c| m.lift(c as u64)).collect();
            // f = Π_{i=1}^{colours} Q_Y(z^i), Q_Y(z) = Σ_j s_j(Y) z^j.
            f.fill(0);
            f[0] = m.one();
            for i in 1..=self.colours {
                for w in 0..=k {
                    let mut acc = 0;
                    for (j, &sj) in s.iter().enumerate() {
                        if i * j > w {
                            break;
                        }
                        acc = m.add(acc, m.mul(sj, f[w - i * j]));
                    }
                    g[w] = acc;
                }
                std::mem::swap(&mut f, &mut g);
            }
            let target = if (self.n as u32 - y.count_ones()) & 1 == 1 { &mut neg } else { &mut pos };
            for w in 0..=k {
                target[w] = m.add(target[w], f[w]);
            }
        }
        (pos, neg)
    }
}

/// The chromatic sum `min_C Σ_v C(v)` over proper colourings `C: V -> {1, 2, ...}`
/// (Proposition 4).
///
/// The paper counts `l`-tuples of stable sets `S_1..S_l` covering `V` with
/// `|S_1| + 2|S_2| + ... + l|S_l| <= k`. We do the same in one go for all
/// weights: with `l = Δ + 1` colours (an optimal sum colouring never gives a
/// vertex a colour above `deg(v) + 1`) and allowing empty classes, the number
/// of covering tuples of weight exactly `w` is the coefficient of `z^w` in
/// `Σ_Y (-1)^{n-|Y|} Π_{i=1}^{l} Q_Y(z^i)`, and the chromatic sum is the least
/// `w` with a nonzero coefficient.
pub fn chromatic_sum(g: &Graph) -> usize {
    chromatic_sum_with(g, Exec::Parallel)
}

/// [`chromatic_sum`] with an explicit execution strategy.
pub fn chromatic_sum_with(g: &Graph, exec: Exec) -> usize {
    let n = g.n();
    if n == 0 {
        return 0;
    }
    let (table, stride) = sized_independent_set_counts(g, exec);
    let colours = g.max_degree() + 1;
    let max_weight = baselines::coloring_sum(&baselines::greedy(g, &(0..n).collect::<Vec<_>>()));
    let total_tuples_bits = colours as u32 * bit_len(table[((1 << n) - 1) * stride..].iter().map(|&c| c as u64).sum());
    let e = WeightedCovers { table: &table, stride, n, colours, max_weight, exec };
    modular::first_nonzero(&e, total_tuples_bits).expect("a greedy colouring achieves the upper bound")
}

// ----- Section 2.5: constructing an optimal colouring --------------------------

/// Removes vertex `x` from a compact adjacency list, renumbering the vertices above it.
fn remove_vertex(adj: &[u64], x: usize) -> Vec<u64> {
    let low = (1u64 << x) - 1;
    adj.iter().enumerate().filter(|&(v, _)| v != x).map(|(_, &row)| (row & low) | ((row >> (x + 1)) << x)).collect()
}

fn graph_from_adjacency(adj: &[u64]) -> Graph {
    let mut g = Graph::new(adj.len());
    for (u, &row) in adj.iter().enumerate() {
        for v in iter_bits(row) {
            if u < v {
                g.add_edge(u, v);
            }
        }
    }
    g
}

/// An optimal colouring (colours `0..χ`) built from the decision procedure, as
/// in Section 2.5 of the paper.
///
/// Repeatedly take a vertex `v` with non-neighbours `u_1, ..., u_m` and let
/// `G_i` be the graph with edges `vu_1, ..., vu_i` added. If `G_m` is still
/// `χ`-colourable, `v` can be a colour class by itself: remove it. Otherwise
/// binary search for the least `i` with `χ(G_i) > χ`; every `χ`-colouring of
/// `G_{i-1}` then gives `v` and `u_i` the same colour, so contract them. Each
/// round removes one vertex using `O(log n)` colourability tests.
pub fn optimal_coloring(g: &Graph) -> Vec<usize> {
    let mut k = chromatic_number(g);
    let mut adj = g.adjacency().to_vec();
    let mut members: Vec<u64> = (0..g.n()).map(|v| 1 << v).collect();
    let mut colour = vec![usize::MAX; g.n()];
    let mut next_colour = 0;
    let with_edges_to = |adj: &[u64], non: &[usize], i: usize| -> Vec<u64> {
        let mut a = adj.to_vec();
        for &u in &non[..i] {
            a[0] |= 1 << u;
            a[u] |= 1;
        }
        a
    };
    let colourable = |adj: &[u64], k: usize| is_k_colorable(&graph_from_adjacency(adj), k);

    while !adj.is_empty() {
        let non: Vec<usize> = iter_bits(full_mask(adj.len()) & !adj[0] & !1).collect();
        if non.is_empty() || colourable(&with_edges_to(&adj, &non, non.len()), k) {
            for v in iter_bits(members[0]) {
                colour[v] = next_colour;
            }
            next_colour += 1;
            k -= 1;
            adj = remove_vertex(&adj, 0);
            members.remove(0);
            continue;
        }
        // G_0 is k-colourable, G_m is not: find the least i with G_i not k-colourable.
        let (mut ok, mut bad) = (0, non.len());
        while bad - ok > 1 {
            let mid = (ok + bad) / 2;
            if colourable(&with_edges_to(&adj, &non, mid), k) {
                ok = mid;
            } else {
                bad = mid;
            }
        }
        let u = non[bad - 1];
        let mut a = with_edges_to(&adj, &non, bad - 1);
        let merged = a[u] & !1;
        a[0] |= merged;
        for w in iter_bits(merged) {
            a[w] |= 1;
        }
        members[0] |= members[u];
        members.remove(u);
        adj = remove_vertex(&a, u);
    }
    debug_assert_eq!(k, 0);
    colour
}

// ----- Proposition 5: approximation --------------------------------------------

/// A maximum stable set of `G[set]` as a bitmask, by simple branching in
/// polynomial space: a vertex of degree at most one is always taken, otherwise
/// branch on a vertex of maximum degree.
pub fn maximum_independent_set(adj: &[u64], set: u64) -> u64 {
    if set == 0 {
        return 0;
    }
    let (mut min_v, mut min_d, mut max_v, mut max_d) = (0, u32::MAX, 0, 0);
    for v in iter_bits(set) {
        let d = (adj[v] & set).count_ones();
        if d < min_d {
            (min_v, min_d) = (v, d);
        }
        if d > max_d {
            (max_v, max_d) = (v, d);
        }
    }
    if min_d <= 1 {
        return 1 << min_v | maximum_independent_set(adj, set & !(adj[min_v] | 1 << min_v));
    }
    let with = 1 << max_v | maximum_independent_set(adj, set & !(adj[max_v] | 1 << max_v));
    let without = maximum_independent_set(adj, set & !(1 << max_v));
    if with.count_ones() >= without.count_ones() { with } else { without }
}

/// A value `χ̄` with `χ <= χ̄ <= ⌈(1 + ε) χ⌉` (Proposition 5): strip maximum
/// stable sets until at most `e^{-ε} n` vertices remain, then colour the rest
/// exactly with the polynomial-space algorithm.
pub fn approximate_chromatic_number(g: &Graph, epsilon: f64) -> usize {
    assert!(epsilon >= 0.0);
    let threshold = ((-epsilon).exp() * g.n() as f64).floor() as u32;
    let adj = g.adjacency();
    let mut rest = g.all();
    let mut stripped = 0;
    while rest.count_ones() > threshold {
        rest &= !maximum_independent_set(adj, rest);
        stripped += 1;
    }
    stripped + chromatic_number_polyspace(&g.induced(rest))
}
