# exact-coloring-inclusion-exclusion

A Rust implementation of

> Andreas Björklund and Thore Husfeldt, *Inclusion–Exclusion Based Algorithms for Graph Colouring*, ECCC TR06-044 (2006).

The paper computes the chromatic number in time `2^n n^O(1)`. It counts the ways `c_k(G)` to **cover** the vertex set with `k` stable sets:

```
c_k(G) = Σ_{X ⊆ V} (-1)^{|X|} s(X)^k          (Lemma 2)
χ(G)   = min { k : c_k(G) > 0 }                (Lemma 1)
```

Here `s(X)` is the number of nonempty stable sets avoiding `X`. A single `2^n` table of `s(X)`, filled by the recursion `s(X) = s(X ∪ {v}) + s(X ∪ N[v]) + 1`, is all that is needed.

## What is implemented

| Paper | Function | Notes |
|---|---|---|
| Prop. 1: χ in time and space `2^n n^O(1)` | `ie::chromatic_number`, `ie::is_k_colorable` | Binary search on `k`, bracketed by the DSATUR heuristic. `ie::chromatic_number_unbracketed` searches `[1, n]` as in the paper. |
| Prop. 2: χ in polynomial space | `ie::chromatic_number_polyspace` | `s(X)` is recomputed by a branching stable-set counter, with closed forms for paths and cycles. One pass accumulates every candidate `k`. |
| Lemma 4 / Prop. 3: number of colourings | `ie::stable_partition_counts`, `ie::count_colorings`, `ie::chromatic_polynomial` | `π_G(r)` from the `A(l, m, X)` DP, then `p_G(k) = Σ π_G(r) (k)_r`. |
| Prop. 4: chromatic sum | `ie::chromatic_sum` | Weight enumerator of covering colour-class tuples with `Δ+1` colours. |
| §2.5: constructing an optimal colouring | `ie::optimal_coloring` | Adds edges and contracts vertices, using `O(log n)` decision calls per vertex. |
| Prop. 5: approximation | `ie::approximate_chromatic_number` | Strips maximum stable sets, then runs Prop. 2 on the rest. Guarantees `χ ≤ χ̄ ≤ ⌈(1+ε)χ⌉`. |

**Exact arithmetic.** The paper's Lemma 3 suggests working modulo random prime powers, which gives a randomised algorithm. This crate keeps the modular speedup but is **deterministic**. It evaluates `c_k` modulo `2^64` and modulo primes just below `2^62` (Montgomery multiplication), see `src/modular.rs`:

- A nonzero residue proves `c_k ≠ 0`.
- Since `0 ≤ c_k ≤ s(∅)^k`, residues of zero modulo enough coprime moduli prove `c_k = 0` by the CRT.

Usually only the first modulus is needed.

**Parallelism.** The `2^n`-term inclusion–exclusion sums (Props. 1–4) and the stable-set tables run in parallel with [rayon](https://crates.io/crates/rayon). Each sum is split into chunks of 4096 subsets. The chunks' partial sums are combined in modular or wrapping arithmetic, which is associative and commutative, so the results do not depend on scheduling. The tables use the recursion on the *highest* vertex of `Y`, so the block of subsets with highest vertex `j` depends only on subsets below `2^j` and is filled in parallel. The tests check that the parallel results are identical to the sequential path. Every entry point has a `*_with(g, Exec::Sequential | Exec::Parallel)` variant, and the plain functions use `Exec::Parallel`.

**Comparison algorithms** (`src/baselines.rs`):

- first-fit greedy;
- largest-first (Welsh–Powell);
- the DSATUR heuristic (Brélaz 1979);
- exact DSATUR branch and bound;
- Lawler's `O(2.4423^n)` maximal-independent-set DP (1976), the first row of the paper's Table 1;
- brute-force backtracking, used in tests.

Graphs are `Graph` (`src/graph.rs`), one `u64` adjacency mask per vertex (≤ 64 vertices). It includes generators for complete graphs, paths, cycles, wheels, complete bipartite graphs, the Petersen graph, Mycielski graphs, queen graphs and `G(n, p)`. The table-based algorithms accept up to 30 vertices (memory `4·2^n` bytes). The counting algorithms accept up to 25.

## Usage

```rust
use ie_coloring::{Graph, ie};

let g = Graph::mycielski(4);                 // Grötzsch graph: triangle-free, χ = 4
assert_eq!(ie::chromatic_number(&g), 4);
let colouring = ie::optimal_coloring(&g);    // Vec<usize>, colours 0..4
assert!(g.is_proper_coloring(&colouring));
let p = ie::chromatic_polynomial(&g);        // coefficients of k^0..k^n
let p4 = ie::count_colorings(&g, 4);         // number of proper 4-colourings
```

## Tests

```
cargo test        # test profile is optimised (opt-level 3)
```

The tests check:

- known chromatic numbers: complete graphs, cycles, wheels, Petersen, Mycielski M2–M5, queen graphs;
- all five exact algorithms against brute force on 117 random graphs (`n ≤ 13`, densities 0.1–0.9);
- Lemma 2 against explicit enumeration of covering `k`-tuples;
- chromatic polynomials against deletion–contraction and the closed forms for `K_n`, trees and cycles;
- `π(E_n, r)` against Stirling numbers;
- colouring counts against enumeration;
- chromatic sums against brute force;
- optimality and properness of `optimal_coloring`;
- the approximation guarantee.

## Benchmarks and comparison

```
cargo bench                                  # Criterion, see benches/coloring.rs
cargo run --release --example compare [max_n]  # Markdown table of times and colour counts
```

### Results

These are from one machine (12 cores, all algorithms single-threaded), Criterion medians.

**Exact χ on `G(n, 1/2)`:**

| n | IE (Prop. 1) | IE, unbracketed | IE polyspace (Prop. 2) | Lawler | DSATUR B&B |
|---:|---:|---:|---:|---:|---:|
| 10 | 5.7 µs | 8.2 µs | 49 µs | 46 µs | 1.5 µs |
| 14 | 94 µs | 176 µs | 1.76 ms | 1.85 ms | 2.0 µs |
| 18 | 1.47 ms | 2.87 ms | 47 ms | 51 ms | 2.6 µs |
| 20 | 5.94 ms | 14.9 ms | 229 ms | 246 ms | 4.6 µs |
| 22 | 23.8 ms | 45.9 ms | – | – | 3.8 µs |

**Density sweep, n = 18:**

| p | IE | IE polyspace | Lawler | DSATUR B&B |
|---:|---:|---:|---:|---:|
| 0.1 | 0.74 ms | 19.8 ms | 46.8 ms | 2.4 µs |
| 0.5 | 1.47 ms | 47.3 ms | 52.8 ms | 2.6 µs |
| 0.9 | 2.44 ms | 33.0 ms | 20.3 ms | 3.7 µs |

**Structured graphs:**

| graph | IE | Lawler | DSATUR B&B |
|---|---:|---:|---:|
| Queen 4×5 (n=20, χ=5) | 5.9 ms | 234 ms | 3.4 µs |
| Mycielski M5 (n=23, χ=5) | 56.7 ms | – | 83 µs |
| Queen 5×5 (n=25, χ=5) | 229 ms | – | 4.9 µs |

**Counting:**

- `chromatic_polynomial` for n = 10 / 14 / 18: 0.16 / 4.8 / 152 ms.
- `chromatic_sum` for n = 10 / 14: 0.67 / 94 ms.

**Heuristics vs. the exact value.** This is from `examples/compare`, `G(n, p)` with `n ≤ 24`. DSATUR found χ on 23 of 24 random instances. Largest-first found it on 18 and plain greedy on 9. Greedy was off by as much as 2.

**Parallel vs. sequential** (12 cores, `G(n, 1/2)`):

| | n | sequential | parallel | speedup |
|---|---:|---:|---:|---:|
| Prop. 1 sum only | 26 | 436 ms | 87 ms | 5.0× |
| `i(Y)` table | 26 | 100 ms | 19 ms | 5.2× |
| Prop. 1 end to end | 26 | 550 ms | 132 ms | 4.2× |
| Prop. 2 (polyspace) | 20 | 230 ms | 39 ms | 6.0× |
| Prop. 3 (`stable_partition_counts`) | 18 | 172 ms | 40 ms | 4.3× |
| Prop. 4 (`chromatic_sum`) | 14 | 58 ms | 30 ms | 1.9× |

For n ≤ 14 the thread overhead outweighs the gain. The single-threaded tables above were measured before parallelisation and correspond to `Exec::Sequential`.

### Takeaways

- **IE scales as predicted.** Prop. 1 grows by about 4× per two extra vertices, i.e. `2^n`. Lawler grows by about 5× (theory: `2.44^n`). The polynomial-space variant is close to Lawler and grows a little faster than the exponential-space IE.
- **IE beats the classical exact DP.** At n = 20 it is roughly 40× faster than Lawler, and the gap widens with `n`.
- **Bracketing helps.** Using DSATUR as the upper bound halves the binary-search cost compared with searching all of `[1, n]` as the paper does.
- **IE's cost is almost independent of structure.** It only grows mildly with χ, because proving `c_k = 0` for larger `k` needs more moduli. Lawler instead gets faster on dense graphs, where there are fewer maximal independent sets.
- **DSATUR branch and bound is far faster on every instance here.** Its clique and heuristic bounds close the gap almost immediately at these sizes. It agreed with IE on all test graphs and on 36 extra random graphs with 24 vertices. Its worst case is still exponential with no `2^n` guarantee: G(60, 1/2) takes about 150 ms and Mycielski M6 (n = 47) about 90 ms, and it degrades sharply beyond that.
- **IE's strengths are elsewhere.** Its running time is predictable and independent of structure, it gives a provable `O*(2^n)` bound, and the same machinery yields exact counts: chromatic polynomial, partitions into stable sets, and chromatic sum. Branch and bound cannot produce those counts.
