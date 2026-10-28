//! A small simple undirected graph on at most 64 vertices, stored as one
//! adjacency bitmask per vertex. Exponential-time colouring algorithms are
//! only practical for a few dozen vertices, so a `u64` per row is plenty and
//! makes the set operations in the algorithms single machine instructions.

/// Maximum number of vertices a [`Graph`] can hold.
pub const MAX_VERTICES: usize = 64;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Graph {
    n: usize,
    adj: Vec<u64>,
}

impl Graph {
    /// The edgeless graph on `n` vertices.
    pub fn new(n: usize) -> Self {
        assert!(n <= MAX_VERTICES, "at most {MAX_VERTICES} vertices supported");
        Graph { n, adj: vec![0; n] }
    }

    pub fn from_edges(n: usize, edges: &[(usize, usize)]) -> Self {
        let mut g = Graph::new(n);
        for &(u, v) in edges {
            g.add_edge(u, v);
        }
        g
    }

    /// Adds the edge `uv`. Self-loops are rejected because a graph with a loop
    /// has no proper colouring.
    pub fn add_edge(&mut self, u: usize, v: usize) {
        assert!(u < self.n && v < self.n, "vertex out of range");
        assert_ne!(u, v, "self-loops are not allowed");
        self.adj[u] |= 1 << v;
        self.adj[v] |= 1 << u;
    }

    pub fn has_edge(&self, u: usize, v: usize) -> bool {
        self.adj[u] >> v & 1 == 1
    }

    pub fn n(&self) -> usize {
        self.n
    }

    /// Bitmask of the neighbours of `v`.
    pub fn neighbors(&self, v: usize) -> u64 {
        self.adj[v]
    }

    pub fn adjacency(&self) -> &[u64] {
        &self.adj
    }

    pub fn degree(&self, v: usize) -> usize {
        self.adj[v].count_ones() as usize
    }

    pub fn max_degree(&self) -> usize {
        (0..self.n).map(|v| self.degree(v)).max().unwrap_or(0)
    }

    pub fn edge_count(&self) -> usize {
        self.adj.iter().map(|a| a.count_ones() as usize).sum::<usize>() / 2
    }

    pub fn edges(&self) -> Vec<(usize, usize)> {
        let mut out = Vec::new();
        for u in 0..self.n {
            for v in iter_bits(self.adj[u] >> u >> 1) {
                out.push((u, u + 1 + v));
            }
        }
        out
    }

    /// Bitmask containing every vertex.
    pub fn all(&self) -> u64 {
        full_mask(self.n)
    }

    pub fn is_independent(&self, set: u64) -> bool {
        iter_bits(set).all(|v| self.adj[v] & set == 0)
    }

    /// Checks that `colours` assigns every vertex a colour and that no edge is
    /// monochromatic.
    pub fn is_proper_coloring(&self, colours: &[usize]) -> bool {
        colours.len() == self.n && self.edges().iter().all(|&(u, v)| colours[u] != colours[v])
    }

    /// The subgraph induced by the vertices in `set`, relabelled `0..|set|` in
    /// increasing order of their original index.
    pub fn induced(&self, set: u64) -> Graph {
        let verts: Vec<usize> = iter_bits(set).collect();
        let mut g = Graph::new(verts.len());
        for (i, &u) in verts.iter().enumerate() {
            for (j, &v) in verts.iter().enumerate().skip(i + 1) {
                if self.has_edge(u, v) {
                    g.add_edge(i, j);
                }
            }
        }
        g
    }

    pub fn complement(&self) -> Graph {
        let all = self.all();
        Graph { n: self.n, adj: (0..self.n).map(|v| all & !self.adj[v] & !(1 << v)).collect() }
    }

    // ----- standard families -------------------------------------------------

    pub fn complete(n: usize) -> Graph {
        let mut g = Graph::new(n);
        for u in 0..n {
            for v in u + 1..n {
                g.add_edge(u, v);
            }
        }
        g
    }

    pub fn path(n: usize) -> Graph {
        let edges: Vec<_> = (1..n).map(|v| (v - 1, v)).collect();
        Graph::from_edges(n, &edges)
    }

    pub fn cycle(n: usize) -> Graph {
        assert!(n >= 3, "a cycle needs at least 3 vertices");
        let mut g = Graph::path(n);
        g.add_edge(n - 1, 0);
        g
    }

    /// The wheel: a cycle on `n - 1` vertices plus a hub adjacent to all of them.
    pub fn wheel(n: usize) -> Graph {
        assert!(n >= 4, "a wheel needs at least 4 vertices");
        let mut g = Graph::new(n);
        for v in 1..n {
            g.add_edge(0, v);
            g.add_edge(v, if v + 1 == n { 1 } else { v + 1 });
        }
        g
    }

    pub fn complete_bipartite(a: usize, b: usize) -> Graph {
        let mut g = Graph::new(a + b);
        for u in 0..a {
            for v in a..a + b {
                g.add_edge(u, v);
            }
        }
        g
    }

    /// The Petersen graph, `GP(5, 2)`.
    pub fn petersen() -> Graph {
        Graph::generalized_petersen(5, 2)
    }

    /// The generalized Petersen graph `GP(n, k)` for `1 <= k < n / 2`: an outer
    /// cycle `u_0..u_{n-1}` (vertices `0..n`), spokes `u_i v_i`, and inner edges
    /// `v_i v_{i+k}` (vertices `n..2n`). It is cubic with `3n` edges.
    pub fn generalized_petersen(n: usize, k: usize) -> Graph {
        assert!(n >= 3 && k >= 1 && 2 * k < n, "GP(n, k) needs n >= 3 and 1 <= k < n / 2");
        let mut g = Graph::new(2 * n);
        for i in 0..n {
            g.add_edge(i, (i + 1) % n); // outer cycle
            g.add_edge(i, n + i); // spokes
            g.add_edge(n + i, n + (i + k) % n); // inner star polygon
        }
        g
    }

    /// Mycielski's construction: triangle-free whenever `self` is, with
    /// chromatic number one higher.
    pub fn mycielskian(&self) -> Graph {
        let n = self.n;
        let mut g = Graph::new(2 * n + 1);
        for (u, v) in self.edges() {
            g.add_edge(u, v);
            g.add_edge(u, n + v);
            g.add_edge(n + u, v);
        }
        for v in 0..n {
            g.add_edge(n + v, 2 * n);
        }
        g
    }

    /// The Mycielski graph `M_k` with chromatic number `k` (`M_2 = K_2`,
    /// `M_3 = C_5`, `M_4` = Grötzsch graph on 11 vertices, `M_5` on 23).
    pub fn mycielski(k: usize) -> Graph {
        assert!(k >= 2);
        let mut g = Graph::complete(2);
        for _ in 2..k {
            g = g.mycielskian();
        }
        g
    }

    /// The `r x c` queen graph: squares of a chessboard, adjacent when a queen
    /// could move between them.
    pub fn queen(r: usize, c: usize) -> Graph {
        let mut g = Graph::new(r * c);
        for a in 0..r * c {
            for b in a + 1..r * c {
                let (ra, ca) = ((a / c) as isize, (a % c) as isize);
                let (rb, cb) = ((b / c) as isize, (b % c) as isize);
                if ra == rb || ca == cb || (ra - rb).abs() == (ca - cb).abs() {
                    g.add_edge(a, b);
                }
            }
        }
        g
    }

    /// Erdős–Rényi `G(n, p)` from a deterministic generator seeded by `seed`.
    pub fn random(n: usize, p: f64, seed: u64) -> Graph {
        let mut rng = SplitMix64(seed);
        let mut g = Graph::new(n);
        for u in 0..n {
            for v in u + 1..n {
                if rng.next_f64() < p {
                    g.add_edge(u, v);
                }
            }
        }
        g
    }
}

/// Bitmask with the lowest `n` bits set.
pub fn full_mask(n: usize) -> u64 {
    if n == 64 { u64::MAX } else { (1u64 << n) - 1 }
}

/// Iterates over the indices of the set bits of `mask`, lowest first.
pub fn iter_bits(mut mask: u64) -> impl Iterator<Item = usize> {
    std::iter::from_fn(move || {
        if mask == 0 {
            None
        } else {
            let v = mask.trailing_zeros() as usize;
            mask &= mask - 1;
            Some(v)
        }
    })
}

/// Tiny deterministic PRNG so the crate needs no dependencies.
pub struct SplitMix64(pub u64);

impl SplitMix64 {
    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
}
