//! Exact graph colouring by inclusion–exclusion, after Björklund and Husfeldt,
//! "Inclusion–Exclusion Based Algorithms for Graph Colouring", ECCC TR06-044.
//!
//! * [`ie::chromatic_number`] – Proposition 1, time and space `2^n n^O(1)`.
//! * [`ie::chromatic_number_polyspace`] – Proposition 2, polynomial space.
//! * [`ie::stable_partition_counts`], [`ie::count_colorings`],
//!   [`ie::chromatic_polynomial`] – Proposition 3.
//! * [`ie::chromatic_sum`] – Proposition 4.
//! * [`ie::optimal_coloring`] – Section 2.5.
//! * [`ie::approximate_chromatic_number`] – Proposition 5.
//!
//! [`baselines`] holds classical algorithms for comparison.

pub mod baselines;
pub mod graph;
pub mod ie;
pub mod modular;

pub use graph::Graph;
