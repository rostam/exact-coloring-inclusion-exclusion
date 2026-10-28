//! Exact zero-testing of huge integers by evaluating them modulo several
//! coprime moduli.
//!
//! The inclusion–exclusion sums of the paper involve numbers like `s(X)^k`
//! with `n * k` bits. Following the idea behind Lemma 3 of the paper we never
//! form these numbers; we evaluate the whole sum modulo `2^64` (free wrapping
//! arithmetic) and then modulo primes just below `2^62`. A nonzero residue
//! proves the sum is nonzero. If the value is known to lie in `[0, B)` and it is
//! zero modulo moduli whose product exceeds `B`, the Chinese remainder theorem
//! proves it is exactly zero. So, unlike the randomised variant in the paper,
//! the answers here are deterministic and always correct.

use std::sync::OnceLock;

/// Arithmetic in `Z / mZ`. Elements are opaque `u64` representatives in which
/// zero is always represented by `0`.
pub trait Modulus: Copy {
    /// `log2` of a lower bound on the modulus.
    fn bits(&self) -> u32;
    /// Maps an ordinary integer into the internal representation.
    fn lift(&self, x: u64) -> u64;
    fn one(&self) -> u64;
    fn add(&self, a: u64, b: u64) -> u64;
    fn sub(&self, a: u64, b: u64) -> u64;
    fn mul(&self, a: u64, b: u64) -> u64;

    fn pow(&self, mut base: u64, mut exp: u64) -> u64 {
        let mut acc = self.one();
        while exp > 0 {
            if exp & 1 == 1 {
                acc = self.mul(acc, base);
            }
            base = self.mul(base, base);
            exp >>= 1;
        }
        acc
    }
}

/// Arithmetic modulo `2^64`.
#[derive(Clone, Copy, Debug)]
pub struct Wrap64;

impl Modulus for Wrap64 {
    fn bits(&self) -> u32 {
        64
    }
    fn lift(&self, x: u64) -> u64 {
        x
    }
    fn one(&self) -> u64 {
        1
    }
    #[inline]
    fn add(&self, a: u64, b: u64) -> u64 {
        a.wrapping_add(b)
    }
    #[inline]
    fn sub(&self, a: u64, b: u64) -> u64 {
        a.wrapping_sub(b)
    }
    #[inline]
    fn mul(&self, a: u64, b: u64) -> u64 {
        a.wrapping_mul(b)
    }
}

/// Montgomery arithmetic modulo an odd prime `p < 2^62`.
#[derive(Clone, Copy, Debug)]
pub struct Montgomery {
    p: u64,
    /// `-p^{-1} mod 2^64`.
    p_neg_inv: u64,
    /// `2^128 mod p`, used to enter Montgomery form.
    r2: u64,
}

impl Montgomery {
    pub fn new(p: u64) -> Self {
        assert!(p % 2 == 1 && p < 1 << 62);
        // Newton iteration for p^{-1} mod 2^64; each step doubles the correct bits.
        let mut inv = p;
        for _ in 0..6 {
            inv = inv.wrapping_mul(2u64.wrapping_sub(p.wrapping_mul(inv)));
        }
        let r = ((1u128 << 64) % p as u128) as u64;
        let r2 = ((r as u128 * r as u128) % p as u128) as u64;
        Montgomery { p, p_neg_inv: inv.wrapping_neg(), r2 }
    }

    #[inline]
    fn redc(&self, t: u128) -> u64 {
        let m = (t as u64).wrapping_mul(self.p_neg_inv);
        let t = ((t + m as u128 * self.p as u128) >> 64) as u64;
        if t >= self.p { t - self.p } else { t }
    }
}

impl Modulus for Montgomery {
    fn bits(&self) -> u32 {
        63 - self.p.leading_zeros()
    }
    fn lift(&self, x: u64) -> u64 {
        self.redc((x % self.p) as u128 * self.r2 as u128)
    }
    fn one(&self) -> u64 {
        self.lift(1)
    }
    #[inline]
    fn add(&self, a: u64, b: u64) -> u64 {
        let s = a + b;
        if s >= self.p { s - self.p } else { s }
    }
    #[inline]
    fn sub(&self, a: u64, b: u64) -> u64 {
        if a >= b { a - b } else { a + self.p - b }
    }
    #[inline]
    fn mul(&self, a: u64, b: u64) -> u64 {
        self.redc(a as u128 * b as u128)
    }
}

/// A computation that can be carried out over any [`Modulus`], producing a
/// vector of residues of integers that are known to lie in `[0, 2^bound_bits)`.
pub trait ModEval {
    fn eval<M: Modulus>(&self, m: &M) -> Vec<u64>;
}

/// Returns the smallest index whose integer value is nonzero, or `None` if all
/// values are zero. All values must lie in `[0, 2^bound_bits)`.
pub fn first_nonzero<E: ModEval>(e: &E, bound_bits: u32) -> Option<usize> {
    let mut best = first_nonzero_residue(&e.eval(&Wrap64));
    let mut covered = Wrap64.bits();
    for &p in primes() {
        if covered >= bound_bits || best == Some(0) {
            break;
        }
        let m = Montgomery::new(p);
        best = min_opt(best, first_nonzero_residue(&e.eval(&m)));
        covered += m.bits();
    }
    assert!(covered >= bound_bits || best == Some(0), "bound of {bound_bits} bits exceeds the prime table");
    best
}

/// Decides whether a single integer in `[0, 2^bound_bits)` is nonzero, stopping
/// at the first modulus that witnesses it.
pub fn is_nonzero<E: ModEval>(e: &E, bound_bits: u32) -> bool {
    if e.eval(&Wrap64)[0] != 0 {
        return true;
    }
    let mut covered = Wrap64.bits();
    for &p in primes() {
        if covered >= bound_bits {
            return false;
        }
        let m = Montgomery::new(p);
        if e.eval(&m)[0] != 0 {
            return true;
        }
        covered += m.bits();
    }
    panic!("bound of {bound_bits} bits exceeds the prime table");
}

fn first_nonzero_residue(v: &[u64]) -> Option<usize> {
    v.iter().position(|&x| x != 0)
}

fn min_opt(a: Option<usize>, b: Option<usize>) -> Option<usize> {
    match (a, b) {
        (Some(x), Some(y)) => Some(x.min(y)),
        (x, None) | (None, x) => x,
    }
}

/// The largest primes below `2^62`, in decreasing order. Together with `2^64`
/// they cover more than 9000 bits, enough for any graph this crate can handle.
pub fn primes() -> &'static [u64] {
    static PRIMES: OnceLock<Vec<u64>> = OnceLock::new();
    PRIMES.get_or_init(|| {
        let mut out = Vec::with_capacity(150);
        let mut c = (1u64 << 62) - 1;
        while out.len() < 150 {
            if is_prime(c) {
                out.push(c);
            }
            c -= 2;
        }
        out
    })
}

/// Deterministic Miller–Rabin for 64-bit integers.
pub fn is_prime(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    const BASES: [u64; 12] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
    for &b in &BASES {
        if n.is_multiple_of(b) {
            return n == b;
        }
    }
    let mulmod = |a: u64, b: u64| (a as u128 * b as u128 % n as u128) as u64;
    let powmod = |mut b: u64, mut e: u64| {
        let mut r = 1;
        while e > 0 {
            if e & 1 == 1 {
                r = mulmod(r, b);
            }
            b = mulmod(b, b);
            e >>= 1;
        }
        r
    };
    let s = (n - 1).trailing_zeros();
    let d = (n - 1) >> s;
    'witness: for &a in &BASES {
        let mut x = powmod(a, d);
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 1..s {
            x = mulmod(x, x);
            if x == n - 1 {
                continue 'witness;
            }
        }
        return false;
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn montgomery_matches_u128_arithmetic() {
        let p = primes()[3];
        let m = Montgomery::new(p);
        let mut rng = crate::graph::SplitMix64(7);
        for _ in 0..1000 {
            let (a, b) = (rng.next_u64(), rng.next_u64());
            let expect = (a as u128 % p as u128) * (b as u128 % p as u128) % p as u128;
            let got = m.mul(m.lift(a), m.lift(b));
            assert_eq!(got, m.lift(expect as u64));
            let e = rng.next_u64() % 100;
            let mut naive = m.one();
            for _ in 0..e {
                naive = m.mul(naive, m.lift(a));
            }
            assert_eq!(m.pow(m.lift(a), e), naive);
        }
    }

    #[test]
    fn prime_table_is_prime_and_large() {
        let ps = primes();
        assert!(ps.iter().all(|&p| is_prime(p) && p > 1 << 61));
        assert!(is_prime((1 << 61) - 1)); // Mersenne prime M61
        assert!(!is_prime(3_215_031_751)); // strong pseudoprime to bases 2, 3, 5, 7
    }
}
