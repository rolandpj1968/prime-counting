# prime-counting

An exploration of **combinatorial approaches to π(N)** — the prime-counting
function — in Zig. Empirical, one knob at a time, theory kept honest by measurement.

## The territory

π(N) counts the primes ≤ N (π(10) = 4, π(10⁹) = 50,847,534). Three families
compute it, in increasing sophistication and decreasing intuition:

- **Sieving / enumeration** — the Sieve of Eratosthenes and its refinements.
  O(N log log N): you find every prime. Feasible to ~10¹⁸–10²⁰.
- **Combinatorial** — Legendre → Meissel → Lehmer → Lagarias–Miller–Odlyzko →
  Deléglise–Rivat. Sub-linear (≈ O(N^(2/3))): counts π(N) *without listing the
  primes*, via a partial-sieve identity π(N) = φ(N,a) + a − 1 − P₂(N,a). This is
  how records reach 10³⁰. **The destination of this project.**
- **Analytic** — Lagarias–Odlyzko. O(N^(1/2+ε)) by contour-integrating ζ(s).
  Asymptotically best, but high-precision-arithmetic constants leave it
  practically dominated; its real use is independent verification.

Meissel — an *astronomer* — computed π(10⁹) by hand in 1885 (off by 56). The
combinatorial method is the intellectual object; the machine is a multiplier.

## Why sieving first

Sieving is the **fundamental operation** everything rests on. Even the
combinatorial methods can't escape it: they need a sieve for the base primes and
for the P₂ term. And sieving is where the *machine* story lives — cache, wheels,
SIMD, word size — which is exactly the part you have to earn empirically rather
than derive. So we characterize sieving deeply before building the combinatorial
layer on top of it.

Everything is one generic sieve, composed from three orthogonal comptime axes
plus a runtime range:

```
Sieve( Wheel,      // coordinate map: {} / {2} / {2,3} / {2,3,5} → all/odds/mod-6/mod-30
       Store,      // backing bits:   flat []bool / bit-packed []uN / std.DynamicBitSet
       seg_bytes ) // segment size = cache-residency knob (≥N ⇒ whole-array)
```

plus `BucketSieve` (two-tier for large N) and `rangesieve` (u128 range counting
past 2⁶⁴).

## Empirical highlights

Full detail and tables in [RESULTS.md](RESULTS.md). Two results stand out.

### Segmentation is about cache, and the cache is legible

Sweeping the segment size at fixed N reads the memory hierarchy straight off the
throughput curve: a peak at the **L1d-sized segment (32 KiB)**, a gentle L1→L2→L3
slope (the prefetcher hides the on-chip latency steps for a strided strike
pattern), and one sharp cliff at **L3→DRAM with the knee exactly at 16 MiB = L3**.
Segmentation makes the sieve *cache-immune to N* — a 10¹⁸ run behaves like a 10⁹
one — buying back ~1.5–3× that a whole-array sieve loses to DRAM, growing with N.

### Wheeling matches Mertens theory, and peaks at mod-30

Strike work ≈ N·Σ_{p≤√N} 1/p, thinned by the wheel density φ(M)/M. Measured
speedups track the theory to two decimals:

| wheel | density | measured step | predicted (Mertens) |
|-------|---------|:-------------:|:-------------------:|
| all → odds | 1 → ½ | 2.39× | 2.48× |
| odds → mod-6 | ½ → ⅓ | 1.80× | 1.78× |
| mod-6 → mod-30 | ⅓ → 4/15 | 1.41× | 1.41× |

Beyond mod-30 it *regresses*: the wheel's own per-prime delta tables spill L2, and
the marginal Mertens gain (p/(p−1) → 1) is too small to pay for it. The champion —
bit-packed `[]u64` + mod-30, ~3.6 G ints/s — is exactly where real sieves sit,
now understood from first principles.

The engine also turns around to *measure the number theory it rests on*: π(N)/N,
the average gap ≈ ln N − 1, and Li(N) tracking π to parts-per-million on the
√N/ln N (Riemann) error scale. And it counts primes **past 2⁶⁴** (u128), validated
not against any table but against the prime number theorem itself.

### Combinatorial π(x) leaves sieving behind

`meissel.zig` implements the Meissel–Lehmer identity π(x) = φ(x,a) + a − 1 −
P₂(x,a), a = π(x^⅓) — counting primes *without enumerating them*. It's
**sub-linear** (empirical exponent ~0.8; sieving is 1.0), so the lead over
sieving *widens* with x:

| x | π(x) | Meissel | vs sieving |
|---|------|--------:|:----------:|
| 10¹² | 37,607,912,018 | 1.18 s | ~236× |
| 10¹³ | 346,065,536,839 | 7.25 s | ~383× |

Sieving 10¹³ would take ~46 minutes; Meissel does it in 7 seconds. This is why
records reach 10³⁰ combinatorially and never by sieving. Meissel computed π(10⁹)
*by hand* in 1885 (off by 56); the same identity is exact here in 4 ms.

## Build & run

Zig **0.16**. No `build.zig` yet (see open tasks):

```
zig build-exe -O ReleaseFast -mcpu=native src/main.zig -femit-bin=./sieve && ./sieve
zig build-exe -O ReleaseSafe ...        # correctness pass: asserts + bounds live
zig test src/stores/bit_packed.zig      # unit tests
```

`src/main.zig` is the current experiment driver (swapped per investigation);
`src/sweep.zig` holds the benchmark harnesses.

## Status

Straight sieving is characterized end to end, and the combinatorial π(x)
(Meissel–Lehmer) works to 10¹³. Next: LMO's special-leaf sieving to remove the
O(x^(2/3)) prefix-π table (the current memory wall), Deléglise–Rivat's log-factor
refinement, and u128 to push past 10¹⁸.

## References

Prime-counting methods:
- A. M. Legendre, *Essai sur la théorie des nombres* (1808) — the Φ(x,a) sieve formula.
- E. Meissel (1870–1885) — the first practical combinatorial method; π(10⁹) by hand.
- D. H. Lehmer, "On the exact number of primes less than a given limit," *Illinois J. Math.* 3 (1959).
- J. C. Lagarias, V. S. Miller, A. M. Odlyzko, "Computing π(x): the Meissel–Lehmer method," *Math. Comp.* 44 (1985).
- M. Deléglise, J. Rivat, "Computing π(x): the Meissel, Lehmer, Lagarias, Miller, Odlyzko method," *Math. Comp.* 65 (1996).
- X. Gourdon, "Computation of π(x): improvements to the … method" (2001).
- J. C. Lagarias, A. M. Odlyzko, "Computing π(x): an analytic method," *J. Algorithms* 8 (1987).

Foundations:
- B. Riemann, "Über die Anzahl der Primzahlen unter einer gegebenen Größe" (1859).
- F. Mertens, "Ein Beitrag zur analytischen Zahlentheorie," *J. reine angew. Math.* 78 (1874).
- Prime Number Theorem — Hadamard and de la Vallée Poussin (1896).

Tools & systems:
- K. Walisch, [primesieve](https://github.com/kimwalisch/primesieve) & [primecount](https://github.com/kimwalisch/primecount) — the modern reference sieve / combinatorial π(x).
- R. Bryant, D. O'Hallaron, *Computer Systems: A Programmer's Perspective* — the "memory mountain."
- [OEIS A006880](https://oeis.org/A006880) — π(10ⁿ). [Zig](https://ziglang.org) 0.16.
