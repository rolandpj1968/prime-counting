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
  how records reach 10³⁰. **Built here through the LMO/DR line** (π(10¹⁹) in ~42
  min, single-threaded).
- **Analytic** — Lagarias–Odlyzko. O(N^(1/2+ε)) by contour-integrating ζ(s).
  Asymptotically best, but high-precision-arithmetic constants leave it
  practically dominated; its real use is independent verification.

Meissel — an *astronomer* — computed π(10⁹) by hand in 1885 (off by 56). The
combinatorial method is the intellectual object; the machine is a multiplier.

## Status

Complete, single-threaded, through the LMO/Deléglise–Rivat combinatorial line and
the full straight-sieving characterization beneath it. Everything below is exact
and verified against every known value 10 → 10¹⁹ (plus exhaustively over [0, 5000]).
Presented most-advanced first.

All timings on one quiet mini-PC: **AVX2 (no AVX-512), L1d 32 KiB / L2 512 KiB /
L3 16 MiB per core, 28 GiB RAM**, `zig 0.16 build-exe -O ReleaseFast -mcpu=native`,
single-threaded.

### Combinatorial π(x): LMO / Deléglise–Rivat (`lmo.zig`)

The current frontier. π(x) = φ(x,a) + a − 1 − P₂(x,a) with y = 4·x^(1/3), a = π(y),
run at the **2/3 exponent in Θ(x^(1/3)) memory**:

| x | π(x) | time | peak RSS |
|---|------|-----:|---------:|
| 10¹⁴ | 3,204,941,750,802 | 1.16 s | 2.7 MB |
| 10¹⁵ | 29,844,570,422,669 | 5.23 s | ~5 MB |
| 10¹⁶ | 279,238,341,033,925 | 23.1 s | 11 MB |
| 10¹⁷ | 2,623,557,157,654,233 | 1.72 min | ~24 MB |
| 10¹⁸ | 24,739,954,287,740,860 | 8.3 min | ~52 MB |
| 10¹⁹ | 234,057,667,276,344,607 | 42.4 min | ~109 MB |
| 10²⁰ | 2,220,819,602,560,918,840 | 3.44 h | ~230 MB |

Least-squares scaling exponent **0.658**, just under the theoretical 2/3. π(10¹⁹)
matches M. Deléglise's 1996 computation — reproduced here single-threaded in less
time than his HP-730 took for a value **four powers of ten smaller**. Sieving 10²⁰ would
take years; the lead over sieving *widens* with x, which is why records reach 10³⁰
combinatorially and never by enumeration. 10²⁰ (2,220,819,602,560,918,840) is the first value
past 2⁶⁴, reached via u128 with a ≤ few-% wide-arithmetic tax, matching the published value.

How it earns the exponent and the footprint:

- **Special-leaf sieve, no table.** φ splits into ordinary leaves (direct) and
  *special* leaves, resolved by a segmented sieve over [1, x/y] with a **3-level
  O(1)-kill alive-counter** — the O(x^(2/3)) prefix-π table is never stored. This
  counter replaces LMO's O(log x) Fenwick-style tree, whose updates are their
  *dominant* cost, so we beat their leading term by a whole log factor.
- **Closed-form leaf classes.** ~52% of leaves (p³ ≥ x) collapse to a single
  binomial C(n₁, 2); a further large class is read from a π-table over [1, y]. Only
  the genuinely hard leaves ever touch the counter.
- **mod-30 wheel fold** and P₂ **fused** onto the same counter (its primes sieved
  on the fly from x^(1/4) base primes, never stored — the last Θ(√x) term gone).

The build was driven **empirically against the two source papers** ([COMBINATORIAL.md](COMBINATORIAL.md)
has the full log). The most interesting results were the negatives: both of
Deléglise–Rivat's headline optimisations — the x^(1/4) sieve bound (§7) and leaf
*clustering* (§6.5) — measured **net-negative here**, because each amortises a cost
(an O(log x) tree; per-leaf evaluation) that the O(1)-kill counter and the π-table
had already removed. We sit at ~1.6× DR's 1996 implementation while using *neither*
of their log-factor tricks: they weren't skipped, they were made unnecessary.

### Meissel–Lehmer baseline (`meissel.zig`)

The classical recursion for φ, kept as the readable reference and as LMO's
correctness oracle (`phiOfXY`). A compact bit-sieve + per-word checkpoints keeps
the π-table cache-resident (memory *and* speed), but the algorithm is unchanged, so
it runs at the historical **~0.8 exponent** and hits an **O(x^(2/3)) memory wall**:
10¹⁴ = 3,204,941,750,802 in 41 s, ~0.5 GB — where LMO does the same in ~1.2 s and
2.7 MB, and reaches four powers of ten further.

### Straight sieving (`sieve.zig` + friends)

The foundation everything rests on — the combinatorial methods still need it for
base primes and P₂ — and where the *machine* story lives (cache, wheels, word
size). One generic sieve, three orthogonal comptime axes plus a runtime range:

```
Sieve( Wheel,      // coordinate map: {} / {2} / {2,3} / {2,3,5} → all/odds/mod-6/mod-30
       Store,      // backing bits:   flat []bool / bit-packed []uN / std.DynamicBitSet
       seg_bytes ) // segment size = cache-residency knob (≥N ⇒ whole-array)
```

plus `BucketSieve` (two-tier for large N) and `rangesieve` (u128 range counting
past 2⁶⁴). Fully characterized — the two headline empirical results are below.

### Parallelism

`piLMOPar` runs the sweep across pinned cores: [1,z] is block-decomposed, each block
swept with a local running φ and stitched by a cross-block reduction, with a
**cost-balanced partition** (blocks of equal width+leaves, since 99.7% of leaves sit
in the bottom of [1,z] and a naïve equal-width split makes one block 36% of the run).
On the 6-core Ryzen 5 6600H it reaches **~4.7× at 10¹⁴**, falling to ~3.8× at 10¹⁶ —
the ceiling is DRAM bandwidth, not cores. Pinned experiments show HT gives a real
~1.3× per-core gain but adds nothing once the cores saturate memory bandwidth.

### Next

Gourdon's algorithm (a better decomposition than DR, ~2.3×) and his per-block
transfer-size reduction (to push past the bandwidth wall); the 3.8× implementation
gap to primecount-DR (fast division).

## Empirical highlights

Full detail and tables in [SIEVING.md](SIEVING.md). Two sieving results stand out.

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

## Build & run

Zig **0.16**, hand-driven (no `build.zig` yet):

```
zig build-exe -O ReleaseFast -mcpu=native src/main.zig -femit-bin=./sieve && ./sieve
zig build-exe -O ReleaseSafe ...        # correctness pass: asserts + bounds live
```

`src/main.zig` is the current experiment driver, swapped per investigation — **the
project is still in flux**, so the tooling (a proper `build.zig` + bench step, a
stable entry point) is itself on the list. More coming.

## References

Prime-counting methods:
- A. M. Legendre, *Essai sur la théorie des nombres* (1808) — the Φ(x,a) sieve formula.
- E. Meissel (1870–1885) — the first practical combinatorial method; π(10⁹) by hand.
- D. H. Lehmer, "On the exact number of primes less than a given limit," *Illinois J. Math.* 3 (1959).
- J. C. Lagarias, V. S. Miller, A. M. Odlyzko, "Computing π(x): the Meissel–Lehmer method," *Math. Comp.* 44 (1985).
- M. Deléglise, J. Rivat, "Computing π(x): the Meissel, Lehmer, Lagarias, Miller, Odlyzko method," *Math. Comp.* 65 (1996).
- T. Oliveira e Silva, "Computing π(x): the combinatorial method," *Revista do DETUA* 4:6 (2006) — the self-contained, C-code-level treatment of LMO/DR; independent π(10ⁿ) tables.
- X. Gourdon, "Computation of π(x): improvements to the … method" (2001).
- J. C. Lagarias, A. M. Odlyzko, "Computing π(x): an analytic method," *J. Algorithms* 8 (1987).

Foundations:
- B. Riemann, "Über die Anzahl der Primzahlen unter einer gegebenen Größe" (1859).
- F. Mertens, "Ein Beitrag zur analytischen Zahlentheorie," *J. reine angew. Math.* 78 (1874).
- Prime Number Theorem — Hadamard and de la Vallée Poussin (1896).

Tools & systems:
- T. Oliveira e Silva, [Fast implementation of the segmented sieve of Eratosthenes](http://sweet.ua.pt/tos/) — the **bucket sieve** for large primes; also prime-gap tables and Goldbach verified to 4×10¹⁸.
- K. Walisch, [primesieve](https://github.com/kimwalisch/primesieve) & [primecount](https://github.com/kimwalisch/primecount) — the modern reference sieve / combinatorial π(x).
- R. Bryant, D. O'Hallaron, *Computer Systems: A Programmer's Perspective* — the "memory mountain."
- [OEIS A006880](https://oeis.org/A006880) — π(10ⁿ). [Zig](https://ziglang.org) 0.16.
