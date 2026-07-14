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

Straight sieving is characterized end to end. Next: the combinatorial π(x)
(Meissel–Lehmer), the algorithm this repo was always heading toward.
