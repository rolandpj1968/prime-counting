# prime-counting

An exploration of **combinatorial approaches to π(N)** — the prime-counting
function — in Zig. Empirical, one knob at a time, theory kept honest by measurement.

## The territory

π(N) counts the primes ≤ N (π(10) = 4, π(10⁹) = 50,847,534). Three families
compute it, in increasing sophistication and decreasing intuition:

- **Sieving / enumeration** — the Sieve of Eratosthenes and its refinements.
  O(N log log N): you find every prime. Feasible to ~10¹⁸–10²⁰.
- **Combinatorial** — Legendre → Meissel → Lehmer → Lagarias–Miller–Odlyzko →
  Deléglise–Rivat → Gourdon. Sub-linear (≈ O(N^(2/3))): counts π(N) *without
  listing the primes*, via a partial-sieve identity π(N) = φ(N,a) + a − 1 − P₂(N,a).
  This is how records reach 10³⁰. **Built here through the LMO/DR line and then
  Gourdon's decomposition** — π(10²²) in 4.74 h on six cores.
- **Analytic** — Lagarias–Odlyzko. O(N^(1/2+ε)) by contour-integrating ζ(s).
  Asymptotically best, but high-precision-arithmetic constants leave it
  practically dominated; its real use is independent verification.

Meissel — an *astronomer* — computed π(10⁹) by hand in 1885 (off by 56). The
combinatorial method is the intellectual object; the machine is a multiplier.

## Status

Complete through **Gourdon's decomposition**, parallel, with the LMO/Deléglise–Rivat
line beneath it as reference and correctness oracle, and the full straight-sieving
characterization beneath that. Everything below is exact and verified against every
known value 10 → 10²² (plus exhaustively over [0, 5000]). Presented most-advanced first.

All timings on one quiet laptop: **AVX2 (no AVX-512), L1d 32 KiB / L2 512 KiB /
L3 16 MiB, 28 GiB RAM** (6-core Ryzen 5 6600H), `zig 0.16 build-exe -O ReleaseFast
-mcpu=native`, six threads pinned one per physical core (SMT siblings left idle).

### Combinatorial π(x): Gourdon (`gourdon.zig`) — the frontier

π(x) = A − B + ω + φ₀ + Σ, with y = α·x^(1/3), z = x/y, and x* = max(x^(1/4), x/y²).

| x | π(x) | 6-core | RSS | ×prev |
|---|------|-----:|----:|----:|
| 10¹³ | 346,065,536,839 | 0.05 s | 2 MB | — |
| 10¹⁴ | 3,204,941,750,802 | 0.16 s | 3 MB | 3.37 |
| 10¹⁵ | 29,844,570,422,669 | 0.61 s | 4 MB | 3.80 |
| 10¹⁶ | 279,238,341,033,925 | 2.52 s | 7 MB | 4.13 |
| 10¹⁷ | 2,623,557,157,654,233 | 9.41 s | 15 MB | 3.73 |
| 10¹⁸ | 24,739,954,287,740,860 | 39.1 s | 36 MB | 4.15 |
| 10¹⁹ | 234,057,667,276,344,607 | 2.82 min | 90 MB | 4.32 |
| 10²⁰ | 2,220,819,602,560,918,840 | **12.3 min** | **215 MB** | 4.37 |
| 10²¹ | 21,127,269,486,018,731,928 | 50.1 min | 433 MB | 4.06 |
| 10²² | 201,467,286,689,315,906,290 | **3.79 h** | 985 MB | 4.54 |

(Timing provenance: cool-start runs; sustained multi-hour load costs up to ~9%
on this laptop — measured, not assumed. Profiles of both top rows show the
sieve's kill loop — `bt`/`btr` — as the hottest instructions, i.e. the machine
spends its time on the algorithm, not on stalls.)
**Memory now scales as O(x^(1/3))**, not O(√x): the largest resident structures
are the y-sized leaf table and prime list. 10²⁴ needs ~4 GB (was ~40 GB) — the
memory wall is gone, leaving runtime (~4 days at 10²⁴ on this laptop) as the
only constraint.

The ×prev column is the empirical growth per power of ten; theory (x^(2/3)) says
**4.64**. It sits at 4.12–4.25 through 10²⁰ and converges *up* to ~4.5 beyond — a
reminder that the plateau was a property of the measured range, not a permanent
state (extrapolating it made the 10²² estimate 21% optimistic). Two unmeasured
candidates for the drift: the π oracle stops fitting any cache (3.5 GB at 10²², and
its probes are random), and Amdahl — oracle build, prime sieve, μ/δ construction and
the cross-block reduction are all serial and all grow with x.

What earns it, beyond the decomposition itself:

- **ω and B fused on one counter.** A single segmented counter is folded to π(√z);
  ω's special leaves query it *during* the fold at their stage, and once fully
  folded it holds φ(·, π(√z)), so B's π(x/p) is read straight off it by Legendre
  (valid since z < y²). Gourdon folds π(√z) primes against LMO's π(y) — one pass
  serves both terms.
- **Segmented π, no resident √x structure at all.** Every π(v) above y is answered
  from a freshly sieved *window*: ω's C-leaves ride the sweep itself (the walk
  guards deliver each query inside the current segment — the sweep *is* the window
  pass), B's cursor walks a descending chunked window, and A's v ≥ y pairs are
  enumerated per v-window off the prime list. The only resident remnants are a
  y-capped oracle (a coprime-30 bitset: one 64-bit word covers exactly 240
  integers, so indexing is a shift) and `bpi` — π at segment boundaries, 8 bytes
  per 262,080 integers. Replacing resident random probes with L1/L2 windows also
  measured **5–9% faster**.
- **A fitted α(x)**, not a constant. α = clamp(−17.727 + 0.5980·ln x, 4, 24), from
  measured optima. Worth 18% at 10¹⁷ and 34% at 10¹⁸ over the hardcoded α = 4.
- **A bucket ring for sparse fold primes.** Once √z > segw, most fold primes hit a
  given segment at most once, and walking all π(√z) of them per segment becomes the
  dominant term — at 10²⁰ it had overtaken the kills themselves. Filing each under
  the segment its next multiple lands in is worth 19% there, and nothing at all
  below 10¹⁹.
- **Two-level counter, uncounted strikes, word-parallel strike.** The kill-heavy
  workload prefers a 2-level counter to a 3-level one; fold stages above π(x*) carry
  no leaves so they skip the bookkeeping entirely; and primes p ≲ 30 are struck a
  whole 64-bit word at a time.

### Combinatorial π(x): LMO / Deléglise–Rivat (`lmo.zig`)

The previous frontier, now the reference implementation and Gourdon's correctness
oracle. π(x) = φ(x,a) + a − 1 − P₂(x,a) with y = 4·x^(1/3), a = π(y), run at the
**2/3 exponent in Θ(x^(1/3)) memory** (times below are 1- and 4-core):

| x | π(x) | 1-core | 4-core | 4-core RSS |
|---|------|-----:|-----:|---------:|
| 10¹⁴ | 3,204,941,750,802 | 1.16 s | 0.38 s | 13 MB |
| 10¹⁵ | 29,844,570,422,669 | 5.23 s | 1.78 s | 27 MB |
| 10¹⁶ | 279,238,341,033,925 | 23.1 s | 7.01 s | 54 MB |
| 10¹⁷ | 2,623,557,157,654,233 | 1.72 min | 31.8 s | 110 MB |
| 10¹⁸ | 24,739,954,287,740,860 | 8.3 min | 2.55 min | 226 MB |
| 10¹⁹ | 234,057,667,276,344,607 | 42.4 min | 13.1 min | 467 MB |
| 10²⁰ | 2,220,819,602,560,918,840 | 3.44 h | 1.08 h | 964 MB |

Least-squares scaling exponent **0.658**, just under the theoretical 2/3. π(10¹⁹)
matches M. Deléglise's 1996 computation — reproduced here single-core in less
time than his HP-730 took for a value **four powers of ten smaller**. Sieving 10²⁰ would
take years; the lead over sieving *widens* with x, which is why records reach 10³⁰
combinatorially and never by enumeration. 10²⁰ (2,220,819,602,560,918,840) is the first value
past 2⁶⁴, reached via u128 with a ≤ few-% wide-arithmetic tax, matching the published value.
The 4-core column is a steady **~3.2×** (DRAM-bandwidth-bound — see Parallelism). The
*single-core* footprint stays Θ(x^(1/3)) — 2.7 MB at 10¹⁴ to ~230 MB at 10²⁰; the 4-core
RSS above is larger because each thread carries its own sweep scratch and the
cross-block reduction arrays scale with block count (O(nb·a)).

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

A full **4-core** sweep (pinned physical cores, k_over=8) holds a steady **~3.0–3.3×**
from 10¹⁴ to 10²⁰ — 4 cores is the bandwidth bang-for-buck point on this box, so it
is the default in the Status table. Peak RSS grows 13 MB → 964 MB across that range;
the dominant term is the O(nb·a) cross-block reduction arrays (block totals + μ-sums
per block). Reading Gourdon (`literature/gourdon.ps`) pinned the natural fix: those
arrays only ever carry non-zero corrections for prime index ≤ π(√(x/y)) (his *M*),
which is ~6.5× below a = π(y) — a ~2× total-memory cut, deferred until it gates.

### What the optimisation log actually taught

The full record is in [COMBINATORIAL.md](COMBINATORIAL.md); the negatives are the
useful part, and they share a cause.

**Every win came from moving less memory. Every miss was an instruction count
already hidden behind memory latency.** The π oracle, α(x) and the bucket ring all
reduce bytes touched, and all paid. Two changes that looked good on an instruction
count did not: hoisting the division out of the ω m-walk (63% of its iterations are
rejected, so removing their divisions *should* have been large — the whole dense
m-walk turned out to be 8% of runtime, and it delivered 1.3%), and batching the
counter's per-kill bookkeeping (the count array is **256 bytes** — permanently
L1-resident — so the update was never expensive, and adding a compare per kill to
avoid it measured **4.7% slower** and was reverted).

The sharpest single lesson cost a night at 10²²: the bucket ring was rebuilt as
a chain-linked arena (exact memory, no slack) and measured *time-neutral* at
every tested scale — then **+33% at 10²²**, where its ~3 MB/thread of chains ×6
threads crossed the 16 MB L3 exactly between 10²¹ (9.6 MB, fits) and 10²²
(18 MB, cliff). A live `perf` attach found **one chain-following load carrying
37.6% of all cycles**. The fix — the same packed 8-byte entries in contiguous
per-slot vectors, so drains stream instead of chase — took 10²² from 22,671 s to
**13,640 s, 20% under the pre-refactor record**. The law it stamped: *every
performance verdict is scale-indexed* — "neutral at tested scales" is a claim
about the tests, and cache boundaries are where it goes to die.

The same lens explains a third observation: **the parallel path is leaf-side
bandwidth-bound.** The fold works on a per-thread 32 KB L1-resident bitset, while
the leaf walk streams the *shared* leaf table (fused μ+δ), and threads on different blocks touch
different regions — so parallelism multiplies leaf traffic while fold cost stays
private. Hence α* is thread-count dependent (4.58 on six threads vs ~5.0 on one at
10¹⁶), and fold-side wins scale poorly: the uncounted-strike change measured 4%
single-threaded and ~0.7% in the six-thread ladder.

### Next

- ~~Pack μ and δ into one u16 array~~ **Done** — fused as `leaf[m]` = μ-sign bit +
  15-bit *spf index* (not value: the leaf test is an order comparison between
  primes, so π(spf) > bi+1 replaces spf > p, and indices compress ~ln p better).
  Measured −6.7% at 10¹⁷ / −6.0% at 10¹⁸ on 6 threads, neutral serial — exactly
  the parallel-specific shape the bandwidth story predicts. Ceiling moved
  10²⁵ → ~10²⁹. Next rung if the leaf stream still binds: a u8 first-level table
  (~92% of entries have spf ≤ p₁₂₆) with a u16 overflow, ~1.15y bytes.
- **`chooseAlpha(x, nthreads)`** — α* demonstrably depends on both; the current fit
  samples one slice. Needs a 2-D sweep.
- **Phase timing at scale** to settle whether the growth-ratio drift past 10²⁰ is
  the oracle leaving cache or the serial phases (Amdahl). `piGourdonV` already laps
  each phase behind a `verbose` flag.
- **A `build.zig`** with a bench step. The CLI (`src/pi.zig`) and the `Config`
  refactor behind it are done — see Build & run.

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

Zig **0.16**, hand-driven (no `build.zig` yet). Start here:

```
zig build-exe -O ReleaseFast -mcpu=native src/pi.zig -femit-bin=./pi
```

`pi` is one binary for every implementation, with the tuning knobs exposed as
options rather than recompiled — in the spirit of primecount's CLI:

```
./pi 1e20                     # fastest algorithm, defaults
./pi 1e20 -t 0 --pin          # one thread per physical core, pinned
./pi 1e18 -a lmo -t 6         # pick algorithm and thread count
./pi 1e17 --alpha 6.5 -v      # override the fitted α, per-phase timing
./pi 1e16 --check             # verify against the known π(10ⁿ) table
```

`--help` lists the rest. x accepts `1e20`, `10^20`, `1_000_000` or plain digits
(`e` is scientific, `^` is exponentiation; both are built by repeated
multiplication, never via f64, so they stay exact at the top of the u128 range).
Options that do not apply to the chosen algorithm say so rather than being silently
dropped, and `--check` exits non-zero on mismatch, so it scripts.

### Experiment drivers

Alongside the CLI, each investigation has a standing driver — built the same way
(`zig build-exe -O ReleaseFast -mcpu=native src/<file> -femit-bin=./<name>`), and
`-O ReleaseSafe` for a correctness pass with asserts and bounds live:

| file | what it does |
|---|---|
| `gourdon.zig` | is itself runnable — the full correctness + benchmark suite |
| `gsweep.zig` | π(10ⁿ) ladder, n = 13…20, parallel, checked against known values |
| `gasweep.zig` | α sweep — times y = α·x^(1/3) over a grid of α and x |
| `g2122.zig` | the 10²¹ / 10²² runs |
| `main.zig` | scratch driver, swapped per investigation (sieve benchmarks) |

**The project is still in flux**, so the tooling — a proper `build.zig` with a bench
step — is itself on the list.

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
