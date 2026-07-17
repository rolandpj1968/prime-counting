# Empirical results

Detailed findings for the sieve and the combinatorial π(x). See [README.md](README.md)
for framing and references; open/future work is tracked as tasks, not here.

Default benchmark: **N = 10⁹**, π(N) = 50,847,534. Best-of-3, sieve timed
separately from count. Machine: AVX2 (no AVX-512), L1d 32 KiB / L2 512 KiB /
L3 16 MiB per core, 28 GiB RAM. Build: `zig build-exe -O ReleaseFast -mcpu=native`.

Everything is now one generic sieve, `Sieve(Wheel, Store, seg_bytes)`:
- **Wheel(wheel_primes)** — coordinate map, comptime: `{}`=all, `{2}`=odds,
  `{2,3}`=mod-6, `{2,3,5}`=mod-30. Density φ(M)/M.
- **Store** — backing bits: flat `[]bool` (1 byte/flag) / `[]u64` / DynamicBitSet.
- **seg_bytes** — segment store size (cache-critical knob); ≥N ⇒ whole-array.

## The two big findings

### 1. Representation optimum flips with strike-boundedness
- **Whole-array (memory-bound):** 954 MiB `[]bool` → DRAM → footprint wins → bit-packing (148 → 202 M/s).
- **Segmented, un-wheeled (strike-bound):** both fit cache, so per-strike cost rules. A byte store beats a bit read-modify-write → `[]bool` ~2.4× faster.
- **Segmented + mod-30 (overhead-bound):** heavy wheeling slashes the strike count, so you stop being strike-bound; now `[]u64`'s 8× fewer segments + tighter footprint win → `[]u64` beats `[]bool` (3489 vs 2449). **The wheel is the knob that moves you across the line.** Champion = bit-packed + mod-30, exactly what real sieves converge on.

### 2. The wheel speedups match Mertens theory
Strike work ≈ N·Σ_{p≤√N} 1/p, thinned by density φ(M)/M. Measured (`[]u64`, 32 KiB seg):

| wheel | density | rate (M ints/s) | measured step | predicted step |
|-------|---------|-----------------|---------------|----------------|
| M=1 (all) | 1 | 577 | — | — |
| M=2 (odds) | 1/2 | 1381 | 2.39× | 2.48× |
| M=6 (mod-6) | 1/3 | 2482 | 1.80× | 1.78× |
| M=30 (mod-30) | 4/15 | 3489 | 1.41× | 1.41× |

Cumulative all→mod-30: **6.05×** (predicted ~6.2×). The marginal wheel gain is
p/(p−1), maximal at 2, dying toward 1 (Mertens' 3rd theorem: full-wheel density
~ e^(−γ)/ln x) — which is why real sieves stop at mod-30 / mod-210.

`[]bool` barely improves mod-6→mod-30 (2381→2449) — it was never strike-bound,
so cutting strikes doesn't help it; `[]u64` takes the full 1.41×.

## Cache-hierarchy sweep (quiet box, core pinned @ 4.54 GHz, 58 °C, no throttle)
Fixed N=1e9, `[]u64`, segment size swept 4 KiB → 256 MiB. Throughput reads the
cache sizes straight off the curve.

| seg | all-wheel M/s | odds M/s | note |
|-----|--------------:|---------:|------|
| 16–32 KiB | 571–573 | 1707–1714 | **peak ≈ L1d (32 KiB)** — optimal segment |
| 64 KiB–8 MiB | ~460–520 | ~1040–1670 | gentle L2→L3 slope (small latency gaps, hidden) |
| 16 MiB | 303 | 764 | **L3 edge (16 MiB): cliff begins** |
| 32 MiB | 228 | 519 | falling into DRAM |
| ≥64 MiB | ~200 | ~490 | **DRAM floor** (= non-segmented []u64 baseline) |

- Optimal segment = **L1d-sized (32 KiB)** — validates the default. Below it,
  per-segment overhead dominates; above it, cache residency erodes.
- L1→L2→L3 steps are gentle (latency gaps small, pipeline hides them); the one
  sharp cliff is **L3→DRAM, knee exactly at 16 MiB = L3 size**, ~2× drop.
- Wheel sets the mountain's *height*, not its *shape* (boundaries are hardware).
- Single-thread runs at max single-core turbo; parallel will share all-core
  turbo, so per-core throughput won't scale ×cores (per-core-scaling caveat).

## Segmentation vs N: the cache tax, decomposed (all-wheel, []u64)
Sweep N = 2^k two ways — whole-array (working set = array, grows with N) vs
segmented 32 KiB (working set pinned in L1). The segmented curve isolates the
*algorithm*; the whole-array curve is *algorithm × cache*; their ratio is the
pure cache penalty segmentation buys back.

| N | array | whole M/s | segmented M/s | seg ÷ whole |
|---|-------|----------:|--------------:|:-----------:|
| 2²³ | 1 MiB (L3) | 365 | 619 | 1.70× |
| 2²⁷ | 16 MiB (L3 edge) | 388 | 594 | 1.53× |
| 2²⁸ | 32 MiB (DRAM) | 259 | 589 | **2.27×** |
| 2³⁰ | 128 MiB (DRAM) | 201 | 576 | **2.87×** |

- **Segmented is cache-immune to N**: flat 619→576 over 2²³→2³⁰ (~7% sag = the
  Mertens strike-density creep, ~algorithm-only, no cliff). That's the whole
  point of segmentation — a 10¹⁸ run behaves like a 10⁹ one.
- **Whole-array crashes at DRAM**; the ratio grows with N (≈1.5× in-L3 → ~2.9× in
  DRAM → keeps rising) — the economic case for segmentation as a function of scale.
- Caveat: as a *cache probe* this is noisier than the fixed-N segment sweep —
  small/mid-N runs are sub-100 ms (overhead- + Mertens-confounded); the 2²⁶ point
  was a huge-page/L3-occupancy artifact. Clean signals here = segmented flatness
  + the DRAM cliff. For intra-cache structure, the segment sweep is the instrument.

## Wheel sweep: diminishing returns → regression (fixed N=1e9, []u64, 32 KiB)
Sweep p_n = largest wheel prime = 1(all),2,3,5,7,11,13 → M = 1,2,6,30,210,2310,30030.

| p_n | M | φ | desc B | M ints/s |
|-----|---|---|--------|---------:|
| 1 | 1 | 1 | 32 | ~500 |
| 2 | 2 | 1 | 32 | 1423 |
| 3 | 6 | 2 | 40 | 2377 |
| **5** | **30** | **8** | 88 | **3637 ← peak** |
| 7 | 210 | 48 | 408 | 1415 |
| 11 | 2310 | 480 | 3864 | 1500 |
| 13 | 30030 | 5760 | 46104 | 1670 |

- **Peak at mod-30; mod-210+ regresses ~2.2–2.6×.** Four factors, all turning at
  mod-210: (1) φ=8 is the last power-of-2 spoke count → `%φ` is a free AND (48
  needs magic-multiply); (2) desc 88 B keeps the ~300 KB of delta tables in L2,
  mod-210's 1.4 MB spills to L3 — the *wheel's own bookkeeping* becomes a
  memory-bound working set (cache story, recursively, on metadata); (3) 8
  residues = 1 byte (packing sweet spot); (4) the Mertens gain past mod-30
  (~1.28×) is too small to pay for 1–3.
- Our **general delta-table** wheel caps at mod-30 on this box; production sieves
  reach mod-210 only by hard-coding unrolled per-residue loops (no per-prime
  delta table). So the sweet spot is part hardware (φ=8, L2), part our honesty
  about bookkeeping cost.
- Perf bug the sweep exposed & fixed: the strike loop copied the whole per-prime
  delta array (`const d = pr.d`, up to 46 KB) every prime every segment → now a
  pointer. mod-30030: 212 → 1670 M/s (8×).

## Bucket sieve: a large-N / small-segment instrument
The bucket sieve is due to [Tomás Oliveira e Silva](http://sweet.ua.pt/tos/)
(and is what `primesieve` uses for large primes). Two-tier PrimeSource
(bucket_sieve.zig): SMALL primes (min delta ≤ seg_slots,
strike a segment ≥ once) stay in the cursor loop; LARGE primes (min delta >
seg_slots, strike ≤ once) are bucketed by next-strike segment (intrusive linked
list, one node/prime; drain-strike-refile). Avoids iterating + cache-streaming
every large prime on every segment.

**Only helps when √N > seg_slots** — else every prime is SMALL and there's
nothing to bucket. Inert (or negative) at our headline N=1e9/L1 config, because
√N=31623 ≪ 262144 slots.

Bucket vs naive, mod-30 []u64, **N=1e10** (√N=1e5), sweeping the segment down:

| seg | seg_slots | naive M/s | bucket M/s | speedup |
|-----|-----------|----------:|-----------:|:-------:|
| 1 KiB | 8192 | 1161 | 1414 | **1.22×** |
| 2 KiB | 16384 | 1527 | 1747 | 1.14× |
| 4 KiB | 32768 | 2086 | 2126 | 1.02× |
| 8 KiB | 65536 | 2588 | 2455 | **0.95×** |
| 16 KiB | 131072 | 2887 | 2825 | 0.98× |

- Win **grows as the segment shrinks** (√N/seg_slots ↑ → more large primes → more
  skip-work + cursor-streaming the bucket cuts). 1.22× at 1 KiB, still climbing.
- 8 KiB is **slower** (0.95×): threshold p>~124k > √N → *zero* large primes, so
  the bucket machinery (alloc + memset heads[] per run) is pure overhead. Honest
  cost of bucketing with nothing to bucket; a `n_large==0` fast path would fix it.
- Extrapolates to the record regime (N=1e18 → √N/seg_slots ~ 1e4, nearly all
  primes large, naive streams MB of cursors/segment from DRAM) where it's
  decisive — but that's a days-long run, not benchmarkable here. The monotonic
  climb is the proof of shape.
- Verified π(1e10)=455,052,511 for both. Current bucket uses a full per-segment
  heads[] array (fine to ~1e12 / ~1e6 segments); record N needs a circular window.

## π(N) scaling — the sieve as a lens on the number theory
Pure-π(N) driver (sweep.piScaling) over the engine; li(x)=Ei(ln x) via Ei's
convergent series. All π exact through 1e11 (=4,118,054,813).

| N | pi(N) | gap=N/pi | lnN | pi·lnN/N | li−pi |
|---|-------|---------:|----:|---------:|------:|
| 1e6 | 78,498 | 12.74 | 13.82 | 1.084 | 130 |
| 1e9 | 50,847,534 | 19.67 | 20.72 | 1.054 | 1,701 |
| 1e11 | 4,118,054,813 | 24.28 | 25.33 | 1.043 | 11,588 |

- **gap ≈ lnN − 1.08** — the second-order PNT (π ~ N/(lnN−1) beats N/lnN); the
  ~1.07 offset is the "−1", drifting toward 1.
- **pi·lnN/N = 1 + ~1/lnN** — crude PNT crawls to 1, always high by ~1/lnN (4.3%
  at 1e11). x/lnx is the bad approximation.
- **li−pi ~ √N/lnN** (RH error scale; 12.5k predicted vs 11.6k at 1e11), always
  positive (Li overshoots; Littlewood sign-flip only past Skewes ~1e316).
  Relative error 2.8 ppm — Li is ~15000× better than x/lnx.

The engineering built the engine; the engine demonstrates the theory we opened
with (PNT, Li, RH). Bucket engine so the harness scales into the large-N regime.

## Beyond u64: Int-parameterized range sieve + density validation
The coordinate/value type is orthogonal to the store word: `rangesieve.countInRange
(comptime Int, lo, hi, base_primes)` sieves [lo,hi) with positions/multiples in
`Int` (u64 or u128), base primes stay u64 (√(1e38) ≪ 2^64). `BitPacked(Word)`
separately parameterizes the store word (u8..u64), verified by unit test.

Validating high ranges where no exact oracle exists: use **theory**. A window
[N, N+Δ) holds ≈ Δ/ln N primes (local density li′=1/ln N), matched to Poisson
noise ~1/√count. Δ=1e8, base primes to 5e9:

| N | type | observed | Δ/ln N | ratio |
|---|------|---------:|-------:|:-----:|
| 1e12 | u64 | 3,618,282 | 3,619,114 | 0.9998 |
| 1e15 | u64 | 2,893,937 | 2,895,297 | 0.9995 |
| 1e18 | u64 | 2,414,886 | 2,412,747 | 1.0009 |
| **2e19** | **u128** | **2,249,954** | 2,250,110 | 0.9999 |

- Every ratio ≈ 1.000 within √count Poisson noise — sieve validated by PNT across
  7 orders of magnitude, no table needed.
- The 2e19 row is **past 2^64** — a tally primesieve can't produce, verified
  against theory alone. u128 works; arithmetic is trusted (provable algorithm +
  no overflow + this code validated at low N).
- Practical wall is base primes to √N (memory), ~1e20; past that needs the
  combinatorial (non-enumerating) methods. u128 buys the strip just past 2^64.

## Combinatorial π(x): Meissel–Lehmer (meissel.zig)
Counts primes *without enumerating them*: π(x) = φ(x,a) + a − 1 − P₂(x,a),
a = π(x^(1/3)). φ(x,a) is the recursive partial sieve φ(x,a)=φ(x,a−1)−φ(x/p_a,a−1)
with a mod-30 wheel base and the **leaf cutoff** φ(x,a)=1+max(0,π(x)−a) when
p_a²≥x (no coprime composites) — which drops it from exponential to O(x^(2/3)).
P₂ and the cutoff share a prefix-π table up to x^(2/3).

| x | π(x) | time | exponent |
|---|------|-----:|:--------:|
| 10¹⁰ | 455,052,511 | 0.03 s | — |
| 10¹¹ | 4,118,054,813 | 0.15 s | 0.71 |
| 10¹² | 37,607,912,018 | 0.89 s | 0.77 |
| 10¹³ | 346,065,536,839 | 5.68 s | 0.81 |
| 10¹⁴ | 3,204,941,750,802 | 41.3 s | 0.86 |

- All exact (matches OEIS A006880). **exponent** = log₁₀ of the time ratio per
  ×10 in x → **~0.8, sub-linear** (sieving is 1.0). Theoretical Meissel–Lehmer is
  2/3; the extra is the O(x^(2/3)) special-leaves / π structure.
- **LMO Stage A (memory):** the π-per-integer table is replaced by a compact
  bit-sieve + per-word count checkpoints (π(y) = `ckpt[y/64] + popcount`, O(1)).
  ~16× less memory (0.5 GB vs 8.6 GB at 10¹⁴) — moving the wall from ~10¹⁴ to
  ~10¹⁶ — *and faster*, because the small table is cache-resident where the fat
  one was DRAM-bound (10¹² 1.18→0.89 s, 10¹³ 7.25→5.68 s). Same algorithm, so the
  ~0.8 exponent is unchanged. The special-leaf sieve that reaches O(x^(1/3)) memory
  and the 2/3 exponent became `lmo.zig` — see the LMO section below.
- vs sieving (~3.6 G ints/s): 10¹³ in 5.7 s where sieving needs ~46 min; the gap
  *widens* with x (O(x) vs O(x^0.8)). This is why records go to 10³⁰
  combinatorially, never by sieving.

## The y knob: why x^(1/3) is a minimax, and how to break it (meissel.zig)

The identity holds for **any** y ≥ x^(1/3) (below it, y-rough n ≤ x can have 3
prime factors and you need P₃ — measured: α=0.7 gives the wrong answer). So
y = α·x^(1/3) is free to tune. Two consumers want the π-table, and they pull
opposite ways:

| consumer | needs π up to | vs α |
|---|---|---|
| φ's cutoff leaf φ(v,b), p_b² ≥ v ⇒ v ≤ p_b² | **y²** | grows |
| P₂'s π(x/p), p ∈ (y, √x] | **x/y** | shrinks |

They cross at exactly α = 1. **The classical x^(1/3) is a minimax point**, not a
coincidence — and raising y is pure loss (10¹¹: 106 ms/5 MB at α=1 → 7631 ms/329 MB
at α=8, 99% of it table build).

Surprise: **φ's leaf count is flat in y** — 5,927,155 at α=1 → 5,934,991 at α=8,
+0.13% while a grows 6.3×. Same reason the identity works: for p > x^(1/3),
φ(x/p, ·) has x/p < x^(2/3) ≤ p², so it cuts *immediately*. Every prime above
x^(1/3) adds exactly one trivial leaf. The recursion is already saturated at α=1.

**Capping the cutoff breaks the minimax.** Add "and v ≤ z" (z = x/y) to the leaf
test. Tightening a cutoff is always *correct* — the recursion is valid at every
node — it only trades tree size for table size. This deletes the y² consumer, so
the table is just z = x^(2/3)/α, which now **shrinks** with α:

```
             alpha=1.0   alpha=1.5   alpha=2.0   alpha=3.0   alpha=4.0   alpha=6.0
10^11  vs best  +2.8%       0.0%       +0.5%      +11.0%      +19.7%      +39.9%
10^12  vs best  +1.7%       0.0%       +3.1%      +10.7%      +20.6%      +39.5%
10^13  vs best  +6.6%       0.0%       +0.7%       +3.3%       +8.1%      +20.8%
```

- **α = 1.5 is the argmin at all three x** → now the default. It beats classical
  α = 1 on *both* axes (10¹¹: 99.6 ms/3.4 MB vs 102.4 ms/5.1 MB) — free, not a trade.
- The optimum **doesn't move** over two decades; the **basin widens**. α=3 costs
  +11% at 10¹¹ but only +3.3% at 10¹³ (3× less memory: 36.9 vs 110.7 MB).
- Mechanism, from the noise-free leaf counts (leaves ~ α^k): k = **0.611** (10¹¹),
  **0.473** (10¹²), **0.363** (10¹³). φ's sensitivity to α decays with x against a
  build cost that is always 1/α — so the free memory α buys **grows with x**.
- Capping alone bounds the table by x^(2/3)/α with α < x^(1/6) (y < √x, where P₂
  vanishes and it degenerates to Legendre) → **Θ(√x) memory**, a whole exponent
  step from one condition. LMO's x^(1/3) needs the leaf restructuring too.
- Caveat: our α_opt should *not* be expected to match the literature's α ~ log³x —
  that tunes LMO's *sieve-range vs special-leaf* trade; ours is *table-build vs
  φ-tree*. Different mechanism. Related: capping is **not** the LMO decomposition —
  capped φ still requires p_b² ≥ v, so it recurses past the m·p frontier into
  3+-factor d. LMO stops at that frontier unconditionally and pays a sieve query:
  ~4× *fewer* leaves (1.4 M vs 5.9 M at 10¹¹), each more expensive.

## LMO end to end: the 2/3 exponent, in Θ(x^(1/3)) memory (lmo.zig)

Leaf structure derived from *our* recursion, not from memory. φ(x,a)=φ(x,a−1)−φ(x/p_a,a−1)
takes primes in **decreasing** index order, so with F(n,b) = μ(n)·φ(x/n,b):
F(n,b) = F(n,b−1) + F(n·p_b,b−1) — the tree sums F exactly. Cut any child with n·p_b > y:

- **Ordinary leaf** — n ≤ y reaches b=0 ⇒ μ(n)·⌊x/n⌋. Every squarefree n ≤ y gets there
  (all prefixes ≤ n ≤ y), so S1 = Σ_{n≤y} μ(n)⌊x/n⌋. O(y), direct — and **c = 0 is exactly
  right**, not a shortcut.
- **Special leaf** — n = m·p, n > y ≥ m = n/p. Primes *descend* along the path, so the last
  one multiplied is the **smallest**: p = P⁻(n), all of m's primes above it. ⇒
  **S2 = Σ −μ(m)·φ(x/(m·p_b), b−1)** over p_b-rough squarefree m ∈ (y/p_b, y].

The b−1 index is what makes the sieve work: at step b the sieve over [1, z], z = x/y, holds
exactly p_1..p_{b−1}, so sweeping b = 1..a is one monotone prime-at-a-time removal and
φ(v, b−1) is a prefix count. Every leaf has v = x/(m·p_b) ≤ z because m·p_b > y.

**Segmentation** rests on two ideas: a *running φ per b* (φ(v,b−1) = phi_run[b−1] + alive in
[lo,v]), so [1,lo) is never re-sieved; and a *descending m-cursor per prime* — segments run lo
ascending and v = x/(m·p) rises as m falls, so every m is touched once in **total**, not once
per segment.

**P₂ needs no Fenwick**: it is the monotone sweep (p ascending ⇔ x/p descending). One linear
walk per segment with a running count answers every π(x/p), and Σ(π(p)−1) collapses to a
closed form since those primes have indices a+1..A.

| x | π(x) | LMO | Meissel | speedup | exponent |
|---|------|----:|--------:|--------:|:--------:|
| 10¹¹ | 4,118,054,813 | 0.09 s | 0.10 s | 1.1× | — |
| 10¹² | 37,607,912,018 | 0.42 s | 0.69 s | 1.6× | 0.655 |
| 10¹³ | 346,065,536,839 | 1.94 s | 5.81 s | 3.0× | 0.661 |
| 10¹⁴ | 3,204,941,750,802 | 1.21 s | 41.3 s | **34.2×** | 0.654 |
| 10¹⁵ | 29,844,570,422,669 | 15.2 s | — (OOM) | — | 0.658 |
| 10¹⁶ | 279,238,341,033,925 | 72.9 s | — (OOM) | — | **0.681** |

10¹⁵ and 10¹⁶ are new ground — Meissel's ceiling was 10¹⁴, and it needed 73.8 MB there. LMO does
10¹⁶ in **~6.9 MB**: μ/lpf tables 4.3 MB + a-sized arrays 1.6 MB + p-range buffer 0.9 MB + counter
0.1 MB.

| 10¹⁷ | 2,623,557,157,654,233 | 369.3 s | — (OOM) | — | 0.705 |
| 10¹⁸ | 24,739,954,287,740,860 | 1686.1 s (28.1 min) | — (OOM) | — | 0.660 |
| 10¹⁹ | **234,057,667,276,344,607** | 8680.3 s (2.41 h) | — (OOM) | — | 0.712 |

10¹⁹ is the **last power of ten reachable in u64** (x < 2^64 = 1.845e19; 10²⁰ is 5.4× over). The
value matches M. Deléglise's 1996 computation — reproduced here on one core with the *predecessor*
algorithm. Everything else has vast headroom (m·p ≈ 7.4e13, 250,000×; z as i64, 8e6×; lpf/primes
u32 good to x ≈ 1e27), so going further needs u128 for **x, z and the divisions only**.

**Lemma 5.1 confirmed to 0.03% at 167 billion leaves**: a = π(y) = 578,568 ⇒ ½π(y)² =
167,370,465,312 vs measured **167,423,661,964**.

**π(10¹⁸) exact, in 28.1 min and ~31 MB** — six decades past Meissel's ceiling, single-threaded.
The leaf count came in at 40,105,355,244 against the a²/2 law's 40,090,000,000 (0.04%), a closing
check on the leaf derivation.

**The exponent: mostly 2/3, with a mild real drift — and a cautionary tale about 2-point deltas.**

Per-decade exponents read 0.665 / 0.655 / 0.658 / 0.681 / **0.705** / 0.660. The middle of that was
taken as a climb away from 2/3 and given two mechanistic explanations (below). Then 10¹⁸ came in at
**0.660** and the story collapsed. Fitting *all seven points* by least squares instead of
differencing adjacent pairs:

| range | slope |
|---|---:|
| 10¹² .. 10¹⁸ | **0.6715** |
| 10¹² .. 10¹⁵ | 0.6586 |
| 10¹⁵ .. 10¹⁸ | 0.6842 |
| 10¹⁶ .. 10¹⁸ | 0.6822 |

With 10¹⁹ added (8 points):

| range | slope | points |
|---|---:|---:|
| 10¹² .. 10¹⁹ | **0.6755** | 8 |
| 10¹² .. 10¹⁵ | 0.6586 | 4 |
| 10¹⁵ .. 10¹⁹ | 0.6879 | 5 |
| 10¹⁶ .. 10¹⁹ | 0.6887 | 4 |

**The drift is real: 0.659 → 0.688, and the overall slope 0.6755 sits above 2/3.** The per-decade
deltas (0.705, 0.660, 0.712) *alternate* — that is run-to-run noise on unrepeated runs, and it
fooled the analysis twice in both directions: first into "the exponent is climbing to 0.705", then
into "0.660 falsifies the drift". Neither 2-point delta was a measurement. Only the multi-point
fits are, and they say the drift is genuine and mild. Roland predicted it; the data agrees.

Candidate mechanisms, still unresolved:

Two candidate mechanisms, **neither established**:

1. *Cache*: the a-sized arrays (`phi_run`, `seg_cnt`, `cur` = a×24 B) are streamed once per
   segment and cross the 512 KB L2 around there — 403 KB at 10¹⁴, 812 KB at 10¹⁵, 1.6 MB at 10¹⁶,
   3.3 MB at 10¹⁷. Same structure that pins the segment size. But the arrays crossed at 10¹⁵ while
   the exponent only turned at 10¹⁶.
2. *Query growth*: the counter query is ~3·nwords^(1/3) with nwords = S/64 and S ≈ y ~ x^(1/3), so
   **the query itself grows as x^(1/9)** — making leaves×query ~ x^(2/3+1/9)/ln²x = x^0.778/ln²x
   against kills at exactly x^(2/3). S cannot simply be frozen to fix this: a·z/S would then grow
   linearly, and S ~ x^(1/3) is precisely what holds that term at x^(2/3). The two constraints fight.

**Three retractions.** (0) The "climbing exponent" itself was over-read from 2-point deltas on
unrepeated runs; least squares over all seven points says 0.6715 ≈ 2/3, with only a mild
0.659 → 0.684 drift surviving. Single-decade deltas at 28 min/run are not a measurement — min-of-N
and a fit over all points should have been the method from the start.
(a) An α sweep at 10¹⁵/10¹⁶ was run as a "discriminating test" — cache would
predict α_opt drops, query-growth that it holds. It discriminates nothing: at 10¹⁶ even α=2 leaves
the arrays at 849 KB, **still outside L2**, so lowering α never restores residency. The hypothesis
was never on trial. (b) The query-growth model was first reported as fitting the measured exponents
across three decades. Recomputed with the actual power-of-two-rounded query costs (nwords
2901/6250/13465/29010 → avg query 22/28.5/39/46.5) against the real leaf counts, it predicts
**0.723 / 0.748 / 0.690** vs measured **0.658 / 0.681 / 0.705**. Direction right, magnitude wrong —
not a fit, and it should not have been called one.

**What is solid:** α_opt = 4.0 measured at 10¹¹, 10¹², 10¹³, 10¹⁴, 10¹⁵ and 10¹⁶ alike — five
decades flat, where the literature's α ~ log³x demands a steady rise.

| α (seg = x^(1/3) fixed) | 10¹⁵ | 10¹⁶ |
|---|---:|---:|
| 2 | +26.2% | +27.9% |
| 3 | +4.6% | +6.0% |
| **4** | **0.0%** | **0.0%** |
| 6 | +11.0% | +9.8% |

Measured: 10¹⁸ took **28.1 min** single-threaded — close to the ~30 min the 0.705 extrapolation
predicted, and to the ~25 min the 0.66 one did, which is itself a hint that a decade of argument
over 0.66-vs-0.705 was arguing about noise.

The 10¹⁴ column, step by step — each is a separate committed measurement:

| change | 10¹⁴ | vs prev |
|---|---:|---:|
| Fenwick, α=1.5 | 16.95 s | — |
| → block Counter (O(1) kill) | 8.91 s | 1.90× |
| → α = 2 | 7.98 s | 1.12× |
| → m-walk split at √y | 6.07 s | 1.31× |
| → P₂ fused onto the counter | 4.97 s | 1.22× |
| → α = 4 | 4.13 s | 1.20× |
| → 3-level counter, power-of-2 blocks | 3.54 s | 1.17× |
| → class-(1) binomial (LMO p.555 / DR §6.1) | 2.83 s | 1.27× |
| → π-table for classes (2)/(3)-easy | 2.42 s | 1.17× |
| → mod-30 wheel fold (DR §9) | 1.59 s | 1.52× |
| → branchless kill (*after* the wheel) | 1.25 s | 1.27× |
| → gate the diagnostics behind comptime | **1.21 s** | 1.03× |

## What the papers actually say (literature/, read after the fact)

LMO 1985 (Math. Comp. 44, 537–560) and Deléglise–Rivat 1996 (Math. Comp. 65, 235–245), read only
after the implementation was measured into its current shape. Scorecard:

**Confirmed, derived independently here:**
- **δ(n) = the SMALLEST prime factor** (DR eq. 10). Our P⁻ derivation from our own recursion was
  right, and the earlier `maxpf` counter was the wrong convention.
- **Our S2 is DR eq. (11)**, term for term.
- **S₀ = Σ_{n≤y} μ(n)⌊x/n⌋ with φ(u,0) = [u]** — c=0 confirmed. LMO's own implementation uses a
  wheel k=5; DR drops it to 0. Roland's instinct ("why not c=0 and be done") is DR's choice.
- **x^(1/3) ≤ y ≤ x^(1/2)** is DR's stated precondition — the invariant the π(2) bug forced on us.
- **DR = clustering**: *"many special leaves could be computed at the same time… saving a log x
  factor… we show it is possible to compute more special leaves at the same time, saving another
  log x factor."* The one claim asserted here without evidence, and it holds.
- **LMO Lemma 5.1: "The number of special leaves is ½·π(y)² + O(y^(3/2)/log y)"** — exactly the
  a²/2 law fitted empirically here (0.03% at 10¹⁹). Its proof *is* our √y enumeration split:
  *"those for which δ(n) < √y and those for which δ(n) ≥ √y… in that case we must have n = pq."*
- **LMO Truncation Rule T′: "x^(2/5) ≥ y ≥ x^(1/3)"** — Roland recalled this bound from memory,
  including his uncertainty about the upper end. It is exact.

**The α mystery, resolved.** LMO p. 556: *"We choose **y = c·x^(1/3)** which balances the sieving
and calculating special leaves… (In the actual implementation, a good value of the constant **c
was determined empirically**.)"* Our α_opt = 4, flat over five decades and found empirically, **is
LMO's own prescription.** The α ~ log³x measured against here is DR's *asymptotic space bound*, and
at 10¹⁸ it is infeasible anyway: log³x = 71,197 while y ≤ √x caps α at x^(1/6) = 1000. Six sweeps
were not contradicting the literature — they were being compared to the wrong paper's constant.

**What we missed — now implemented, 1.27×.** LMO p. 555 / DR §6.1: for **p > x^(1/3)**, x/(pq) < x^(1/3) < p
so **φ(x/(pq), π(p)−1) = 1 identically** — not "a cheap lookup", *literally 1*. The whole class
collapses to a binomial: S₁ = C(π(y) − π(x^(1/3)), 2), **constant time**. At α=4 that is **52% of
our leaves** — 20.9 × 10⁹ counter queries at 10¹⁸ computing what is one formula. The "easy leaves
cannot be cheaper" argument recorded above is true and irrelevant: for this class you do not
evaluate φ, you *count* it. Our `easy` predicate (p² > v ⟺ m·p³ > x) is the weak version; the
collapsing condition is **v < p ⟺ m·p² > x**. LMO's own "easy" is different again — x/(pq) ≤ y,
readable from a π table over [1, y], which is O(x^(1/3)) and affordable; the prefix-array failure
here assumed the table had to span [1, z].

Implemented as S₁ = C(n₁, 2) with n₁ = #{primes p : p³ ≥ x, p ≤ y}. The threshold comes from icbrt,
not `p*p*p`, which would overflow (p ≤ 4x^(1/3) ⇒ p³ ≤ 64x = 6.4×10²⁰ at 10¹⁹). Measured **1.27×**
at 10¹²–10¹⁴ (3577 → 2827 ms at 10¹⁴). Less than the 52% leaf share implies — back out the ratio
and leaf work is ~41% of runtime, the rest being kills, folds and per-segment overhead.

Verified differentially, which is the point: the fused path never enumerates class (1), yet its S₂
is bit-identical to `specialS2Segmented`, which enumerates every leaf — and the leaf *counts* match,
so the binomial counts exactly the set it skips.

**Classes (2) and (3)-easy: one rule, not four cases.** Working the conditions through, every
cheap class reduces to the same test — **v ≤ y and p² > v** ⇒ φ(v,b−1) = 1 + max(0, π(v) − (b−1))
straight from a π table over **[1, y]** (O(y) space, 0.74 MB at 10¹⁴, ~34 MB at 10¹⁹):

- **class (2)**, p > √(x/y): v < x/p² < y, and p > x^(1/3)/2 > x^(1/4) ⇒ p⁴ > x ⇒ v ≤ x/p² < p².
- **class (3)-easy**, q ≥ x/(yp): v ≤ y, and p² > (x/y²)² > y ≥ v — which requires **x² > y⁵, i.e.
  y ≤ x^(2/5)**. *That is why LMO's Truncation Rule T′ caps y at x^(2/5)*: it is exactly the
  condition that keeps this class's m prime and p² above v. α=4 satisfies it for x ≥ 4¹⁵ ≈ 10⁹.

The table only spans [1, y], never [1, z] — every leaf it serves has v ≤ y by construction. That
distinction is the whole reason the earlier flat-prefix-array attempt failed (it built an O(z)
structure to serve O(leaves) queries at z/leaves ≈ 27).

Measured share of the leaves surviving the binomial that this rule catches — **it grows with x**:

| x | 10⁹ | 10¹⁰ | 10¹¹ | 10¹² | 10¹³ | 10¹⁴ |
|---|---:|---:|---:|---:|---:|---:|
| π-table evaluable | 43.1% | 47.8% | 52.5% | 55.9% | 60.2% | **64.2%** |

Worth **1.17×** at 10¹⁴ (2827 → 2426 ms), also growing with x (1.11× at 10⁹). Cumulative with the
binomial: **1.46×** (3.54 → 2.42 s). Remaining for the sieve: class (3)-hard (v > y) and class (4).

**Still on the table:** DR §6.2+ *clustering* — computing many special leaves at once, which is the
second log factor and DR's actual novelty. The local PDF is 6 of 11 pages, so §6.5 onward and the
complexity analysis are missing.

**Where we are ahead.** LMO's array a(i,j) in (3.9)–(3.10) — a binary hierarchy of counts, queried
by decomposing l into powers of two — **is a Fenwick tree**, never named as such. LMO's own analysis
says *"the number of updating operations of the array {a(i,j)} during each block takes O(x^(1/3)
log x)… adds up to O(x^(2/3+ε))"*: the updates are the **dominant** term at x^(2/3)·log x. The O(1)
kill Counter here removes that log from the paper's leading cost — found not by cleverness but by
measuring the 60:1 kill:query traffic. Our √y prime-m split also beats DR's x^(1/4) threshold
(2·x^(1/6) ≈ 2000 vs 31,623 at 10¹⁸), though it is LMO Lemma 5.1's split rediscovered.

## The profile, and the mod-30 wheel (DR §9)

Nine cost models into this file, we finally ran `perf`. It disagreed with all of them.

```
IPC = 2.65 (measured)          L1-dcache miss 0.78%       branch-miss 2.94%

CYCLES by source line          BRANCH MISSES by source line
 7.44%  lmo.zig:843  fold      13.72%  lmo.zig:370  if (bits[w] & b != 0)
 6.40%  lmo.zig:370  kill      12.01%  lmo.zig:375
 4.83%  lmo.zig:370            6.37%   lmo.zig:843
 4.30%  lmo.zig:375            4.36%   lmo.zig:371
```

**The fold and kill are ~50% of all cycles.** The leaves — the thing every optimisation above
attacked — are the other half. And L1 misses at 0.78% are evidence *against* the a-array/L2 story
for the residual exponent drift.

**A falsified inference, from misreading our own profile.** "Line 370 = 13.72% of branch misses"
looked like *unpredictable*, so the branch was removed (`alive = @intFromBool(...)`, subtract 0 or
1 unconditionally). Measured **0.75× — a 33% regression**. The arithmetic: 13.72% of 137.8M is
18.9M misses over ~1.6×10⁹ fold iterations = a **1.18% miss rate**. That branch is *highly*
predictable — the already-dead pattern within one prime's strike loop is periodic, and TAGE
handles it. It tops the chart because it **executes most**, not because it misses most. *Share of
misses ≠ miss rate.* Reverted.

**So the fold's 50% is VOLUME.** It runs Σ_p z/p ≈ 2.76z but only z kills land: **66% of visits
hit an already-dead element**, and the predictor sees them coming. The only lever is fewer visits
— which is DR §9: *"Precomputing the sieving by the first primes 2, 3, 5 saves some more time."*

**The decomposition that made it cheap.** The wheel's *stepping* (which multiples to visit) is
independent of the array's *indexing* (bit-per-integer vs 8-per-30). All the fold win is in the
stepping, so take that and leave the indexing — and every leaf/prefix path — untouched:

- `reset()` starts from a mod-30 mask instead of all-ones (lcm(30,64) = 960, so the u64 pattern
  repeats every 15 words; segments align to 960). That state *is* φ(·,3).
- Fold p ≥ 7 by stepping j = p·m over m coprime to 30 (gaps 6,4,2,4,2,4,6,2). 2/3/5 are never
  folded at all.
- Σ_{p≥7} (8/30)·z/p ≈ **0.46z visits vs 2.76z — 6× fewer**.
- Leaves with p ∈ {2,3,5} need φ(v,0..2), which a 2·3·5-wheel cannot represent — but those are
  meissel.zig's closed-form base cases, and there are only ~0.5y of them (~93k of 141M at 10¹⁴).

**1.52× at 10¹⁴** (2422 → 1591 ms), growing with x (1.30× at 10⁹). π(10¹⁴) = **1.59 s, 26.0×
Meissel**. Verified differentially: `specialS2Segmented` still uses an all-ones reset and a plain
fold over every prime, and the wheel path matches it bit for bit.

Two bugs the small-x sweep caught, both from `lo` now starting at 0 (mask alignment): `x / lo`
trapped with SIGFPE, and `for (3..a+1)` had start > end when a = π(y) < 2.

**Re-profiling after the wheel: the balance inverted.**

```
                        before wheel      after wheel
instructions            29.6e9            11.5e9        2.6x fewer
IPC                      2.65              1.74         WORSE
branch-misses          137.8M            128.8M        barely moved
branch-miss rate         2.94%             6.66%
```

2.6× fewer instructions bought only 1.5× the time, because the fold was the
*ILP-friendly* half (independent strikes, predictable branches) and what remains is
dependent loads and data-dependent branches. Branch misses barely moved — they were never
the fold's, they were the leaves' — and at ~17 cycles they are now ~⅓ of runtime.

**The same patch, reverted an hour earlier, is now a 27% win.** The branchless kill measured
**0.75×** before the wheel and **1.27×** after. Nothing about the code changed; the *context* did:

| | pre-wheel | post-wheel |
|---|---:|---:|
| kill alive-check miss rate | 1.18% | **9.0%** |
| fold visits | 1.6×10⁹ | 2.5×10⁸ |

The wheel deleted p = 2/3/5 — exactly the primes whose already-dead patterns are short and
periodic, and which the predictor nailed. What is left is p ≥ 7, whose aliveness follows irregular
factorisations. So the branch got 8× *less* predictable while the cost of removing it fell 6.4×.
Both terms moved and the sign flipped. **An optimisation's sign can depend on another
optimisation** — no amount of reasoning about "is this branch predictable" survives that, only
re-measuring after each change.

**DR §7's x^(1/4) bound: measured, and NOT worth it for us.** DR sieves [1, x/y] "by all primes
less than x^(1/4)" — for us that is π(x^(1/4)) = 446 primes instead of π(y) = 16,801 at 10¹⁴, an
83× smaller `a` at 10¹⁸. It looked like the biggest remaining algorithmic win. Measured:

| x | fold visits | share with p ≥ x^(1/4) |
|---|---:|---:|
| 10⁹ | 89,674 | 34.2% |
| 10¹¹ | 2,175,571 | 28.9% |
| 10¹³ | 51,332,456 | 24.9% |
| 10¹⁴ | 247,641,412 | **23.6%** |

The share **shrinks** with x (Σ_{x^(1/4)≤p≤y} 1/p → ln(4/3), constant, while Σ_{7≤p≤y} 1/p grows as
ln ln y). So it saves 38% × 0.236 ≈ **9%**, plus ~2% from smaller a-arrays — against **~13.5%** to
un-fuse P₂ (the bound leaves the counter at φ(·,π(x^(1/4))), so π(v) = φ(v,a) − 1 + a breaks, and
P₂ needs DR §4's own [1,z] sieve: 5.4×10⁸ at ~3.6 G/s ≈ 0.15 s of a 1.11 s run). **Net ≈ −2%.**

DR §8.3 explains why it pays for *them*: *"Owing to the necessity of quickly retrieving the values
φ(u,b), we have to maintain a data structure such that each access costs **O(log x) instead of O(1)
in a normal sieve**."* **The x^(1/4) bound exists to ration expensive tree operations.** Our O(1)-kill
counter is a whole log x cheaper than their S3 sieve, which devalues the optimisation the bound
provides — and makes the P₂ fusion (criticised above as "optimising the 13% and inflating the term
that pins everything else") the right trade after all. Both halves of that judgement were wrong,
and identically so: pricing a paper's design decision without pricing the data structure it was
designed around.

**The instrumentation was costing 3.3%.** Nobody had ever measured the measurement. The
`walk`/`leaves`/`easy` counters sit in the hot leaf loop, and `easy` recomputes the exact
`v ≤ y and p*p > v` test that `leafPhi` does again one line later:

| x | 10⁹ | 10¹¹ | 10¹² | 10¹³ | 10¹⁴ |
|---|---:|---:|---:|---:|---:|
| cost of diagnostics | 9.3% | 7.9% | 6.7% | 6.5% | **3.3%** |

The share falls with x as the class-(1) binomial absorbs a growing fraction of leaves (they are
counted arithmetically, not enumerated). Gated behind a comptime `INST` flag rather than deleted —
`leaves` is what confirmed Lemma 5.1 to 0.03% at 10¹⁹, `walk` is the regression guard on the √y
split, and `easy` sized the π-table class. Both paths verified to produce identical S₂ and P₂.

**Fixed-trip prefix: theory confirmed, trade still loses.** `prefix()`'s three loops stop at
sblk / blk / w — all functions of the query's position, so every exit branch is unpredictable by
construction: 1.36 mispredicts per query, 39% of all branch misses, ~23 of ~62 cycles/query. Made
them fixed-trip (walk a constant nsuper / 2^s2 / 2^s1 slots, mask the out-of-range ones to zero,
pad the tree so the over-read is always zeros). Same trade as the branchless kill.

```
                    variable-trip      fixed-trip
cycles                  5.01e9           5.42e9      +8%
instructions           13.26e9          17.33e9      +31%
branch-misses          56.7M            17.6M        -69%
branch-miss rate         3.67%            1.01%
IPC                      2.65             3.19
```

It removed **69% of every branch miss in the program** and IPC jumped 2.65 → 3.19 — the mechanism
worked exactly as predicted. And it is **4.3% slower**. Backing the numbers out: the misses were
worth ~0.87×10⁹ cycles (~17%, ≈22 cycles each), the extra 4.07×10⁹ instructions cost ~1.28×10⁹
(~25%). **So prefix's 62 cycles/query are WORK-dominated, not mispredict-dominated** — the opposite
of the branchless kill, and the reason that one won and this one loses.

Also: `s += if (cond) v else 0` compiled to a **branch, not a cmov** (39M misses vs the AND-mask's
17.6M, same instruction count, worse time). The `-@intFromBool(...)` AND-mask is the better form,
and still not enough.

**So the prefix opportunity is smaller than it looked.** SIMD is the only remaining path, and it is
unpromising: the `cnt1`/`cnt2` u32 sums vectorise well, but the third loop needs popcount over 16
u64s, and this box has AVX2 without **VPOPCNTDQ** — the `vpshufb` nibble-LUT is ~1.5 ops/u64
against scalar `popcnt`'s 1/cycle throughput. Reverted.

## DR §6.5 clustering: measured reach, and why y = x^(1/3)·log³x exists

DR's actual novelty: *"for each p we can split the summation over q into sums over intervals where
the function q ↦ π(x/pq) is constant."* q ascending ⇒ v = x/(pq) descending ⇒ π(v) non-increasing,
so the runs are contiguous. Measured the reach **before** building it (the last two "biggest
remaining wins" both measured negative):

```
x = 10^14
 alpha         y  pi-tab leaves           runs     ctr leaves  leaf/run    reach
     2     92830       13683218        7311759       17274212      1.87    44.2%
     4    185660       53037184       21780399       14674931      2.44    78.3%
     8    371320      127796929       36082736        8867918      3.54    93.5%
    16    742640      262535484       45657832        3541600      5.75    98.7%
    32   1485280      507399194       49087615        1266051     10.34    99.8%
```

**At our α=4 it is worth ~8–12%** — reach 78.3%, leaf/run 2.44, and finding each run boundary costs
about what evaluating a leaf costs. Not the log x the paper advertises.

**But the reach grows with α, and that is the point.** Three things move together:
- reach 44% → **99.8%**: at large y essentially every leaf is clusterable;
- leaf/run 1.87 → **10.34**: bigger y ⇒ v smaller and flatter in q ⇒ longer runs;
- **the runs saturate** — 21.8M → 49.1M while the leaves go 53M → **507M**. Leaves ride α²; runs
  ride ~α^0.4. *Clustering converts the α² leaf term into a nearly-flat one.*

Meanwhile the two costs of large α collapse: **counter leaves 14.7M → 1.3M** (11× fewer — the ones
clustering cannot reach), and z = x/y shrinks 8× (5.4×10⁸ → 6.7×10⁷), taking the fold with it. At
α=32 that is 2.25× more runs to buy 8× less fold and 11× fewer counter queries.

**So DR's y = x^(1/3)·log³x is not a separate choice from clustering — it is what clustering pays
for.** Predicted on that basis that α_opt (measured at 4 on an implementation where every leaf costs
a lookup) should move a long way up once the α² term is flattened.

**Built it. Correct — and a net loss.** S₂ bit-identical to `specialS2Segmented` (which evaluates
every leaf individually), leaf counts matching, 0 failures; run counts matched the standalone
prediction exactly (5,510,828 at 10¹³).

```
x = 10^14                      (alpha=4 without clustering: 1209.2 ms)
 alpha         y       z = x/y          ms    vs best
     4    185660     538618980      1346.8       4.1%
     8    371320     269309490      1294.2       0.0%     <- alpha_opt DID move 4 -> 8
    16    742640     134654745      1515.1      17.1%
    32   1485280      67327372      2064.6      59.5%
    64   2970560      33663686      3521.8     172.1%
```

**α_opt moved 4 → 8, exactly as predicted** — the α² leaf term really did flatten. And the best
clustered configuration is still **7% worse than no clustering at α=4**. Every step of the reasoning
held and the answer is still no. Two reasons: a **run costs ~2.5 leaves** (three random lookups —
`pi_tab[v]`, `primes[c]`, `pi_tab[qmin−1]` — against a leaf's one, and the leaf's `primes[qi]` walk
is sequential), so leaf/run must clear ~2.5 to break even; and the **dense m-walk rides α^1.5**
(flagged before building), which clustering cannot reach and which swamps z's 1/α shrink past α≈8.
Hoisting the two p-invariant divisions out of the run finder moved it 0.844 → 0.854 — the divisions
were never the problem.

## The pattern: DR's optimisations ration costs we had already removed

Both of DR's headline optimisations measure **negative** here, and for the same reason:

| DR optimisation | what it rations | what we did | measured |
|---|---|---|---|
| §7, x^(1/4) bound on the φ-sieve | tree ops at O(log x) each (§8.3) | O(1)-kill counter | **−2%** |
| §6.5, clustering | leaf evaluation | π table: one lookup | **−7%** |

**They are amortisation schemes for costs we had already made cheap.** That is why we sit at ~1.6×
their 1996 implementation while implementing *neither* of their log factors — we did not skip them,
we made them unnecessary. It also retires the α_opt puzzle: six sweeps said α_opt = 4 flat across
five decades against a literature that says α ~ log³x, and that was never a contradiction — DR's α
is large *because* their leaves are expensive. Ours are not.

## The α knob, where it is finally real (lmo.zig)

## The α knob, where it is finally real (lmo.zig)

In capped Meissel α was nearly free — φ's leaf count barely moved, so the basin was flat. In
LMO it is a **sharp interior optimum**, because the two terms genuinely fight:

- z = x/y = x^(2/3)/α → sieve kills fall as **1/α**
- leaves ≈ π(y)²/2 → rise as **α²**

Segment size held fixed at x^(1/3) across the sweep, so this is one knob, not two:

| α | z = x/y (10¹³) | leaves | φ ms | P₂ ms | total | vs best |
|--:|---:|---:|---:|---:|---:|---:|
| 1.0 | 464,166,357 | 2,939,588 | 2042 | 528 | 2570 | +35.8% |
| 1.5 | 309,444,238 | 6,035,944 | 1667 | 341 | 2008 | +6.1% |
| **2.0** | 232,083,178 | 10,182,481 | 1640 | 252 | **1892** | **0.0%** |
| 3.0 | 154,722,119 | 20,976,262 | 1908 | 163 | 2072 | +9.5% |
| 4.0 | 116,041,589 | 35,211,802 | 2475 | 122 | 2597 | +37.2% |
| 16.0 | 29,010,397 | 437,536,611 | 22003 | 35 | 22038 | **+1064.7%** |

- **α_opt = 2.0 at 10¹¹, 10¹² and 10¹³ alike** → now the default. Worth 6–7% in the sweep,
  ~11% end-to-end (π(10¹⁴) 8.91 → 7.98 s).
- The curve is **steeply asymmetric** — gentle below (+36% at α=1), catastrophic above
  (+1065% at α=16). That is the α² leaf term, exactly as predicted.
- **P₂ monotonically wants large α** (528 → 35 ms as z shrinks), but it is too small a share
  to move the optimum.
- **α_opt does not scale with x** — 2.0/2.0/2.0 over two decades, where α ~ log³x predicts a
  ~1.65× rise. Same non-result as capped Meissel's flat 1.5. Two decades is not a law, but
  there is no visible trend to extrapolate.

**The confound, and α_opt = 4 after fixing it.** `w/leaf` ran 15–23, matching the ~2·ln y
prediction: we rescanned (y/p_b, y] per b and rejected most of it, so the m-walk *also* rode α².
Fixed by splitting the enumeration at p = √y (see below): w/leaf → **1.12**. That flattened the
high-α side — α=8 cost +250%, now +41% — and the optimum drifted **2 → 4**, worth a further 23%:

|  | α=1 | α=2 | α=3 | **α=4** | α=5 | α=6 | α=8 | α=12 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 10¹² | +98.5% | +24.7% | +4.8% | **0.0%** | +2.1% | +10.0% | +36.3% | +109.6% |
| 10¹³ | +105.2% | +24.0% | +5.6% | **0.0%** | +3.4% | +12.6% | +40.9% | +121.1% |
| 10¹⁴ | +104.8% | +23.2% | +4.7% | **0.0%** | +4.1% | +13.1% | +45.5% | +136.6% |

**α_opt still does not scale with x** — flat 1.5 in capped Meissel, then flat 2, now flat 4 over
four decades, where α ~ log³x demands a steady rise. Three independent sweeps, same non-result.

**Enumeration: split at √y.** A per-lpf linked list does NOT work — the valid-m set depends on b
and we re-walk b inside every segment, so it would be rebuilt per segment at y·lnln y × nsegs ≈
1.4×10⁹, *worse* than the 1.9×10⁸ walk it replaces. Instead:
- **p > √y** — any p-rough m ≤ y is 1 or a prime (a composite would be ≥ lpf² > y), and m=1 fails
  m·p > y. So m = q prime in (p, y], with m·p > y automatic. Enumerable straight from primes[] —
  no walk, no rejects. This is the a²/2 bulk.
- **p ≤ √y** — keep the walk; it now costs π(√y)·y ≈ 2×10⁶ instead of 1.9×10⁸.

**P₂ fuses onto the same counter, for free.** After folding all a primes, the alive set in [1,z]
is exactly {1} ∪ primes in (y, z] — a composite with lpf > y is ≥ y², and y² > z whenever α > 1.
So π(v) = φ(v,a) − 1 + a, and P₂'s π(x/p) is a counter query at the end of each segment's b loop.
P₂'s whole separate sweep disappears: 1.24×, more than its 13% share, because it deletes a prime
sieve rather than just bookkeeping.

**Falsified: "easy leaves are cheap".** Predicted that resolving easy leaves (p² > v) via π(v)
from a flat O(1)-query prefix array would defuse the α² term. Measured **5% slower**, and the
analysis says it never could work: (1) the build is O(z) but the queries are only the leaves, and
z/leaves ≈ 27 — the Fenwick mistake inverted, optimise the query and pay on the build; (2) more
fundamentally, for an easy leaf the survivors ≤ v with lpf ≥ p_b are *already* just 1 plus the
primes in [p_b, v], so `ctr.prefix(v)` **is** φ(v,b−1) — the π-formula is a different route to the
identical number, not a cheaper one. **So DR's gain cannot be cheap easy-leaf evaluation; it must
be clustering** — consecutive q sharing a π value, i.e. *fewer* queries.

**Segment size: already right, for a reason worth knowing.** S rode y from the start and was
never swept. Predicted the counter's O(√S) query dominated, so S ≈ 2570 (~72× smaller) would win
~2.8×. Measured: **S=y is within 0.2–1.4% of optimal**, and shrinking S is catastrophic.

| seg (10¹⁴) | ~query | a·nsegs | total ms | vs best |
|---:|---:|---:|---:|---:|
| 2048 | 12 | 4.4×10⁹ | 13456 | +249.5% |
| 8192 | 23 | 1.1×10⁹ | 6323 | +64.2% |
| 32768 | 46 | 2.8×10⁸ | 4298 | +11.6% |
| **131072** | 91 | 6.9×10⁷ | **3850** | **0.0%** |
| 185660 (=y) | 108 | 4.9×10⁷ | 3904 | +1.4% |

Optima are 32768/65536/131072 at 10¹²/10¹³/10¹⁴ ≈ **0.75y** — not worth a knob; S=y self-scales.
The model was wrong because it priced the `a·z/S` fold at ~1 op. Backing it out of the data, each
(segment, prime) fold step costs **~26 cycles**. First guess: the `((lo+p−1)/p)·p` division, so
"the division pins S high" — replace it with a next-multiple cursor and S should un-pin (~2×).

**Falsified. The cursor is a no-op — 1.00× at essentially every (x, S):**

| seg (10¹⁴) | before | with cursor | ratio | vs best |
|---:|---:|---:|---:|---:|
| 2048 | 13309.8 | 13205.6 | 1.01× | +243.9% |
| 8192 | 6331.9 | 6321.4 | 1.00× | +64.6% |
| 32768 | 4301.2 | 4364.4 | 0.99× | +13.6% |
| **131072** | 3852.3 | **3840.4** | 1.00× | **0.0%** |
| 185660 (=y) | 3913.0 | 3892.4 | 1.01× | +1.4% |

Identical curve, identical optimum. The division's ~30-cycle **latency was already hidden by ILP**
— the prime loop's iterations are independent, so the divider pipelines against everything else.
Latency ≠ cost when there is that much parallelism to hide it under. Cursor reverted (neutral on
time, and it adds another a-sized array to the exact pressure point below).

**What actually pins S is structural.** The ~26 cycles are memory traffic on the *a-sized* arrays:
every segment streams `seg_cnt`, `cur`, `next`, `phi_run` — 4 × 16700 × 8 B ≈ **536 KB at 10¹⁴,
right at this box's 512 KB L2 cliff**. And it is irreducible given the algorithm:
`phi_run[bi] += seg_cnt[bi]` must run for *every* b at *every* segment boundary, because the
running φ per b is exactly the mechanism that lets us never re-sieve [1,lo). So a·nsegs is a hard
floor, and **S ≈ 0.75y is right for a real reason** rather than a tuning accident. S=y stays.

**The counter, revisited — and the two divisions that disagree.** The 2-level Counter's own
docstring justified itself with "kills outnumber queries ~60:1", measured at α=1.5. At α=4 that
**inverted**: 5.4×10⁸ kills vs 1.4×10⁸ leaves is ~4:1 by count, and in *cost* the queries dominate
~12:1 the other way. So the balance was retested with three variants, separating the two effects:

| x | Counter (2-level, ÷) | Counter2P (2-level, shift) | Counter3P (3-level, shift) |
|---|---:|---:|---:|
| 10¹³ | 854.2 ms | 745.8 (1.15×) | 740.6 (**1.15×**) |
| 10¹⁴ | 3989.0 ms | 3667.5 (1.09×) | 3345.2 (**1.19×**) |

Predicted ~20% from the third level (query 3·nwords^(1/3) vs 2·√nwords — 43 vs 108 at 10¹⁴). The
**total** is ~19% at 10¹⁴, but the attribution was wrong: most of it is removing the **division**,
and at 10¹³ the third level adds *nothing* over a 2-level shift. Adopted anyway — 3P is ≥ 2P at
both, and its margin grows with x (1.15 → 1.19), which is the right direction for 10¹⁸.

**Why this division cost when the fold's didn't.** In `kill()`, `cnt[w / wpb] -= 1`: wpb is a
runtime divisor and the quotient is an **address**, so the dependent load-modify-store sits on the
critical path. The fold's `((lo+p−1)/p)·p` measured *free* because its quotient only seeded a loop
whose iterations were independent across primes — ILP swallowed the latency. **Division latency is
free when it feeds independent work, and expensive when it feeds an address you must immediately
touch.** Two opposite results, one rule.

**Parallelism: the leaves are wildly non-uniform.** Kills are uniform by construction (each
element of [1,z] dies once), but leaves are not — measured at 10¹⁴:

| decile of [1,z] | 1 | 2 | 5 | 10 |
|---|---:|---:|---:|---:|
| share of leaves | **99.72%** | 0.11% | 0.02% | 0.01% |

Leaves are dominated by n = p·q with both primes large, so pq clusters near y² and v = x/(pq)
clusters near x/y² — the very bottom of [1,z]. Since leaf queries dominate the work
(~1.4×10⁸ × ~108 ops vs 5.4×10⁸ kills), **decile 1 holds ~90% of total work**: an equal-width
block split would give one thread nearly everything. Blocks must be sized by leaf count, which is
cheap to precompute.

The prefix dependency is *not* an obstacle, and the P₂ fusion does not worsen it — both halves
always had one (S2 needs phi_run[bi], standalone P₂ needed π(lo−1)); fusing merged two identical
prefix problems into one. φ enters every leaf linearly, so a block-and-scan works: each thread
takes a contiguous block, computes relative to phi_run = 0 at its start, and returns local S2 plus
`block_total[bi]` and `mu_sum[bi] = Σ−μ(m)`; then
S2 = Σ_t [local_t + Σ_bi phi_run_at_block_start[t][bi]·mu_sum[t][bi]], with
phi_run_at_block_start[t] = Σ_{t'<t} block_total[t'] — O(nthreads × a) storage (~267k values at
10¹⁴). **P₂ is easier fused than standalone**: π(v) = φ(v,a) − 1 + a makes its per-block
correction a scalar (the count of p handled), not a per-bi vector.

**Bug found by the small-x spot checks.** π(2) returned 2. y = 4·icbrt(2) = 4 **> x**, so
a = π(4) = 2 counts the prime 3 > x, and π(x) = φ + a − 1 − P₂ over-counts. The invariant is
**x^(1/3) ≤ y ≤ √x** — floor is the P₃ bound, ceiling keeps every p_a ≤ √x ≤ x. Always
satisfiable (icbrt ≤ isqrt) and only binds below x = 4096, which is why α=2 masked it. Now
verified exhaustively for every x in [0, 5000].

- Exact at all 14 known values through 10¹⁴ plus 12 small-x / non-power-of-ten spot checks.
- **Exponent 0.68 vs Meissel's 0.86** — the theoretical 2/3, achieved. Crossover ~10¹²; the
  gap only widens.
- **φ footprint 410 KB at 10¹³** vs capped Meissel's 73.8 MB — 180×, and Θ(x^(1/3)).
- Segmenting made S2 *faster* as well as smaller: 2.5× vs the flat [1,z] version at 10¹¹
  (132.9 → 336.4 ms), because the Fenwick went cache-resident.

**Leaf law — and a correction.** An earlier counter used `maxpf` (p = P⁺(m), the *opposite*
convention) and reported Θ(x^(2/3)/ln x). Against the real P⁻ set the ratio to x^(2/3)/ln x
**drifts** (0.805 → 0.474), while the ratio to x^(2/3)/**ln²**x is **flat at ~12.0**, and a²/2
predicts the count almost exactly:

| x | a = π(y) | leaves | a²/2 | leaves/(x^(2/3)/ln²x) |
|---|---:|---:|---:|---:|
| 10⁸ | 125 | 7,815 | 7,813 | 12.31 |
| 10¹¹ | 894 | 403,134 | 399,618 | 12.00 |

So **leaves ≈ π(y)²/2 = Θ(x^(2/3)/ln²x)**, dominated by two-prime leaves n = p·q. At 10¹⁸ that
is ~6.5×10⁹ leaves, not the 42.6×10⁹ the wrong convention predicted.

**Fenwick was the wrong data structure — the traffic is lopsided.** Folding the primes
kills every element of [1,z] exactly once (z ≈ 1.4×10⁹ at 10¹⁴), while the leaves only query
a²/2 ≈ 2.4×10⁷ times — **~60:1**. Fenwick charges O(log S) for *both*, taxing the hot side to
subsidise the cold one. Replaced with a two-level counter (a bit per element + an alive-count
per block of √nwords words): **O(1) kill**, O(√S) query, O(1) segment total. That last one
matters more than it looks — the old code spent a full `prefix(len)` per (segment, b), which
is 1.4×10⁸ tree walks at 10¹⁴, now a single register read.

| x | Fenwick | Counter | speedup |
|---|--------:|--------:|--------:|
| 10¹¹ | 155.2 ms | 93.8 ms | 1.65× |
| 10¹² | 714.2 ms | 424.0 ms | 1.68× |
| 10¹³ | 3408 ms | 1943 ms | 1.75× |
| 10¹⁴ | 16947 ms | 8912 ms | **1.90×** |

The gain *grows* with x, as it must when you drop a log factor off the dominant term — and it
pulled the exponent from 0.68 to **0.66**, essentially the theoretical 2/3. The Fenwick is kept
in the flat [1,z] path as a cross-check: flat (Fenwick) == segmented (Counter) == oracle.

**Where the cost actually is.** Not the leaves — the *sieve*: z = 3.1×10⁸ at 10¹³ against
6.0×10⁶ leaves. So α is a **real lever** here, unlike capped Meissel: z ~ 1/α but leaves ~ α²,
a genuine interior optimum. Balancing them buys only ~ln^(2/3)x, though — extrapolating 16.9 s
at 10¹⁴ by the measured 0.68 gives **~2 h single-threaded at 10¹⁸**, so plain LMO is nowhere
near single-digit seconds. That gap is what Deléglise–Rivat and parallelism must cover.

**P₂'s memory: fixed, and free.** It stored every prime ≤ √x — Θ(√x/ln x), 5 MB at 10¹⁴ and
386 MB at 10¹⁸ — which dominated everything else and gave back part of what segmenting won. But
for segment [lo, hi) the p with x/p ∈ [lo, hi) are exactly **p ∈ (⌊x/hi⌋, ⌊x/lo⌋]**, and as lo
sweeps upward those ranges **tile (y, √x]**: disjoint, contiguous, each provably ≤ seg wide. So
each is sieved on the fly from base primes ≤ **x^(1/4)**, and π(√x) is just counted as we go
(never stored) for the Σ(π(p)−1) closed form.

| | 10¹⁴ | 10¹⁸ |
|---|---:|---:|
| primes ≤ √x (was) | 5.0 MB | 386.0 MB |
| primes ≤ x^(1/4) (now) | 3.1 KB | 24.4 KB |
| **total footprint** | **6.5 → 1.5 MB** | **416.9 → 30.8 MB** |

Costs nothing in time (4101 vs 4134 ms at 10¹⁴ — the on-the-fly sieving is free), and total
memory is now genuinely **Θ(x^(1/3))**, dominated by the μ/lpf tables at O(y).

## Sieve of Atkin vs Eratosthenes: op-count is the wrong metric (atkin.zig)
Atkin toggles a bit per quadratic-form solution (4x²+y², 3x²+y², 3x²−y², mod-12
residues) → O(N/log log N) ops, *fewer* than Eratosthenes' O(N log log N).
Full-array bit sieve, both counting π(N):

| N | Atkin M/s | Eratosthenes M/s |
|---|----------:|-----------------:|
| 10⁶ (in-cache) | 1131 | 630 |
| 10⁸ | 642 | 484 |
| 10⁹ (DRAM) | 200 | 201 |

- **In-cache (10⁶): Atkin wins ~1.8×** — its fewer ops show when compute-bound.
- **Memory-bound (10⁹): they converge (~200)** — DRAM bandwidth erases the
  op-count advantage; both just push bits to RAM.
- But the champion **segmented + wheeled Eratosthenes is ~3.6 G/s — 18× faster**
  than full-array Atkin, and Atkin *can't* segment/wheel to follow (scattered
  quadratic-form writes, not dense resumable strides). The memory optimizations
  win 18×; Atkin structurally can't have them.
- Verdict: Atkin's better op-count is real but visible only when compute-bound;
  at scale it's memory-bound (moot) and crushed by memory-optimized Eratosthenes.
  For π(x) at scale it's doubly moot — the cost is combinatorial (φ, P₂), not
  enumeration. Op-count is the wrong metric on this hardware.

## Architecture notes
- Three orthogonal comptime axes (wheel × store × traversal-via-seg-size) +
  a runtime interval. Segmentation = repeated **range sieve** over [lo,hi);
  whole-array is the single-interval degenerate case.
- Sieving primes computed **once** (base sieve to ceil(√N)), reused across
  segments; completeness limit L carried and asserted **L² ≥ top** (active in
  ReleaseSafe/Debug). A partial/P_n source would skip the check.
- Counting folds into the segmented sieve (segments discarded); count() ~0.
- Coordinate/value type parameterized (`Int` = u64/u128) for ranges past 2⁶⁴;
  store word parameterized (`BitPacked(Word)`); both orthogonal.
