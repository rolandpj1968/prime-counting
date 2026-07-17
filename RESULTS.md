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
  ~0.8 exponent is unchanged; **Stage B** (special-leaf sieve) targets O(x^(1/3))
  memory and the 2/3 exponent.
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
| 10¹⁴ | 3,204,941,750,802 | 4.13 s | 41.3 s | **10.0×** | 0.657 |

The 10¹⁴ column, step by step — each is a separate committed measurement:

| change | 10¹⁴ | vs prev |
|---|---:|---:|
| Fenwick, α=1.5 | 16.95 s | — |
| → block Counter (O(1) kill) | 8.91 s | 1.90× |
| → α = 2 | 7.98 s | 1.12× |
| → m-walk split at √y | 6.07 s | 1.31× |
| → P₂ fused onto the counter | 4.97 s | 1.22× |
| → α = 4 | **4.13 s** | 1.20× |

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
The model was wrong because it priced the `a·z/S` fold at ~1 op. Back it out of the data: at
S=32768 the query cost halves (91→46), which should save ~6×10⁹ ops, yet it runs *slower*. That
only balances if each (segment, prime) fold step costs **~30 cycles** — i.e. the
`((lo+p−1)/p)·p` division. **That division is what pins S high**, and it is still ~12% of runtime
at the optimum. Replacing it with a per-prime next-multiple cursor should un-pin S (≈32768 viable,
queries halve) — model says ~2×. A lock on a knob, not a 12% cleanup.

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

**Open: P₂ is now the memory bottleneck.** It stores the primes ≤ √x — 5 MB at 10¹⁴ vs φ's
890 KB — so total memory is Θ(√x/ln x), not Θ(x^(1/3)). Fix: each v-segment needs p ∈
(x/hi, x/lo]; those ranges are disjoint, each ~S wide, and their union is exactly (y, √x], so
sieve them on the fly instead of storing.

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
