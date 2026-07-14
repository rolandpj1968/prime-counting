# Sieve benchmark ladder

Fixed benchmark: **N = 10⁹**, π(N) = 50,847,534. Best-of-3, sieve timed
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
Two-tier PrimeSource (bucket_sieve.zig): SMALL primes (min delta ≤ seg_slots,
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

## Architecture notes
- Three orthogonal comptime axes (wheel × store × traversal-via-seg-size) +
  a runtime interval. Segmentation = repeated **range sieve** over [lo,hi);
  whole-array is the single-interval degenerate case.
- Sieving primes computed **once** (base sieve to ceil(√N)), reused across
  segments; completeness limit L carried and asserted **L² ≥ top** (active in
  ReleaseSafe/Debug). A partial/P_n source would skip the check.
- Counting folds into the segmented sieve (segments discarded); count() ~0.

## TODO / open threads
- Runtime N (scaling study: π(N)/N, rate vs N).
- Segment-byte sweep on an idle box (find L1/L2/L3 crossover per store).
- Bucket `PrimeSource` for large primes (hybrid with cursor for small) — the
  co-habitation payoff; slots into the range-sieve seam without touching the core.
- Numbers still contention-affected; re-verify absolutes on a quiet machine.
