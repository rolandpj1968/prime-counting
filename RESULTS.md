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
