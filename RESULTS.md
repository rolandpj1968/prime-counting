# Sieve benchmark ladder

Fixed benchmark: **N = 10⁹**, π(N) = 50,847,534. Best-of-3, sieve timed
separately from count. Machine: AVX2 (no AVX-512), 28 GiB RAM.
Build: `zig build-exe -O ReleaseFast -mcpu=native`.

Each rung flips exactly one knob vs. the previous, so every delta is attributable.

Axes: **traversal** (whole-array / segmented) × **store** (bool / []u64 / bitset).
Segment span fixed at 262144 integers (store decides its own byte size).

| traversal | store | footprint / seg | rate (M ints/s) |
|-----------|-------|-----------------|-----------------|
| whole-array | flat `[]bool` | 954 MiB | 148 |
| whole-array | `std.DynamicBitSet` | 119 MiB | 201 |
| whole-array | hand-rolled `[]u64` | 119 MiB | 202 |
| segmented | hand-rolled `[]u64` | 32 KiB (L1d) | 581 |
| segmented | `std.DynamicBitSet` | 32 KiB (L1d) | 578 |
| segmented | flat `[]bool` | 256 KiB (L2) | **1063** |

(counting folded into the segmented sieves; count() ~0. Contention eased —
these are stable across reps, but still re-verify on a fully idle box.)

## The headline finding: the best representation FLIPS with the traversal
- **Whole-array = memory-bound.** 954 MiB bool blows all caches → DRAM-bound →
  footprint dominates → **bit-packing wins** (148 → 202).
- **Segmented = cache-bound.** Both segments fit cache (no DRAM), so footprint
  stops mattering and **per-strike cost** takes over. `[]bool` strike is one
  byte store; `[]u64` strike is a load-modify-write (shift/mask/or). Byte store
  wins → **`[]bool` segment is ~1.8× faster** despite being 8× bigger and
  spilling L1→L2. The L1-vs-L2 latency gap is too small to overcome the RMW.
- **Implication for the SIMD rung:** bit-packing only earns its keep in-cache if
  the per-strike RMW is amortized over many bits — i.e. SIMD/wheel word-strikes.
  Naive bit-by-bit striking in cache loses to bytes. This is *why* real sieves
  (primesieve) are bit-packed AND wheeled/SIMD, never bit-packed + scalar.

(Representative run; numbers wobble ±5% run-to-run. ReleaseFast.)

> ⚠️ **Provisional** — these were taken while the box had other load (large Ruby
> processes contending). Re-run the whole ladder on a quiet machine before
> trusting absolute numbers or cross-impl deltas. An isolated microbench showed
> hand-rolled `[]u64` and `DynamicBitSet` strike loops are actually equal
> (~5100 ms); the 181-vs-200 split above was contention noise, not a real gap.

## Notes
- **flat `[]bool`**: memory-bound. 954 MiB buffer blows all caches; every strike
  is effectively a DRAM access. First-run page-fault tax discarded by best-of-3.
- **std.DynamicBitSet**: 8× smaller footprint → **1.35× faster sieve** (149→200)
  despite per-bit bounds-checked `set()`/`isSet()`. Confirms we're memory-bound:
  shrinking the working set helps even when we *add* per-strike compute.
  Moves TWO variables vs flat (byte→bit AND raw→abstraction), so not a clean
  single-knob delta — the hand-rolled `[]u64` rung isolates the bit-packing win.
- **The count blew up 24×**: 151 ms → 6.2 ms. `count()` is a streaming
  `@popCount` over 119 MiB of words; flat's count is a 1e9-element per-byte
  branch over 954 MiB. Packing helps striking a little, counting enormously.
