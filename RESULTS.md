# Sieve benchmark ladder

Fixed benchmark: **N = 10⁹**, π(N) = 50,847,534. Best-of-3, sieve timed
separately from count. Machine: AVX2 (no AVX-512), L1d 32 KiB / L2 512 KiB /
L3 16 MiB per core, 28 GiB RAM. Build: `zig build-exe -O ReleaseFast -mcpu=native`.

Three orthogonal axes: **traversal** (whole-array / segmented / segmented-odds)
× **store** (flat []bool / hand-rolled []u64 / std.DynamicBitSet) × coordinate
map (all / odds / …). Segment knob = **store BYTES** (the cache-critical
quantity — see methodology note), held at 32 KiB (L1d) for all segmented rows.

| traversal | store | footprint | ints/seg | rate (M ints/s) |
|-----------|-------|-----------|----------|-----------------|
| whole-array | flat `[]bool` | 954 MiB | — | 141 |
| whole-array | `std.DynamicBitSet` | 119 MiB | — | 199 |
| whole-array | hand-rolled `[]u64` | 119 MiB | — | 201 |
| segmented (all) | hand-rolled `[]u64` | 32 KiB | 262144 | 471 |
| segmented (all) | `std.DynamicBitSet` | 32 KiB | 262144 | 510 |
| segmented (all) | flat `[]bool` | 32 KiB | 32768 | **1140** |
| segmented (odds) | hand-rolled `[]u64` | 32 KiB | 524288 | 1435 |
| segmented (odds) | flat `[]bool` | 32 KiB | 65536 | **2618** |

(counting folded into segmented sieves; stable across reps this run.)

## Finding 1 — the representation optimum FLIPS with the traversal
- **Whole-array = memory-bound.** 954 MiB `[]bool` blows all caches → every
  strike is a DRAM access → **footprint dominates** → bit-packing wins (141→201).
- **Segmented = cache-bound, now byte-matched at 32 KiB L1.** Cache footprint is
  EQUAL across stores, so this is the clean representation comparison: a `[]bool`
  strike is one **byte store**; a `[]u64` strike is a **read-modify-write**
  (load word, shift/mask/or). Byte store wins → **`[]bool` is ~2.4× faster**
  (1140 vs 471) — winning *harder* than in the earlier (byte-unmatched) run, and
  despite covering 8× fewer integers/segment (8× more segment overhead). Strike
  cost dominates once you're out of DRAM.
- **Implication for SIMD:** bit-packing only earns its keep in-cache once the RMW
  is amortized over many bits (SIMD / wheel word-strikes). Scalar bit-by-bit
  striking loses to bytes. This is *why* real sieves are bit-packed AND wheeled.

## Finding 2 — odds-only (wheel-2) adds ~2–3×
- Drops prime 2 (the densest striker) and halves the array; each segment covers
  2× the integers. `[]bool` odds reaches **2618 M ints/s**.

## Methodology note (corrected mid-session)
- The segment knob must be **store bytes**, not flags/slots or integers-per-
  segment. Bytes is the cache-critical quantity, and *only* bytes holds cache
  footprint constant across stores: `[]bool` is 1 byte/flag, the bit stores are
  1 bit/flag (8× denser). Tuning by integers or flags silently gives `[]bool` an
  8× larger store, confounding the cache variable with the representation.
- Aside: the naive `[]u64` vs `DynamicBitSet` strike loops are equal (~5100 ms
  whole-array); earlier apparent gaps were CPU contention, since confirmed by an
  isolated microbench.

## TODO / open threads
- Runtime-parameterize N (scaling study: π(N)/N and rate vs N).
- Segment-BYTE sweep (16/32/64/128/256 KiB…) on an idle box — find the L1/L2/L3
  crossover per store.
- Wheel-30 (3rd coordinate-map point) → then extract the coordinate abstraction.
