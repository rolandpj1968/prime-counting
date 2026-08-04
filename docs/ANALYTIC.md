# Analytic π(x): the explicit formula

The third arc. Below the combinatorial line ([COMBINATORIAL.md](COMBINATORIAL.md))
and the sieving foundation ([SIEVING.md](SIEVING.md)) sits a different animal
entirely: computing π(x) from the **nontrivial zeros of ζ**, via the
Riemann–von Mangoldt explicit formula — the method of Lagarias–Odlyzko (1987),
practiced at scale by Platt (π(10²⁴), 2012) and Franke–Kleinjung–Büthe–Jost
(π(10²⁵)). Zeros in, prime counts out; the primes are never enumerated, sieved,
or φ-recursed — they interfere into existence out of a sum of cosines.

As with every arc in this repo, the build is **empirical and naive-first**: start
with the textbook formula, measure where it breaks, and let each failure dictate
the next structure. The referee at every rung is the combinatorial ladder — every
known value 10 → 10²⁴ computed and verified by `gourdon.zig`.

Code: `src/analytic/` (built like the other drivers:
`zig build-exe -O ReleaseFast -mcpu=native --dep rs -Mroot=explicit.zig
-Mrs=../rangesieve.zig`).

## Inputs: the zeros

Two Odlyzko tables (`zeros/`, gitignored — fetch from
[dtc.umn.edu/~odlyzko/zeta_tables](https://www.dtc.umn.edu/~odlyzko/zeta_tables/)):
the first 10⁵ zeros (T ≈ 7.5×10⁴) and the first 2,001,052 zeros
(T ≈ 1.13×10⁶), imaginary parts to **±4×10⁻⁹**.

For later rungs: LMFDB carries Platt's ~1.038×10¹¹ zeros to height
T ≈ 3.06×10¹⁰ at ±2⁻¹⁰² — rigorous and Turing-verified complete. The bulk
format is decoded and validated (`lmfdb2txt.py`): files `zeros_<t>.dat` at
beta.lmfdb.org/riemann-zeta-zeros/, blocks of
`[t0:f64][t1:f64][N(t0):u64][N(t1):u64]` + 13-byte zeros (u64+u32+u8 =
104-bit **delta** to the previous zero in units 2⁻¹⁰¹; accumulate in exact
integer arithmetic — a float cumulative sum random-walks ~10⁻⁶ over millions
of zeros). Our first prefix: 4,826,908 zeros to T = 2,546,000 (~62 MB);
the 2,001,052-zero overlap with Odlyzko agrees to worst 3.0×10⁻⁹, mutually
validating both tables and both decoders.

**Zero precision is not the binding constraint.** The sensitivity of the zero sum
to a perturbation δγ is ~√x·ln x·δγ per zero; over the 2M-zero table the RMS
contribution at x = 10¹² is ~0.017 against the ±0.5 needed to pin an integer, and
even the paranoid all-coherent worst case only reaches 0.77. Double precision and
9-digit zeros carry the early rungs comfortably. The binding constraint is zero
**height** — see below.

## Rung 1: the naive formula, watched failing (`explicit.zig`)

The unsmoothed formula for the Chebyshev function ψ(x) = Σ_{pᵏ≤x} ln p:

```
ψ(x) = x − Σ_ρ x^ρ/ρ − ln 2π − ½ln(1 − x⁻²)
```

with each conjugate zero pair ρ = ½ ± iγ contributing (L = ln x):

```
2√x · (½·cos γL + γ·sin γL) / (¼ + γ²)
```

Truncate the sum at the first N zeros and compare against the **direct referee**:
a sieve summing ln p over primes and prime powers (Neumaier–Kahan on both sides —
the naive-summation error bound n·ε·Σ|terms| is ~11 at x = 10⁹, larger than the
thing being measured). Measured, err = ψ_T − ψ_direct, and `rms_pred` is the
dropped-tail model √x·√((ln(T/2π)+1)/(2πT)):

| x | N zeros | T | err | rms_pred |
|---|--------:|--:|----:|---------:|
| 10⁶ | 100 | 236 | +51.1 | 55.8 |
| 10⁶ | 10⁵ | 7.5×10⁴ | +3.2 | 4.7 |
| 10⁸ | 100 | 236 | −860.4 | 558.1 |
| 10⁸ | 10⁵ | 7.5×10⁴ | +29.3 | 47.0 |
| 10⁸ | 2.0×10⁶ | 1.13×10⁶ | +8.6 | 13.6 |
| 10⁹ | 10⁵ | 7.5×10⁴ | +40.5 | 148.5 |
| 10⁹ | 2.0×10⁶ | 1.13×10⁶ | +51.3 | 42.9 |

Three empirical facts:

1. **Convergence in T is √-slow.** At x = 10⁸, a 20,000× increase in zero count
   (100 → 2M) buys a factor-100 error reduction — exactly the 1/√T law. The error
   is Gibbs oscillation around the jumps of ψ, and the tail-RMS model tracks the
   measured error within a factor ~2–3 at every checkpoint.
2. **At fixed T the error grows like √x.** Every decade of x demands a decade of
   T just to stand still.
3. **The kill shot.** Setting err ≈ 0.5 requires T of order x·ln x. At x = 10¹²
   that is ~6×10¹² zeros — 60× more than exist in any database on Earth. The
   sharp cutoff is not slow; it is *structurally incapable* of pinning an integer
   at any interesting x. This is why every serious analytic computation smooths.

Positive result banked alongside the failure: the arithmetic side (a sieve
summing logs of primes) and the analytic side (a sum of cosines over ζ zeros)
agree to within the predicted oscillation at three scales — formula, pair-term
algebra, constants, and trivial-zero tail all certified. And f64 holds: at
x = 10⁹ the phase γL reaches 2×10⁷ rad, whose argument-reduction error (~5×10⁻⁹)
is still below the tables' own ±4×10⁻⁹.

## Rung 2: Gaussian smoothing (`smooth.zig`)

Kernel provenance, verified against the sources (both in `literature/`):
Lagarias–Odlyzko 1987 leave the Mellin pair φ̂, φ "suitable" but unspecified;
**Galway's thesis** (UIUC 2004) proposed the Gaussian pair
φ̂(s) = (x^s/s)·e^(λ²s²/2), φ(t) = ½erfc(log(t/x)/(√2λ)) and showed it
near-optimal by the uncertainty principle (the Gaussian is its own Fourier
transform — fastest simultaneous decay in γ and in n); **Platt** (arXiv:
1203.5712) implemented it rigorously for unconditional π(10²⁴). FKBJ used a
different kernel family (compactly-supported Fourier transform), announced
π(10²⁴) in 2010 *contingent on RH* — the kernel choice is a live fork, not a
convention, and the conditionality lives in how the tail beyond the last zero
is bounded. (Our β = ½ pair term is verified fact at tabulated heights, not
hypothesis.)

Applied to ψ, the transform is exact and simple: each pair term gains damping
e^(λ²(¼−γ²)/2) *and a phase shift* θ = γ(L + λ²/2), and the x-term becomes
x·e^(λ²/2) (a 0.02 correction at 10⁹ — not optional). Direct side: the same
sieve referee, now summing Λ(n)·½erfc(log(n/x)/(√2λ)) — the weight differs
from the sharp step only on a window of x·2√2·c·λ integers (λ = c/T, c = 7.5,
so the dropped-tail factor is e^(−c²/2) ≈ 10⁻¹²). log(n/x) via `log1p` on the
exact integer offset — no cancellation; erfc from libc.

| x | zeros | naive err (rung 1) | smoothed err | window |
|---|------:|-------------------:|-------------:|--------|
| 10⁶ | 10⁵ | +3.2 | −3×10⁻⁷ | 2,128 ints, 152 primes |
| 10⁸ | 2×10⁶ | +8.6 | **+9×10⁻⁶** | 14,053 ints, 776 primes |
| 10⁹ | 2×10⁶ | +51.3 | **+5×10⁻⁵** | 140,490 ints, 6,763 primes |

The Gibbs oscillation is annihilated: four to five orders below the ±0.5
integer threshold, at the price of erfc-weighting a few thousand primes near x.
The checkpoint tables show the mechanism — at N = 10⁶ (T = 6×10⁵) the error is
already ~10⁻³; the final million zeros arrive pre-damped to ~10⁻⁴ amplitude.
The λ knob works as the theory says: zero count needed ~ 8/λ, window width
~ x·λ, product fixed — the uncertainty tradeoff made tangible.

## Rung 3: an integer π(x) from the zeros (`pistar.zig`)

A one-point ψ(x) cannot yield π(x) (∫dψ/ln t needs the whole range), so this
rung switches central object to Riemann's **π*(x) = Σ_{pᵐ≤x} 1/m**, from
log ζ rather than ζ′/ζ (Platt eq. 3.1):

```
π*(x) = li(x) − Σ_ρ li(x^ρ) − ln 2 + ∫ₓ^∞ dt/(t(t²−1)ln t)
```

smoothed with the same Gaussian pair. The per-zero term becomes a li-type
integral F(ρ) = ∫ x^z e^(λ²z²/2)/z dz along a horizontal path — evaluated as
e^(λ²ρ²/2)·x^ρ·[S/(ρL) − λ²/L²], where S is the asymptotic li series
Σ k!/(ρL)^k (|ρL| ≥ 195 everywhere here, so K = 13 terms reach ~10⁻²⁵), and
−λ²/L² is the first-order kernel correction, written in its cancelled form
(the naive (1 − λ²ρ²)·li form loses precision at the damped tail). The pole
term gains λ²x(L−1)/(2L²) — 0.77 at 10¹², not optional. li(x) = Ei(ln x) by
the convergent all-positive series; the trivial-zero integral (~1/(2x²L)) is
dropped. Sharp π* = smoothed zero sum + Σ_window (1/m)(χ − φ) over an exactly
sieved window; then Möbius with *exact* tiny π*(x^(1/m)) from direct sieves:
π(x) = Σ_m (μ(m)/m)·π*(x^(1/m)). Round; report the margin; compare published.

| x | π(x) from zeros | pre-round error | margin | window | time |
|---|---|---:|---:|---|---:|
| 10⁶ | 78,498 | ~10⁻⁴ | 0.5000 | 2.1×10³ ints | 0.0 s |
| 10⁸ | 5,761,455 | ~10⁻⁶ | 0.5000 | 1.4×10⁴ | 0.4 s |
| 10⁹ | 50,847,534 | +2×10⁻⁶ | 0.5000 | 1.4×10⁵ | 0.4 s |
| 10¹⁰ | 455,052,511 | +3×10⁻⁶ | 0.5000 | 1.4×10⁶ | 0.3 s |
| 10¹¹ | 4,118,054,813 | −1.5×10⁻⁵ | 0.5000 | 1.4×10⁷ | 0.5 s |
| 10¹² | 37,607,912,018 | −4.6×10⁻⁵ | 0.5000 | 1.4×10⁸ | 1.8 s |
| 10¹³ | 346,065,536,839 | −5×10⁻⁷ | 0.5000 | 1.4×10⁹ | 16.7 s |
| 10¹⁴ | 3,204,941,750,802 | +3×10⁻³ | 0.4971 | 6.2×10⁹ (LMFDB prefix) | 83 s |

**MATCH against the published value at every row.** Beyond 10⁹ no full-sieve
referee exists — the zeros carry everything except the ~10⁻⁵-relative-width
window.

### The 10¹³ bug: a margin-column success story

The first run at 10¹³ showed error +0.0714 (margin 0.4286), which we
initially — plausibly, wrongly — attributed to the Odlyzko tables' ±4×10⁻⁹
zero precision (budgeted RMS ≈ 0.06 there). The margin column then earned its
keep. Re-running with the LMFDB prefix — different table, different T,
different λ, different window — reproduced the error to 10⁻⁴: **two
independent zero tables agreeing on the same error means machinery, not
noise** (and its λ-independence acquitted the kernel analytics wholesale,
since the kernel corrections differed 5× between runs). A 45-digit reference
acquitted li(x) (f64 off by only 1.7×10⁻⁴). The pipeline was then split at
π*: computing exact π*(10¹³) = Σ π(x^(1/m))/m via `gourdon.zig` at every root
matched the analytic π* to 2×10⁻⁵ — the entire zero side was *correct*, so
the bug had to be in the Möbius unwind. It was: the hardcoded μ table stopped
at m = 40, and 2⁴¹⁻⁴³ ≤ 10¹³ — three dropped terms summing to
1/41 + 1/42 + 1/43 = **0.071456**. Observed error: +0.0714. Threshold
behavior explained too: 2⁴¹ = 2.2×10¹², so 10¹² and below were clean. Fix:
compute μ(m) by trial factorization to m ≤ 64. Margins after: 0.5000 at 10¹³
(the Odlyzko run lands within 5×10⁻⁷ of the integer).

Lessons banked: (1) margins catch bugs that MATCH alone would forgive at
smaller x; (2) table-independence of an error is the tell for a shared-
machinery systematic; (3) exact-π*-split is the right bisection knife, and
`gourdon.zig` referees it at any x we can reach combinatorially.

## Rung 4: parallel, segmented, f128 aggregates — to 10¹⁷ and beyond

Three walls fell in one evening, each diagnosed by the same epistemics.

**The window array is gone.** The window sum is a streaming reduction, so it
segments: an atomic dispenser hands 16M-integer strips to workers, each
sieved against shared base primes and folded into a per-worker Kahan.
Memory is O(threads·strip + π(√x)) *regardless of x* — 10¹⁴ went from
6.2 GB / 83 s to 151 MB / 24.9 s (t = 6), identical output. The parallel
geometry deliberately mirrors the eventual distribution seam: zero ranges
and window segments are both additive fragments (a future fragment is a
(partial, error-radius) pair — the bounds ride the same seam as the work).

**The f64 li(x) systematic.** At 10¹⁵ both zero tables agreed on err −0.027
— the table-independence tell again. A 45-digit referee convicted li(x):
ε·e^L/L series noise plus ln(x)-in-f64 argument error (x/L·4×10⁻¹⁵ ≈ 0.1
at 10¹⁵), plus the f64 ulp of the ~3×10¹³ aggregates themselves (~4×10⁻³,
reaching ±0.5 by ~10¹⁷ — a hard representation wall). Fix: Ei and all large
aggregates in f128, rounded once at the end.

**Zig's `@log(f128)` is silently 53-bit.** The first f128 attempt landed
−0.034: `compiler_rt`'s `logq` is literally `return log(@floatCast(a))` — an
f64 stub (as are `expq`/`exp2q`, with TODO comments; `sinq`/`cosq`/`tanq`
are real 113-bit ports). Upstream archaeology: ziglang/zig#4026 ("implement
all the math functions for all the floating point types") was closed in 2022
against PR #11532 — a *reorg* whose own description says "functions with
missing implementations call other functions and have TODO comments" — so
the accuracy work lost its tracking issue and still ships stubbed in 0.16.
Fix here: `ln128()` from correctly-rounded f128 *arithmetic* only (exact
power-of-two reduction + atanh series, ln 2 as a 34-digit constant) — the
softfloat ± × ÷ are trustworthy; the transcendentals are not.

Results on the 251M-zero LMFDB table (T = 1.01×10⁸, converted with N(t)
continuity validated block-by-block; one md5-canonical trailer block per
some files skipped — decodes non-monotone, not zero data):

| x | π(x) from zeros | err | margin | window | time (t=6) |
|---|---|---:|---:|---|---:|
| 10¹⁵ | 29,844,570,422,669 | +2.0×10⁻⁵ | 0.5000 | 1.6×10⁹ | 27 s |
| 10¹⁶ | 279,238,341,033,925 | −2.5×10⁻⁵ | 0.5000 | 1.6×10¹⁰ | 94 s |
| 10¹⁷ | 2,623,557,157,654,233 | −1.2×10⁻⁴ | 0.4999 | 1.6×10¹¹ | 747 s |
| 10¹⁸ | 24,739,954,287,740,860 | +5.8×10⁻⁴ | 0.4994 | 1.6×10¹² | 2.17 h |

**MATCH at every row.** The noise floor sits at ~10⁻⁴⁻⁵ (zero-sum f64
texture); the window sieve is ≥97% of runtime and ∝ x — exactly the
complexity law's prediction (T = 10⁸ balances at x ≈ T² = 10¹⁶), and
exactly the side that distributes with ~zero data dependency. Remaining
walls, in expected order: window cost (fleet-shaped), then f64 phase γL and
text-f64 zeros (~γ·2⁻⁵² — the double-double phase + hi/lo zeros rung, with
the LMFDB 2⁻¹⁰¹ format supplying the pairs natively).

## Rung 5: the kernel seam, and Gaussian vs Logan measured

The Mellin pair is now a comptime plug (`init / lo,hi / phi / zeroTerm /
cutoff / poleCorr`), and plug #2 is **Büthe's Logan-function kernel**
(arXiv:1410.7008, Thms 3.1 + 4.1) — the FKBJ band-limited line. Structure:
the prime-side window is *exactly* [xe^(−ε), xe^(ε)] (ε = c/T; the
convolution bump is a compactly-supported Bessel I₀, its μ/ν integrals
tabulated once per run), while the zero side uses every tabulated zero with
weight ℓ_c(εγ) and per-zero Ei₁/Ei₂/Ei₃ asymptotic series; the tail beyond T
is ~e^(−c), so c tracks ½ln x + O(loglog) — a **linear** price where the
Gaussian pays **quadratically** (window ∝ c/T vs ∝ c²/T).

Head-to-head on the 251M-zero table, both kernels MATCHing published π(x):

| x | kernel | window ints | window time | total |
|---|---|---:|---:|---:|
| 10¹⁶ | Gaussian (c=7.5) | 1.57×10¹⁰ | 68.8 s | 93.8 s |
| 10¹⁶ | Logan (c=27.4) | 5.4×10⁹ | 23.6 s | 51.1 s |

**~3× narrower window at equal table and accuracy — Büthe's efficiency claim
([4]) confirmed empirically.** 10¹³ drops to 1.2 s; projected 10¹⁸ ~45 min
(was 2.17 h). The shared 16 s text-parse of the zeros table is now a visible
overhead (binary format queued). Next plugs: Beurling–Selberg
majorant/minorant sandwich (certified interval output — rigor's kernel), and
the Slepian/prolate experiment (the true extremal for fixed-T concentration;
apparently unexplored for π(x)).

### The 10¹⁷ bug: a prime hiding in a half-ulp

Logan's first run at 10¹⁷ landed at **published − 1**, with the analytic
residual a clean +6×10⁻⁵ — integer-exact wrongness, the signature of a
counting bug, not an analytic one. The diagnosis chain, each step killing a
hypothesis class:

1. **Exact π\*(10¹⁷) referee** (combinatorial `pi` at every root, rational
   arithmetic): π\* = …114.798822, so the smooth+window side was short by
   exactly 1.000000 and the Möbius unwind was acquitted digit-for-digit.
2. **c = 31 rerun**: same deficit at a different ε, window, and zero
   weighting — kills tail/truncation and window-edge theories.
3. **Post-seam Gaussian control at 10¹⁷**: MATCH — acquits the entire shared
   pipeline (segmented sieve, χ, unwind, f128 aggregates). Logan-specific.
4. **x-bisection** at 3, 5.5, 7.3, 8.6×10¹⁶: all MATCH. Not monotone in any
   natural threshold (7.3×10¹⁶ > 2⁵⁶ already) — which was the tell.

Cause: `bigM` selected its μ-branch by `sign(y)` (f64) while χ cuts on the
integer `n <= x`. Above 2⁵³, every integer in (x − ulp/2, x] *rounds to
y == 0* and took the μ(0⁺) = +½ branch — φ off by exactly 1/λ ≈ 1 for any
**prime** in that band. The bands at all five MATCHing x's are prime-free;
the band at 10¹⁷ contains exactly one prime, **99999999999999997** — the
missing unit, found by Miller–Rabin, not by staring at the sieve. (The 10¹⁸
band holds two primes; the pre-fix code would have been low by exactly 2.)
The Gaussian is immune: its φ is continuous at u = 0, so the same rounding
costs ~10⁻¹⁰. The fix is one line — branch on `n <= x`, same as χ.

The lesson generalizes: **when a smoothed weight is split into a
discontinuous-exact piece plus a smooth-float piece, both must cut on the
same comparison** — any float-sign branch owns a half-ulp band of integers
that the exact side assigns to the other shore. Third kill for the margin
column: the error was invisible to MATCH at four nearby probe points and
pinned to one prime by the fractional residue.

## The ladder ahead

1. **Precision.** LMFDB/Platt zeros (±2⁻¹⁰²) to kill the table noise;
   f64 → f128 (native in Zig, softfloat) or Dekker double-double when term
   counts × cancellation outgrow 53 bits.
2. **Scale + distribute.** Zero-range partial sums are *purely additive*
   fragments — the fleet architecture (plan / controller / merge / tiling proof,
   see [COMBINATORIAL.md](COMBINATORIAL.md)) transfers verbatim. Bulk zeros from
   LMFDB/Platt.
3. **Rigor last.** Ball arithmetic (midpoint + radius), *not* directed rounding —
   LLVM reorders float ops assuming round-to-nearest, so fesetround-style
   interval arithmetic is unreliable; balls need no rounding modes. Certified
   tail bounds close the argument.

Ceiling of the method with existing zeros: Platt's T ≈ 3×10¹⁰ supports
x ≈ 10²⁴⁻²⁵ (correction window ~x/T·polylog). Beyond that, computing new zeros
(Odlyzko–Schönhage) becomes its own sub-project — also fleet-friendly.

See [README.md](../README.md) for framing and references.
