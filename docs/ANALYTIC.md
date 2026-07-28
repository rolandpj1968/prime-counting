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
T ≈ 3.06×10¹⁰ at ±2⁻¹⁰² — rigorous and Turing-verified complete.

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

## The ladder ahead

1. **Unwind ψ → θ → π** (Möbius/partial summation — elementary bookkeeping,
   verified at every rung against the combinatorial table).
2. **Precision ladder.** f64 → f128 (native in Zig, softfloat) or Dekker
   double-double when term counts × cancellation outgrow 53 bits.
3. **Scale + distribute.** Zero-range partial sums are *purely additive*
   fragments — the fleet architecture (plan / controller / merge / tiling proof,
   see [COMBINATORIAL.md](COMBINATORIAL.md)) transfers verbatim. Bulk zeros from
   LMFDB/Platt.
4. **Rigor last.** Ball arithmetic (midpoint + radius), *not* directed rounding —
   LLVM reorders float ops assuming round-to-nearest, so fesetround-style
   interval arithmetic is unreliable; balls need no rounding modes. Certified
   tail bounds close the argument.

Ceiling of the method with existing zeros: Platt's T ≈ 3×10¹⁰ supports
x ≈ 10²⁴⁻²⁵ (correction window ~x/T·polylog). Beyond that, computing new zeros
(Odlyzko–Schönhage) becomes its own sub-project — also fleet-friendly.

See [README.md](../README.md) for framing and references.
