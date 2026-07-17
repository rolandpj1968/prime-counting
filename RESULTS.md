# Empirical results — split into two docs

This file has been split by topic:

- **[SIEVING.md](SIEVING.md)** — the Sieve of Eratosthenes (and Atkin): wheels,
  segmentation, the cache hierarchy, the bucket sieve, u128 range counting, and
  the sieve as a lens on the number theory (PNT / Li / RH).
- **[COMBINATORIAL.md](COMBINATORIAL.md)** — combinatorial π(x): the Meissel–Lehmer
  baseline and the LMO / Deléglise–Rivat implementation that reaches π(10¹⁹) at the
  2/3 exponent in Θ(x^(1/3)) memory. Optimisation ladder, dead ends (including both
  DR headline optimisations), and detailed per-optimisation analysis.

See [README.md](README.md) for framing and references.
