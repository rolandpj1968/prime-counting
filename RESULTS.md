# Empirical results — split into two docs

This file has been split by topic:

- **[SIEVING.md](SIEVING.md)** — the Sieve of Eratosthenes (and Atkin): wheels,
  segmentation, the cache hierarchy, the bucket sieve, u128 range counting, and
  the sieve as a lens on the number theory (PNT / Li / RH).
- **[COMBINATORIAL.md](COMBINATORIAL.md)** — combinatorial π(x): the Meissel–Lehmer
  baseline, the LMO / Deléglise–Rivat implementation at the 2/3 exponent in
  Θ(x^(1/3)) memory, and **Gourdon's decomposition** on the same O(1)-kill counter,
  which reaches π(10²²) in 4.74 h on six cores. Optimisation ladders, dead ends
  (both DR headline optimisations; batching the counter bookkeeping), and why the
  parallel path is leaf-side bandwidth-bound.

See [README.md](README.md) for framing and references.
