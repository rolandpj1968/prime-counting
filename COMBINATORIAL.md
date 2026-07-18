# Combinatorial π(x): Meissel–Lehmer → LMO / Deléglise–Rivat

Counting primes *without enumerating them*, via π(x) = φ(x,a) + a − 1 − P₂(x,a),
a = π(y). Two implementations: `meissel.zig` (the classical recursion, baseline and
correctness oracle) and `lmo.zig` (Lagarias–Miller–Odlyzko / Deléglise–Rivat, the
2/3 exponent in Θ(x^(1/3)) memory). See [SIEVING.md](SIEVING.md) for the sieve
foundation and [README.md](README.md) for framing.

Machine: AVX2 (no AVX-512), L1d 32 KiB / L2 512 KiB / L3 16 MiB per core, 28 GiB RAM,
single-threaded, `zig 0.16 -O ReleaseFast -mcpu=native`. Everything exact and verified
against every known value 10 → 10¹⁹ plus exhaustively over [0, 5000], with S₂ and P₂
each cross-checked differentially against an independent reference (`specialS2Segmented`
/ `p2Segmented`) that evaluates every leaf off the sieve.

## Results

| x | π(x) | time | exponent | peak RSS |
|---|------|-----:|:--------:|---------:|
| 10¹⁴ | 3,204,941,750,802 | 1.16 s | — | 2.7 MB |
| 10¹⁵ | 29,844,570,422,669 | 5.23 s | 0.65 | ~5 MB |
| 10¹⁶ | 279,238,341,033,925 | 23.1 s | 0.66 | 11 MB |
| 10¹⁷ | 2,623,557,157,654,233 | 1.72 min | 0.66 | ~24 MB |
| 10¹⁸ | 24,739,954,287,740,860 | 8.3 min | 0.66 | ~52 MB |
| 10¹⁹ | 234,057,667,276,344,607 | 42.4 min | 0.66 | ~109 MB |
| 10²⁰ | 2,220,819,602,560,918,840 | 3.44 h | 0.69 | ~230 MB |

Least-squares scaling exponent **0.658** over 10¹²–10¹⁹, just under the theoretical 2/3.
π(10²⁰) is the first value **past 2⁶⁴** (via u128; π(10¹⁹) is the last under it) and matches the published value (Oliveira e Silva 2006 Table IV / Gourdon 2001). π(10¹⁹) matches M. Deléglise's 1996 computation —
reproduced here single-threaded in less time than his HP-730 took for a value four
powers of ten smaller. We implement *neither* of DR's two log-factor optimisations
(see [dead ends](#what-didnt-work-dead-ends-and-unhelpful-literature)) — they were made
unnecessary on this hardware, not skipped.

A cross-era *speed* comparison is not really recoverable, and the honest lean is that
hardware-adjusted we are at par or behind DR 1996, with the modern advantage dominated by
**native 64-bit width**. This computation lives entirely in 64-bit arithmetic (products
m·p ≤ y² ≈ 3×10¹⁴, divides x/(m·p)); their HP-730 *emulated* 64-bit in software (§10),
and DR report the native-64-bit DEC Alpha as ">3× faster" at similar SPEC — that factor is
precisely the emulation penalty, invisible to a clock×IPC adjustment. We measured the same
phenomenon one level up: the u128-divide tax at 10²⁰ (2.5–2.9× per op, compiler-rt emulating
128-bit on 64-bit hardware) is the 2026 version of their 1996 penalty. What *is* solid, and
independent of any cross-era guess, is that §7 and §6.5 measure net-negative **here** — a
direct measurement on this machine, not an adjusted comparison.

**π(10²⁰) = 2,220,819,602,560,918,840** (u128, past 2⁶⁴) — computed here in 3.44 h single-threaded and matching the published value (Gourdon 2001, in Oliveira e Silva 2006, Table IV), an independent check not a self-check. The u128-division tax landed ≤ a few % (10²⁰ sits on the existing 0.688 local drift) — exactly the modest, throughput-bound cost predicted, since divisions are OoO-overlapped and not the bottleneck. Oliveira's leaf taxonomy (trivial / easy / clustered-easy / hard),
π-table over [1,y] for easy leaves, and empirical flat-basin α all match what was reconstructed
here independently; clustered easy leaves are the one thing he implements and we measured
net-negative, for the reason above.

## Calibration against primecount (same machine, single thread)

primecount (Kim Walisch, ~13 yr) is the reference combinatorial π(x). Its README benchmarks
(π(10¹⁸) 1.58 s Gourdon / 3.72 s DR) are on a **32-core EPYC Zen4** — so the raw ~133–386×
over our single-threaded mini-PC number is mostly cores, not code. Built primecount 8.6 from
source and ran it **single-threaded on this box** to remove the machine and parallelism guesses:

| x | ours 1T | primecount DR 1T | primecount Gourdon 1T | vs DR | vs Gourdon |
|---|---:|---:|---:|---:|---:|
| 10¹⁴ | 1.18 s | 0.351 s | 0.168 s | 3.4× | 7.0× |
| 10¹⁶ | 23.16 s | 6.039 s | 2.595 s | 3.8× | 8.9× |

So the ~133× raw gap vs their 32-thread EPYC decomposes as **~26× parallelism × ~1.3× machine ×
3.8× implementation**. Same box, same thread, same algorithm class (DR), we are **3.8× off the
reference** — the honest implementation-maturity gap (libdivide fast division, hand-tuned inner
loops). Of the further gap to their *best*, **2.3× is the Gourdon algorithm** (their DR→Gourdon
speedup), which we have not implemented. Priorities in order of payoff: parallelism (the ~26×,
[scoped](#parallelism-measured-not-yet-built)), then Gourdon, then the 3.8× of division/micro-opt.

## What worked: the optimisation ladder

Each row is a separate committed measurement at 10¹⁴, in build order. 16.95 s → 1.21 s
is **14×**; the mod-30 wheel and the O(1)-kill counter are the two biggest single steps.

| change | 10¹⁴ | vs prev | section |
|---|---:|---:|---|
| Fenwick tree, α=1.5 (starting point) | 16.95 s | — | [counter](#the-o1-kill-counter-fenwick--2-level--3-level) |
| → block Counter (O(1) kill) | 8.91 s | 1.90× | [counter](#the-o1-kill-counter-fenwick--2-level--3-level) |
| → α = 2 | 7.98 s | 1.12× | [α knob](#the-α-knob) |
| → m-walk split at √y | 6.07 s | 1.31× | [√y split](#the-y-enumeration-split) |
| → P₂ fused onto the counter | 4.97 s | 1.22× | [P₂ fusion](#p₂-fusion-and-on-the-fly-p-range-sieving) |
| → α = 4 | 4.13 s | 1.20× | [α knob](#the-α-knob) |
| → 3-level counter, power-of-2 blocks | 3.54 s | 1.17× | [counter](#the-o1-kill-counter-fenwick--2-level--3-level) |
| → class-(1) binomial (LMO p.555 / DR §6.1) | 2.83 s | 1.27× | [leaf classes](#closed-form-leaf-classes-binomial--π-table) |
| → π-table for classes (2)/(3)-easy | 2.42 s | 1.17× | [leaf classes](#closed-form-leaf-classes-binomial--π-table) |
| → mod-30 wheel fold (DR §9) | 1.59 s | 1.52× | [wheel fold](#the-mod-30-wheel-fold-dr-9) |
| → branchless kill (*after* the wheel) | 1.25 s | 1.27× | [branchless kill](#the-branchless-kill-a-sign-that-flipped) |
| → gate diagnostics behind comptime | **1.21 s** | 1.03× | [verification](#verification-methodology) |

(The final sweep numbers above are lower still — 1.16 s at 10¹⁴ — from best-of-3 on a
quiet box.)

## What didn't work: dead ends and unhelpful literature

Roughly **ten cost models were falsified by measurement** over this build, several of
them the "obvious biggest remaining win". Recorded because the negatives were the most
informative results — including both of Deléglise–Rivat's headline optimisations.

| attempt | predicted | measured | why |
|---|---|---|---|
| **DR §7** x^(1/4) φ-sieve bound | biggest algorithmic win (83× smaller `a`) | **−2%** | rations O(log x) tree ops we don't have; un-fusing P₂ costs more |
| **DR §6.5** leaf clustering | the paper's log factor | **−7%** | a run costs ~2.5 leaf-evals; α_opt *did* move 4→8 as predicted, still lost |
| branchless kill (pre-wheel) | remove a mispredictor | **−33%** | branch was 98.8% predictable; became a win only *after* the wheel |
| fixed-trip prefix loops | delete 39% of branch misses | **−4.3%** | removed 69% of misses (worked!) but query is work-dominated, not mispredict-dominated |
| flat prefix-array for easy leaves | O(1) leaf eval | **−5%** | O(z) build for O(leaves) queries; and `ctr.prefix(v)` *is* φ already |
| next-multiple cursor (no division) | un-pin the segment size | **no-op** | the division's latency was already hidden by ILP |
| per-lpf linked list for m-walk | O(leaves) enumeration | worse | must be rebuilt per segment; √y split is the fix |
| SIMD prefix | vectorise the query | unpromising | AVX2 has no VPOPCNTDQ; nibble-LUT barely beats scalar popcnt |

**The pattern.** Both DR headline optimisations measure negative for the *same* reason —
they are amortisation schemes for costs our O(1)-kill counter and π-table had already
made cheap:

| DR optimisation | what it rations | what we did | measured |
|---|---|---|---|
| §7, x^(1/4) bound on the φ-sieve | tree ops at O(log x) each (§8.3) | O(1)-kill counter | **−2%** |
| §6.5, clustering | leaf evaluation | π-table: one lookup | **−7%** |

That is why using neither of their log factors is not a handicap here: we did not skip them,
we made them unnecessary on this hardware (the raw speed gap over their 1996 run is largely
native 64-bit width, not algorithm — see Results). It also retires the α_opt puzzle (below):
six sweeps said α_opt = 4 flat over five powers of ten against a literature that says α ~ log³x —
never a contradiction, because DR's α is large *because* their leaves are expensive.

**Meta-lesson.** Cost models were unreliable at nearly every step. Two rules survived:
(1) an optimisation's *sign* can flip when another lands (branchless kill: −33% before the
wheel, +27% after); (2) *always profile before AND after* — nine cost models preceded the
first `perf` run, which contradicted all of them. `perf` works here
(`perf_event_paranoid=1`); callgrind is the unprivileged fallback.

---

# Detailed sections

## The LMO decomposition (S1 + S2, δ = smallest prime factor)

Leaf structure derived from *our* recursion, not from memory (later confirmed against the
papers — see [scorecard](#what-the-papers-actually-say-scorecard)).
φ(x,a)=φ(x,a−1)−φ(x/p_a,a−1) takes primes in **decreasing** index order, so with
F(n,b) = μ(n)·φ(x/n,b): F(n,b) = F(n,b−1) + F(n·p_b,b−1) — the tree sums F exactly. Cut any
child with n·p_b > y:

- **Ordinary leaf** — n ≤ y reaches b=0 ⇒ μ(n)·⌊x/n⌋. Every squarefree n ≤ y gets there
  (all prefixes ≤ n ≤ y), so S1 = Σ_{n≤y} μ(n)⌊x/n⌋. O(y), direct — and **c = 0 is exactly
  right**, not a shortcut.
- **Special leaf** — n = m·p, n > y ≥ m = n/p. Primes *descend* along the path, so the last
  one multiplied is the **smallest**: p = P⁻(n), all of m's primes above it. ⇒
  **S2 = Σ −μ(m)·φ(x/(m·p_b), b−1)** over p_b-rough squarefree m ∈ (y/p_b, y].

The b−1 index is what makes the sieve work: at step b the sieve over [1, z], z = x/y, holds
exactly p_1..p_{b−1}, so sweeping b = 1..a is one monotone prime-at-a-time removal and
φ(v, b−1) is a prefix count. Every leaf has v = x/(m·p_b) ≤ z because m·p_b > y.

**Leaf law.** An earlier counter used `maxpf` (p = P⁺(m), the *opposite* convention) and
reported Θ(x^(2/3)/ln x). Against the real P⁻ set the ratio to x^(2/3)/ln x **drifts**
(0.805 → 0.474), while the ratio to x^(2/3)/**ln²**x is **flat at ~12.0**, and a²/2 predicts
the count almost exactly:

| x | a = π(y) | leaves | a²/2 | leaves/(x^(2/3)/ln²x) |
|---|---:|---:|---:|---:|
| 10⁸ | 125 | 7,815 | 7,813 | 12.31 |
| 10¹¹ | 894 | 403,134 | 399,618 | 12.00 |

So **leaves ≈ π(y)²/2 = Θ(x^(2/3)/ln²x)**, dominated by two-prime leaves n = p·q. This is
exactly LMO Lemma 5.1, confirmed to 0.03% at 10¹⁹: a = 578,568 ⇒ ½π(y)² = 167,370,465,312
vs measured 167,423,661,964.

**Segmentation** rests on two ideas: a *running φ per b* (φ(v,b−1) = phi_run[b−1] + alive in
[lo,v]), so [1,lo) is never re-sieved; and a *descending cursor per prime* — segments run lo
ascending and v = x/(m·p) rises as m falls, so every candidate is touched once in **total**,
not once per segment. Segmenting made S2 *faster* as well as smaller: 2.5× vs the flat [1,z]
version at 10¹¹ (336.4 → 132.9 ms), because the counter went cache-resident.

## The O(1)-kill counter (Fenwick → 2-level → 3-level)

The core structural win, and the reason we beat DR's leading term by a log factor.

**Fenwick was the wrong data structure — the traffic is lopsided.** Folding the primes
kills every element of [1,z] exactly once (z ≈ 1.4×10⁹ at 10¹⁴), while the leaves only
query a²/2 ≈ 2.4×10⁷ times — **~60:1**. Fenwick charges O(log S) for *both*, taxing the hot
side to subsidise the cold one. Replaced with a two-level counter (a bit per element + an
alive-count per block of √nwords words): **O(1) kill**, O(√S) query, O(1) segment total.
That last one matters more than it looks — the old code spent a full `prefix(len)` per
(segment, b), which is 1.4×10⁸ tree walks at 10¹⁴, now a single register read.

| x | Fenwick | Counter | speedup |
|---|--------:|--------:|--------:|
| 10¹¹ | 155.2 ms | 93.8 ms | 1.65× |
| 10¹² | 714.2 ms | 424.0 ms | 1.68× |
| 10¹³ | 3408 ms | 1943 ms | 1.75× |
| 10¹⁴ | 16947 ms | 8912 ms | **1.90×** |

The gain *grows* with x, as it must when you drop a log factor off the dominant term, and
pulled the exponent from 0.68 to 0.66. LMO's array a(i,j) in (3.9)–(3.10) — a binary
hierarchy of counts queried by decomposing an index into powers of two — **is a Fenwick
tree, never named as such**, and their own analysis makes its updates the *dominant* term at
x^(2/3)·log x. The Fenwick is kept as the flat [1,z] path's reference: flat (Fenwick) ==
segmented (Counter) == oracle.

**Two levels → three, and the two divisions that disagree.** By α=4 the traffic had
*inverted*: 5.4×10⁸ kills vs 1.4×10⁸ leaves is ~4:1 by count, but in *cost* the queries
dominate ~12:1 the other way. Retested with three variants separating the two effects:

| x | Counter (2-level, ÷) | Counter2P (2-level, shift) | Counter3P (3-level, shift) |
|---|---:|---:|---:|
| 10¹³ | 854.2 ms | 745.8 (1.15×) | 740.6 (**1.15×**) |
| 10¹⁴ | 3989.0 ms | 3667.5 (1.09×) | 3345.2 (**1.19×**) |

Predicted ~20% from the third level (query 3·nwords^(1/3) vs 2·√nwords). The total is ~19%,
but most of it is removing the **division** — at 10¹³ the third level adds nothing over a
2-level shift. Adopted 3P anyway: ≥ 2P at both, margin growing with x.

**Why this division cost when the [fold's](#the-mod-30-wheel-fold-dr-9) didn't.** In
`kill()`, `cnt[w / wpb] -= 1`: wpb is a runtime divisor and the quotient is an **address**,
so the dependent load-modify-store sits on the critical path. The fold's `((lo+p−1)/p)·p`
measured *free* because its quotient only seeded a loop whose iterations were independent
across primes — ILP swallowed the latency. **Division latency is free when it feeds
independent work, and expensive when it feeds an address you must immediately touch.**

## The √y enumeration split

The special-leaf sum enumerates m ∈ (y/p, y] per prime p. Done naïvely it rescans and
rejects most of the range: `w/leaf` ran 15–23 (≈ 2·ln y), so the m-walk *also* rode α²,
worth ~⅓ of φ. A per-lpf linked list does NOT fix it — the valid-m set depends on b and we
re-walk b inside every segment, so it would be rebuilt per segment at y·lnln y × nsegs ≈
1.4×10⁹, *worse* than the 1.9×10⁸ walk it replaces. Instead, split at p = √y:

- **p > √y** — any p-rough m ≤ y is 1 or a prime (a composite would be ≥ lpf² > y), and m=1
  fails m·p > y. So m = q prime in (p, y], with m·p > y automatic. Enumerable straight from
  primes[] — no walk, no rejects. This is the a²/2 bulk.
- **p ≤ √y** — keep the walk; it now costs π(√y)·y ≈ 2×10⁶ instead of 1.9×10⁸.

`w/leaf` → **1.12** (near-zero waste). Worth 1.31× at 10¹⁴. This is also LMO Lemma 5.1's own
proof split (δ(n) < √y vs ≥ √y), rediscovered as an optimisation.

## Closed-form leaf classes (binomial + π-table)

Most leaves never need the counter. Two cheap classes, from LMO p.555 / DR §6:

**Class (1) — the binomial.** For **p³ ≥ x**, x/(pq) < x^(1/3) < p, so φ(x/(pq), π(p)−1) = 1
*identically* — not a cheap lookup, literally 1. With m = q prime, μ(q) = −1 so −μ(m) = +1,
and the entire class collapses to one pair-count:

    S₁ = C(n₁, 2),   n₁ = #{primes p : p³ ≥ x, p ≤ y}

That is **~52% of all leaves at α=4** — 20.9×10⁹ counter queries at 10¹⁸ computing what is
one binomial. (The threshold comes from icbrt, not `p*p*p`, which overflows: p ≤ 4x^(1/3) ⇒
p³ ≤ 64x = 6.4×10²⁰ at 10¹⁹.) Worth 1.27×. The earlier "easy leaves cannot be cheaper"
[dead end](#what-didnt-work-dead-ends-and-unhelpful-literature) was true and irrelevant: for
this class you do not *evaluate* φ, you *count* it.

**Classes (2) and (3)-easy — one rule, not four cases.** Every remaining cheap class reduces
to the same test: **v ≤ y and p² > v** ⇒ φ(v,b−1) = 1 + max(0, π(v) − (b−1)) straight from a
π-table over **[1, y]** (O(y) space, 0.74 MB at 10¹⁴):

- **class (2)**, p > √(x/y): v < x/p² < y, and p > x^(1/3)/2 > x^(1/4) ⇒ p⁴ > x ⇒ v ≤ x/p² < p².
- **class (3)-easy**, q ≥ x/(yp): v ≤ y, and p² > (x/y²)² > y ≥ v — which requires **x² > y⁵,
  i.e. y ≤ x^(2/5)**. *That is why LMO's Truncation Rule T′ caps y at x^(2/5).*

The table spans [1, y], **never [1, z]** — every leaf it serves has v ≤ y by construction.
That distinction is the whole reason the earlier flat-prefix-array attempt failed (an O(z)
build for O(leaves) queries at z/leaves ≈ 27). Share caught grows with x — 43.1% at 10⁹ →
64.2% at 10¹⁴ — worth 1.17×, cumulative 1.46× with the binomial. Only class (3)-hard (v > y)
and class (4) ever touch the counter (~16M queries at 10¹⁴).

## P₂ fusion and on-the-fly p-range sieving

P₂(x,a) = Σ_{y<p≤√x} (π(x/p) − π(p) + 1). No Fenwick needed — it is the monotone sweep
(p ascending ⇔ x/p descending). But it can be **fused onto the S2 counter for free**: after
folding all a primes, the alive set in [1,z] is exactly {1} ∪ primes in (y, z] (a composite
with lpf > y is ≥ y², and y² > z whenever α > 1). So π(v) = φ(v,a) − 1 + a, and P₂'s π(x/p)
is a counter query at the end of each segment's b loop. Worth 1.24× — more than its 13%
share, because it deletes a whole prime sieve, not just bookkeeping.

**P₂'s memory, fixed and free.** Storing every prime ≤ √x is Θ(√x/ln x) — 5 MB at 10¹⁴,
386 MB at 10¹⁸ — which dominated the footprint. But for segment [lo, hi) the p with
x/p ∈ [lo, hi) are exactly **p ∈ (⌊x/hi⌋, ⌊x/lo⌋]**, and as lo sweeps those ranges **tile
(y, √x]**: disjoint, contiguous, each provably ≤ seg wide. So each is sieved on the fly from
base primes ≤ **x^(1/4)** (24 KB at 10¹⁸), and π(√x) counted as we go for the closed-form
Σ(π(p)−1).

| | 10¹⁴ | 10¹⁸ |
|---|---:|---:|
| primes ≤ √x (was) | 5.0 MB | 386.0 MB |
| primes ≤ x^(1/4) (now) | 3.1 KB | 24.4 KB |

Costs nothing in time; total memory now genuinely **Θ(x^(1/3))**, dominated by the μ/lpf
tables at O(y).

## The mod-30 wheel fold (DR §9)

The single biggest step (1.52×), and the first profile-driven one. `perf` said the fold
(`while (j < hi) ctr.kill(...)`) and kill were **~50% of all cycles** — the leaves, which
every prior optimisation attacked, were the other half. The fold's 50% is *volume*: it runs
Σ_p z/p ≈ 2.76z but only z kills land, so 66% of visits hit an already-dead element. The
only lever is fewer visits — DR §9: *"Precomputing the sieving by the first primes 2, 3, 5."*

**The decomposition that made it cheap.** The wheel's *stepping* (which multiples to visit)
is independent of the array's *indexing* (bit-per-integer vs 8-per-30). All the win is in the
stepping, so take that and leave the indexing — and every leaf/prefix path — untouched:

- `reset()` starts from a mod-30 mask instead of all-ones (lcm(30,64) = 960, so the u64
  pattern repeats every 15 words; segments align to 960). That state *is* φ(·,3).
- Fold p ≥ 7 by stepping j = p·m over m coprime to 30 (gaps 6,4,2,4,2,4,6,2). 2/3/5 are
  never folded at all.
- Σ_{p≥7} (8/30)·z/p ≈ **0.46z visits vs 2.76z — 6× fewer**.
- Leaves with p ∈ {2,3,5} need φ(v,0..2) — meissel.zig's closed-form base cases, and only
  ~0.5y of them (~93k of 141M at 10¹⁴).

Two bugs the small-x sweep caught, both from `lo` now starting at 0: `x / lo` trapped SIGFPE,
and `for (3..a+1)` had start > end when a = π(y) < 2.

## The branchless kill (a sign that flipped)

The clearest example of an optimisation whose sign depends on context. `kill()`'s
`if (bits[w] & b != 0)` guard was made unconditional (`alive = @intFromBool(...)`, subtract
0 or 1). Measured **0.75× — a 33% regression** when first tried, because the branch missed at
just 1.18% (the already-dead pattern within one prime's strike is periodic; TAGE nails it).
"Line 370 = 13.72% of branch misses" had *looked* unpredictable — but that is *share of
misses*, not *miss rate*; the line topped the chart because it executed most.

Re-profiling after the wheel, the balance had inverted:

| | pre-wheel | post-wheel |
|---|---:|---:|
| kill alive-check miss rate | 1.18% | **9.0%** |
| fold visits | 1.6×10⁹ | 2.5×10⁸ |

The wheel deleted p = 2/3/5 — exactly the primes whose already-dead patterns are periodic.
What remains is p ≥ 7, whose aliveness follows irregular factorisations: the branch got 8×
*less* predictable while the cost of removing it fell 6.4×. Both terms moved and the sign
flipped — the identical patch now measures **+27%**. IPC fell 2.65 → 1.74 across the wheel
(the fold was the ILP-friendly half; what remains is dependent loads), so branch misses,
barely moved in absolute terms, became ~⅓ of runtime and worth removing.

## The α knob

y = α·x^(1/3). In capped Meissel (see [minimax](#the-meissel-baseline-and-the-y-minimax-meisselzig)) α
was nearly free; in LMO it is a **sharp interior optimum**, because two terms fight:
z = x/y ~ 1/α (sieve kills) against leaves ≈ π(y)²/2 ~ α². Measured at 10¹³ (seg fixed at
x^(1/3), so one knob):

| α | z = x/y | leaves | total ms | vs best |
|--:|---:|---:|---:|---:|
| 1.0 | 464,166,357 | 2,939,588 | 2570 | +35.8% |
| **4.0** (after √y split) | 116,041,589 | 35,211,802 | (opt) | **0.0%** |
| 16.0 | 29,010,397 | 437,536,611 | 22038 | +1064.7% |

The curve is steeply asymmetric — gentle below, catastrophic above (the α² leaf term).
α_opt moved 1.5 → 2 → **4** as the m-walk waste was removed (each fix flattened the high-α
side), and then held at **4 across five powers of ten (10¹¹–10¹⁶)**. That flatness looked like a
contradiction of the literature's α ~ log³x for a long time; it is not — see the
[scorecard](#what-the-papers-actually-say-scorecard). P₂ monotonically wants large α (its
cost ~1/α) but is too small a share to move the optimum.

## Segment size (why S = y is right)

S rode y from the start and was never swept. Predicted the O(√S) query dominated, so
S ≈ 2570 (72× smaller) would win ~2.8×. Measured: **S = y is within 0.2–1.4% of optimal**,
and shrinking S is catastrophic (+249% at S=2048). Each (segment, prime) fold step costs
~26 cycles — first guessed to be the `((lo+p−1)/p)·p` division, but replacing it with a
next-multiple cursor was a **no-op** (the division's latency was already hidden by ILP; see
[the counter](#the-o1-kill-counter-fenwick--2-level--3-level)).

**What actually pins S is structural.** The ~26 cycles are memory traffic on the *a-sized*
arrays: every segment streams `seg_cnt`, `cur`, `next`, `phi_run` ≈ 536 KB at 10¹⁴, right at
the 512 KB L2 cliff. And it is irreducible: `phi_run[bi] += seg_cnt[bi]` must run for every b
at every segment boundary, because the running φ per b is exactly what lets us never re-sieve
[1,lo). So a·nsegs is a hard floor, and S ≈ 0.75y is right for a real reason. S = y stays.

## The scaling exponent and the drift

A cautionary tale about 2-point deltas. The per-step exponents (each from one ×10 in x) read
0.665 / 0.655 / 0.658 / 0.681 / **0.705** / 0.660 / 0.712 — they *alternate*, which is
run-to-run noise on single unrepeated runs (28 min each at the top), and it fooled the
analysis twice in *both* directions: first into "the exponent is climbing to 0.705", then
into "0.660 falsifies the drift". Only the least-squares fit over all points is a
measurement:

| range | before DR classes | after |
|---|---:|---:|
| 10¹² .. 10¹⁹ | 0.6755 | **0.6582** |
| 10¹² .. 10¹⁵ | 0.6586 | 0.6471 |
| 10¹⁵ .. 10¹⁹ | 0.6879 | 0.6705 |

The overall exponent is **0.658**, just under 2/3, with a real but mild drift (+0.023). The
DR leaf classes *lowered* the exponent but barely touched the drift — informative, because
stripping ~half the leaf work should have weakened a *leaf-driven* drift and did not. So the
residual lives on the kill/fold/counter side, consistent with the L1-miss climb
(0.78% → 2.15%) after the wheel shrank the working set. ~3% over seven powers of ten; not chased.

Two candidate mechanisms for the residual, neither established: (1) *cache* — the a-sized
arrays cross L2 around 10¹⁵; (2) *query growth* — the counter query ~3·nwords^(1/3) grows as
x^(1/9). An α sweep run as a "discriminating test" discriminated nothing (at 10¹⁶ even α=2
leaves the arrays outside L2), and a query-growth model reported as "fitting" did not on
recomputation. Both are noted here as retracted over-claims.

## What the papers actually say (scorecard)

LMO 1985 (Math. Comp. 44, 537–560) and Deléglise–Rivat 1996 (Math. Comp. 65, 235–245), read
only *after* the implementation was measured into its current shape.

**Confirmed, derived independently here:**
- **δ(n) = the SMALLEST prime factor** (DR eq. 10). The P⁻ derivation was right; the earlier
  `maxpf` counter was the wrong convention.
- **Our S2 is DR eq. (11)**, term for term.
- **S₀ = Σ_{n≤y} μ(n)⌊x/n⌋ with φ(u,0) = [u]** — c=0 confirmed. LMO's own code uses a wheel
  k=5; DR drops it to 0.
- **x^(1/3) ≤ y ≤ x^(1/2)** is DR's stated precondition — the invariant the π(2) bug forced
  on us (y = 4·icbrt(2) = 4 > x made a = π(4) count the prime 3 > x). Now clamped and verified
  exhaustively over [0, 5000].
- **LMO Lemma 5.1: "½·π(y)² + O(y^(3/2)/log y)"** — exactly the a²/2 law (0.03% at 10¹⁹), and
  its proof *is* the √y split.
- **LMO Truncation Rule T′: "x^(2/5) ≥ y ≥ x^(1/3)"** — recalled from memory, exact.

**The α_opt puzzle, resolved.** LMO p. 556: *"We choose **y = c·x^(1/3)** … a good value of
the constant **c was determined empirically**."* Our α_opt = 4, flat over five powers of ten and
found empirically, **is LMO's own prescription.** The α ~ log³x measured against is DR's
*asymptotic space bound*, infeasible at 10¹⁸ anyway (log³x = 71,197 vs the y ≤ √x cap of
α ≤ 1000). Six sweeps were compared to the wrong paper's constant.

**Where we are ahead of all three papers.** LMO's a(i,j) array is a Fenwick tree — Oliveira e
Silva 2006 (the most implementation-focused treatment) names it as such and cites Fenwick [9]
directly. All three (LMO, DR, Oliveira) keep **O(log z) on both update and query**; LMO §3.2
makes the *updates* the dominant term at x^(2/3)·log x. Our O(1)-kill counter (O(√S) query,
because kills outnumber queries ~60:1) removes that log from the leading cost — found by
measuring the traffic, and not present in even the most practical prior implementation.

## DR §7's x^(1/4) bound: measured, and not worth it

DR sieves [1, x/y] "by all primes less than x^(1/4)" — for us that is π(x^(1/4)) = 446 primes
vs π(y) = 16,801 at 10¹⁴, an 83× smaller `a` at 10¹⁸. It looked like the biggest remaining
algorithmic win. Measured, the fold-visit share with p ≥ x^(1/4) **shrinks** with x
(Σ_{x^(1/4)≤p≤y} 1/p → ln(4/3) constant, while Σ_{7≤p≤y} 1/p grows as ln ln y): 34.2% at 10⁹
→ 23.6% at 10¹⁴. So it saves ~9%, plus ~2% from smaller a-arrays — against **~13.5%** to
un-fuse P₂ (the bound leaves the counter at φ(·,π(x^(1/4))), breaking π(v) = φ(v,a) − 1 + a).
**Net ≈ −2%.**

DR §8.3 says why it pays for *them*: *"each access costs O(log x) instead of O(1) in a normal
sieve."* The x^(1/4) bound rations expensive tree operations. Our O(1)-kill counter is a whole
log x cheaper, which devalues the bound — and makes the P₂ fusion (which I had criticised) the
right trade after all. Both halves of that judgement were wrong, identically: pricing a paper's
design decision without pricing the data structure it was designed around.

## DR §6.5 clustering: measured reach, and why y = x^(1/3)·log³x exists

DR's actual novelty: *"for each p, split the sum over q into intervals where q ↦ π(x/pq) is
constant."* q ascending ⇒ v = x/(pq) descending ⇒ π(v) non-increasing, so runs are contiguous.
Measured the reach *before* building it:

```
x = 10^14
 alpha         y  pi-tab leaves           runs     ctr leaves  leaf/run    reach
     4    185660       53037184       21780399       14674931      2.44    78.3%
    16    742640      262535484       45657832        3541600      5.75    98.7%
    32   1485280      507399194       49087615        1266051     10.34    99.8%
```

At α=4 it is worth ~8–12%. **But the reach grows with α, and that is the point:** the runs
*saturate* (21.8M → 49.1M) while the leaves explode (53M → 507M), so clustering converts the
α² leaf term into a nearly-flat one. Meanwhile large-α's two costs collapse: counter leaves
14.7M → 1.3M and z shrinks 8×. So DR's y = x^(1/3)·log³x is not a separate choice from
clustering — it is what clustering *pays for*.

**Built it. Correct — and a net loss.** S₂ bit-identical to the reference; run counts matched
the standalone prediction exactly.

```
x = 10^14                      (alpha=4 without clustering: 1209.2 ms)
 alpha         y       z = x/y          ms    vs best
     4    185660     538618980      1346.8       4.1%
     8    371320     269309490      1294.2       0.0%     <- alpha_opt DID move 4 -> 8
    16    742640     134654745      1515.1      17.1%
```

**α_opt moved 4 → 8, exactly as predicted** — and the best clustered configuration is still
7% worse than no clustering at α=4. Every step of the reasoning held; the answer is still no.
A run costs ~2.5 leaves (three random lookups — `pi_tab[v]`, `primes[c]`, `pi_tab[qmin−1]` —
against a leaf's one sequential walk), so leaf/run must clear ~2.5 to break even; and the
dense m-walk (p ≤ √y) rides α^1.5, which clustering cannot reach and which swamps z's 1/α
shrink past α≈8. Hoisting the p-invariant divisions out of the run finder moved it
0.844 → 0.854 — the divisions were never the problem.

## The Meissel baseline, and the y-minimax (meissel.zig)

`meissel.zig` is the classical Meissel–Lehmer recursion φ(x,a)=φ(x,a−1)−φ(x/p_a,a−1) with a
mod-30 wheel base and the **leaf cutoff** φ(x,a)=1+max(0,π(x)−a) when p_a²≥x — which drops it
from exponential to O(x^(2/3)). A compact bit-sieve + per-word checkpoints (`PiTable`) keeps
π(y) at O(1) and the table cache-resident (~16× less memory: 0.5 GB vs 8.6 GB at 10¹⁴), but
the algorithm is unchanged:

| x | π(x) | time | exponent |
|---|------|-----:|:--------:|
| 10¹¹ | 4,118,054,813 | 0.15 s | 0.71 |
| 10¹² | 37,607,912,018 | 0.89 s | 0.77 |
| 10¹³ | 346,065,536,839 | 5.68 s | 0.81 |
| 10¹⁴ | 3,204,941,750,802 | 41.3 s | 0.86 |

Sub-linear (~0.8) but with an **O(x^(2/3)) memory wall** at ~10¹⁴. Kept as LMO's φ oracle
(`phiOfXY`).

**The y-minimax.** The identity holds for any y ≥ x^(1/3) (below it you need P₃). Two consumers
want the π-table: φ's cutoff needs it up to y² (grows with α), P₂ up to x/y (shrinks). They
cross at α = 1 — **x^(1/3) is a minimax**, and raising y un-capped is pure loss (99% table
build). *Capping* the cutoff at v ≤ z = x/y deletes the y² consumer, so the table is x^(2/3)/α
and shrinks with α — measured α=1.5 as the argmin, beating classical α=1 on *both* axes. This
is the Meissel precursor to LMO's α knob; capping alone takes Meissel from Θ(x^(2/3)) toward
Θ(√x), but reaching Θ(x^(1/3)) needs the full LMO leaf restructuring (lmo.zig).

## Parallelism (measured, not yet built)

The leaves are **wildly non-uniform**: leaves are dominated by n = p·q with pq near y², so
v = x/(pq) clusters near x/y² — the very bottom of [1,z].

| decile of [1,z] | 1 | 2 | 5 | 10 |
|---|---:|---:|---:|---:|
| share of leaves | **99.72%** | 0.11% | 0.02% | 0.01% |

Since leaf queries dominate the work, decile 1 holds ~90% of it: an equal-width block split
would starve every thread but one. **Blocks must be sized by leaf count**, which is cheap to
precompute. The prefix dependency (phi_run per b) is not an obstacle — φ enters every leaf
linearly, so a block-and-scan works: each thread takes a contiguous block, computes relative
to phi_run = 0 at its start, returns local S2 plus `block_total[bi]` and `mu_sum[bi]`, and a
serial O(nthreads × a) scan corrects. P₂ is *easier* fused than standalone: π(v) = φ(v,a) − 1 + a
makes its per-block correction a scalar.

## Verification methodology

Correctness held through ~40 commits of aggressive optimisation via three independent oracles:
- **Differential**: S₂ and P₂ each checked against `specialS2Segmented` / `p2Segmented` — an
  independent path (all-ones reset, plain fold over every prime, branchy kill, Fenwick) that
  evaluates every leaf individually. Every closed-form class and the wheel path produce
  *bit-identical* S₂, and the leaf *counts* match too, so a shortcut is proven to account for
  exactly the set it skips.
- **Known values**: exact at all of OEIS A006880's π(10ⁿ) through 10¹⁹.
- **Exhaustive small-x**: every x in [0, 5000] vs a plain sieve — this is what caught the π(2)
  = 2 bug (y > x) that the α=2 default had masked.

The diagnostics (`walk`/`leaves`/`easy` counters) that power the analysis above cost 3.3% at
10¹⁴ and are gated behind a comptime `INST` flag — kept because `leaves` confirmed Lemma 5.1,
`walk` guards the √y split, and `easy` sized the π-table class, but not paid for on production
runs.
