//! HKM rung 0 — the segmentation geometry, refereed by brute force.
//!
//! Hirsch, Kessler & Mendlovic, arXiv:2212.09857, "Counting primes in
//! O(sqrt N log^3 N) time". This rung builds NONE of the algorithm: it
//! establishes the two objects everything else is defined against, and
//! measures the raw segmentation error before any correction exists.
//!
//! Two things get checked, both against independent brute force:
//!
//!   A. Claim 2 (§2.3.1). The paper builds mu-hat as a CONVOLUTION of
//!      per-prime arrays, mu-hat = *_{p <= sqrt N} (delta_0 - delta_kbar(p)),
//!      and asserts it equals the direct bucket sum
//!      mu-hat[k] = sum_{n : khat(n) = k} mu_{<=sqrt N}(n).
//!      Note khat(n) = sum_i kbar(p_i) over the factorization, NOT kbar(n) —
//!      this is the "alternative world" where each prime p is replaced by
//!      2^{kbar(p) Delta}, so 42 = 2*3*7 lands at 2^1 * 2^1 * 2^2 = 16 even
//!      though 42 is between 32 and 64. Getting this wrong (bucketing the
//!      product instead of summing the buckets) silently builds a different
//!      and wrong array, so it is checked by subset enumeration.
//!
//!   B. Lemma 1 (§2.1). sum_{n<=N} (1 * mu_{<=sqrt N})(n) = pi(N) - pi(sqrt N) + 1.
//!      This is Legendre's phi(N, pi(sqrt N)) — the identity the whole
//!      algorithm evaluates — refereed against a sieve.
//!
//! Then C reports the number R1/R2 exist to kill: the segmented estimate
//! sum_{k1+k2 <= kbar(N)} 1bar[k1] * mu-hat[k2] against B's exact answer.
//! That gap is the segmentation error, repaired later by sieving the
//! critical interval (N, N+S]. Measured here as a function of Delta, with
//! no correction, so the correction has something to be judged against.

const std = @import("std");
const rs = @import("rs");

/// kbar(n) = floor(log2 n / Delta): the cell of the geometric segmentation
/// with thresholds 2^{k Delta}. n = 1 gives cell 0, which is exactly why the
/// per-prime Moebius factor is delta_0 - delta_{kbar(p)}.
///
/// FLOAT SEAM, flagged not hidden: the cell walls 2^{k Delta} are irrational,
/// so this comparison is inherently inexact and a prime sitting a hair from a
/// wall can land either side. At rung-0 sizes f64 is far more precision than
/// the question needs; monotonicity is asserted below. When Delta shrinks to
/// Theta(log2 N / sqrt N) this needs revisiting — it is the same class of bug
/// as pistar's half-ulp band.
pub fn kbar(n: u64, delta: f64) usize {
    const l = @log2(@as(f64, @floatFromInt(n)));
    return @intFromFloat(@floor(l / delta));
}

/// mu-hat by the paper's route: convolve (delta_0 - delta_{kbar(p)}) over all
/// primes <= sqrt N. Truncating at `cap` cells is free — 1bar's indices are
/// >= 0, so a mu-hat cell beyond kbar(N) can never reach the final sum. That
/// truncation is what makes the array O(sqrt N) once Delta = Theta(log2 N / sqrt N).
fn muHatConv(gpa: std.mem.Allocator, primes: []const u64, delta: f64, cap: usize) ![]i64 {
    const a = try gpa.alloc(i64, cap + 1);
    @memset(a, 0);
    a[0] = 1;
    for (primes) |p| {
        const j = kbar(p, delta);
        if (j > cap) continue; // shifts entirely off the end
        // in place, descending: a[i] -= a[i - j] with a[] the previous array
        var i: usize = cap + 1;
        while (i > j) {
            i -= 1;
            a[i] -= a[i - j];
        }
    }
    return a;
}

/// mu-hat by direct enumeration of squarefree sqrt-N-smooth n: sign (-1)^r
/// into cell sum_i kbar(p_i). Independent of the convolution by construction
/// (no shared code path), which is the point. 2^|primes| subsets, so this
/// referee only runs at rung-0 sizes.
fn muHatEnum(gpa: std.mem.Allocator, primes: []const u64, delta: f64, cap: usize) ![]i64 {
    const a = try gpa.alloc(i64, cap + 1);
    @memset(a, 0);
    const R = struct {
        fn go(out: []i64, ks: []const usize, i: usize, k: usize, sign: i64, capa: usize) void {
            out[k] += sign;
            var j = i;
            while (j < ks.len) : (j += 1) {
                if (k + ks[j] <= capa) go(out, ks, j + 1, k + ks[j], -sign, capa);
            }
        }
    };
    const ks = try gpa.alloc(usize, primes.len);
    defer gpa.free(ks);
    for (primes, 0..) |p, i| ks[i] = kbar(p, delta);
    R.go(a, ks, 0, 0, 1, cap);
    return a;
}

/// Legendre's phi(N, pi(sqrt N)) = sum over squarefree sqrt-N-smooth d of
/// mu(d) * floor(N/d), by pruned DFS. This is the left-hand side of Lemma 1
/// evaluated EXACTLY — the referee for everything the segmented world later
/// approximates.
pub fn legendre(N: u64, primes: []const u64, i: usize, d: u64, sign: i64) i64 {
    var s: i64 = sign * @as(i64, @intCast(N / d));
    var j = i;
    while (j < primes.len) : (j += 1) {
        if (d > N / primes[j]) break; // primes ascend, so nothing beyond fits
        s += legendre(N, primes, j + 1, d * primes[j], -sign);
    }
    return s;
}

pub fn main(init: std.process.Init) !void {
    var args = try init.minimal.args.iterateAllocator(init.gpa);
    _ = args.skip();
    const gpa = std.heap.page_allocator;

    const usage = "Usage: hkm0 <N> [delta ...]\n";
    const n_str = args.next() orelse return std.debug.print(usage, .{});
    const N = try std.fmt.parseInt(u64, n_str, 10);

    var deltas = try std.ArrayList(f64).initCapacity(gpa, 8);
    defer deltas.deinit(gpa);
    while (args.next()) |a| try deltas.append(gpa, try std.fmt.parseFloat(f64, a));
    if (deltas.items.len == 0) {
        for ([_]f64{ 1.0, 0.5, 0.25, 0.1, 0.05, 0.02, 0.01 }) |d| try deltas.append(gpa, d);
    }

    const sq: u64 = std.math.sqrt(N);
    const primes = try rs.basePrimes(gpa, sq);
    defer gpa.free(primes);

    // ---- kbar sanity: cell 0 holds n = 1, and the map is non-decreasing.
    std.debug.assert(kbar(1, 1.0) == 0);
    {
        var n: u64 = 1;
        var prev: usize = 0;
        while (n <= @min(N, 1 << 20)) : (n += 1) {
            const k = kbar(n, 0.1);
            std.debug.assert(k >= prev);
            prev = k;
        }
    }

    // ---- B: Lemma 1, exact, refereed by a sieve.
    const pi_N = try rs.countInRange(u64, gpa, 2, N + 1, primes);
    const pi_sq = try rs.countInRange(u64, gpa, 2, sq + 1, primes);
    const lhs = legendre(N, primes, 0, 1, 1);
    const rhs: i64 = @as(i64, @intCast(pi_N)) - @as(i64, @intCast(pi_sq)) + 1;
    std.debug.print("N = {d}  sqrt N = {d}  pi(N) = {d}  pi(sqrt N) = {d}\n", .{ N, sq, pi_N, pi_sq });
    std.debug.print("B  Lemma 1: legendre = {d}   pi(N)-pi(sqrt N)+1 = {d}   {s}\n\n", .{ lhs, rhs, if (lhs == rhs) "MATCH" else "MISMATCH" });
    if (lhs != rhs) return error.Lemma1Failed;

    std.debug.print("{s:>7} {s:>7} {s:>8} {s:>14} {s:>14} {s:>12}\n", .{ "delta", "cells", "ClaimA", "segmented", "exact", "seg-error" });
    for (deltas.items) |delta| {
        const cap = kbar(N, delta);

        const conv = try muHatConv(gpa, primes, delta, cap);
        defer gpa.free(conv);

        // A: the enumeration is pruned by `cap` (k + kbar(p) <= cap is exactly
        // "the product stays <= N in the alternative world"), so it costs
        // Psi-many nodes, not 2^|primes| — the same order as the legendre DFS.
        var claim: []const u8 = "skipped";
        if (primes.len <= 2000) {
            const enu = try muHatEnum(gpa, primes, delta, cap);
            defer gpa.free(enu);
            claim = if (std.mem.eql(i64, conv, enu)) "MATCH" else "MISMATCH";
            if (!std.mem.eql(i64, conv, enu)) return error.Claim2Failed;
        }

        // C: 1bar[k] = #{ m >= 1 : kbar(m) = k }, exact integer counts.
        const one = try gpa.alloc(i64, cap + 1);
        defer gpa.free(one);
        @memset(one, 0);
        {
            var k: usize = 0;
            while (k <= cap) : (k += 1) {
                const lo = std.math.pow(f64, 2.0, @as(f64, @floatFromInt(k)) * delta);
                const hi = std.math.pow(f64, 2.0, @as(f64, @floatFromInt(k + 1)) * delta);
                const a: u64 = @intFromFloat(@ceil(lo));
                const b: u64 = @intFromFloat(@ceil(hi));
                one[k] = @as(i64, @intCast(b)) - @as(i64, @intCast(a));
            }
        }

        // the segmented estimate of Lemma 1's LHS: every (m, d) pair whose
        // approximate index k1 + k2 lands at or below kbar(N)
        var seg: i64 = 0;
        for (one, 0..) |c1, k1| {
            if (c1 == 0) continue;
            var k2: usize = 0;
            while (k1 + k2 <= cap) : (k2 += 1) seg += c1 * conv[k2];
        }

        std.debug.print("{d:>7.3} {d:>7} {s:>8} {d:>14} {d:>14} {d:>12}\n", .{ delta, cap + 1, claim, seg, rhs, seg - rhs });
    }
}
