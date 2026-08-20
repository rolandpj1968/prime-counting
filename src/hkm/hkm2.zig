//! HKM rung 2 — the critical-interval correction. First end-to-end pi(N).
//!
//! R1 proved the error is one-sided and localised: every unit of it is a
//! SPURIOUS pair (m, d) with m*d > N that the approximation nevertheless
//! admits, and all of them satisfy m*d < N * 2^(D(1+omega(d))). So
//!
//!     exact LHS = segmented - C,
//!     C = sum_{n = N+1}^{N+S} sum_{d | n, d sqfree sqrt-N-smooth,
//!                              kbar(n/d) + khat(d) <= K}  mu(d)
//!
//! and pi(N) = exact LHS + pi(sqrt N) - 1 by Lemma 1.
//!
//! Two things R1 measured drive the design:
//!
//!   * The correction MUST carry the mu signs. At 1e8 and the designed Delta,
//!     2.1M spurious pairs cancel to a net -15240 (139:1). Counting anything
//!     over the interval — pairs, integers, primes — cannot produce that.
//!
//!   * S is taken from the PROVEN bound, not R1's measured maximum. The
//!     measured spurious extent runs ~2.2x inside the bound; both are printed,
//!     so the slack stays visible instead of being quietly assumed away.
//!
//! The sieve is simpler than it looks: d is squarefree, so only the SET of
//! distinct primes <= sqrt N dividing n matters. No multiplicities, no
//! dividing out, no rough-part bookkeeping — a prime just marks the n it
//! divides. Primes in (sqrt N, sqrt(N+S)] are correctly absent: they are not
//! sqrt-N-smooth, so they can never appear in d.
//!
//! CRITICAL DIVISORS (HKM §3.6). Not every squarefree smooth divisor of n can
//! contribute. Because kbar(x) >= log2(x)/D - 1 and khat is a sum of omega(d)
//! such floors,
//!
//!     kbar(n/d) + khat(d) >= kbar(n) - 1 - omega(d)
//!
//! so the admission test kbar(n/d) + khat(d) <= K can only pass when
//!
//!     omega(d) >= kbar(n) - K - 1  =:  t(n)
//!
//! The critical interval spans only ~12 cells at 1e12, and t(n) is one less
//! than the cell index, so the first cell prunes nothing and the last admits
//! only divisors with omega(d) >= 11 — of which there are almost none, since
//! omega_max is 11. The DFS therefore carries omega and prunes any branch that
//! can no longer reach t(n) even by taking every remaining prime. n with
//! omega(n) < t(n) are skipped outright, which is the cheap half of §3.7 (in
//! the k-th segment only n with omega(n) >= k can contribute); §3.7's other
//! half — a restricted sieve that avoids FACTORING those n at all — only pays
//! past k ~ 18 log2 log2 N, far beyond any N reachable here.
//!
//! The bound is stated for ideal walls 2^(kD), and it applies to seg.zig's
//! integer table UNCHANGED, because that table computes the ideal floor
//! exactly. wall[k] = ceil(2^(kD)) >= 2^(kD) gives kbar(x) <= floor(log2 x/D);
//! and at k = floor(log2 x / D) we have wall[k] < 2^(kD) + 1 <= x + 1, so
//! wall[k] <= x (both integers) and kbar(x) >= k. Hence kbar(x) =
//! floor(log2 x / D) on the nose — which is what the double-double wall
//! recurrence buys, beyond merely not corrupting the answer. --crit-off
//! exposes the constant: 1 is HKM's, 0 is the default here and is one
//! sharper. The sharpening is proved, not measured. Writing A = kbar(n/d) +
//! khat(d) and summing omega(d)+1 STRICT inequalities floor(y) > y - 1,
//!     A > log2(n)/D - (omega(d) + 1) >= kbar(n) - omega(d) - 1
//! and A, kbar(n), omega(d) are all integers, so A >= kbar(n) - omega(d).
//! Admission A <= K then forces omega(d) >= kbar(n) - K. All three settings
//! produce a bit-identical correction at every N tested, which is the check
//! that the prune changes work and not the answer.
//!
//! Storage is CSR (count, prefix, fill) rather than a fixed [k]u32 per n:
//! omega is small on average (~ln ln N) but its max grows, and a fixed row
//! width would pay the max everywhere. The interval is processed in BLOCKS —
//! S itself is Otilde(sqrt N) but its CSR is a few times larger, and the
//! whole point of the method is Otilde(sqrt N) SPACE, so materialising the
//! entire interval at once would give the space bound away for nothing.

const std = @import("std");
const rs = @import("rs");
const r0 = @import("hkm0.zig");
const sg = @import("seg.zig");
const common = @import("common");

fn nowNs() u64 {
    var ts: std.os.linux.timespec = undefined;
    _ = std.os.linux.clock_gettime(.MONOTONIC, &ts);
    return @as(u64, @intCast(ts.sec)) * 1_000_000_000 + @as(u64, @intCast(ts.nsec));
}
fn secs(a: u64, b: u64) f64 {
    return @as(f64, @floatFromInt(b - a)) / 1e9;
}

/// Walk the squarefree sqrt-N-smooth divisors d of n, accumulating mu(d) over
/// those the approximation admits. Pruned on khat(d) <= K: past that no m at
/// all survives, so the whole subtree is dead.
fn corr(g: *const sg.Geom, n: u64, plist: []const u32, primes: []const u64, ks: []const usize, i: usize, d: u64, khat: usize, w: usize, t: usize, sign: i64, acc: *i64) void {
    // even taking every prime left, can this branch still reach omega >= t?
    if (w + (plist.len - i) < t) return;
    if (w >= t and g.kbar(n / d) + khat <= g.K) acc.* += sign;
    var j = i;
    while (j < plist.len) : (j += 1) {
        const pi = plist[j];
        if (khat + ks[pi] > g.K) continue;
        corr(g, n, plist, primes, ks, j + 1, d * primes[pi], khat + ks[pi], w + 1, t, -sign, acc);
    }
}

pub fn main(init: std.process.Init) !void {
    var args = try init.minimal.args.iterateAllocator(init.gpa);
    _ = args.skip();
    const gpa = std.heap.page_allocator;

    const usage = "Usage: hkm2 <N> [delta] [--no-critical] [--crit-off n]\n";
    const n_str = args.next() orelse return std.debug.print(usage, .{});
    const N = try std.fmt.parseInt(u64, n_str, 10);
    const sq: u64 = std.math.sqrt(N);
    const nf: f64 = @floatFromInt(N);
    var delta: f64 = @log2(nf) / @as(f64, @floatFromInt(sq));
    var critical = true;
    var crit_off: i64 = 0; // t(n) = kbar(n) - K - crit_off; 0 is proved below
    while (args.next()) |a| {
        if (std.mem.eql(u8, a, "--no-critical")) {
            critical = false;
        } else if (std.mem.eql(u8, a, "--crit-off")) {
            const v = args.next() orelse return std.debug.print(usage, .{});
            crit_off = try std.fmt.parseInt(i64, v, 10);
        } else delta = try std.fmt.parseFloat(f64, a);
    }
    std.debug.assert(delta <= 1.0);

    const t0 = nowNs();
    const primes = try rs.basePrimes(gpa, sq);
    defer gpa.free(primes);
    var g = try sg.build(gpa, N, delta);
    defer g.deinit(gpa);
    const K = g.K;
    const ks = try gpa.alloc(usize, primes.len);
    defer gpa.free(ks);
    for (primes, 0..) |p, i| ks[i] = g.kbar(p);

    // ---- omega_max: the largest number of distinct primes any admissible d
    // can carry, = greedily take the cheapest kbar(p) (i.e. smallest primes)
    // while khat stays <= K. This is what sets S.
    var omega_max: usize = 0;
    {
        var acc: usize = 0;
        for (ks) |k| {
            if (acc + k > K) break;
            acc += k;
            omega_max += 1;
        }
    }
    const s_bound_f = nf * (std.math.pow(f64, 2.0, delta * @as(f64, @floatFromInt(1 + omega_max))) - 1.0);
    const S: u64 = @as(u64, @intFromFloat(@ceil(s_bound_f))) + 2; // strict ineq + float cushion

    std.debug.print("N = {d}  sqrt N = {d}  delta = {d:.6}  K = {d}  omega_max = {d}\n", .{ N, sq, delta, K, omega_max });
    std.debug.print("critical interval (N, N+S], S = {d}  ({d:.2} x sqrt N)\n", .{ S, @as(f64, @floatFromInt(S)) / @sqrt(nf) });
    const t_setup = nowNs();

    // ---- segmented estimate: 1bar convolved with mu-hat, truncated at K.
    // Prefix-summing mu-hat makes this O(K) instead of R0's O(K^2).
    var segmented: i128 = 0;
    {
        const conv = try gpa.alloc(i64, K + 1);
        defer gpa.free(conv);
        @memset(conv, 0);
        conv[0] = 1;
        for (ks) |j| {
            if (j > K) continue;
            var t: usize = K + 1;
            while (t > j) {
                t -= 1;
                conv[t] -= conv[t - j];
            }
        }
        const pre = try gpa.alloc(i128, K + 1);
        defer gpa.free(pre);
        var run: i128 = 0;
        for (conv, 0..) |v, i| {
            run += v;
            pre[i] = run;
        }
        var k: usize = 0;
        while (k <= K) : (k += 1) {
            const c1: i128 = @intCast(g.cell(k)); // the wall table, not a fresh exp2
            if (c1 != 0) segmented += c1 * pre[K - k];
        }
    }
    const t_seg = nowNs();

    // ---- sieve the critical interval in blocks, CSR per block: the distinct
    // smooth primes of each n. Only the SET matters (d is squarefree), so a
    // prime just marks the n it divides — no multiplicities, no dividing out.
    const BLK: usize = 1 << 22;
    var C: i64 = 0;
    var max_hit: u64 = 0;
    var max_omega_seen: usize = 0;
    var ns_sieve: u64 = 0;
    var ns_corr: u64 = 0;
    var skipped: u64 = 0;
    {
        const off = try gpa.alloc(u32, BLK + 1);
        defer gpa.free(off);
        const fill = try gpa.alloc(u32, BLK);
        defer gpa.free(fill);
        var ent = try gpa.alloc(u32, BLK * 4);
        defer gpa.free(ent);

        var base: u64 = N + 1;
        while (base <= N + S) {
            const hi = @min(base + BLK - 1, N + S);
            const w: usize = @intCast(hi - base + 1);
            const ta = nowNs();
            @memset(off[0 .. w + 1], 0);
            for (primes) |p| {
                var q = ((base + p - 1) / p) * p;
                while (q <= hi) : (q += p) off[@intCast(q - base)] += 1;
            }
            var total: u32 = 0;
            for (off[0..w]) |*c| {
                const v = c.*;
                c.* = total;
                total += v;
            }
            off[w] = total;
            if (total > ent.len) {
                gpa.free(ent);
                ent = try gpa.alloc(u32, total);
            }
            @memset(fill[0..w], 0);
            for (primes, 0..) |p, pi| {
                var q = ((base + p - 1) / p) * p;
                while (q <= hi) : (q += p) {
                    const s: usize = @intCast(q - base);
                    ent[off[s] + fill[s]] = @intCast(pi);
                    fill[s] += 1;
                }
            }
            const tb = nowNs();
            var s: usize = 0;
            while (s < w) : (s += 1) {
                const n = base + @as(u64, s);
                const plist = ent[off[s]..off[s + 1]];
                if (plist.len > max_omega_seen) max_omega_seen = plist.len;
                // t(n) = kbar(n) - K - 1, floored at 0, minus the slack
                var t: usize = 0;
                if (critical) {
                    const v = @as(i64, @intCast(g.kbar(n))) - @as(i64, @intCast(g.K)) - crit_off;
                    if (v > 0) t = @intCast(v);
                }
                if (plist.len < t) {
                    skipped += 1; // cheap half of §3.7
                    continue;
                }
                var acc: i64 = 0;
                corr(&g, n, plist, primes, ks, 0, 1, 0, 0, t, 1, &acc);
                if (acc != 0) max_hit = n;
                C += acc;
            }
            ns_sieve += tb - ta;
            ns_corr += nowNs() - tb;
            base = hi + 1;
        }
    }
    const t_corr = nowNs();

    const lhs: i128 = segmented - C;

    // ---- referee: the repo's published-pi table where it has N, else a full
    // sieve. The sieve dominates the whole run past 1e9 — the algorithm is
    // already faster than the thing checking it.
    const pi_sq = try rs.countInRange(u64, gpa, 2, sq + 1, primes);
    const from_table = common.expectedPi(N);
    const pi_N = from_table orelse try rs.countInRange(u64, gpa, 2, N + 1, primes);
    const exact: i128 = @as(i128, pi_N) - @as(i128, pi_sq) + 1;
    const t_ref = nowNs();

    std.debug.print("\nsegmented   = {d}\ncorrection  = {d}\nLHS         = {d}\nexact LHS   = {d}   {s}\n", .{ segmented, C, lhs, exact, if (lhs == exact) "MATCH" else "MISMATCH" });
    const pi_computed: i128 = lhs + @as(i128, pi_sq) - 1;
    std.debug.print("\npi({d}) = {d}\nreferee  = {d} ({s})   {s}\n", .{ N, pi_computed, pi_N, if (from_table != null) "published table" else "sieve", if (pi_computed == pi_N) "MATCH" else "MISMATCH" });
    std.debug.print("\nspurious extent: largest contributing n - N = {d}  (S bound {d}, slack {d:.2}x)\n", .{ max_hit - N, S, @as(f64, @floatFromInt(S)) / @as(f64, @floatFromInt(@max(1, max_hit - N))) });
    std.debug.print("max omega in interval = {d}  (omega_max bound {d})\n", .{ max_omega_seen, omega_max });
    std.debug.print("critical divisors: {s}", .{if (critical) "on" else "OFF"});
    if (critical) std.debug.print(" (offset {d}) — {d} of {d} interval integers skipped outright ({d:.1}%)", .{ crit_off, skipped, S, 100.0 * @as(f64, @floatFromInt(skipped)) / @as(f64, @floatFromInt(S)) });
    std.debug.print("\n", .{});
    std.debug.print("\ntime: setup {d:.2}s  segmented {d:.2}s  sieve {d:.2}s  correction {d:.2}s  = {d:.2}s  | referee {d:.2}s\n", .{ secs(t0, t_setup), secs(t_setup, t_seg), @as(f64, @floatFromInt(ns_sieve)) / 1e9, @as(f64, @floatFromInt(ns_corr)) / 1e9, secs(t0, t_corr), secs(t_corr, t_ref) });

    if (lhs != exact) return error.CorrectionWrong;
}
