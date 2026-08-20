//! HKM rung 1 — where the segmentation error lives.
//!
//! R0 measured the error. R1 localises it *exactly*, before any correction
//! exists, so that R2's correction is refereed by something measured rather
//! than by the paper's O(·). Two structural facts, both checked rather than
//! assumed, and one closed form that ties them to R0's number.
//!
//! ONE-SIDEDNESS.  m·d <= N  =>  kbar(m) + khat(d) <= K,  K = kbar(N).
//!   kbar(x) <= log2(x)/D for every x, and khat is a sum of such floors, so
//!   kbar(m) + khat(d) <= log2(m·d)/D <= log2(N)/D. The left side is an
//!   integer, hence <= floor(log2 N / D) = K.
//!   CONSEQUENCE: the approximation never DROPS a true pair. Every unit of
//!   error is a SPURIOUS pair with m·d > N — so the correction is a pure
//!   subtraction, and the sign R0 measured (always negative) is a statement
//!   about mu on the spurious set, not about anything being lost.
//!
//! CRITICAL INTERVAL.  kbar(x) > log2(x)/D - 1, so a surviving pair obeys
//!   K >= kbar(m) + khat(d) > log2(m·d)/D - (1 + omega(d)), i.e.
//!   m·d < N · 2^(D(1+omega(d))).
//!   omega(d) is capped by the largest primorial that survives the khat
//!   prune (~log N / log log N), which is where the paper's S = O(D N log N)
//!   comes from. Both the bound and the ACTUAL maximum are reported below —
//!   the gap between them is what R2 gets to exploit.
//!
//! THE CLOSED FORM. No pair enumeration is needed. For fixed d the
//! approximate predicate is kbar(m) <= K - khat(d), i.e. m <= M(d), while the
//! true predicate is m <= floor(N/d). So
//!
//!     segmented - exact = sum_d mu(d) * ( M(d) - floor(N/d) )
//!
//! with every bracket >= 0 by one-sidedness. Reproducing R0's measured error
//! from this form is the proof that the localisation is COMPLETE — that
//! nothing outside the critical interval contributes.
//!
//! One trap: d here ranges over squarefree sqrt-N-smooth d with khat(d) <= K,
//! which includes some d > N (where floor(N/d) = 0 but M(d) >= 1). Pruning the
//! DFS on d <= N — the natural instinct, and what R0's Legendre walk does —
//! drops them and silently changes the answer.

const std = @import("std");
const rs = @import("rs");
const r0 = @import("hkm0.zig");

const Row = struct {
    err_closed: i64 = 0,
    spurious: u64 = 0,
    max_md: u128 = 0,
    max_omega: usize = 0,
    neg_bracket: bool = false,
    d_over_n: u64 = 0,
};

/// DFS over squarefree sqrt-N-smooth d, pruned by khat(d) <= K (NOT by d <= N).
fn walk(N: u64, primes: []const u64, ks: []const usize, mtab: []const u64, i: usize, d: u64, khat: usize, omega: usize, sign: i64, row: *Row) void {
    const m_approx = mtab[mtab.len - 1 - khat];
    const m_true: u64 = N / d;
    if (m_approx < m_true) {
        row.neg_bracket = true;
    } else {
        const bracket = m_approx - m_true;
        row.err_closed += sign * @as(i64, @intCast(bracket));
        row.spurious += bracket;
        if (m_approx > 0) {
            const md = @as(u128, d) * @as(u128, m_approx);
            if (md > row.max_md) row.max_md = md;
        }
    }
    if (d > N) row.d_over_n += 1;
    if (omega > row.max_omega) row.max_omega = omega;

    var j = i;
    while (j < primes.len) : (j += 1) {
        if (khat + ks[j] > mtab.len - 1) continue;
        walk(N, primes, ks, mtab, j + 1, d * primes[j], khat + ks[j], omega + 1, -sign, row);
    }
}

pub fn main(init: std.process.Init) !void {
    var args = try init.minimal.args.iterateAllocator(init.gpa);
    _ = args.skip();
    const gpa = std.heap.page_allocator;

    const usage = "Usage: hkm1 <N> [delta ...]\n";
    const n_str = args.next() orelse return std.debug.print(usage, .{});
    const N = try std.fmt.parseInt(u64, n_str, 10);

    var deltas = try std.ArrayList(f64).initCapacity(gpa, 8);
    defer deltas.deinit(gpa);
    while (args.next()) |a| try deltas.append(gpa, try std.fmt.parseFloat(f64, a));
    if (deltas.items.len == 0) {
        for ([_]f64{ 0.1, 0.05, 0.02, 0.01, 0.005 }) |d| try deltas.append(gpa, d);
    }

    const sq: u64 = std.math.sqrt(N);
    const primes = try rs.basePrimes(gpa, sq);
    defer gpa.free(primes);

    const pi_N = try rs.countInRange(u64, gpa, 2, N + 1, primes);
    const pi_sq = try rs.countInRange(u64, gpa, 2, sq + 1, primes);
    const exact: i64 = @as(i64, @intCast(pi_N)) - @as(i64, @intCast(pi_sq)) + 1;
    const leg = r0.legendre(N, primes, 0, 1, 1);
    std.debug.print("N = {d}  sqrt N = {d}  exact LHS = {d}  (legendre {s})\n", .{ N, sq, exact, if (leg == exact) "MATCH" else "MISMATCH" });
    if (leg != exact) return error.Lemma1Failed;
    std.debug.print("designed delta = log2(N)/sqrt(N) = {d:.5}\n\n", .{@log2(@as(f64, @floatFromInt(N))) / @as(f64, @floatFromInt(sq))});

    std.debug.print("{s:>7} {s:>11} {s:>11} {s:>7} {s:>7} {s:>13} {s:>13} {s:>6} {s:>15}\n", .{ "delta", "err(closed)", "err(direct)", "agree", "1-sided", "S measured", "S bound", "max w", "spurious pairs" });

    for (deltas.items) |delta| {
        std.debug.assert(delta <= 1.0);
        const K = r0.kbar(N, delta);

        const ks = try gpa.alloc(usize, primes.len);
        defer gpa.free(ks);
        for (primes, 0..) |p, i| ks[i] = r0.kbar(p, delta);

        const mtab = try gpa.alloc(u64, K + 1);
        defer gpa.free(mtab);
        {
            var acc: u64 = 0;
            var k: usize = 0;
            while (k <= K) : (k += 1) {
                const lo = std.math.pow(f64, 2.0, @as(f64, @floatFromInt(k)) * delta);
                const hi = std.math.pow(f64, 2.0, @as(f64, @floatFromInt(k + 1)) * delta);
                acc += @as(u64, @intFromFloat(@ceil(hi))) - @as(u64, @intFromFloat(@ceil(lo)));
                mtab[k] = acc;
            }
        }

        var row = Row{};
        walk(N, primes, ks, mtab, 0, 1, 0, 0, 1, &row);

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
        var seg: i64 = 0;
        {
            var k1: usize = 0;
            while (k1 <= K) : (k1 += 1) {
                const c1: i64 = @intCast(mtab[k1] - (if (k1 == 0) @as(u64, 0) else mtab[k1 - 1]));
                if (c1 == 0) continue;
                var k2: usize = 0;
                while (k1 + k2 <= K) : (k2 += 1) seg += c1 * conv[k2];
            }
        }
        const err_direct = seg - exact;

        const s_meas: i128 = @as(i128, @intCast(row.max_md)) - @as(i128, N);
        const nf: f64 = @floatFromInt(N);
        const s_bound = nf * (std.math.pow(f64, 2.0, delta * @as(f64, @floatFromInt(1 + row.max_omega))) - 1.0);

        std.debug.print("{d:>7.4} {d:>11} {d:>11} {s:>7} {s:>7} {d:>13} {d:>13.0} {d:>6} {d:>15}\n", .{ delta, row.err_closed, err_direct, if (row.err_closed == err_direct) "MATCH" else "DIFFER", if (row.neg_bracket) "BROKEN" else "ok", @as(i64, @intCast(s_meas)), s_bound, row.max_omega, row.spurious });
        if (row.err_closed != err_direct) return error.LocalisationIncomplete;
        if (row.neg_bracket) return error.OneSidednessViolated;
    }
}
