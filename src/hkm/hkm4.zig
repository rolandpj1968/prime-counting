//! HKM rung 4 — mu-hat by Newton's identities, in O(log^2 N) convolutions.
//!
//! R2 measured the naive build at O(K pi(sqrt N)) = O(N/log N), growing ~9.7x
//! per decade against the correction's ~4.3x and overtaking it near 1e15.
//! This rung replaces it.
//!
//! As a polynomial, mu-hat(x) = prod_p (1 - x^j_p) with j_p = kbar(p). Setting
//! t_p = x^j_p, that is prod (1 - t_p) = sum_r (-1)^r e_r(t), so with
//!
//!     C_r = e_r(t) = sum over r-subsets of x^(j_{p1} + ... + j_{pr})
//!     E_m = power sum = sum_p x^(m j_p)
//!
//! Newton's identities give  r C_r = sum_{m=1..r} (-1)^(m-1) C_{r-m} * E_m,
//! and mu-hat = sum_{r=0..R} (-1)^r C_r. R = omega_max: past it every
//! r-subset already exceeds degree K, so C_r truncates to zero.
//!
//! WHY IT IS WRITTEN IN THE FREQUENCY DOMAIN. R3 measured, at K = 1e6, one
//! convolve at 0.951 s and one bare transform at 0.296 s. The identities need
//! R(R+1)/2 = 66 products at R = 11; as 66 convolve calls that is 63 s
//! against the naive build's 10.19 s — six times WORSE than the thing being
//! replaced. But those 66 products involve only 3R distinct transforms if
//! each E_m and each C_j is transformed once and reused, which is ~10 s. The
//! batched form is not an optimisation here, it is the difference between the
//! rung working and not working.
//!
//! Truncation to degree K still happens in the time domain once per r: higher
//! degrees would wrap around the cyclic transform and corrupt low
//! coefficients. That is why the count is 3R (E forward, C inverse, C
//! forward) and not 2R.
//!
//! ONE MODULUS IS PROVABLY ENOUGH, which R3 left open. Every C_r[k] is a
//! count of r-subsets, all of one sign, so there is no cancellation to lean
//! on. But sum over r and k of C_r[k] counts each admissible squarefree
//! sqrt-N-smooth n exactly once, and R1's bound says such an n satisfies
//! n < N 2^(D(1+omega)) = N + S. So every C_r[k] <= N + S < 2N, and a single
//! Goldilocks prime carries the whole computation for any N below ~4.6e18 —
//! asserted at startup rather than hoped for. No CRT layer is built.

const std = @import("std");
const rs = @import("rs");
const sg = @import("seg.zig");
const ntt = @import("ntt.zig");

fn nowNs() u64 {
    var ts: std.os.linux.timespec = undefined;
    _ = std.os.linux.clock_gettime(.MONOTONIC, &ts);
    return @as(u64, @intCast(ts.sec)) * 1_000_000_000 + @as(u64, @intCast(ts.nsec));
}
fn sec(a: u64, b: u64) f64 {
    return @as(f64, @floatFromInt(b - a)) / 1e9;
}

/// R2's route, kept only as the referee: fold in (1 - x^j) one prime at a time.
fn muHatSparse(gpa: std.mem.Allocator, ks: []const usize, K: usize) ![]i64 {
    const a = try gpa.alloc(i64, K + 1);
    @memset(a, 0);
    a[0] = 1;
    for (ks) |j| {
        if (j > K) continue;
        var t: usize = K + 1;
        while (t > j) {
            t -= 1;
            a[t] -= a[t - j];
        }
    }
    return a;
}

pub fn main(init: std.process.Init) !void {
    var args = try init.minimal.args.iterateAllocator(init.gpa);
    _ = args.skip();
    const gpa = std.heap.page_allocator;

    const usage = "Usage: hkm4 <N> [--no-referee]\n";
    const n_str = args.next() orelse return std.debug.print(usage, .{});
    const N = try std.fmt.parseInt(u64, n_str, 10);
    var referee = true;
    while (args.next()) |a| {
        if (std.mem.eql(u8, a, "--no-referee")) referee = false;
    }

    const sq = std.math.sqrt(N);
    const delta = @log2(@as(f64, @floatFromInt(N))) / @as(f64, @floatFromInt(sq));
    const primes = try rs.basePrimes(gpa, sq);
    defer gpa.free(primes);
    var g = try sg.build(gpa, N, delta);
    defer g.deinit(gpa);
    const K = g.K;
    const ks = try gpa.alloc(usize, primes.len);
    defer gpa.free(ks);
    for (primes, 0..) |p, i| {
        ks[i] = g.kbar(p);
        std.debug.assert(ks[i] >= 1); // else khat never grows and R is unbounded
    }

    // R = omega_max, and the modulus argument that goes with it
    var R: usize = 0;
    {
        var acc: usize = 0;
        for (ks) |k| {
            if (acc + k > K) break;
            acc += k;
            R += 1;
        }
    }
    const s_bound = @as(f64, @floatFromInt(N)) * (std.math.pow(f64, 2.0, delta * @as(f64, @floatFromInt(1 + R))) - 1.0);
    const cmax_bound: f64 = @as(f64, @floatFromInt(N)) + s_bound;
    if (cmax_bound >= @as(f64, @floatFromInt(ntt.P / 2))) return error.ModulusTooSmall;

    var L: usize = 1;
    while (L < 2 * (K + 1)) L <<= 1;
    std.debug.print("N = {d}  sqrt N = {d}  K = {d}  pi(sqrt N) = {d}  R = {d}  transform 2^{d}\n", .{ N, sq, K, primes.len, R, std.math.log2_int(usize, L) });
    std.debug.print("modulus: every C_r[k] <= N+S = {e:.3} < P/2 = {e:.3}  ({e:.3}x headroom, one prime, no CRT)\n", .{ cmax_bound, @as(f64, @floatFromInt(ntt.P / 2)), @as(f64, @floatFromInt(ntt.P / 2)) / cmax_bound });
    std.debug.print("frequency-domain arrays: {d} x {d:.1} MB = {d:.1} MB\n\n", .{ 2 * R, @as(f64, @floatFromInt(L * 8)) / 1e6, @as(f64, @floatFromInt(2 * R * L * 8)) / 1e6 });

    const t0 = nowNs();

    // ---- E_m in the time domain, then transformed once each.
    // E_m[k] = #{ p : m*kbar(p) = k }. For m >= 3 only primes below N^(1/m)
    // can land at all, so these get sparse fast — but they are still
    // transformed, because R3 measured sparse-times-dense losing to the NTT
    // right down to m = 3.
    const Ehat = try gpa.alloc([]u64, R + 1);
    defer gpa.free(Ehat);
    for (1..R + 1) |m| {
        const e = try gpa.alloc(u64, L);
        @memset(e, 0);
        for (ks) |j| {
            const k = m * j;
            if (k <= K) e[k] = ntt.add(e[k], 1);
        }
        ntt.transform(e, false);
        Ehat[m] = e;
    }
    defer for (1..R + 1) |m| gpa.free(Ehat[m]);
    const t_e = nowNs();

    // ---- Newton. Chat[0] is the transform of delta_0, i.e. all ones, so it
    // is never materialised: the m = r term is just Ehat[r].
    const Chat = try gpa.alloc([]u64, R + 1);
    defer gpa.free(Chat);
    const muhat = try gpa.alloc(i64, K + 1);
    defer gpa.free(muhat);
    @memset(muhat, 0);
    muhat[0] = 1; // r = 0 term

    const acc = try gpa.alloc(u64, L);
    defer gpa.free(acc);
    var ns_point: u64 = 0;
    var ns_trans: u64 = 0;
    var cmax_seen: i64 = 0;

    for (1..R + 1) |r| {
        const ta = nowNs();
        @memcpy(acc, Ehat[r]); // m = r, sign (-1)^(r-1) applied below
        if (r % 2 == 0) for (acc) |*v| {
            v.* = ntt.sub(0, v.*);
        };
        var m: usize = 1;
        while (m < r) : (m += 1) {
            const c = Chat[r - m];
            const e = Ehat[m];
            if (m % 2 == 1) {
                for (acc, c, e) |*v, cv, ev| v.* = ntt.add(v.*, ntt.mul(cv, ev));
            } else {
                for (acc, c, e) |*v, cv, ev| v.* = ntt.sub(v.*, ntt.mul(cv, ev));
            }
        }
        const tb = nowNs();
        ntt.transform(acc, true);
        const tc = nowNs();

        const rinv = ntt.inv(@intCast(r));
        const cr = try gpa.alloc(u64, L);
        @memset(cr[K + 1 ..], 0); // truncate: degrees > K would wrap next round
        for (0..K + 1) |k| cr[k] = ntt.mul(acc[k], rinv);
        for (0..K + 1) |k| {
            const v = ntt.toSigned(cr[k]);
            cmax_seen = @max(cmax_seen, v);
            if (r % 2 == 1) muhat[k] -= v else muhat[k] += v;
        }
        const td = nowNs();
        ntt.transform(cr, false);
        Chat[r] = cr;
        const te = nowNs();
        ns_point += (tb - ta) + (td - tc);
        ns_trans += (tc - tb) + (te - td);
    }
    for (1..R + 1) |r| gpa.free(Chat[r]);
    const t_newton = nowNs();

    std.debug.print("time: E build+transform {d:.2}s   Newton {d:.2}s  (pointwise {d:.2}s, transforms {d:.2}s)   total {d:.2}s\n", .{ sec(t0, t_e), sec(t_e, t_newton), @as(f64, @floatFromInt(ns_point)) / 1e9, @as(f64, @floatFromInt(ns_trans)) / 1e9, sec(t0, t_newton) });
    std.debug.print("transforms: {d} for E + {d} inverse + {d} forward = {d}   (vs {d} if written as convolve calls)\n", .{ R, R, R, 3 * R, 3 * R * (R + 1) / 2 });
    std.debug.print("max C_r[k] seen = {d}  ({e:.3}x under the N+S bound)\n", .{ cmax_seen, cmax_bound / @as(f64, @floatFromInt(@max(1, cmax_seen))) });

    if (referee) {
        const tr = nowNs();
        const want = try muHatSparse(gpa, ks, K);
        defer gpa.free(want);
        const ts = nowNs();
        var bad: usize = 0;
        var first: usize = 0;
        for (want, muhat, 0..) |w, got, i| {
            if (w != got) {
                if (bad == 0) first = i;
                bad += 1;
            }
        }
        std.debug.print("referee (sequential sparse build, {d:.2}s): {d} differing cells   {s}", .{ sec(tr, ts), bad, if (bad == 0) "MATCH\n" else "MISMATCH" });
        if (bad != 0) {
            std.debug.print(" first at k={d}: want {d} got {d}\n", .{ first, want[first], muhat[first] });
            return error.NewtonMismatch;
        }
        std.debug.print("speedup vs naive build: {d:.2}x\n", .{sec(tr, ts) / sec(t0, t_newton)});
    }
}
