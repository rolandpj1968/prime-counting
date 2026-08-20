//! HKM rung 5 — §3.2 (Newton pointwise in Fourier space, no truncation) and
//! §3.3 (partitioning primes to reduce the padding that costs).
//!
//! R4 implements §3.1: Newton in the time domain, truncating to degree K once
//! per r, which costs 2R round-trip transforms on top of the R for the E_m —
//! 3R = 36 at 1e13, and 16.1 s of the rung's 18.2 s.
//!
//! §3.2 removes the round trips. At a FIXED Fourier coordinate the recurrence
//!     r c_r = sum_{m=1..r} (-1)^(m-1) c_{r-m} e_m,  c_0 = 1
//! is a scalar recurrence in r alone: every coordinate is independent. So run
//! it pointwise, accumulate mu-hat-hat = sum_r (-1)^r c_r as you go, and take
//! ONE inverse transform at the end. R+1 transforms instead of 3R.
//! (HKM go further and observe x c'(x) = e(x) c(x), so c = exp(int e/t), which
//! computes the c_r in O(r_max log r_max) rather than O(r_max^2). At
//! r_max = 12 that is not worth a length-13 FFT, so the scalar recurrence is
//! kept — and note that c evaluated at x = -1 is exactly
//! mu-hat = exp(-sum_m E_m/m), the same identity along the other axis.)
//!
//! THE PRICE, AND THE PREDICTION. Truncation is what kept degrees bounded by
//! K. Without it, deg C_r <= r * max_p j_p, so the arrays must be padded to
//! R * max_j ~ 6K instead of 2K — a 4x longer transform to save a factor 3 in
//! transform count, and 4x as many coordinates paying the same O(R^2/2)
//! pointwise recurrence. Prediction before running: §3.2 on the whole prime
//! set LOSES to R4 at these N. The paper agrees in advance — it says only
//! that this "recovered almost the same complexity", and that its value is in
//! enabling §3.3.
//!
//! §3.3 claws the padding back by running §3.2 separately on each prime
//! interval [N^(1/2^(m+1)), N^(1/2^m)] and convolving the pieces at the end.
//! On such an interval at most ~2^(m+1) primes can multiply to <= N while
//! max_j is only K/2^m, so the required padding is ~2K for EVERY interval,
//! independent of m. Also r_max is capped by how many primes the interval
//! actually holds, which for the small-prime intervals is tiny.
//!
//! --parts 1 forces a single interval, i.e. §3.2 as written; the default
//! partitions, i.e. §3.3. Same code path, one knob.

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

/// §3.2 on one prime interval. `ks` must be ascending. Returns the interval's
/// mu-hat factor truncated to degree K, in the time domain.
fn factorViaFourier(gpa: std.mem.Allocator, ks: []const usize, K: usize, stat: *Stats) ![]i64 {
    // r_max: capped both by how many primes are here and by how many of the
    // smallest can multiply while staying at degree <= K.
    var rmax: usize = 0;
    {
        var acc: usize = 0;
        for (ks) |j| {
            if (acc + j > K) break;
            acc += j;
            rmax += 1;
        }
    }
    const out = try gpa.alloc(i64, K + 1);
    @memset(out, 0);
    if (rmax == 0) {
        out[0] = 1;
        return out;
    }
    const maxj = ks[ks.len - 1];
    const D = rmax * maxj; // strict bound on any degree that can arise
    var L: usize = 1;
    while (L < D + 1) L <<= 1;
    stat.rmax_sum += rmax;
    stat.L_max = @max(stat.L_max, L);
    stat.work += L * rmax * (rmax + 1) / 2;

    const E = try gpa.alloc([]u64, rmax + 1);
    defer gpa.free(E);
    for (1..rmax + 1) |m| {
        const e = try gpa.alloc(u64, L);
        @memset(e, 0);
        for (ks) |j| {
            const k = m * j;
            if (k < L) e[k] = ntt.add(e[k], 1);
        }
        ntt.transform(e, false);
        E[m] = e;
    }
    defer for (1..rmax + 1) |m| gpa.free(E[m]);
    stat.transforms += rmax;

    const rinv = try gpa.alloc(u64, rmax + 1);
    defer gpa.free(rinv);
    for (1..rmax + 1) |r| rinv[r] = ntt.inv(@intCast(r));

    const acc = try gpa.alloc(u64, L);
    defer gpa.free(acc);
    var c: [64]u64 = undefined;
    std.debug.assert(rmax < c.len);
    for (0..L) |w| {
        c[0] = 1;
        var sum: u64 = 1; // r = 0 term, sign +
        for (1..rmax + 1) |r| {
            var s: u64 = 0;
            for (1..r + 1) |m| {
                const t = ntt.mul(c[r - m], E[m][w]);
                if (m % 2 == 1) s = ntt.add(s, t) else s = ntt.sub(s, t);
            }
            c[r] = ntt.mul(s, rinv[r]);
            if (r % 2 == 1) sum = ntt.sub(sum, c[r]) else sum = ntt.add(sum, c[r]);
        }
        acc[w] = sum;
    }
    ntt.transform(acc, true);
    stat.transforms += 1;
    for (0..@min(K + 1, L)) |k| out[k] = ntt.toSigned(acc[k]);
    return out;
}

const Stats = struct {
    transforms: usize = 0,
    rmax_sum: usize = 0,
    L_max: usize = 0,
    work: usize = 0,
    parts: usize = 0,
};

pub fn main(init: std.process.Init) !void {
    var args = try init.minimal.args.iterateAllocator(init.gpa);
    _ = args.skip();
    const gpa = std.heap.page_allocator;

    const usage = "Usage: hkm5 <N> [--parts 1] [--no-referee]\n";
    const n_str = args.next() orelse return std.debug.print(usage, .{});
    const N = try std.fmt.parseInt(u64, n_str, 10);
    var one_part = false;
    var referee = true;
    while (args.next()) |a| {
        if (std.mem.eql(u8, a, "--parts")) {
            const v = args.next() orelse return std.debug.print(usage, .{});
            one_part = std.mem.eql(u8, v, "1");
        } else if (std.mem.eql(u8, a, "--no-referee")) referee = false;
    }

    const sq = std.math.sqrt(N);
    const B = @log2(@as(f64, @floatFromInt(N)));
    const delta = B / @as(f64, @floatFromInt(sq));
    const primes = try rs.basePrimes(gpa, sq);
    defer gpa.free(primes);
    var g = try sg.build(gpa, N, delta);
    defer g.deinit(gpa);
    const K = g.K;
    const ks = try gpa.alloc(usize, primes.len);
    defer gpa.free(ks);
    for (primes, 0..) |p, i| ks[i] = g.kbar(p);

    // level m: p in [N^(1/2^(m+1)), N^(1/2^m)), i.e. m = floor(log2(B/log2 p)) - 1
    const cuts = try gpa.alloc(usize, 40);
    defer gpa.free(cuts);
    var nparts: usize = 0;
    if (one_part) {
        cuts[0] = 0;
        nparts = 1;
    } else {
        var prev_level: usize = std.math.maxInt(usize);
        for (primes, 0..) |p, i| {
            const u = B / @log2(@as(f64, @floatFromInt(p)));
            const lvl: usize = @intFromFloat(@floor(@log2(u)));
            if (lvl != prev_level) {
                cuts[nparts] = i;
                nparts += 1;
                prev_level = lvl;
            }
        }
    }

    std.debug.print("N = {d}  K = {d}  pi(sqrt N) = {d}  parts = {d}{s}\n", .{ N, K, primes.len, nparts, if (one_part) "  (section 3.2 as written)" else "  (section 3.3 partition)" });

    const t0 = nowNs();
    var stat = Stats{ .parts = nparts };
    var acc = try gpa.alloc(i64, K + 1);
    defer gpa.free(acc);
    @memset(acc, 0);
    acc[0] = 1;
    var first = true;
    for (0..nparts) |i| {
        const lo = cuts[i];
        const hi = if (i + 1 < nparts) cuts[i + 1] else primes.len;
        const piece = try factorViaFourier(gpa, ks[lo..hi], K, &stat);
        defer gpa.free(piece);
        if (first) {
            @memcpy(acc, piece);
            first = false;
            continue;
        }
        // combine: acc = acc * piece  mod x^(K+1)
        const af = try gpa.alloc(u64, K + 1);
        defer gpa.free(af);
        const bf = try gpa.alloc(u64, K + 1);
        defer gpa.free(bf);
        for (acc, af) |x, *y| y.* = ntt.fromSigned(x);
        for (piece, bf) |x, *y| y.* = ntt.fromSigned(x);
        const prod = try ntt.convolve(gpa, af, bf, K);
        defer gpa.free(prod);
        stat.transforms += 3;
        for (prod, 0..) |v, k| acc[k] = ntt.toSigned(v);
    }
    const t_done = nowNs();

    std.debug.print("transforms {d}   sum of r_max {d}   largest transform 2^{d}   pointwise modmuls {e:.2}\n", .{ stat.transforms, stat.rmax_sum, std.math.log2_int(usize, @max(2, stat.L_max)), @as(f64, @floatFromInt(stat.work)) });
    std.debug.print("time {d:.2}s\n", .{sec(t0, t_done)});

    if (referee) {
        const tr = nowNs();
        const want = try muHatSparse(gpa, ks, K);
        defer gpa.free(want);
        const ts = nowNs();
        var bad: usize = 0;
        var firstbad: usize = 0;
        for (want, acc, 0..) |w, got, i| {
            if (w != got) {
                if (bad == 0) firstbad = i;
                bad += 1;
            }
        }
        std.debug.print("referee (sparse build, {d:.2}s): {d} differing cells   {s}\n", .{ sec(tr, ts), bad, if (bad == 0) "MATCH" else "MISMATCH" });
        if (bad != 0) {
            std.debug.print("  first at k={d}: want {d} got {d}\n", .{ firstbad, want[firstbad], acc[firstbad] });
            return error.Rung5Mismatch;
        }
    }
}
