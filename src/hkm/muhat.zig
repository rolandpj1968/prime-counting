//! mu-hat construction, shared by the rungs so they cannot drift apart.
//!
//! Two routes, both exact and both verified against each other at every N the
//! ladder has run:
//!
//!   sparse  — fold in (1 - x^j) one prime at a time, truncating at K.
//!             O(K pi(sqrt N)), i.e. O(N/log N). Wins below ~1e11.
//!   fourier — HKM §3.1 Newton's identities, run pointwise in Fourier space
//!             per §3.2, partitioned by prime size per §3.3. O(R K log K).
//!             Wins above, and uses 2.8x less memory at 1e13.
//!
//! `auto` picks by K, and the crossover was measured rather than estimated —
//! a first guess of K ~ 3e5 was wrong by 2x. Segmented-stage seconds,
//! fourier against sparse: K = 316227 (1e11) 1.77 / 1.07; K = 547722 (3e11)
//! 2.98 / 3.09; K = 774596 (6e11) 4.45 / 5.99; K = 1000000 (1e12) 4.56 / 9.81.
//! The lines cross at K ~ 5.5e5, i.e. N ~ 3e11.
//!
//! Only `ks` (the per-prime shifts kbar(p)) is needed, never the primes
//! themselves: the §3.3 partition level of a prime is floor(log2(K/ks[i])),
//! since log2 p ~ ks*Delta and log2 N ~ K*Delta make N^(1/u) with u = K/ks.

const std = @import("std");
const ntt = @import("ntt.zig");

pub const Mode = enum { auto, sparse, fourier };

pub const Stats = struct {
    part_rmax: [16]u32 = [_]u32{0} ** 16,
    part_L: [16]u32 = [_]u32{0} ** 16,
    part_n: usize = 0,
    transforms: usize = 0,
    rmax_sum: usize = 0,
    L_max: usize = 0,
    work: usize = 0,
    parts: usize = 0,
    route: []const u8 = "",
};

/// Fold in (1 - x^j) one prime at a time, truncated at K.
pub fn sparse(gpa: std.mem.Allocator, ks: []const usize, K: usize) ![]i64 {
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

/// §3.2 on one prime interval: Newton pointwise in Fourier space, no
/// truncation, one inverse transform. Returns the interval's factor
/// truncated to degree K.
fn factorViaFourier(gpa: std.mem.Allocator, ks: []const usize, K: usize, stat: *Stats) ![]i64 {
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
    if (stat.part_n < 16) {
        stat.part_rmax[stat.part_n] = @intCast(rmax);
        stat.part_L[stat.part_n] = @intCast(@min(L, std.math.maxInt(u32)));
        stat.part_n += 1;
    }
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
        var sum: u64 = 1;
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

/// §3.3: split the primes by size level and convolve the per-level factors.
/// `one_part` forces a single interval, i.e. §3.2 as written.
/// `ratio` is the size ratio of a §3.3 prime interval: 2.0 is HKM's
/// [N^(1/2^(m+1)), N^(1/2^m)]. It is the space/time knob. A part needs
/// L >= D = r_max * max_j ~ K * (log p_hi / log p_lo) = K * ratio, so peak
/// memory per part falls with the ratio while the number of parts — and hence
/// the transform count — rises as log(log2 N)/log(ratio).
pub fn fourierR(gpa: std.mem.Allocator, ks: []const usize, K: usize, one_part: bool, ratio: f64, stat: *Stats) ![]i64 {
    var cuts: [64]usize = undefined;
    var nparts: usize = 0;
    if (one_part or ks.len == 0) {
        cuts[0] = 0;
        nparts = 1;
    } else {
        var prev: usize = std.math.maxInt(usize);
        for (ks, 0..) |j, i| {
            const u = @as(f64, @floatFromInt(K)) / @as(f64, @floatFromInt(@max(1, j)));
            const lvl: usize = @intFromFloat(@max(0.0, @floor(@log(@max(1.0, u)) / @log(ratio))));
            if (lvl != prev) {
                cuts[nparts] = i;
                nparts += 1;
                prev = lvl;
            }
        }
    }
    stat.parts = nparts;

    var acc: ?[]i64 = null;
    for (0..nparts) |i| {
        const lo = cuts[i];
        const hi = if (i + 1 < nparts) cuts[i + 1] else ks.len;
        const piece = try factorViaFourier(gpa, ks[lo..hi], K, stat);
        if (acc == null) {
            acc = piece;
            continue;
        }
        defer gpa.free(piece);
        const a = acc.?;
        const af = try gpa.alloc(u64, K + 1);
        defer gpa.free(af);
        const bf = try gpa.alloc(u64, K + 1);
        defer gpa.free(bf);
        for (a, af) |x, *y| y.* = ntt.fromSigned(x);
        for (piece, bf) |x, *y| y.* = ntt.fromSigned(x);
        const prod = try ntt.convolve(gpa, af, bf, K);
        defer gpa.free(prod);
        stat.transforms += 3;
        for (prod, 0..) |v, k| a[k] = ntt.toSigned(v);
    }
    return acc.?;
}

/// Measured crossover: sparse below, fourier above. See the table above.
pub const AUTO_K = 550_000;

pub fn fourier(gpa: std.mem.Allocator, ks: []const usize, K: usize, one_part: bool, stat: *Stats) ![]i64 {
    return fourierR(gpa, ks, K, one_part, 2.0, stat);
}

pub fn build(gpa: std.mem.Allocator, ks: []const usize, K: usize, mode: Mode, stat: *Stats) ![]i64 {
    const use_fourier = switch (mode) {
        .sparse => false,
        .fourier => true,
        .auto => K >= AUTO_K,
    };
    if (!use_fourier) {
        stat.route = "sparse";
        return sparse(gpa, ks, K);
    }
    stat.route = "fourier(3.3)";
    return fourier(gpa, ks, K, false, stat);
}
