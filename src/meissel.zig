//! Meissel–Lehmer combinatorial π(x): count primes WITHOUT enumerating them.
//!
//!   π(x) = φ(x, a) + a − 1 − P₂(x, a),   a = π(x^(1/3))
//!
//! φ(x,a) counts n ≤ x with no prime factor among the first a primes; P₂(x,a)
//! counts n ≤ x that are a product of two primes both > p_a. With a = π(x^(1/3))
//! such n has ≤ 2 prime factors, which is what makes the identity hold — the
//! astronomer's method Meissel ran by hand to π(10⁹) in 1885.
//!
//! φ is the recursion φ(x,a)=φ(x,a−1)−φ(x/p_a,a−1) with a mod-30 wheel base and
//! the leaf cutoff φ=1+max(0,π(x)−a) when p_a²≥x (LMO Stage A: π comes from a
//! COMPACT bit-sieve+checkpoints table, ~16× smaller than a π-per-integer array).
//!
//! The identity holds for any y ≥ x^(1/3), so y = α·x^(1/3) is a free knob. Two
//! consumers want the π-table: φ's cutoff (v ≤ y², grows with α) and P₂ (x/y,
//! shrinks). Un-capped they cross at α=1 — the classical exponent is a minimax.
//! CAPPING the cutoff at v ≤ z = x/y deletes the y² consumer, so the table is just
//! z = x^(2/3)/α and α becomes a real memory knob. See piWithY / the y-sweep.

const std = @import("std");
const common = @import("common.zig");
const rs = @import("rangesieve.zig");

/// floor(x^(1/3)).
fn icbrt(x: u64) u64 {
    if (x == 0) return 0;
    var r: u64 = @intFromFloat(std.math.pow(f64, @floatFromInt(x), 1.0 / 3.0));
    if (r == 0) r = 1;
    while (r * r * r > x) r -= 1;
    while ((r + 1) * (r + 1) * (r + 1) <= x) r += 1;
    return r;
}

fn gcd(a: u64, b: u64) u64 {
    var x = a;
    var y = b;
    while (y != 0) {
        const t = x % y;
        x = y;
        y = t;
    }
    return x;
}

// mod-30 wheel base for φ(x, 3): count of n ≤ x coprime to 2,3,5.
// smallphi[r] = # coprime-to-30 residues in [1, r].
const smallphi: [30]u64 = blk: {
    var t: [30]u64 = undefined;
    var r: u64 = 0;
    var c: u64 = 0;
    while (r < 30) : (r += 1) {
        if (r >= 1 and gcd(r, 30) == 1) c += 1;
        t[r] = c;
    }
    break :blk t;
};

/// Compact prime-counting table over [0, limit]: a bit per integer (prime = 1)
/// plus a per-word prefix-count checkpoint, so π(y) is one popcount. ~16× less
/// memory than a π-per-integer array; O(1) queries.
const PiTable = struct {
    bits: []u64, // bit i set ⇒ i is prime
    ckpt: []u64, // ckpt[w] = # primes in words [0, w) = π(64·w − 1)

    fn init(gpa: std.mem.Allocator, limit: u64) !PiTable {
        const nwords: usize = @intCast(limit / 64 + 1);
        const bits = try gpa.alloc(u64, nwords);
        @memset(bits, ~@as(u64, 0)); // all candidate-prime
        bits[0] &= ~@as(u64, 0b11); // 0, 1 not prime
        // clear padding past limit
        var b = limit + 1;
        while (b < nwords * 64) : (b += 1) bits[@intCast(b >> 6)] &= ~(@as(u64, 1) << @as(u6, @intCast(b & 63)));
        // sieve: clear composites
        var i: u64 = 2;
        while (i * i <= limit) : (i += 1) {
            if ((bits[@intCast(i >> 6)] >> @as(u6, @intCast(i & 63))) & 1 == 1) {
                var j = i * i;
                while (j <= limit) : (j += i) bits[@intCast(j >> 6)] &= ~(@as(u64, 1) << @as(u6, @intCast(j & 63)));
            }
        }
        const ckpt = try gpa.alloc(u64, nwords + 1);
        var c: u64 = 0;
        for (bits, 0..) |w, wi| {
            ckpt[wi] = c;
            c += @popCount(w);
        }
        ckpt[nwords] = c;
        return .{ .bits = bits, .ckpt = ckpt };
    }

    fn deinit(self: *PiTable, gpa: std.mem.Allocator) void {
        gpa.free(self.bits);
        gpa.free(self.ckpt);
    }

    /// π(y) = number of primes ≤ y.
    fn pi(self: *const PiTable, y: u64) u64 {
        const wy: usize = @intCast(y >> 6);
        const r: u6 = @intCast(y & 63);
        const mask: u64 = if (r == 63) ~@as(u64, 0) else (@as(u64, 1) << (r + 1)) - 1;
        return self.ckpt[wy] + @popCount(self.bits[wy] & mask);
    }
};

pub const Counters = struct { leaves: u64 = 0, nodes: u64 = 0 };

/// φ(x, a): count of n in [1, x] with no prime factor among the first a primes.
///
/// `zmax` caps the cutoff: we may only take the π-leaf when x ≤ zmax, since that
/// is all the table holds. Tightening a cutoff is always CORRECT — the recursion
/// is valid at every node — it just trades tree size for table size.
fn phi(x: u64, a: usize, primes: []const u64, tbl: *const PiTable, zmax: u64, c: ?*Counters) u64 {
    if (c) |cc| cc.nodes += 1;
    if (x == 0) return 0;
    if (a == 0) return x;
    if (a == 1) return x - x / 2;
    if (a == 2) return x - x / 2 - x / 3 + x / 6;
    if (a == 3) return (x / 30) * 8 + smallphi[@intCast(x % 30)]; // coprime to 2,3,5
    const pa = primes[a - 1];
    if (pa * pa >= x and x <= zmax) {
        // no coprime composites ≤ x → φ = 1 + primes in (p_a, x] = 1 + max(0, π(x)−a)
        if (c) |cc| cc.leaves += 1;
        const pix = tbl.pi(x);
        const au: u64 = @intCast(a);
        return 1 + (if (pix > au) pix - au else 0);
    }
    return phi(x, a - 1, primes, tbl, zmax, c) - phi(x / pa, a - 1, primes, tbl, zmax, c);
}

/// φ(x, π(x^(1/3))) — the oracle LMO's S1+S2 must reproduce.
pub fn phiOfX(gpa: std.mem.Allocator, x: u64) !u64 {
    if (x < 2) return 0;
    const cbrt_x = icbrt(x);
    const sqrt_x = common.isqrt(x);
    const primes = try rs.basePrimes(gpa, sqrt_x);
    defer gpa.free(primes);
    var a: usize = 0;
    for (primes) |p| {
        if (p <= cbrt_x) a += 1 else break;
    }
    var tbl = try PiTable.init(gpa, x / cbrt_x);
    defer tbl.deinit(gpa);
    return phi(x, a, primes, &tbl, std.math.maxInt(u64), null);
}

/// Meissel–Lehmer with y as a free knob (a = π(y)), for the y-sweep experiment.
///
/// The identity needs every y-rough n ≤ x to have ≤ 2 prime factors, i.e. y ≥ x^(1/3);
/// below that you'd need P₃ (Lehmer). Above it the identity still holds, but the
/// π-table must serve BOTH consumers:
///   • φ's cutoff leaf φ(v,b), p_b² ≥ v  ⇒  v ≤ p_b² ≤ y²      (grows with y)
///   • P₂'s π(x/p) for p ∈ (y, √x]       ⇒  argument ≤ x/y      (shrinks with y)
/// so the bound is max(y², x/y) — minimised exactly at y = x^(1/3), where both are x^(2/3).
///
/// `capped` adds "and v ≤ z" to the cutoff (z = x/y). Tightening a cutoff is always
/// correct — the recursion is valid at every node — so this only trades tree for table.
/// It deletes the y² consumer, leaving the bound at z = x/y, which now SHRINKS with α.
pub fn piWithY(gpa: std.mem.Allocator, x: u64, y: u64, capped: bool) !YResult {
    const sqrt_x = common.isqrt(x);
    const primes = try rs.basePrimes(gpa, sqrt_x);
    defer gpa.free(primes);

    var a: usize = 0;
    for (primes) |p| {
        if (p <= y) a += 1 else break;
    }

    // z = x/y serves BOTH consumers once the cutoff is capped at z: φ's leaves
    // (by construction) and P₂'s π(x/p), p > y. No y² term — it shrinks with α.
    const limit = @min(x, if (capped) x / y else @max(y *| y, x / y));
    const zmax = if (capped) limit else std.math.maxInt(u64);
    var t0 = common.nowNs();
    var tbl = try PiTable.init(gpa, limit);
    defer tbl.deinit(gpa);
    const build_ns = common.nowNs() - t0;

    var ctr = Counters{};
    t0 = common.nowNs();
    const phi_val = phi(x, a, primes, &tbl, zmax, &ctr);
    const phi_ns = common.nowNs() - t0;

    t0 = common.nowNs();
    var p2: i64 = 0;
    for (primes, 0..) |p, j| {
        if (p <= y) continue;
        if (p > sqrt_x) break;
        p2 += @as(i64, @intCast(tbl.pi(x / p))) - @as(i64, @intCast(j + 1)) + 1;
    }
    const p2_ns = common.nowNs() - t0;

    const result = @as(i64, @intCast(phi_val)) + @as(i64, @intCast(a)) - 1 - p2;
    return .{
        .pi = @intCast(result),
        .y = y,
        .a = a,
        .limit = limit,
        .bytes = (limit / 64 + 1) * 16, // bits[] + ckpt[], 8 bytes per word each
        .leaves = ctr.leaves,
        .nodes = ctr.nodes,
        .build_ns = build_ns,
        .phi_ns = phi_ns,
        .p2_ns = p2_ns,
    };
}

pub const YResult = struct {
    pi: u64,
    y: u64,
    a: usize,
    limit: u64,
    bytes: usize,
    leaves: u64,
    nodes: u64,
    build_ns: u64,
    phi_ns: u64,
    p2_ns: u64,
    pub fn totalNs(self: YResult) u64 {
        return self.build_ns + self.phi_ns + self.p2_ns;
    }
};

/// Tuning constant: y = α·x^(1/3), as the 3/2 numerator over ALPHA_DEN.
/// α = 3/2 with the capped cutoff measured strictly better than the classical
/// α = 1 on BOTH axes at x = 10^11 (100.0 ms / 3.4 MB vs 105.5 ms / 5.1 MB) —
/// the capped table is x^(2/3)/α, so raising α shrinks it while φ barely moves.
const ALPHA_NUM = 3;
const ALPHA_DEN = 2;

pub fn pi(gpa: std.mem.Allocator, x: u64) !u64 {
    if (x < 2) return 0;
    // Integer y keeps the y ≥ x^(1/3) correctness floor exact (below it the
    // identity needs P₃): icbrt(x) already satisfies it, and we only go up.
    const y = icbrt(x) * ALPHA_NUM / ALPHA_DEN;
    const r = try piWithY(gpa, x, y, true);
    return r.pi;
}
