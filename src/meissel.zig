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

/// φ(x, a): count of n in [1, x] with no prime factor among the first a primes.
fn phi(x: u64, a: usize, primes: []const u64, tbl: *const PiTable) u64 {
    if (x == 0) return 0;
    if (a == 0) return x;
    if (a == 1) return x - x / 2;
    if (a == 2) return x - x / 2 - x / 3 + x / 6;
    if (a == 3) return (x / 30) * 8 + smallphi[@intCast(x % 30)]; // coprime to 2,3,5
    const pa = primes[a - 1];
    if (pa * pa >= x) {
        // no coprime composites ≤ x → φ = 1 + primes in (p_a, x] = 1 + max(0, π(x)−a)
        const pix = tbl.pi(x);
        const au: u64 = @intCast(a);
        return 1 + (if (pix > au) pix - au else 0);
    }
    return phi(x, a - 1, primes, tbl) - phi(x / pa, a - 1, primes, tbl);
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
    return phi(x, a, primes, &tbl);
}

pub fn pi(gpa: std.mem.Allocator, x: u64) !u64 {
    if (x < 2) return 0;
    const cbrt_x = icbrt(x); // ⌊x^(1/3)⌋
    const sqrt_x = common.isqrt(x); // ⌊x^(1/2)⌋
    const primes = try rs.basePrimes(gpa, sqrt_x); // primes ≤ √x (need up to √x for P₂)
    defer gpa.free(primes);

    // a = π(x^(1/3))
    var a: usize = 0;
    for (primes) |p| {
        if (p <= cbrt_x) a += 1 else break;
    }

    // compact π-table up to x^(2/3): used by both the φ leaf cutoff and P₂.
    var tbl = try PiTable.init(gpa, x / cbrt_x);
    defer tbl.deinit(gpa);

    const phi_val = phi(x, a, primes, &tbl);

    // P₂(x,a) = Σ_{cbrt_x < p ≤ sqrt_x} (π(x/p) − π(p) + 1)
    var p2: i64 = 0;
    for (primes, 0..) |p, j| {
        if (p <= cbrt_x) continue;
        if (p > sqrt_x) break;
        const pi_xp: i64 = @intCast(tbl.pi(x / p));
        const pi_p: i64 = @intCast(j + 1); // p is the (j+1)-th prime
        p2 += pi_xp - pi_p + 1;
    }

    const result = @as(i64, @intCast(phi_val)) + @as(i64, @intCast(a)) - 1 - p2;
    return @intCast(result);
}
