//! Meissel–Lehmer combinatorial π(x): count primes WITHOUT enumerating them.
//!
//!   π(x) = φ(x, a) + a − 1 − P₂(x, a),   a = π(x^(1/3))
//!
//! where φ(x, a) counts n ≤ x with no prime factor among the first a primes, and
//! P₂(x, a) counts n ≤ x that are a product of exactly two primes both > p_a.
//! With a = π(x^(1/3)), any such n has ≤ 2 prime factors (three would exceed x),
//! which is what makes the identity hold. This is the astronomer's method — the
//! same one Meissel ran by hand to π(10⁹) in 1885.
//!
//! φ is the recursive partial sieve φ(x,a) = φ(x,a−1) − φ(x/p_a, a−1), bottoming
//! at a mod-30 wheel base. P₂ needs π(x/p), got from a prefix-π table up to x^(2/3).

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

// mod-30 wheel base for φ(x, 3): count of n ≤ x coprime to 2,3,5.
// smallphi[r] = # coprime-to-30 residues in [1, r].
const smallphi: [30]u64 = blk: {
    var t: [30]u64 = undefined;
    var r: u64 = 0;
    var c: u64 = 0;
    while (r < 30) : (r += 1) {
        if (gcd(r, 30) == 1 and r >= 1) c += 1;
        t[r] = c;
    }
    break :blk t;
};
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

/// φ(x, a): count of n in [1, x] with no prime factor among the first a primes.
/// `pref` is the prefix-π table (π(i) for i ≤ x^(2/3)), used by the leaf cutoff.
fn phi(x: u64, a: usize, primes: []const u64, pref: []const u32) u64 {
    if (x == 0) return 0;
    if (a == 0) return x;
    if (a == 1) return x - x / 2;
    if (a == 2) return x - x / 2 - x / 3 + x / 6;
    if (a == 3) return (x / 30) * 8 + smallphi[@intCast(x % 30)]; // coprime to 2,3,5
    const pa = primes[a - 1];
    if (pa * pa >= x) {
        // no coprime composites ≤ x → φ = 1 + primes in (p_a, x] = 1 + max(0, π(x)−a)
        const pix: u64 = pref[@intCast(x)];
        const au: u64 = @intCast(a);
        return 1 + (if (pix > au) pix - au else 0);
    }
    return phi(x, a - 1, primes, pref) - phi(x / pa, a - 1, primes, pref);
}

/// π(i) for all i ≤ limit (prefix counts).
fn prefixPi(gpa: std.mem.Allocator, limit: u64) ![]u32 {
    const m: usize = @intCast(limit + 1);
    const c = try gpa.alloc(bool, m);
    defer gpa.free(c);
    @memset(c, false);
    var i: usize = 2;
    while (i * i < m) : (i += 1) {
        if (!c[i]) {
            var j = i * i;
            while (j < m) : (j += i) c[j] = true;
        }
    }
    const pref = try gpa.alloc(u32, m);
    var cnt: u32 = 0;
    var k: usize = 0;
    while (k < m) : (k += 1) {
        if (k >= 2 and !c[k]) cnt += 1;
        pref[k] = cnt;
    }
    return pref;
}

pub fn pi(gpa: std.mem.Allocator, x: u64) !u64 {
    if (x < 2) return 0;
    const x13 = icbrt(x);
    const x12 = common.isqrt(x);
    const primes = try rs.basePrimes(gpa, x12); // primes ≤ √x (need up to √x for P₂)
    defer gpa.free(primes);

    // a = π(x^(1/3))
    var a: usize = 0;
    for (primes) |p| {
        if (p <= x13) a += 1 else break;
    }

    // prefix-π up to x^(2/3): used by both the φ leaf cutoff and P₂.
    const pref = try prefixPi(gpa, x / x13);
    defer gpa.free(pref);

    const phi_val = phi(x, a, primes, pref);

    // P₂(x,a) = Σ_{x13 < p ≤ x12} (π(x/p) − π(p) + 1)
    var p2: i64 = 0;
    for (primes, 0..) |p, j| {
        if (p <= x13) continue;
        if (p > x12) break;
        const pi_xp: i64 = pref[@intCast(x / p)];
        const pi_p: i64 = @intCast(j + 1); // p is the (j+1)-th prime
        p2 += pi_xp - pi_p + 1;
    }

    const result = @as(i64, @intCast(phi_val)) + @as(i64, @intCast(a)) - 1 - p2;
    return @intCast(result);
}
