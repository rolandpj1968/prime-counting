//! Int-parameterized range sieve: count primes in [lo, hi) with the coordinate
//! type a comptime parameter (u64 for the native-fast case, u128 to address
//! past 2^64 — the whole point). Base primes stay u64 (√(10^38) still < 2^64...
//! well past anything feasible); only positions/multiples need the wide type.
//!
//! Correctness at height comes from within: provable algorithm + no overflow +
//! validate this exact code at low N against known π + cross-impl agreement.

const std = @import("std");

/// Sorted primes ≤ limit (simple sieve). Base primes are always u64.
pub fn basePrimes(gpa: std.mem.Allocator, limit: u64) ![]u64 {
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
    var n: usize = 0;
    for (c[2..]) |x| {
        if (!x) n += 1;
    }
    const p = try gpa.alloc(u64, n);
    var k: usize = 0;
    var v: usize = 2;
    while (v < m) : (v += 1) {
        if (!c[v]) {
            p[k] = @intCast(v);
            k += 1;
        }
    }
    return p;
}

/// Count primes in [lo, hi). `primes` must include all primes ≤ √(hi-1).
/// Int = u64 or u128 (the coordinate/value type that goes beyond 2^64).
pub fn countInRange(comptime Int: type, gpa: std.mem.Allocator, lo: Int, hi: Int, primes: []const u64) !u64 {
    std.debug.assert(hi > lo);
    const width: usize = @intCast(hi - lo);
    const nwords = (width + 63) / 64;
    const w = try gpa.alloc(u64, nwords);
    defer gpa.free(w);
    @memset(w, 0);

    for (primes) |p| {
        const P: Int = p;
        if (P * P >= hi) break; // p > √(hi-1): nothing more to strike
        var pos: Int = ((lo + P - 1) / P) * P; // first multiple of p ≥ lo
        if (pos < P * P) pos = P * P; // ...but not below p²
        while (pos < hi) : (pos += P) {
            const b: usize = @intCast(pos - lo);
            w[b >> 6] |= @as(u64, 1) << @as(u6, @intCast(b & 63));
        }
    }

    // unmarked bits in [0, width) are prime; exclude 0 and 1 if in range
    var comp: u64 = 0;
    const full = width / 64;
    for (w[0..full]) |word| comp += @popCount(word);
    const rem = width % 64;
    if (rem != 0) {
        const mask = (@as(u64, 1) << @as(u6, @intCast(rem))) - 1;
        comp += @popCount(w[full] & mask);
    }
    var count: u64 = width - comp;
    if (lo == 0) {
        count -= 2; // 0 and 1 are unmarked but not prime
    } else if (lo == 1) {
        count -= 1;
    }
    return count;
}
