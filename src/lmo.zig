//! LMO (Lagarias–Miller–Odlyzko) φ(x,a), building toward sub-linear π(x) with
//! O(x^(1/3)) memory — the path to π(10^18) in seconds.
//!
//! Split the Legendre sum φ(x,a) = Σ_{d | P_a, μ²(d)=1} μ(d)⌊x/d⌋ at d = y = x^(1/3):
//!   φ(x,a) = S1 + S2
//!   S1 (ordinary leaves) = Σ_{n≤y} μ(n)⌊x/n⌋            ← here; O(x^(1/3)), direct
//!   S2 (special leaves)  = Σ_{d>y, sqfree, y-smooth} μ(d)⌊x/d⌋  ← TODO: segmented
//!                          sieve + Fenwick tree, so the O(x^(2/3)) leaves are
//!                          counted incrementally, never stored.
//! meissel.phiOfX is the oracle: S1 + S2 must equal it exactly.

const std = @import("std");

fn icbrt(x: u64) u64 {
    if (x == 0) return 0;
    var r: u64 = @intFromFloat(std.math.pow(f64, @floatFromInt(x), 1.0 / 3.0));
    if (r == 0) r = 1;
    while (r * r * r > x) r -= 1;
    while ((r + 1) * (r + 1) * (r + 1) <= x) r += 1;
    return r;
}

pub const Foundation = struct { s1: i128, a: usize, y: u64 };

/// Ordinary leaves S1 = Σ_{n≤y} μ(n)⌊x/n⌋, with y = ⌊x^(1/3)⌋ and a = π(y).
pub fn ordinaryS1(gpa: std.mem.Allocator, x: u64) !Foundation {
    const y = icbrt(x);
    const m: usize = @intCast(y + 1);
    const comp = try gpa.alloc(bool, m);
    defer gpa.free(comp);
    @memset(comp, false);
    const mu = try gpa.alloc(i8, m);
    defer gpa.free(mu);
    @memset(mu, 1);

    // μ sieve: flip μ for every prime's multiples, zero on square factors.
    var i: u64 = 2;
    while (i <= y) : (i += 1) {
        if (!comp[@intCast(i)]) {
            var j = i;
            while (j <= y) : (j += i) {
                if (j > i) comp[@intCast(j)] = true;
                mu[@intCast(j)] = -mu[@intCast(j)];
            }
            var k = i * i;
            while (k <= y) : (k += i * i) mu[@intCast(k)] = 0;
        }
    }

    var a: usize = 0;
    var n2: u64 = 2;
    while (n2 <= y) : (n2 += 1) {
        if (!comp[@intCast(n2)]) a += 1;
    }

    var s1: i128 = 0;
    var n: u64 = 1;
    while (n <= y) : (n += 1) {
        const mn = mu[@intCast(n)];
        if (mn != 0) s1 += @as(i128, mn) * @as(i128, @intCast(x / n));
    }
    return .{ .s1 = s1, .a = a, .y = y };
}
