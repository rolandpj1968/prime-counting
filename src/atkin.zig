//! Sieve of Atkin (Atkin–Bernstein) — for the head-to-head vs Eratosthenes.
//! Toggle a bit per quadratic-form solution; odd count ⇒ candidate; then remove
//! squareful numbers. Full-array (both sides full-array for a fair comparison).

const std = @import("std");

inline fn toggle(b: []u64, n: u64) void {
    b[@intCast(n >> 6)] ^= @as(u64, 1) << @as(u6, @intCast(n & 63));
}
inline fn get(b: []const u64, n: u64) bool {
    return (b[@intCast(n >> 6)] >> @as(u6, @intCast(n & 63))) & 1 == 1;
}
inline fn clear(b: []u64, n: u64) void {
    b[@intCast(n >> 6)] &= ~(@as(u64, 1) << @as(u6, @intCast(n & 63)));
}
inline fn set(b: []u64, n: u64) void {
    b[@intCast(n >> 6)] |= @as(u64, 1) << @as(u6, @intCast(n & 63));
}

pub fn countAtkin(gpa: std.mem.Allocator, limit: u64) !u64 {
    if (limit < 5) return switch (limit) {
        0, 1 => 0,
        2 => 1,
        else => 2, // 3, 4 → {2,3}
    };
    const b = try gpa.alloc(u64, @intCast(limit / 64 + 1));
    defer gpa.free(b);
    @memset(b, 0);

    // Form 1: 4x²+y², n ≡ 1 or 5 (mod 12)
    var x: u64 = 1;
    while (4 * x * x < limit) : (x += 1) {
        const fx = 4 * x * x;
        var y: u64 = 1;
        while (fx + y * y <= limit) : (y += 1) {
            const n = fx + y * y;
            const r = n % 12;
            if (r == 1 or r == 5) toggle(b, n);
        }
    }
    // Form 2: 3x²+y², n ≡ 7 (mod 12)
    x = 1;
    while (3 * x * x < limit) : (x += 1) {
        const fx = 3 * x * x;
        var y: u64 = 1;
        while (fx + y * y <= limit) : (y += 1) {
            const n = fx + y * y;
            if (n % 12 == 7) toggle(b, n);
        }
    }
    // Form 3: 3x²−y² with x>y, n ≡ 11 (mod 12). For fixed x, n grows as y shrinks.
    x = 2;
    while (2 * x * x + 2 * x - 1 <= limit) : (x += 1) {
        const fx = 3 * x * x;
        var y: u64 = x - 1;
        while (y >= 1) : (y -= 1) {
            const n = fx - y * y;
            if (n > limit) break;
            if (n % 12 == 11) toggle(b, n);
        }
    }
    // Remove squareful numbers: clear multiples of p² for each candidate p.
    var n: u64 = 5;
    while (n * n <= limit) : (n += 1) {
        if (get(b, n)) {
            var k = n * n;
            while (k <= limit) : (k += n * n) clear(b, k);
        }
    }
    // Candidates are exactly the primes ≥ 5; add 2 and 3.
    var c: u64 = 2;
    for (b) |w| c += @popCount(w);
    return c;
}

/// Plain full-array bit Sieve of Eratosthenes, for comparison. bit = composite.
pub fn countEratosthenes(gpa: std.mem.Allocator, limit: u64) !u64 {
    if (limit < 2) return 0;
    const b = try gpa.alloc(u64, @intCast(limit / 64 + 1));
    defer gpa.free(b);
    @memset(b, 0);
    var i: u64 = 2;
    while (i * i <= limit) : (i += 1) {
        if (!get(b, i)) {
            var j = i * i;
            while (j <= limit) : (j += i) set(b, j); // composite = bit set (OR, not XOR!)
        }
    }
    var comp: u64 = 0;
    for (b) |w| comp += @popCount(w);
    return limit - 1 - comp; // primes in [2,limit] = (limit-1) − composites
}
