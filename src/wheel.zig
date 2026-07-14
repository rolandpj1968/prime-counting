//! Wheel — the coordinate map, computed at comptime from a set of wheel primes.
//!   Wheel(&.{})      → all numbers (M=1)
//!   Wheel(&.{2})     → odds        (M=2)
//!   Wheel(&.{2,3})   → mod-6       (M=6)
//!   Wheel(&.{2,3,5}) → mod-30      (M=30)
//!
//! It supplies everything a traversal needs to be coordinate-agnostic: the
//! index<->value map, slot counting, the π formula, and per-prime strike
//! descriptors — a start index (of p²) plus φ cyclic index-deltas, computed by
//! direct index-differencing so it is correct for any wheel with no hand-derived
//! closed forms.

const std = @import("std");

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

pub fn Wheel(comptime wheel_primes: []const u64) type {
    return struct {
        pub const M: u64 = blk: {
            var m: u64 = 1;
            for (wheel_primes) |p| m *= p;
            break :blk m;
        };

        /// φ(M): number of residues coprime to M = spokes of the wheel.
        pub const spokes: usize = blk: {
            var c: usize = 0;
            var r: u64 = 0;
            while (r < M) : (r += 1) if (gcd(r, M) == 1) {
                c += 1;
            };
            break :blk c;
        };

        const S: u64 = spokes; // φ as u64, for index arithmetic

        /// Sorted residues in [0, M) coprime to M (R[0] == 1 for M≥2, == 0 for M==1).
        const R: [spokes]u64 = blk: {
            var arr: [spokes]u64 = undefined;
            var n: usize = 0;
            var r: u64 = 0;
            while (r < M) : (r += 1) if (gcd(r, M) == 1) {
                arr[n] = r;
                n += 1;
            };
            break :blk arr;
        };

        /// gap[i] = distance from residue R[i] to the next coprime residue (wrapping +M).
        const gap: [spokes]u64 = blk: {
            var g: [spokes]u64 = undefined;
            for (0..spokes) |i| {
                const nxt = if (i + 1 < spokes) R[i + 1] else R[0] + M;
                g[i] = nxt - R[i];
            }
            break :blk g;
        };

        /// rank[r] = spoke index of residue r, or `spokes` (sentinel) if not coprime.
        const rank: [M]usize = blk: {
            var t: [M]usize = undefined;
            for (&t) |*x| x.* = spokes;
            for (R, 0..) |res, i| t[@intCast(res)] = i;
            break :blk t;
        };

        const removed: usize = wheel_primes.len; // primes dropped by the wheel (counted separately)
        // π(N) = slotCount(N) - composites + correction, where
        //   correction = removed - (#represented values in {0,1})
        //   {0,1} represented: M==1 → both coprime → 2;  M≥2 → only 1 → 1
        const correction: i64 = @as(i64, @intCast(removed)) - @as(i64, if (M == 1) 2 else 1);

        pub const Prime = struct {
            start: u64, // index of p²
            d: [spokes]u64, // cyclic index-deltas
            idx: u64 = 0, // strike cursor (reset per run)
            spoke: usize = 0, // strike cursor (reset per run)
        };

        pub inline fn index(m: u64) u64 {
            return (m / M) * S + @as(u64, @intCast(rank[@intCast(m % M)]));
        }

        pub fn slotCount(n: u64) u64 {
            var c = (n / M) * S;
            const r = n % M;
            inline for (R) |res| {
                if (res <= r) c += 1;
            }
            return c;
        }

        pub fn piFromComposites(n: u64, composites: u64) u64 {
            const s: i64 = @intCast(slotCount(n));
            const cc: i64 = @intCast(composites);
            return @intCast(s - cc + correction);
        }

        pub fn isWheelPrime(p: u64) bool {
            for (wheel_primes) |wp| if (wp == p) return true;
            return false;
        }

        /// Strike descriptor for prime p: start index (of p²) and φ cyclic
        /// index-deltas, obtained by walking φ+1 coprime multiples and differencing.
        pub fn makePrime(p: u64) Prime {
            var w = p;
            var iprev = index(p * w); // index(p²)
            const start = iprev;
            var d: [spokes]u64 = undefined;
            var j: usize = 0;
            while (j < spokes) : (j += 1) {
                const sp = rank[@intCast(w % M)];
                w += gap[sp];
                const inext = index(p * w);
                d[j] = inext - iprev;
                iprev = inext;
            }
            return .{ .start = start, .d = d };
        }
    };
}
