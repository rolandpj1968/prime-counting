//! SegmentedSieve — the segmented traversal, generic over a backing Store AND a
//! segment span. A segment is just a small bit-array (a Store sized to `span`
//! integers, reused per block), so segmentation composes with the same stores
//! the naive sieve uses. Three orthogonal axes: coordinate map × traversal ×
//! store.
//!
//! The knob is `span` = integers per segment (NOT bytes) — so the store decides
//! its own footprint. That makes the []u64-vs-[]bool comparison fair and turns
//! it into a direct measurement of cache residency:
//!   span 262144 integers  =>  []u64 segment 32 KiB (fits L1d)
//!                             []bool segment 256 KiB (spills to L2)
//! Same algorithm, same span; the gap is purely "does the segment fit in L1".
//!
//! Counting is folded in (segments are discarded as we go). The segment Store
//! only needs clearAll/set/countComposites — `get` is unused here (base primes
//! come from a separate small sieve).

const std = @import("std");
const common = @import("common.zig");

const Prime = struct { p: u64, next: u64 };

pub fn SegmentedSieve(comptime Store: type, comptime span: u64) type {
    return struct {
        const cap: u64 = span - 1; // store "max index" for a span-slot store
        const seg_bytes = Store.footprintBytes(cap); // comptime-known segment size

        pub const name = std.fmt.comptimePrint("segmented[{s}], span {d}, seg {d} KiB", .{ Store.name, span, seg_bytes / 1024 });

        pub const State = struct {
            seg: Store.State,
            primes: []Prime,
            composites: u64 = 0,
        };

        pub fn footprintBytes(n: u64) usize {
            const lim = common.isqrt(n);
            const est_primes: usize = if (lim < 3) 2 else @intFromFloat(1.2 * @as(f64, @floatFromInt(lim)) / @log(@as(f64, @floatFromInt(lim))));
            return seg_bytes + est_primes * @sizeOf(Prime);
        }

        pub fn init(gpa: std.mem.Allocator, n: u64) !State {
            const primes = try collectBasePrimes(gpa, common.isqrt(n));
            const seg = try Store.init(gpa, cap);
            return .{ .seg = seg, .primes = primes };
        }

        pub fn deinit(st: *State, gpa: std.mem.Allocator) void {
            Store.deinit(&st.seg, gpa);
            gpa.free(st.primes);
        }

        pub fn sieve(st: *State, n: u64) void {
            for (st.primes) |*pr| pr.next = pr.p * pr.p; // reset cursors for this run
            var comp: u64 = 0;
            var lo: u64 = 0;
            while (lo <= n) : (lo += span) {
                const hi = @min(lo + span, n + 1); // exclusive integer bound
                Store.clearAll(&st.seg, cap);
                for (st.primes) |*pr| {
                    const p = pr.p;
                    var j = pr.next;
                    while (j < hi) : (j += p) Store.set(&st.seg, j - lo); // local index
                    pr.next = j; // carry cursor into the next segment
                }
                comp += Store.countComposites(&st.seg, cap); // fold the count in
            }
            st.composites = comp;
        }

        pub fn count(st: *const State, n: u64) u64 {
            return n - 1 - st.composites; // set bits are composites; 0,1 unset
        }
    };
}

/// Base primes up to `limit` (== isqrt(n)), each with its cursor at p*p.
/// Evens included, so this starts at 2.
fn collectBasePrimes(gpa: std.mem.Allocator, limit: u64) ![]Prime {
    const m: usize = @intCast(limit + 1);
    const is_comp = try gpa.alloc(bool, m);
    defer gpa.free(is_comp);
    @memset(is_comp, false);

    var i: usize = 2;
    while (i * i < m) : (i += 1) {
        if (!is_comp[i]) {
            var j = i * i;
            while (j < m) : (j += i) is_comp[j] = true;
        }
    }

    var cnt: usize = 0;
    for (is_comp[2..]) |c| {
        if (!c) cnt += 1;
    }

    const primes = try gpa.alloc(Prime, cnt);
    var k: usize = 0;
    var v: usize = 2;
    while (v < m) : (v += 1) {
        if (!is_comp[v]) {
            primes[k] = .{ .p = @intCast(v), .next = @as(u64, @intCast(v)) * @as(u64, @intCast(v)) };
            k += 1;
        }
    }
    return primes;
}
