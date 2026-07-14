//! SegmentedSieve — the segmented traversal, generic over a backing Store and a
//! segment size. A segment is just a small bit-array (a Store sized to `slots`
//! flags, reused per block), so segmentation composes with the same stores the
//! naive sieve uses. Three orthogonal axes: coordinate map × traversal × store.
//!
//! The knob is `slots` = the store's flag count (NOT integers, NOT bytes). This
//! is the coordinate-independent, cache-critical quantity: equal `slots` means
//! equal store bytes for a given store type, so comparisons hold cache footprint
//! constant. For all-numbers, 1 flag/integer, so slots == integers per segment.
//!   slots 262144  =>  []u64 store 32 KiB (fits L1d), covers 262144 integers
//!                     []bool store 256 KiB (spills L2), covers 262144 integers
//! (odds-only, at the same `slots`, has the same store bytes but covers 2×.)
//!
//! Counting is folded in (segments are discarded as we go). The segment Store
//! only needs clearAll/set/countComposites — `get` is unused here (base primes
//! come from a separate small sieve).

const std = @import("std");
const common = @import("common.zig");

const Prime = struct { p: u64, next: u64 };

pub fn SegmentedSieve(comptime Store: type, comptime seg_bytes: u64) type {
    return struct {
        const slots: u64 = seg_bytes * Store.flags_per_byte; // flags that fit in the target bytes
        const cap: u64 = slots - 1;
        const store_bytes = Store.footprintBytes(cap); // actual footprint (~ seg_bytes)

        pub const name = std.fmt.comptimePrint("segmented[{s}], {d} KiB store, {d} slots", .{ Store.name, store_bytes / 1024, slots });

        pub const State = struct {
            seg: Store.State,
            primes: []Prime,
            composites: u64 = 0,
        };

        pub fn footprintBytes(n: u64) usize {
            const lim = common.isqrt(n);
            const est_primes: usize = if (lim < 3) 2 else @intFromFloat(1.2 * @as(f64, @floatFromInt(lim)) / @log(@as(f64, @floatFromInt(lim))));
            return store_bytes + est_primes * @sizeOf(Prime);
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
            while (lo <= n) : (lo += slots) {
                const hi = @min(lo + slots, n + 1); // exclusive integer bound
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
