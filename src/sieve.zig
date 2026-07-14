//! Sieve — the unified segmented driver. ONE traversal, generic over a Wheel
//! (coordinate map) and a Store (backing bits), with segment size (store bytes)
//! as the knob. A segment ≥ N degenerates to a whole-array sieve, so this
//! subsumes every rung: naive/segmented × all/odds/mod-6/mod-30 × any store.
//!
//! The range-sieve core (clear → strike cursors → count) is invariant; the Wheel
//! supplies per-prime (start, φ deltas) so the strike loop is coordinate-blind.
//! Sieving primes are computed ONCE (a base sieve to ceil(√N)) and reused across
//! all segments; the completeness limit L is carried and asserted (L² ≥ top) —
//! the sanity check we designed (default-on; a partial/P_n source would skip it).

const std = @import("std");
const common = @import("common.zig");

pub fn Sieve(comptime W: type, comptime Store: type, comptime seg_bytes: u64) type {
    return struct {
        const slots: u64 = seg_bytes * Store.flags_per_byte;
        const cap: u64 = slots - 1;
        const store_bytes = Store.footprintBytes(cap);

        pub const name = std.fmt.comptimePrint("sieve(M={d})[{s}], {d} KiB", .{ W.M, Store.name, store_bytes / 1024 });

        pub const State = struct {
            seg: Store.State,
            primes: []W.Prime,
            limit: u64, // completeness: contains ALL primes ≤ limit
            composites: u64 = 0,
        };

        pub fn footprintBytes(n: u64) usize {
            const lim = common.isqrt(n);
            const est: usize = if (lim < 8) 8 else @intFromFloat(1.2 * @as(f64, @floatFromInt(lim)) / @log(@as(f64, @floatFromInt(lim))));
            return store_bytes + est * @sizeOf(W.Prime);
        }

        pub fn init(gpa: std.mem.Allocator, n: u64) !State {
            var limit = common.isqrt(n);
            if (limit * limit < n) limit += 1; // ceil(√n): guarantee limit² ≥ n
            const primes = try collectPrimes(gpa, limit);
            const seg = try Store.init(gpa, cap);
            return .{ .seg = seg, .primes = primes, .limit = limit };
        }

        pub fn deinit(st: *State, gpa: std.mem.Allocator) void {
            Store.deinit(&st.seg, gpa);
            gpa.free(st.primes);
        }

        /// Base sieve: all primes ≤ limit that are coprime to the wheel (i.e. not
        /// wheel primes), each turned into a strike descriptor. Computed once.
        fn collectPrimes(gpa: std.mem.Allocator, limit: u64) ![]W.Prime {
            const mm: usize = @intCast(limit + 1);
            const is_comp = try gpa.alloc(bool, mm);
            defer gpa.free(is_comp);
            @memset(is_comp, false);
            var i: usize = 2;
            while (i * i < mm) : (i += 1) {
                if (!is_comp[i]) {
                    var j = i * i;
                    while (j < mm) : (j += i) is_comp[j] = true;
                }
            }
            var cnt: usize = 0;
            var v: usize = 2;
            while (v < mm) : (v += 1) {
                if (!is_comp[v] and !W.isWheelPrime(@intCast(v))) cnt += 1;
            }
            const primes = try gpa.alloc(W.Prime, cnt);
            var k: usize = 0;
            v = 2;
            while (v < mm) : (v += 1) {
                if (!is_comp[v] and !W.isWheelPrime(@intCast(v))) {
                    primes[k] = W.makePrime(@intCast(v));
                    k += 1;
                }
            }
            return primes;
        }

        pub fn sieve(st: *State, n: u64) void {
            std.debug.assert(st.limit * st.limit >= n); // primes cover √(top): the sanity check
            const total = W.slotCount(n);
            for (st.primes) |*pr| {
                pr.idx = pr.start;
                pr.spoke = 0;
            }
            var comp: u64 = 0;
            var lo: u64 = 0;
            while (lo < total) : (lo += slots) {
                const hi = @min(lo + slots, total);
                Store.clearAll(&st.seg, cap);
                for (st.primes) |*pr| {
                    var idx = pr.idx;
                    var spoke = pr.spoke;
                    const d = &pr.d; // pointer, not a copy — the delta table can be 46 KB (mod-30030)
                    while (idx < hi) {
                        Store.set(&st.seg, idx - lo);
                        idx += d[spoke];
                        spoke = (spoke + 1) % W.spokes;
                    }
                    pr.idx = idx;
                    pr.spoke = spoke;
                }
                comp += Store.countComposites(&st.seg, cap);
            }
            st.composites = comp;
        }

        pub fn count(st: *const State, n: u64) u64 {
            return W.piFromComposites(n, st.composites);
        }
    };
}
