//! SegmentedOdds — segmented traversal on the ODD numbers only (wheel-2 on the
//! coordinate-map axis), generic over a backing Store.
//!
//! Coordinate map: bit index t <-> odd value v = 2t + 1  (t=0 is value 1).
//! Only odd numbers get a slot, so the footprint and the strike count both
//! roughly halve, and prime 2 — the densest striker — is gone entirely.
//!
//! Striking prime p (odd): its odd multiples are p*p, p*p+2p, ... (all odd,
//! step 2p), which in index space is start (p*p-1)/2, step p.
//!
//! Counting: set bits are the odd composites. π(N) = (N+1)/2 - odd_composites.
//! The value-1 slot (index 0, never struck, not prime) stands in for the prime
//! 2 we dropped, so the formula needs no special case.
//!
//! Kept concrete (not a Coord abstraction yet): once the mod-30 wheel exists we
//! extract the shared coordinate map from {all, odds, wheel} — three points.

const std = @import("std");
const common = @import("common.zig");

const Prime = struct { p: u64, next: u64 };

pub fn SegmentedOdds(comptime Store: type, comptime seg_bytes: u64) type {
    return struct {
        const slots: u64 = seg_bytes * Store.flags_per_byte; // odd-flags that fit in the target bytes
        const cap: u64 = slots - 1;
        const store_bytes = Store.footprintBytes(cap); // actual footprint (~ seg_bytes)

        pub const name = std.fmt.comptimePrint("seg-odds[{s}], {d} KiB store, {d} slots", .{ Store.name, store_bytes / 1024, slots });

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
            const primes = try collectOddBasePrimes(gpa, common.isqrt(n));
            const seg = try Store.init(gpa, cap);
            return .{ .seg = seg, .primes = primes };
        }

        pub fn deinit(st: *State, gpa: std.mem.Allocator) void {
            Store.deinit(&st.seg, gpa);
            gpa.free(st.primes);
        }

        pub fn sieve(st: *State, n: u64) void {
            const total: u64 = (n + 1) / 2; // number of odd slots (values 1,3,5,...)
            for (st.primes) |*pr| pr.next = (pr.p * pr.p - 1) / 2; // reset cursor to odd-index of p*p
            var comp: u64 = 0;
            var tlo: u64 = 0;
            while (tlo < total) : (tlo += slots) {
                const thi = @min(tlo + slots, total); // exclusive odd-index bound
                Store.clearAll(&st.seg, cap);
                for (st.primes) |*pr| {
                    const p = pr.p;
                    var t = pr.next;
                    while (t < thi) : (t += p) Store.set(&st.seg, t - tlo);
                    pr.next = t;
                }
                comp += Store.countComposites(&st.seg, cap);
            }
            st.composites = comp;
        }

        pub fn count(st: *const State, n: u64) u64 {
            return (n + 1) / 2 - st.composites; // T - odd_composites (value-1 slot covers prime 2)
        }
    };
}

/// Odd primes in [3, limit], each with its cursor at the odd-index of p*p.
fn collectOddBasePrimes(gpa: std.mem.Allocator, limit: u64) ![]Prime {
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
    var v: usize = 3;
    while (v < m) : (v += 2) {
        if (!is_comp[v]) cnt += 1;
    }

    const primes = try gpa.alloc(Prime, cnt);
    var k: usize = 0;
    v = 3;
    while (v < m) : (v += 2) {
        if (!is_comp[v]) {
            const p: u64 = @intCast(v);
            primes[k] = .{ .p = p, .next = (p * p - 1) / 2 };
            k += 1;
        }
    }
    return primes;
}
