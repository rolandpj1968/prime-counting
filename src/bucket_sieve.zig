//! BucketSieve — segmented sieve with a bucketed PrimeSource for large primes.
//!
//! Two tiers, split by whether a prime strikes a segment more than once:
//!   - SMALL (min delta ≤ seg_slots): strikes a segment ≥ once → cursor loop,
//!     visited every segment (same as the plain Sieve).
//!   - LARGE (min delta > seg_slots): strikes a segment ≤ once → BUCKETED. Each
//!     large prime is one node in an intrusive singly-linked list, filed under
//!     the segment where it next strikes. Processing segment S drains its list
//!     (one strike each), then re-files each node under its next-strike segment.
//!
//! This only helps when √N > seg_slots (else every prime is SMALL and there's
//! nothing to bucket) — i.e. large N or small segments, the record regime. It
//! avoids iterating (and cache-streaming) every large prime on every segment.
//!
//! Bucket structure is a full per-segment head array (heads[num_segments]) — fine
//! for demo N (≤1e12 → ≤~1e6 segments); record N would need a circular window.

const std = @import("std");
const common = @import("common.zig");

pub fn BucketSieve(comptime W: type, comptime Store: type, comptime seg_bytes: u64) type {
    return struct {
        const seg_slots: u64 = seg_bytes * Store.flags_per_byte;
        const cap: u64 = seg_slots - 1;
        const store_bytes = Store.footprintBytes(cap);
        const phi = W.spokes;
        const SENT: u32 = 0xFFFF_FFFF;

        pub const name = std.fmt.comptimePrint("bucket(M={d})[{s}], {d} KiB", .{ W.M, Store.name, store_bytes / 1024 });

        const Entry = struct { idx: u64, spoke: usize, next: u32 };

        pub const State = struct {
            seg: Store.State,
            small: []W.Prime,
            large: []W.Prime,
            entries: []Entry, // one per large prime (list node; entry i ↔ large[i])
            heads: []u32, // per-segment list head, SENT = empty
            composites: u64 = 0,
            n_small: usize = 0,
            n_large: usize = 0,
        };

        fn minDelta(d: [phi]u64) u64 {
            var m: u64 = d[0];
            for (d[1..]) |x| m = @min(m, x);
            return m;
        }

        pub fn footprintBytes(n: u64) usize {
            const lim = common.isqrt(n);
            const est: usize = if (lim < 8) 8 else @intFromFloat(1.2 * @as(f64, @floatFromInt(lim)) / @log(@as(f64, @floatFromInt(lim))));
            const num_seg = W.slotCount(n) / seg_slots + 1;
            return store_bytes + est * @sizeOf(W.Prime) + @as(usize, @intCast(num_seg)) * 4;
        }

        pub fn init(gpa: std.mem.Allocator, n: u64) !State {
            var limit = common.isqrt(n);
            if (limit * limit < n) limit += 1;

            // base sieve → all wheel-coprime primes as strike descriptors
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
            var np: usize = 0;
            var v: usize = 2;
            while (v < mm) : (v += 1) {
                if (!is_comp[v] and !W.isWheelPrime(@intCast(v))) np += 1;
            }
            const all = try gpa.alloc(W.Prime, np);
            defer gpa.free(all);
            var k: usize = 0;
            v = 2;
            while (v < mm) : (v += 1) {
                if (!is_comp[v] and !W.isWheelPrime(@intCast(v))) {
                    all[k] = W.makePrime(@intCast(v));
                    k += 1;
                }
            }

            // classify small vs large
            var nl: usize = 0;
            for (all) |pr| {
                if (minDelta(pr.d) > seg_slots) nl += 1;
            }
            const ns = np - nl;
            const small = try gpa.alloc(W.Prime, ns);
            const large = try gpa.alloc(W.Prime, nl);
            var si: usize = 0;
            var li: usize = 0;
            for (all) |pr| {
                if (minDelta(pr.d) > seg_slots) {
                    large[li] = pr;
                    li += 1;
                } else {
                    small[si] = pr;
                    si += 1;
                }
            }

            const total = W.slotCount(n);
            const num_seg: usize = @intCast(total / seg_slots + 1);
            const heads = try gpa.alloc(u32, num_seg);
            const entries = try gpa.alloc(Entry, nl);
            const seg = try Store.init(gpa, cap);
            return .{ .seg = seg, .small = small, .large = large, .entries = entries, .heads = heads, .n_small = ns, .n_large = nl };
        }

        pub fn deinit(st: *State, gpa: std.mem.Allocator) void {
            Store.deinit(&st.seg, gpa);
            gpa.free(st.small);
            gpa.free(st.large);
            gpa.free(st.entries);
            gpa.free(st.heads);
        }

        pub fn sieve(st: *State, n: u64) void {
            const total = W.slotCount(n);
            // reset small cursors
            for (st.small) |*pr| {
                pr.idx = pr.start;
                pr.spoke = 0;
            }
            // reset + seed buckets
            @memset(st.heads, SENT);
            for (st.large, 0..) |*pr, i| {
                st.entries[i] = .{ .idx = pr.start, .spoke = 0, .next = SENT };
                if (pr.start < total) {
                    const s: usize = @intCast(pr.start / seg_slots);
                    st.entries[i].next = st.heads[s];
                    st.heads[s] = @intCast(i);
                }
            }

            var comp: u64 = 0;
            var seg_idx: usize = 0;
            var lo: u64 = 0;
            while (lo < total) : (lo += seg_slots) {
                const hi = @min(lo + seg_slots, total);
                Store.clearAll(&st.seg, cap);

                // SMALL: cursor loop, every segment
                for (st.small) |*pr| {
                    var idx = pr.idx;
                    var spoke = pr.spoke;
                    const d = &pr.d;
                    while (idx < hi) {
                        Store.set(&st.seg, idx - lo);
                        idx += d[spoke];
                        spoke = (spoke + 1) % phi;
                    }
                    pr.idx = idx;
                    pr.spoke = spoke;
                }

                // LARGE: drain this segment's bucket, one strike each, refile
                var i = st.heads[seg_idx];
                st.heads[seg_idx] = SENT;
                while (i != SENT) {
                    const e = &st.entries[i];
                    const ni = e.next;
                    Store.set(&st.seg, e.idx - lo);
                    const d = &st.large[i].d;
                    e.idx += d[e.spoke];
                    e.spoke = (e.spoke + 1) % phi;
                    if (e.idx < total) {
                        const ns2: usize = @intCast(e.idx / seg_slots);
                        e.next = st.heads[ns2];
                        st.heads[ns2] = i;
                    }
                    i = ni;
                }

                comp += Store.countComposites(&st.seg, cap);
                seg_idx += 1;
            }
            st.composites = comp;
        }

        pub fn count(st: *const State, n: u64) u64 {
            return W.piFromComposites(n, st.composites);
        }
    };
}
