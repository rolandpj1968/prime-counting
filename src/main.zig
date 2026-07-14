//! Prime sieve benchmark ladder — entry point.
//!
//! Three orthogonal axes, all composable:
//!   - traversal: NaiveSieve (whole-array) / SegmentedSieve (cache blocks)
//!   - backing store: flat []bool / hand-rolled []u64 / std.DynamicBitSet
//!   - (coming) coordinate map: all / odds / wheel
//! Both traversals are generic over the same stores.

const std = @import("std");
const bench = @import("bench.zig");

const naive = @import("naive_sieve.zig");
const segmented = @import("segmented_sieve.zig");
const segmented_odds = @import("segmented_odds.zig");

const store_flat_bool = @import("stores/flat_bool.zig");
const store_bit_packed = @import("stores/bit_packed.zig");
const store_bit_set_std = @import("stores/bit_set_std.zig");

const N: u64 = 1_000_000_000;
const REPEATS: usize = 3;
const SEG_BYTES: u64 = 32 * 1024; // store byte size per segment — the cache-critical knob (L1d)

pub fn main() !void {
    const gpa = std.heap.page_allocator;

    // Whole-array traversal, one logic, three stores:
    try bench.run(naive.NaiveSieve(store_flat_bool), gpa, N, REPEATS);
    try bench.run(naive.NaiveSieve(store_bit_set_std), gpa, N, REPEATS);
    try bench.run(naive.NaiveSieve(store_bit_packed), gpa, N, REPEATS);

    // Segmented, all numbers — knob is store BYTES, so every store has the SAME
    // 32 KiB cache footprint (bit stores cover 8× more integers/segment).
    try bench.run(segmented.SegmentedSieve(store_bit_packed, SEG_BYTES), gpa, N, REPEATS);
    try bench.run(segmented.SegmentedSieve(store_bit_set_std, SEG_BYTES), gpa, N, REPEATS);
    try bench.run(segmented.SegmentedSieve(store_flat_bool, SEG_BYTES), gpa, N, REPEATS);

    // Segmented, odds only (wheel-2) — same 32 KiB store footprint.
    try bench.run(segmented_odds.SegmentedOdds(store_bit_packed, SEG_BYTES), gpa, N, REPEATS);
    try bench.run(segmented_odds.SegmentedOdds(store_flat_bool, SEG_BYTES), gpa, N, REPEATS);
}
