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

const store_flat_bool = @import("stores/flat_bool.zig");
const store_bit_packed = @import("stores/bit_packed.zig");
const store_bit_set_std = @import("stores/bit_set_std.zig");

const N: u64 = 1_000_000_000;
const REPEATS: usize = 3;
const SPAN: u64 = 262144; // integers per segment (2^18); []u64 => 32 KiB (L1d)

pub fn main() !void {
    const gpa = std.heap.page_allocator;

    // Whole-array traversal, one logic, three stores:
    try bench.run(naive.NaiveSieve(store_flat_bool), gpa, N, REPEATS);
    try bench.run(naive.NaiveSieve(store_bit_set_std), gpa, N, REPEATS);
    try bench.run(naive.NaiveSieve(store_bit_packed), gpa, N, REPEATS);

    // Segmented traversal, SAME span, three stores — the L1-residency experiment.
    // []u64 / bitset segments are 32 KiB (L1d); []bool segment is 256 KiB (L2).
    try bench.run(segmented.SegmentedSieve(store_bit_packed, SPAN), gpa, N, REPEATS);
    try bench.run(segmented.SegmentedSieve(store_bit_set_std, SPAN), gpa, N, REPEATS);
    try bench.run(segmented.SegmentedSieve(store_flat_bool, SPAN), gpa, N, REPEATS);
}
