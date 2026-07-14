//! Prime sieve benchmark ladder — entry point.
//!
//! Everything is one generic sieve, `Sieve(Wheel, Store, seg_bytes)`. This entry
//! runs the cache-hierarchy sweep: fixed N, segment size swept across the cache
//! levels, throughput cliffs = cache sizes.

const std = @import("std");
const wheel = @import("wheel.zig");
const sweep = @import("sweep.zig");

const store_bit_packed = @import("stores/bit_packed.zig");

const N: u64 = 1_000_000_000;
const REPEATS: usize = 3;

const W_all = wheel.Wheel(&[_]u64{});
const W_odds = wheel.Wheel(&[_]u64{2});

// Dense sampling around L1 (32 KiB) and L2 (512 KiB); out to whole-array (DRAM).
const SEGS = [_]u64{
    4 * 1024,         8 * 1024,          16 * 1024,         24 * 1024,
    32 * 1024,        48 * 1024,         64 * 1024,         96 * 1024,
    128 * 1024,       192 * 1024,        256 * 1024,        384 * 1024,
    512 * 1024,       768 * 1024,        1024 * 1024,       2 * 1024 * 1024,
    4 * 1024 * 1024,  8 * 1024 * 1024,   16 * 1024 * 1024,  32 * 1024 * 1024,
    64 * 1024 * 1024, 128 * 1024 * 1024, 256 * 1024 * 1024,
};

// Powers of 2 spanning L1 (arr=32K @ 2^18) → L2 (2^22) → L3 (2^27) → DRAM.
const NS = [_]u64{
    1 << 16, 1 << 17, 1 << 18, 1 << 19, 1 << 20, 1 << 21, 1 << 22, 1 << 23,
    1 << 24, 1 << 25, 1 << 26, 1 << 27, 1 << 28, 1 << 29, 1 << 30,
};

pub fn main() !void {
    const gpa = std.heap.page_allocator;
    // Two curves: whole-array (algorithm × cache) vs segmented (algorithm only).
    try sweep.nSweep(W_all, store_bit_packed, &NS, true, gpa, REPEATS);
    try sweep.nSweep(W_all, store_bit_packed, &NS, false, gpa, REPEATS);
    _ = SEGS;
    _ = W_odds;
    _ = N;
}
