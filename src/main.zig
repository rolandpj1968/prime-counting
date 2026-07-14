//! Prime sieve benchmark ladder — entry point.
//!
//! Wheel sweep: p_n = 1 (all) → 13, one generic sieve, fixed store/segment/N.

const std = @import("std");
const sweep = @import("sweep.zig");

const store_bit_packed = @import("stores/bit_packed.zig");

const N: u64 = 1_000_000_000;
const REPEATS: usize = 3;
const SEG: u64 = 32 * 1024;

// Primorial wheels: p_n = 1(all), 2, 3, 5, 7, 11, 13 → M = 1,2,6,30,210,2310,30030.
const WHEELS = .{
    &[_]u64{},
    &[_]u64{2},
    &[_]u64{ 2, 3 },
    &[_]u64{ 2, 3, 5 },
    &[_]u64{ 2, 3, 5, 7 },
    &[_]u64{ 2, 3, 5, 7, 11 },
    &[_]u64{ 2, 3, 5, 7, 11, 13 },
};

pub fn main() !void {
    const gpa = std.heap.page_allocator;
    try sweep.wheelSweep(WHEELS, store_bit_packed, SEG, gpa, N, REPEATS);
}
