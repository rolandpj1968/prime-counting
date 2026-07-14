//! Prime sieve benchmark ladder — entry point.
//!
//! Bucket vs naive across segment sizes at N=1e10. As the segment shrinks,
//! more base primes strike each segment ≤ once (become "large"), so the naive
//! all-cursor sieve wastes more skip-work — and the bucket sieve's advantage
//! grows from ~1× (no large primes) toward a real win.

const std = @import("std");
const wheel = @import("wheel.zig");
const sweep = @import("sweep.zig");

const store_bit_packed = @import("stores/bit_packed.zig");

const N: u64 = 10_000_000_000;
const REPEATS: usize = 2;
const W = wheel.Wheel(&[_]u64{ 2, 3, 5 }); // mod-30

const SEGS = [_]u64{ 1024, 2048, 4096, 8192, 16384, 32768 };

pub fn main() !void {
    const gpa = std.heap.page_allocator;
    try sweep.bucketCompare(W, store_bit_packed, &SEGS, gpa, N, REPEATS);
}
