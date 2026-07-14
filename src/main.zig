//! Prime sieve benchmark ladder — entry point.
//!
//! π(N) scaling study. The sieve is now just the engine; this tabulates the
//! number theory — prime density, average gap vs ln N, the PNT ratio → 1, and
//! Li(N) − π(N) staying tiny where x/ln x is ~10% off. BucketSieve as engine so
//! the same harness scales into the large-N (bucketing) regime.

const std = @import("std");
const wheel = @import("wheel.zig");
const sweep = @import("sweep.zig");
const bucket_sieve = @import("bucket_sieve.zig");

const store_bit_packed = @import("stores/bit_packed.zig");

const W = wheel.Wheel(&[_]u64{ 2, 3, 5 }); // mod-30
const Engine = bucket_sieve.BucketSieve(W, store_bit_packed, 32 * 1024);

const NS = [_]u64{
    10_000,      100_000,       1_000_000,      10_000_000,
    100_000_000, 1_000_000_000, 10_000_000_000, 100_000_000_000,
};

pub fn main() !void {
    const gpa = std.heap.page_allocator;
    try sweep.piScaling(Engine, &NS, gpa);
}
