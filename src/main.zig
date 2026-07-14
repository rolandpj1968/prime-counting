//! Prime sieve benchmark ladder — entry point.
//!
//! Everything is now ONE generic sieve, `Sieve(Wheel, Store, seg_bytes)`:
//!   - Wheel(wheel_primes): coordinate map (all / odds / mod-6 / mod-30 / …)
//!   - Store: backing bits (flat []bool / hand-rolled []u64 / std.DynamicBitSet)
//!   - seg_bytes: segment size; ≥ N ⇒ whole-array
//! Every earlier rung is a point in this space.

const std = @import("std");
const bench = @import("bench.zig");

const wheel = @import("wheel.zig");
const sieve = @import("sieve.zig");

const store_flat_bool = @import("stores/flat_bool.zig");
const store_bit_packed = @import("stores/bit_packed.zig");
const store_bit_set_std = @import("stores/bit_set_std.zig");

const N: u64 = 1_000_000_000;
const REPEATS: usize = 3;
const SEG: u64 = 32 * 1024; // store bytes per segment (L1d)

const W_all = wheel.Wheel(&[_]u64{});
const W_odds = wheel.Wheel(&[_]u64{2});
const W_mod6 = wheel.Wheel(&[_]u64{ 2, 3 });
const W_mod30 = wheel.Wheel(&[_]u64{ 2, 3, 5 });

pub fn main() !void {
    const gpa = std.heap.page_allocator;

    // All four wheels from one sieve — first, verify π is correct for each.
    try bench.run(sieve.Sieve(W_all, store_bit_packed, SEG), gpa, N, REPEATS);
    try bench.run(sieve.Sieve(W_odds, store_bit_packed, SEG), gpa, N, REPEATS);
    try bench.run(sieve.Sieve(W_mod6, store_bit_packed, SEG), gpa, N, REPEATS);
    try bench.run(sieve.Sieve(W_mod30, store_bit_packed, SEG), gpa, N, REPEATS);

    // mod-30 across two stores (the byte-vs-bit strike comparison, wheeled).
    try bench.run(sieve.Sieve(W_mod30, store_flat_bool, SEG), gpa, N, REPEATS);
}
