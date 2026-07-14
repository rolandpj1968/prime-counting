//! Prime-density validation of high ranges. No exact oracle exists up here, so
//! we validate against THEORY: a window [N, N+Δ) should hold ≈ Δ/ln N primes
//! (local density li′(N)=1/ln N), matched to Poisson noise ~1/√count. Ratio ≈
//! 1.000 validates the sieve at that range. Top point is u128, past 2^64 —
//! a tally no public tool will hand you.

const std = @import("std");
const rs = @import("rangesieve.zig");

const DELTA: u64 = 100_000_000; // window width (count ~Δ/ln N → ~0.1% Poisson noise)
const BASE_LIMIT: u64 = 5_000_000_000; // √(2.5e19): covers windows just past 2^64

/// Expected count = ∫[N,N+Δ] dt/ln t ≈ Δ/ln(N+Δ/2)  (Δ ≪ N, midpoint is exact enough).
fn expected(n: f64) f64 {
    return @as(f64, @floatFromInt(DELTA)) / @log(n + @as(f64, @floatFromInt(DELTA)) / 2.0);
}

fn row(comptime Int: type, gpa: std.mem.Allocator, tag: []const u8, n: Int, primes: []const u64) !void {
    const obs = try rs.countInRange(Int, gpa, n, n + DELTA, primes);
    const exp = expected(@floatFromInt(n));
    std.debug.print("{d:>22}  {s:>4}  {d:>11}  {d:>11.0}  {d:>8.4}\n", .{ n, tag, obs, exp, @as(f64, @floatFromInt(obs)) / exp });
}

pub fn main() !void {
    const gpa = std.heap.page_allocator;
    std.debug.print("generating base primes ≤ {d} (~1 min)...\n", .{BASE_LIMIT});
    const bp = try rs.basePrimes(gpa, BASE_LIMIT);
    defer gpa.free(bp);
    std.debug.print("  {d} base primes\n\n", .{bp.len});

    std.debug.print("prime density validation — count [N, N+{d}) vs Δ/ln N\n", .{DELTA});
    std.debug.print("{s:>22}  {s:>4}  {s:>11}  {s:>11}  {s:>8}\n", .{ "N", "type", "observed", "expected", "ratio" });

    try row(u64, gpa, "u64", 1_000_000_000_000, bp); // 1e12
    try row(u64, gpa, "u64", 1_000_000_000_000_000, bp); // 1e15
    try row(u64, gpa, "u64", 1_000_000_000_000_000_000, bp); // 1e18
    try row(u128, gpa, "u128", 20_000_000_000_000_000_000, bp); // 2e19 > 2^64
}
