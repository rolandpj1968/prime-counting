//! Shared utilities used by the harness and by every sieve implementation.

const std = @import("std");
const linux = std.os.linux;

/// floor(sqrt(n)), exact for u64.
pub fn isqrt(n: u64) u64 {
    if (n < 2) return n;
    var x: u64 = @intFromFloat(@sqrt(@as(f64, @floatFromInt(n))));
    while (x * x > n) x -= 1;
    while ((x + 1) * (x + 1) <= n) x += 1;
    return x;
}

/// Monotonic nanoseconds. std.time.Timer was removed in 0.16's Io rework, so we
/// read the clock directly (VDSO fast path, no syscall in the common case).
pub fn nowNs() u64 {
    var ts: linux.timespec = undefined;
    _ = linux.clock_gettime(.MONOTONIC, &ts);
    return @as(u64, @intCast(ts.sec)) * 1_000_000_000 + @as(u64, @intCast(ts.nsec));
}

pub fn ms(ns: u64) f64 {
    return @as(f64, @floatFromInt(ns)) / 1e6;
}

pub const KnownPi = struct { n: u64, pi: u64 };

pub const known = [_]KnownPi{
    .{ .n = 10, .pi = 4 },
    .{ .n = 100, .pi = 25 },
    .{ .n = 1000, .pi = 168 },
    .{ .n = 1_000_000, .pi = 78498 },
    .{ .n = 1_000_000_000, .pi = 50_847_534 },
    .{ .n = 10_000_000_000, .pi = 455_052_511 },
    .{ .n = 100_000_000_000, .pi = 4_118_054_813 },
    .{ .n = 1_000_000_000_000, .pi = 37_607_912_018 },
};

pub fn expectedPi(n: u64) ?u64 {
    for (known) |k| if (k.n == n) return k.pi;
    return null;
}
