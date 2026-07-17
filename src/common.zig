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

/// Logarithmic integral li(x) = ∫₀ˣ dt/ln t = Ei(ln x), via Ei's convergent
/// series Ei(z) = γ + ln z + Σ_{k≥1} z^k/(k·k!). Accurate to f64 for our range.
pub fn li(x: f64) f64 {
    const gamma = 0.5772156649015329;
    const z = @log(x);
    var sum = gamma + @log(z);
    var t: f64 = 1.0; // running z^k / k!
    var k: f64 = 1.0;
    while (k < 400) : (k += 1) {
        t *= z / k;
        const add = t / k;
        sum += add;
        if (k > z and add < 1e-15 * sum) break;
    }
    return sum;
}

pub const KnownPi = struct { n: u64, pi: u64 };

pub const known = [_]KnownPi{
    .{ .n = 10, .pi = 4 },
    .{ .n = 100, .pi = 25 },
    .{ .n = 1000, .pi = 168 },
    .{ .n = 10_000, .pi = 1229 },
    .{ .n = 100_000, .pi = 9592 },
    .{ .n = 1_000_000, .pi = 78498 },
    .{ .n = 10_000_000, .pi = 664579 },
    .{ .n = 100_000_000, .pi = 5761455 },
    .{ .n = 1_000_000_000, .pi = 50_847_534 },
    .{ .n = 10_000_000_000, .pi = 455_052_511 },
    .{ .n = 100_000_000_000, .pi = 4_118_054_813 },
    .{ .n = 1_000_000_000_000, .pi = 37_607_912_018 },
    .{ .n = 10_000_000_000_000, .pi = 346_065_536_839 },
    .{ .n = 100_000_000_000_000, .pi = 3_204_941_750_802 },
    .{ .n = 1_000_000_000_000_000, .pi = 29_844_570_422_669 },
    .{ .n = 10_000_000_000_000_000, .pi = 279_238_341_033_925 },
    .{ .n = 100_000_000_000_000_000, .pi = 2_623_557_157_654_233 },
    .{ .n = 1_000_000_000_000_000_000, .pi = 24_739_954_287_740_860 },
    .{ .n = 10_000_000_000_000_000_000, .pi = 234_057_667_276_344_607 },
};

pub fn expectedPi(n: u64) ?u64 {
    for (known) |k| if (k.n == n) return k.pi;
    return null;
}
