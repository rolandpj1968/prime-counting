//! Meissel–Lehmer π(x) scaling: x = 10^k, timed, with the empirical exponent
//! (log₁₀ of the time ratio per ×10 in x). Sub-linear ⇒ exponent well below 1.

const std = @import("std");
const meissel = @import("meissel.zig");
const common = @import("common.zig");

pub fn main() !void {
    const gpa = std.heap.page_allocator;
    const xs = [_]struct { x: u64, k: u32 }{
        .{ .x = 1_000_000, .k = 6 },
        .{ .x = 10_000_000, .k = 7 },
        .{ .x = 100_000_000, .k = 8 },
        .{ .x = 1_000_000_000, .k = 9 },
        .{ .x = 10_000_000_000, .k = 10 },
        .{ .x = 100_000_000_000, .k = 11 },
        .{ .x = 1_000_000_000_000, .k = 12 },
        .{ .x = 10_000_000_000_000, .k = 13 },
    };
    std.debug.print("Meissel π(x) scaling\n", .{});
    std.debug.print("{s:>5}  {s:>18}  {s:>10}  {s:>8}\n", .{ "x", "pi(x)", "time", "exponent" });
    var prev_ns: f64 = 0;
    for (xs) |it| {
        const t0 = common.nowNs();
        const p = try meissel.pi(gpa, it.x);
        const dt: f64 = @floatFromInt(common.nowNs() - t0);
        const ok = if (common.expectedPi(it.x)) |e| (p == e) else true;
        const expo = if (prev_ns > 0) @log(dt / prev_ns) / @log(10.0) else 0.0;
        const secs = dt / 1e9;
        std.debug.print("10^{d:<2}  {d:>18}  {d:>8.3}s  {d:>8.3}  {s}\n", .{ it.k, p, secs, expo, if (ok) "OK" else "FAIL" });
        prev_ns = dt;
    }
}
