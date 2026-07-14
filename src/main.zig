//! Meissel–Lehmer combinatorial π(x) — verify against known values, timed.

const std = @import("std");
const meissel = @import("meissel.zig");
const common = @import("common.zig");

pub fn main() !void {
    const gpa = std.heap.page_allocator;
    const xs = [_]u64{
        1_000_000,      100_000_000,     1_000_000_000,
        10_000_000_000, 100_000_000_000, 1_000_000_000_000,
    };
    std.debug.print("Meissel π(x):\n", .{});
    for (xs) |x| {
        const t0 = common.nowNs();
        const p = try meissel.pi(gpa, x);
        const ms = @as(f64, @floatFromInt(common.nowNs() - t0)) / 1e6;
        const chk = if (common.expectedPi(x)) |e| (if (p == e) "OK" else "FAIL!") else "?";
        std.debug.print("  π({d:>12}) = {d:>10}   {s:<6} [{d:.1} ms]\n", .{ x, p, chk, ms });
    }
}
