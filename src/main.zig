//! Meissel π(x) with the compact π-table (LMO Stage A) — verify + time, now
//! reaching 10^14 (the fat π-per-integer table would need ~8.6 GB there).

const std = @import("std");
const meissel = @import("meissel.zig");
const common = @import("common.zig");

pub fn main() !void {
    const gpa = std.heap.page_allocator;
    const xs = [_]u64{
        1_000_000,           1_000_000_000,     10_000_000_000,
        100_000_000_000,     1_000_000_000_000, 10_000_000_000_000,
        100_000_000_000_000,
    };
    std.debug.print("Meissel π(x), compact table:\n", .{});
    for (xs) |x| {
        const t0 = common.nowNs();
        const p = try meissel.pi(gpa, x);
        const secs = @as(f64, @floatFromInt(common.nowNs() - t0)) / 1e9;
        const chk = if (common.expectedPi(x)) |e| (if (p == e) "OK" else "FAIL!") else "?";
        std.debug.print("  π({d:>15}) = {d:>14}   {s:<6} [{d:.2} s]\n", .{ x, p, chk, secs });
    }
}
