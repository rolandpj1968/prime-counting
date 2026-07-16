//! How many LMO special leaves? Measure vs the x^(2/3)/ln x asymptotic.

const std = @import("std");
const lmo = @import("lmo.zig");

pub fn main() !void {
    const gpa = std.heap.page_allocator;
    const xs = [_]u64{
        1_000_000_000,         1_000_000_000_000,
        1_000_000_000_000_000, 1_000_000_000_000_000_000,
    };
    std.debug.print("LMO special-leaf count\n", .{});
    std.debug.print("{s:>22} {s:>9} {s:>18} {s:>16} {s:>7}\n", .{ "x", "y=x^1/3", "special leaves", "x^2/3 / ln x", "ratio" });
    for (xs) |x| {
        const r = try lmo.countSpecialLeaves(gpa, x);
        const x23 = std.math.pow(f64, @floatFromInt(x), 2.0 / 3.0);
        const est = x23 / @log(@as(f64, @floatFromInt(x)));
        const cf: f64 = @floatFromInt(r.count);
        std.debug.print("{d:>22} {d:>9} {d:>18} {d:>16.0} {d:>7.3}\n", .{ x, r.y, r.count, est, cf / est });
    }
}
