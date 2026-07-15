//! LMO Stage B foundation: φ(x,a) = S1 (ordinary, done) + S2 (special, TODO).
//! Show S1 and the oracle φ so S2 has an exact target to reproduce.

const std = @import("std");
const lmo = @import("lmo.zig");
const meissel = @import("meissel.zig");

pub fn main() !void {
    const gpa = std.heap.page_allocator;
    const xs = [_]u64{ 1_000_000, 100_000_000, 10_000_000_000, 1_000_000_000_000 };
    std.debug.print("LMO φ(x,a) = S1 (ordinary) + S2 (special leaves, TODO)\n", .{});
    std.debug.print("{s:>14} {s:>4} {s:>20} {s:>16} {s:>20}\n", .{ "x", "a", "S1 = Σμ(n)⌊x/n⌋", "φ(x,a) oracle", "target S2 = φ−S1" });
    for (xs) |x| {
        const f = try lmo.ordinaryS1(gpa, x);
        const phi: i128 = @intCast(try meissel.phiOfX(gpa, x));
        std.debug.print("{d:>14} {d:>4} {d:>20} {d:>16} {d:>20}\n", .{ x, f.a, f.s1, phi, phi - f.s1 });
    }
}
