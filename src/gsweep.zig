const std = @import("std");
const common = @import("common.zig");
const g = @import("gourdon.zig");

/// π(10^n) ladder, n = 13..20, parallel. Known values: Gourdon 2001 / Oliveira e
/// Silva Table IV. Prints time, growth ratio and peak RSS, so the memory scaling
/// law (oracle 1.125·√x/30 + 9·y bytes) can be checked against reality and the
/// time exponent read off the ratio column (theory: ~10^(2/3) ≈ 4.6 per step).
pub fn main() !void {
    const gpa = std.heap.page_allocator;
    const pins = [_]u32{ 0, 2, 4, 6, 8, 10 }; // one thread per physical core

    const want = [_]i128{
        346065536839, // 10^13
        3204941750802, // 10^14
        29844570422669, // 10^15
        279238341033925, // 10^16
        2623557157654233, // 10^17
        24739954287740860, // 10^18
        234057667276344607, // 10^19
        2220819602560918840, // 10^20
    };

    std.debug.print("{s:>3} {s:>22} {s:>7} {s:>10} {s:>8} {s:>8}\n", .{ "n", "pi(10^n)", "ok", "secs", "ratio", "RSS_MB" });
    var prev: f64 = 0;
    var x: u128 = 10_000_000_000_000; // 10^13
    for (want) |w| {
        const t0 = common.nowNs();
        const r = try g.piGourdonPar(gpa, x, pins.len, &pins);
        const secs = @as(f64, @floatFromInt(common.nowNs() - t0)) / 1e9;
        const ru = std.posix.getrusage(std.posix.rusage.SELF);
        std.debug.print("{d:>3} {d:>22} {s:>7} {d:>10.2} {d:>8.2} {d:>8}\n", .{
            std.math.log10_int(@as(u128, x)), r.pi,
            if (r.pi == w) "y" else "**NO**",  secs,
            if (prev > 0) secs / prev else 0.0, @divTrunc(ru.maxrss, 1024),
        });
        prev = secs;
        x *= 10;
    }
}
