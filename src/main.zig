//! Parallel sweep 10^14..10^21 at 4 pinned cores, with peak-RSS memory.
//! Ascending x => getrusage.maxrss after each x is that x's peak (new high-water).
//! 10^20, 10^21 need u128 x (past 2^64) and u128 pi (pi(10^21) > u64). Published
//! values (Gourdon/Deleglise, in Oliveira e Silva 2006 Table IV) for the check.
const std = @import("std");
const common = @import("common.zig");
const lmo = @import("lmo.zig");

fn maxrssMB() f64 {
    const ru = std.posix.getrusage(std.posix.rusage.SELF);
    return @as(f64, @floatFromInt(ru.maxrss)) / 1024.0; // maxrss is KB on Linux
}

pub fn main() !void {
    const gpa = std.heap.page_allocator;
    const cores = [_]u32{ 0, 2, 4, 6 }; // 4 physical cores (bang/buck)
    const K = 8;
    const cases = [_]struct { x: u128, want: u128 }{
        .{ .x = 100_000_000_000_000, .want = 3_204_941_750_802 },
        .{ .x = 1_000_000_000_000_000, .want = 29_844_570_422_669 },
        .{ .x = 10_000_000_000_000_000, .want = 279_238_341_033_925 },
        .{ .x = 100_000_000_000_000_000, .want = 2_623_557_157_654_233 },
        .{ .x = 1_000_000_000_000_000_000, .want = 24_739_954_287_740_860 },
        .{ .x = 10_000_000_000_000_000_000, .want = 234_057_667_276_344_607 },
        .{ .x = 100_000_000_000_000_000_000, .want = 2_220_819_602_560_918_840 }, // 10^20
        .{ .x = 1_000_000_000_000_000_000_000, .want = 21_127_269_486_018_731_928 }, // 10^21
    };
    std.debug.print("Parallel sweep, 4 pinned cores, k_over={d}\n", .{K});
    std.debug.print("{s:>4} {s:>22} {s:>11} {s:>10} {s:>6}\n", .{ "10^N", "pi(x)", "time s", "peak MB", "ok" });
    var N: u64 = 14;
    for (cases) |c| {
        const t0 = common.nowNs();
        const r = try lmo.piLMOPar(gpa, c.x, cores.len, K, &cores);
        const sec = common.ms(common.nowNs() - t0) / 1000.0;
        const ok = r.pi == c.want;
        std.debug.print("{d:>4} {d:>22} {d:>11.2} {d:>9.0} {s:>6}\n", .{ N, r.pi, sec, maxrssMB(), if (ok) "y" else "WRONG" });
        N += 1;
    }
}
