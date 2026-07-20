const std = @import("std");
const common = @import("common.zig");
const g = @import("gourdon.zig");

pub fn main() !void {
    const gpa = std.heap.page_allocator;
    const x: u128 = 100_000_000_000_000_000_000; // 10^20, past 2^64
    const want: i128 = 2_220_819_602_560_918_840; // Gourdon 2001 / Oliveira Table IV
    const pins = [_]u32{ 0, 2, 4, 6, 8, 10 }; // one thread per physical core (SMT siblings left free)
    const t0 = common.nowNs();
    const r = try g.piGourdonPar(gpa, x, pins.len, &pins);
    const secs = @as(f64, @floatFromInt(common.nowNs() - t0)) / 1e9;
    const ru = std.posix.getrusage(std.posix.rusage.SELF);
    std.debug.print("pi(10^20) = {d}\n  want    = {d}\n  {s}   {d:.1} s ({d:.2} h)   {d} threads   peakRSS = {d} MB\n", .{
        r.pi, want, if (r.pi == want) "MATCH" else "MISMATCH",
        secs,             secs / 3600.0,
        pins.len,         @divTrunc(ru.maxrss, 1024),
    });
}
