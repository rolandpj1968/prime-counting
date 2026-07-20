const std = @import("std");
const common = @import("common.zig");
const g = @import("gourdon.zig");

/// π(10^21) and π(10^22), parallel. Long-running (hours) — prints each result as it
/// lands rather than at the end.
///
/// The reference values are Oliveira e Silva / Walisch's published π(10^n). A
/// MISMATCH here means EITHER a bug in this code OR a mistyped reference constant,
/// so the computed value is printed prominently either way — check it against the
/// literature before believing a failure.
pub fn main() !void {
    const gpa = std.heap.page_allocator;
    const pins = [_]u32{ 0, 2, 4, 6, 8, 10 }; // one thread per physical core

    const cases = [_]struct { n: usize, x: u128, want: i128 }{
        .{ .n = 21, .x = 1_000_000_000_000_000_000_000, .want = 21_127_269_486_018_731_928 },
        .{ .n = 22, .x = 10_000_000_000_000_000_000_000, .want = 201_467_286_689_315_906_290 },
    };

    for (cases) |c| {
        const t0 = common.nowNs();
        const r = try g.piGourdonPar(gpa, c.x, pins.len, &pins);
        const secs = @as(f64, @floatFromInt(common.nowNs() - t0)) / 1e9;
        const ru = std.posix.getrusage(std.posix.rusage.SELF);
        std.debug.print(
            \\
            \\pi(10^{d}) = {d}
            \\  ref      = {d}   [{s}]
            \\  {d:.1} s ({d:.2} h)   {d} threads   peakRSS = {d} MB
            \\
        , .{
            c.n,   r.pi,
            c.want, if (r.pi == c.want) "MATCH" else "DIFFERS - verify the reference before assuming a bug",
            secs,  secs / 3600.0,
            pins.len, @divTrunc(ru.maxrss, 1024),
        });
    }
}
