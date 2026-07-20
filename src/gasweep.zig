const std = @import("std");
const common = @import("common.zig");
const g = @import("gourdon.zig");

/// α sweep: y = α·x^(1/3). chooseY hardcodes α = 4 and has never been swept, yet α
/// sets the whole fold/leaf balance — fold work ∝ z = x/y, leaf work ∝ y. Raising α
/// shrinks z, which cuts kills AND nseg AND π(√z) at once, so the optimum is not
/// obvious a priori. Constraint: z < y² (the Legendre read-off for B) ⟺ α > 1.
///
/// Reports time per α at several x, so the optimum can be seen to move with x — if it
/// does, chooseY needs an x-dependent α, which is what primecount does. Every result
/// is checked against the α=4 run at the same x: α must not change π(x).
fn icbrt(x: u128) u64 {
    var r: u64 = @intFromFloat(std.math.cbrt(@as(f64, @floatFromInt(x))));
    while (r > 1 and @as(u128, r) * r * r > x) r -= 1;
    while (@as(u128, r + 1) * (r + 1) * (r + 1) <= x) r += 1;
    return r;
}

/// Same u64/u128 dispatch as piGourdon, so the sweep is not biased by u128 division
/// on x that the real entry point would have run on the u64 path.
fn run(gpa: std.mem.Allocator, x: u128, y: u64, nthreads: usize, pins: []const u32) !i128 {
    const r = if (x <= std.math.maxInt(u64))
        try g.piGourdonV(u64, gpa, @intCast(x), y, false, nthreads, pins)
    else
        try g.piGourdonV(u128, gpa, x, y, false, nthreads, pins);
    return r.pi;
}

pub fn main() !void {
    const gpa = std.heap.page_allocator;
    const pins = [_]u32{ 0, 2, 4, 6, 8, 10 };
    const NLO: usize = 13;
    const NHI: usize = 18;
    const alphas = [_]f64{ 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 10, 12 };

    const nlo: usize = NLO;
    const nhi: usize = NHI;

    std.debug.print("α sweep (y = α·x^(1/3)), {d} threads. α=4 is chooseY's hardcoded value.\n\n", .{pins.len});
    std.debug.print("{s:>5}", .{"α"});
    for (nlo..nhi + 1) |n| std.debug.print("   10^{d:<8}", .{n});
    std.debug.print("\n", .{});

    var ref: [40]i128 = @splat(0);
    var best: [40]f64 = @splat(1e30);
    var best_a: [40]f64 = @splat(0);

    for (alphas) |al| {
        std.debug.print("{d:>5.1}", .{al});
        for (nlo..nhi + 1) |n| {
            var x: u128 = 1;
            for (0..n) |_| x *= 10;
            const cr = icbrt(x);
            var y: u64 = @intFromFloat(al * @as(f64, @floatFromInt(cr)));
            if (y <= cr) y = cr + 1; // z < y² needs α > 1
            const t0 = common.nowNs();
            const pi = try run(gpa, x, y, pins.len, &pins);
            const secs = @as(f64, @floatFromInt(common.nowNs() - t0)) / 1e9;
            if (ref[n] == 0) ref[n] = pi;
            if (secs < best[n]) {
                best[n] = secs;
                best_a[n] = al;
            }
            std.debug.print("  {d:>8.2}s{s}", .{ secs, if (pi == ref[n]) " " else "!" });
        }
        std.debug.print("\n", .{});
    }
    std.debug.print("\nbest: ", .{});
    for (nlo..nhi + 1) |n| std.debug.print("10^{d}: α={d:.1} ({d:.2}s)   ", .{ n, best_a[n], best[n] });
    std.debug.print("\n('!' = π(x) differed from the α=4 run — a correctness bug, not a tuning result)\n", .{});
}
