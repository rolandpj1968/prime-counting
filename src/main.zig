//! Atkin vs Eratosthenes — full-array bit sieve, same N, count π(N).

const std = @import("std");
const atkin = @import("atkin.zig");
const common = @import("common.zig");

pub fn main() !void {
    const gpa = std.heap.page_allocator;
    const ns = [_]u64{ 1_000_000, 100_000_000, 1_000_000_000 };
    std.debug.print("Atkin vs Eratosthenes (full-array bit sieve)\n", .{});
    std.debug.print("{s:>13}  {s:>9}  {s:>9}  {s:>9}  {s:>9}  {s}\n", .{ "N", "atkin ms", "atk M/s", "era ms", "era M/s", "check" });
    for (ns) |N| {
        var t = common.nowNs();
        const pa = try atkin.countAtkin(gpa, N);
        const a_ms = @as(f64, @floatFromInt(common.nowNs() - t)) / 1e6;
        t = common.nowNs();
        const pe = try atkin.countEratosthenes(gpa, N);
        const e_ms = @as(f64, @floatFromInt(common.nowNs() - t)) / 1e6;
        const nf: f64 = @floatFromInt(N);
        const ok = pa == pe and (if (common.expectedPi(N)) |e| pa == e else true);
        std.debug.print("{d:>13}  {d:>9.1}  {d:>9.0}  {d:>9.1}  {d:>9.0}  {s} (π={d})\n", .{ N, a_ms, nf / (a_ms / 1e3) / 1e6, e_ms, nf / (e_ms / 1e3) / 1e6, if (ok) "OK" else "FAIL", pa });
    }
}
