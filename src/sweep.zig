//! Cache-hierarchy sweep. Fix N, sweep the segment size (= the store's working
//! set), and watch throughput: it rises as per-segment overhead amortizes, then
//! falls off a cliff each time the segment outgrows a cache level. The cliffs
//! ARE the cache sizes, read straight off the curve.
//!
//! Use a strike-heavy wheel (all / odds) + []u64 so the working set is hammered
//! by memory ops → sharp cliffs. mod-30 is too efficient (barely touches memory)
//! and washes the steps out.

const std = @import("std");
const common = @import("common.zig");
const sieve = @import("sieve.zig");

fn level(seg: u64) []const u8 {
    if (seg <= 32 * 1024) return "L1d";
    if (seg <= 512 * 1024) return "L2";
    if (seg <= 16 * 1024 * 1024) return "L3";
    return "DRAM";
}

/// Sweep N (powers of 2). `whole=true`: one segment = the whole array, so the
/// working set grows with N and crosses cache levels (memory mountain in N).
/// `whole=false`: fixed 32 KiB segment → working set pinned in L1, so no cache
/// signal — only the Mertens strike-density decay. The two curves' ratio is the
/// pure cache penalty (algorithm×cache ÷ algorithm).
pub fn nSweep(
    comptime W: type,
    comptime Store: type,
    comptime ns: []const u64,
    comptime whole: bool,
    gpa: std.mem.Allocator,
    repeats: usize,
) !void {
    const mode = if (whole) "WHOLE-ARRAY (working set grows with N)" else "SEGMENTED 32 KiB (working set pinned in L1)";
    std.debug.print("\n# N sweep — {s}, wheel M={d}, store {s}\n", .{ mode, W.M, Store.name });
    std.debug.print("{s:>5}  {s:>10}  {s:>5}  {s:>10}\n", .{ "N", "workset", "level", "M ints/s" });
    inline for (ns) |n| {
        const seg: u64 = if (whole) n / Store.flags_per_byte + 64 else 32 * 1024;
        const S = sieve.Sieve(W, Store, seg);
        var st = try S.init(gpa, n);
        defer S.deinit(&st, gpa);
        var best: u64 = std.math.maxInt(u64);
        var r: usize = 0;
        while (r < repeats) : (r += 1) {
            const t0 = common.nowNs();
            S.sieve(&st, n);
            best = @min(best, common.nowNs() - t0);
        }
        _ = S.count(&st, n);
        const rate = @as(f64, @floatFromInt(n)) / (@as(f64, @floatFromInt(best)) / 1e9) / 1e6;
        const ws: u64 = if (whole) n / Store.flags_per_byte else 32 * 1024;
        const k = comptime @ctz(n);
        std.debug.print("  2^{d:<2}  {d:>8} KiB  {s:>5}  {d:>10.0}\n", .{ k, ws / 1024, level(ws), rate });
    }
}

pub fn segSweep(
    comptime W: type,
    comptime Store: type,
    comptime segs: []const u64,
    gpa: std.mem.Allocator,
    n: u64,
    repeats: usize,
) !void {
    std.debug.print("\n# segment sweep — wheel M={d}, store {s}, N={d}\n", .{ W.M, Store.name, n });
    std.debug.print("{s:>10}  {s:>5}  {s:>10}  {s:>9}\n", .{ "seg", "level", "M ints/s", "best ms" });
    inline for (segs) |seg| {
        const S = sieve.Sieve(W, Store, seg);
        var st = try S.init(gpa, n);
        defer S.deinit(&st, gpa);
        var best: u64 = std.math.maxInt(u64);
        var r: usize = 0;
        while (r < repeats) : (r += 1) {
            const t0 = common.nowNs();
            S.sieve(&st, n);
            best = @min(best, common.nowNs() - t0);
        }
        const pi = S.count(&st, n);
        const ok = if (common.expectedPi(n)) |e| pi == e else true;
        const rate = @as(f64, @floatFromInt(n)) / (@as(f64, @floatFromInt(best)) / 1e9) / 1e6;
        const kib = seg / 1024;
        std.debug.print("{d:>7} KiB  {s:>5}  {d:>10.0}  {d:>9.1}  {s}\n", .{ kib, level(seg), rate, common.ms(best), if (ok) "" else "PI FAIL" });
    }
}
