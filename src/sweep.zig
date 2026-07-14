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
const wheel = @import("wheel.zig");
const bucket_sieve = @import("bucket_sieve.zig");

fn timeSieve(comptime S: type, gpa: std.mem.Allocator, n: u64, repeats: usize) !f64 {
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
    if (common.expectedPi(n)) |e| {
        if (pi != e) std.debug.print("  !! PI FAIL: {d} != {d}\n", .{ pi, e });
    }
    return @as(f64, @floatFromInt(n)) / (@as(f64, @floatFromInt(best)) / 1e9) / 1e6;
}

/// Naive (all-cursor) vs bucketed, across segment sizes. Bucket win grows as the
/// segment shrinks (more primes become "large" → more skip-work the bucket cuts).
pub fn bucketCompare(
    comptime W: type,
    comptime Store: type,
    comptime segs: []const u64,
    gpa: std.mem.Allocator,
    n: u64,
    repeats: usize,
) !void {
    std.debug.print("\n# bucket vs naive — wheel M={d}, store {s}, N={d}\n", .{ W.M, Store.name, n });
    std.debug.print("{s:>8}  {s:>9}  {s:>10}  {s:>9}  {s:>8}\n", .{ "seg", "seg_slots", "naive M/s", "bkt M/s", "speedup" });
    inline for (segs) |seg| {
        const naive = try timeSieve(sieve.Sieve(W, Store, seg), gpa, n, repeats);
        const bkt = try timeSieve(bucket_sieve.BucketSieve(W, Store, seg), gpa, n, repeats);
        std.debug.print("{d:>6} B  {d:>9}  {d:>10.0}  {d:>9.0}  {d:>7.2}x\n", .{ seg, seg * Store.flags_per_byte, naive, bkt, bkt / naive });
    }
}

/// Sweep the wheel (p_n = largest wheel prime), fixed store/segment/N. Shows the
/// diminishing Mertens returns and where the wheel's own bookkeeping (φ spokes,
/// per-prime delta table) starts to overwhelm the shrinking gain.
pub fn wheelSweep(
    comptime wheels: anytype,
    comptime Store: type,
    comptime seg: u64,
    gpa: std.mem.Allocator,
    n: u64,
    repeats: usize,
) !void {
    std.debug.print("\n# wheel sweep — store {s}, seg {d} KiB, N={d}\n", .{ Store.name, seg / 1024, n });
    std.debug.print("{s:>4}  {s:>7}  {s:>6}  {s:>8}  {s:>9}  {s:>10}\n", .{ "p_n", "M", "spokes", "density", "desc B", "M ints/s" });
    inline for (wheels) |wp| {
        const W = wheel.Wheel(wp);
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
        const pn: u64 = if (wp.len == 0) 1 else wp[wp.len - 1];
        const density = @as(f64, @floatFromInt(W.spokes)) / @as(f64, @floatFromInt(W.M));
        std.debug.print("{d:>4}  {d:>7}  {d:>6}  {d:>8.4}  {d:>9}  {d:>10.0}  {s}\n", .{ pn, W.M, W.spokes, density, @sizeOf(W.Prime), rate, if (ok) "OK" else "PI FAIL" });
    }
}

fn level(seg: u64) []const u8 {
    if (seg <= 32 * 1024) return "L1d";
    if (seg <= 512 * 1024) return "L2";
    if (seg <= 16 * 1024 * 1024) return "L3";
    return "DRAM";
}

/// Pure-π(N) scaling study: compute π(N) across N and tabulate the asymptotics.
/// Duck-typed over any Impl with init/sieve/count/deinit.
pub fn piScaling(comptime Impl: type, ns: []const u64, gpa: std.mem.Allocator) !void {
    std.debug.print("\n# π(N) scaling — {s}\n", .{Impl.name});
    std.debug.print("{s:>14}  {s:>15}  {s:>9}  {s:>7}  {s:>7}  {s:>9}  {s:>12}\n", .{ "N", "pi(N)", "pi/N", "gap", "lnN", "pi*lnN/N", "li(N)-pi(N)" });
    for (ns) |n| {
        var st = try Impl.init(gpa, n);
        Impl.sieve(&st, n);
        const pi = Impl.count(&st, n);
        Impl.deinit(&st, gpa);

        if (common.expectedPi(n)) |e| {
            if (pi != e) std.debug.print("  !! PI FAIL {d} != {d}\n", .{ pi, e });
        }
        const nf: f64 = @floatFromInt(n);
        const pf: f64 = @floatFromInt(pi);
        const lnN = @log(nf);
        std.debug.print("{d:>14}  {d:>15}  {d:>9.6}  {d:>7.3}  {d:>7.3}  {d:>9.5}  {d:>12.1}\n", .{
            n, pi, pf / nf, nf / pf, lnN, pf * lnN / nf, common.li(nf) - pf,
        });
    }
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
