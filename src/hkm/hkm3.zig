//! HKM rung 3 — the exact-convolution layer, refereed.
//!
//! R2's timings showed the naive mu-hat build is O(K pi(sqrt N)) = O(N/log N)
//! and overtakes the correction at ~1e15. R4 fixes that with Newton's
//! identities, which need O(log^2 N) DENSE convolutions of length K. This
//! rung builds that layer and does nothing else — no algorithm, just the
//! primitive, checked five ways.
//!
//! A correction to the plan as originally stated: there is no O(K^2)
//! convolution here for an NTT to replace. R2's segmented sum is already O(K)
//! by prefix sum, and mu-hat is built from pi(sqrt N) SPARSE two-term
//! convolutions at O(K) each — an NTT would make each of those worse. Nor
//! does a product tree help: mean kbar(p) is about K/2 (the primes near
//! sqrt N dominate the sum), so truncation bites at the very first level and
//! the tree costs pi nodes at O(K log K) instead of pi nodes at O(K). The NTT
//! earns its place only once Newton collapses pi(sqrt N) convolutions into
//! O(log^2 N) of them. So this rung is infrastructure, and its job is to be
//! correct, not yet to be fast.
//!
//! Test 4 is the one that matters: mu-hat is built for a real N by two
//! genuinely different routes — R2's sequential sparse product, and a
//! split-halves recombination through one dense NTT convolution — and the two
//! must agree cell for cell. It also reports how close max |mu-hat| runs to
//! P/2, which is the point at which a single Goldilocks prime stops being
//! enough and CRT over two moduli becomes necessary. Measured, not assumed.

const std = @import("std");
const rs = @import("rs");
const sg = @import("seg.zig");
const ntt = @import("ntt.zig");

fn nowNs() u64 {
    var ts: std.os.linux.timespec = undefined;
    _ = std.os.linux.clock_gettime(.MONOTONIC, &ts);
    return @as(u64, @intCast(ts.sec)) * 1_000_000_000 + @as(u64, @intCast(ts.nsec));
}

/// R2's route: fold in (1 - x^j) one prime at a time, truncating at cap.
fn muHatSparse(gpa: std.mem.Allocator, ks: []const usize, cap: usize) ![]i64 {
    const a = try gpa.alloc(i64, cap + 1);
    @memset(a, 0);
    a[0] = 1;
    for (ks) |j| {
        if (j > cap) continue;
        var t: usize = cap + 1;
        while (t > j) {
            t -= 1;
            a[t] -= a[t - j];
        }
    }
    return a;
}

pub fn main(init: std.process.Init) !void {
    var args = try init.minimal.args.iterateAllocator(init.gpa);
    _ = args.skip();
    const gpa = std.heap.page_allocator;
    var prng = std.Random.DefaultPrng.init(0x5eed_1234_abcd_0001);
    const rnd = prng.random();
    var fails: usize = 0;

    // ---- 1. fast reduction vs the u128-remainder reference
    {
        const edges = [_]u64{ 0, 1, 2, ntt.P - 1, ntt.P - 2, ntt.P / 2, 0xFFFF_FFFF, 0x1_0000_0000, 0xFFFF_FFFF_FFFF_FFFF % ntt.P };
        var bad: usize = 0;
        for (edges) |a| for (edges) |b| {
            if (ntt.mul(a, b) != ntt.mulRef(a, b)) bad += 1;
        };
        var i: usize = 0;
        while (i < 20_000_000) : (i += 1) {
            const a = rnd.uintLessThan(u64, ntt.P);
            const b = rnd.uintLessThan(u64, ntt.P);
            if (ntt.mul(a, b) != ntt.mulRef(a, b)) bad += 1;
        }
        std.debug.print("1. mul vs mulRef ({d} edge pairs + 20M random): {d} mismatches  {s}\n", .{ edges.len * edges.len, bad, if (bad == 0) "OK" else "FAIL" });
        fails += @intFromBool(bad != 0);
    }

    // ---- 2. transform round-trip
    {
        var bad: usize = 0;
        var m: usize = 3;
        while (m <= 20) : (m += 1) {
            const n = @as(usize, 1) << @intCast(m);
            const a = try gpa.alloc(u64, n);
            defer gpa.free(a);
            const c = try gpa.alloc(u64, n);
            defer gpa.free(c);
            for (a) |*x| x.* = rnd.uintLessThan(u64, ntt.P);
            @memcpy(c, a);
            ntt.transform(a, false);
            ntt.transform(a, true);
            if (!std.mem.eql(u64, a, c)) bad += 1;
        }
        std.debug.print("2. round-trip, lengths 2^3..2^20: {d} failures  {s}\n", .{ bad, if (bad == 0) "OK" else "FAIL" });
        fails += @intFromBool(bad != 0);
    }

    // ---- 3. NTT convolution vs schoolbook, signed data
    {
        var bad: usize = 0;
        var trial: usize = 0;
        while (trial < 200) : (trial += 1) {
            const la = 1 + rnd.uintLessThan(usize, 300);
            const lb = 1 + rnd.uintLessThan(usize, 300);
            const ai = try gpa.alloc(i64, la);
            defer gpa.free(ai);
            const bi = try gpa.alloc(i64, lb);
            defer gpa.free(bi);
            for (ai) |*x| x.* = @as(i64, @intCast(rnd.uintLessThan(u64, 2_000_001))) - 1_000_000;
            for (bi) |*x| x.* = @as(i64, @intCast(rnd.uintLessThan(u64, 2_000_001))) - 1_000_000;
            const af = try gpa.alloc(u64, la);
            defer gpa.free(af);
            const bf = try gpa.alloc(u64, lb);
            defer gpa.free(bf);
            for (ai, af) |x, *y| y.* = ntt.fromSigned(x);
            for (bi, bf) |x, *y| y.* = ntt.fromSigned(x);
            const maxdeg = la + lb - 2;
            const got = try ntt.convolve(gpa, af, bf, maxdeg);
            defer gpa.free(got);
            const want = try gpa.alloc(i64, maxdeg + 1);
            defer gpa.free(want);
            @memset(want, 0);
            for (ai, 0..) |x, ia| for (bi, 0..) |y, ib| {
                want[ia + ib] += x * y;
            };
            for (want, got) |w, g| {
                if (ntt.toSigned(g) != w) {
                    bad += 1;
                    break;
                }
            }
        }
        std.debug.print("3. convolve vs schoolbook (200 random signed trials): {d} failures  {s}\n", .{ bad, if (bad == 0) "OK" else "FAIL" });
        fails += @intFromBool(bad != 0);
    }

    // ---- 4. mu-hat two ways, on real geometry
    std.debug.print("\n4. mu-hat: sequential sparse product vs split-halves through one NTT\n", .{});
    std.debug.print("{s:>14} {s:>9} {s:>9} {s:>7} {s:>16} {s:>12}\n", .{ "N", "K", "pi(sqrtN)", "agree", "max |mu-hat|", "P/2 room" });
    for ([_]u64{ 1_000_000, 100_000_000, 10_000_000_000 }) |N| {
        const sq = std.math.sqrt(N);
        const primes = try rs.basePrimes(gpa, sq);
        defer gpa.free(primes);
        var g = try sg.build(gpa, N, @log2(@as(f64, @floatFromInt(N))) / @as(f64, @floatFromInt(sq)));
        defer g.deinit(gpa);
        const K = g.K;
        const ks = try gpa.alloc(usize, primes.len);
        defer gpa.free(ks);
        for (primes, 0..) |p, i| ks[i] = g.kbar(p);

        const seq = try muHatSparse(gpa, ks, K);
        defer gpa.free(seq);

        const half = ks.len / 2;
        const lo = try muHatSparse(gpa, ks[0..half], K);
        defer gpa.free(lo);
        const hi = try muHatSparse(gpa, ks[half..], K);
        defer gpa.free(hi);
        const lof = try gpa.alloc(u64, K + 1);
        defer gpa.free(lof);
        const hif = try gpa.alloc(u64, K + 1);
        defer gpa.free(hif);
        for (lo, lof) |x, *y| y.* = ntt.fromSigned(x);
        for (hi, hif) |x, *y| y.* = ntt.fromSigned(x);
        const got = try ntt.convolve(gpa, lof, hif, K);
        defer gpa.free(got);

        var agree = got.len == seq.len;
        if (agree) for (seq, got) |w, v| {
            if (ntt.toSigned(v) != w) {
                agree = false;
                break;
            }
        };
        var mx: i64 = 0;
        for (seq) |v| mx = @max(mx, if (v < 0) -v else v);
        std.debug.print("{d:>14} {d:>9} {d:>9} {s:>7} {d:>16} {e:>12.3}\n", .{ N, K, primes.len, if (agree) "MATCH" else "DIFFER", mx, @as(f64, @floatFromInt(ntt.P / 2)) / @as(f64, @floatFromInt(@max(1, mx))) });
        fails += @intFromBool(!agree);
    }

    // ---- 5. what one convolution costs at the sizes R4 will want
    // A convolve() call is 2 forward transforms + 1 inverse. Newton needs
    // omega(omega+1)/2 products but only ~3 omega DISTINCT transforms if the
    // operands are transformed once and reused, so the single-transform cost
    // — not the convolve cost — is what sets R4's budget. Both are timed.
    std.debug.print("\n5. dense convolution cost (one call, single-threaded)\n", .{});
    for ([_]usize{ 20, 21, 22, 24 }) |m| {
        const K = (@as(usize, 1) << @intCast(m)) / 2;
        const a = try gpa.alloc(u64, K + 1);
        defer gpa.free(a);
        const b = try gpa.alloc(u64, K + 1);
        defer gpa.free(b);
        for (a) |*x| x.* = rnd.uintLessThan(u64, ntt.P);
        for (b) |*x| x.* = rnd.uintLessThan(u64, ntt.P);
        const t0 = nowNs();
        const c = try ntt.convolve(gpa, a, b, K);
        const t1 = nowNs();
        gpa.free(c);
        // one bare forward transform at the same length
        var L: usize = 1;
        while (L < 2 * K + 1) L <<= 1;
        const t = try gpa.alloc(u64, L);
        defer gpa.free(t);
        @memset(t, 0);
        @memcpy(t[0 .. K + 1], a);
        const t2 = nowNs();
        ntt.transform(t, false);
        const t3 = nowNs();
        std.debug.print("   K = {d:>9} (transform 2^{d})  convolve {d:>7.3} s   one transform {d:>7.3} s   ~= N = {e:.1}\n", .{ K, m + 1, @as(f64, @floatFromInt(t1 - t0)) / 1e9, @as(f64, @floatFromInt(t3 - t2)) / 1e9, @as(f64, @floatFromInt(K)) * @as(f64, @floatFromInt(K)) });
    }

    std.debug.print("\n{d} failing checks\n", .{fails});
    if (fails != 0) return error.NttChecksFailed;
}
