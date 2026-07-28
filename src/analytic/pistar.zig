//! Rung 3+: an integer pi(x) from the zeros, via Riemann's pi* formula —
//! parallel, with a segmented window (no resident window array).
//!
//!   pi*(x) = sum_{p^m<=x} 1/m
//!          = li(x) - sum_rho li(x^rho) - ln 2 + integral_x^inf dt/(t(t^2-1)ln t)
//!
//! smoothed with Galway's Gaussian pair (Platt eq 3.3). Each zero's term:
//!   F(rho) = e^(l^2 rho^2/2) x^rho [ S/(rho L) - l^2/L^2 ],
//! S = sum_{k<=13} k!/(rho L)^k (|rho L| >= 195 at x >= 1e6 => ~1e-25), and
//! -l^2/L^2 is the first-order kernel correction in its cancelled form. Pole
//! term li(x) + l^2 x(L-1)/(2L^2); trivial-zero integral (~1/(2x^2 L)) dropped.
//! Sharp pi* = zero sum + sum_window (1/m)(chi - phi) over an exactly sieved
//! window; Moebius with exact tiny pi*(x^(1/m)) gives pi(x); round and check.
//!
//! Parallel geometry mirrors the eventual distribution seam:
//!  - zero sum: static contiguous zero ranges, one Kahan per worker (a zero
//!    range is an additive fragment; the table streams);
//!  - window: atomic dispenser over fixed-width segments, each sieved in a
//!    per-worker strip (a segment is an additive fragment; needs no zeros).
//! Memory is O(nt * SEGW + pi(sqrt(hi))) regardless of window length.

const std = @import("std");
// zig build-exe -O ReleaseFast -mcpu=native -lc --dep rs -Mroot=pistar.zig -Mrs=../rangesieve.zig -femit-bin=pistar
const rs = @import("rs");

extern "c" fn erfc(x: f64) f64;

fn nowNs() u64 {
    var ts: std.os.linux.timespec = undefined;
    _ = std.os.linux.clock_gettime(.MONOTONIC, &ts);
    return @as(u64, @intCast(ts.sec)) * 1_000_000_000 + @as(u64, @intCast(ts.nsec));
}

const Kahan = struct {
    s: f64 = 0,
    c: f64 = 0,
    fn add(k: *Kahan, v: f64) void {
        const t = k.s + v;
        if (@abs(k.s) >= @abs(v)) {
            k.c += (k.s - t) + v;
        } else {
            k.c += (v - t) + k.s;
        }
        k.s = t;
    }
    fn val(k: *const Kahan) f64 {
        return k.s + k.c;
    }
};

const C = struct {
    re: f64,
    im: f64,
    fn mul(a: C, b: C) C {
        return .{ .re = a.re * b.re - a.im * b.im, .im = a.re * b.im + a.im * b.re };
    }
    fn inv(a: C) C {
        const d = a.re * a.re + a.im * a.im;
        return .{ .re = a.re / d, .im = -a.im / d };
    }
    fn scale(a: C, s: f64) C {
        return .{ .re = a.re * s, .im = a.im * s };
    }
};

/// Ei(t) for real t > 0 via the convergent series gamma + ln t + sum t^k/(k k!)
/// (all terms positive -- no cancellation; ~120 terms at t = 30).
fn realEi(t: f64) f64 {
    const euler = 0.577215664901532861;
    var k = Kahan{};
    k.add(euler);
    k.add(@log(t));
    var p: f64 = 1;
    var i: f64 = 1;
    while (i < 200) : (i += 1) {
        p *= t / i;
        const term = p / i;
        k.add(term);
        if (term < 1e-18 * @abs(k.val())) break;
    }
    return k.val();
}

const ERFC_CUT = 7.5; // erfc(7.5) ~ 4e-26: outside this, phi is exactly 0 or 1
const SEGW: u64 = 1 << 24; // window segment width (16M integers per strip)

/// floor(x^(1/m)) exactly, by float guess + integer correction.
fn iroot(x: u64, m: u32) u64 {
    if (m == 1) return x;
    var v: u64 = @intFromFloat(std.math.pow(f64, @floatFromInt(x), 1.0 / @as(f64, @floatFromInt(m))));
    while (ipow(v + 1, m) <= x) v += 1;
    while (v > 1 and ipow(v, m) > x) v -= 1;
    return v;
}
fn ipow(v: u64, m: u32) u64 {
    var r: u128 = 1;
    var i: u32 = 0;
    while (i < m) : (i += 1) {
        r *= v;
        if (r > 1 << 70) return 1 << 63; // saturate: only compared against x
    }
    return @intCast(@min(r, 1 << 63));
}

/// Exact pi*(v) = sum_{p^k <= v} 1/k for small v, off a fresh sieve.
fn pistarExact(gpa: std.mem.Allocator, v: u64) !f64 {
    if (v < 2) return 0;
    const primes = try rs.basePrimes(gpa, v);
    defer gpa.free(primes);
    var k = Kahan{};
    for (primes) |p| {
        k.add(1.0);
        var pk: u64 = p;
        var m: f64 = 1;
        while (pk <= v / p) {
            pk *= p;
            m += 1;
            k.add(1.0 / m);
        }
    }
    return k.val();
}

/// mu(m) by trial factorization — m never exceeds log2(x) <= 64 here.
fn mobius(m0: u32) i32 {
    var m = m0;
    var mu: i32 = 1;
    var p: u32 = 2;
    while (p * p <= m) : (p += 1) {
        if (m % p == 0) {
            m /= p;
            if (m % p == 0) return 0;
            mu = -mu;
        }
    }
    if (m > 1) mu = -mu;
    return mu;
}

const KNOWN = [_]struct { x: u64, pi: u64 }{
    .{ .x = 1_000_000, .pi = 78_498 },
    .{ .x = 10_000_000, .pi = 664_579 },
    .{ .x = 100_000_000, .pi = 5_761_455 },
    .{ .x = 1_000_000_000, .pi = 50_847_534 },
    .{ .x = 10_000_000_000, .pi = 455_052_511 },
    .{ .x = 100_000_000_000, .pi = 4_118_054_813 },
    .{ .x = 1_000_000_000_000, .pi = 37_607_912_018 },
    .{ .x = 10_000_000_000_000, .pi = 346_065_536_839 },
    .{ .x = 100_000_000_000_000, .pi = 3_204_941_750_802 },
    .{ .x = 1_000_000_000_000_000, .pi = 29_844_570_422_669 },
    .{ .x = 10_000_000_000_000_000, .pi = 279_238_341_033_925 },
};

const ZeroCtx = struct {
    zeros: []const f64,
    nt: usize,
    sx: f64,
    l2: f64,
    theta_l: f64,
    L: f64,
    sums: []f64,
};

fn zeroWorker(ctx: *ZeroCtx, wid: usize) void {
    const n = ctx.zeros.len;
    const chunk = (n + ctx.nt - 1) / ctx.nt;
    const a = @min(wid * chunk, n);
    const b = @min(a + chunk, n);
    var k = Kahan{};
    for (ctx.zeros[a..b]) |g| {
        const damp = @exp(ctx.l2 * (0.25 - g * g) / 2.0);
        if (damp < 1e-18) break; // gammas ascend within a range
        const rl = C{ .re = 0.5 * ctx.L, .im = g * ctx.L };
        const irl = rl.inv();
        var S = C{ .re = 1, .im = 0 };
        var kk: f64 = 13;
        while (kk >= 1) : (kk -= 1) {
            S = irl.scale(kk).mul(S);
            S.re += 1;
        }
        var br = S.mul(irl); // S/(rho L)
        br.re -= ctx.l2 / (ctx.L * ctx.L); // first-order kernel correction
        const th = g * ctx.theta_l;
        const e = C{ .re = ctx.sx * damp * @cos(th), .im = ctx.sx * damp * @sin(th) };
        k.add(-2.0 * e.mul(br).re);
    }
    ctx.sums[wid] = k.val();
}

const WinCtx = struct {
    lo: u64, // window is (lo, hi]
    hi: u64,
    x: u64,
    xf: f64,
    inv_s2l: f64,
    wsq: u64,
    base: []const u64,
    next: std.atomic.Value(u64) = std.atomic.Value(u64).init(0),
    sums: []f64,
    counts: []u64,
    gpa: std.mem.Allocator,
};

fn winWorker(ctx: *WinCtx, wid: usize) void {
    const seg = ctx.gpa.alloc(bool, SEGW) catch unreachable;
    defer ctx.gpa.free(seg);
    var k = Kahan{};
    var cnt: u64 = 0;
    while (true) {
        const si = ctx.next.fetchAdd(1, .monotonic);
        const s_lo = ctx.lo + 1 + si * SEGW;
        if (s_lo > ctx.hi) break;
        const s_hi = @min(s_lo + SEGW - 1, ctx.hi);
        const len: usize = @intCast(s_hi - s_lo + 1);
        @memset(seg[0..len], false);
        for (ctx.base) |p| {
            var q = ((s_lo + p - 1) / p) * p;
            while (q <= s_hi) : (q += p) seg[@intCast(q - s_lo)] = true;
        }
        for (seg[0..len], 0..) |comp, i| {
            if (comp) continue;
            const n = s_lo + i;
            cnt += 1;
            const chi: f64 = if (n <= ctx.x) 1.0 else 0.0;
            const u = std.math.log1p((@as(f64, @floatFromInt(n)) - ctx.xf) / ctx.xf) * ctx.inv_s2l;
            k.add(chi - 0.5 * erfc(u));
        }
    }
    ctx.sums[wid] = k.val();
    ctx.counts[wid] = cnt;
}

pub fn main(init: std.process.Init) !void {
    var args = try init.minimal.args.iterateAllocator(init.gpa);
    _ = args.skip();
    const usage = "Usage: pistar <x> <zeros-file> [c: lambda = c/T, default 7.5] [-t nthreads]\n";
    const x_str = args.next() orelse return std.debug.print(usage, .{});
    const x = try std.fmt.parseInt(u64, x_str, 10);
    const zpath = args.next() orelse return std.debug.print(usage, .{});
    var cpar: f64 = 7.5;
    var nt: usize = 6;
    while (args.next()) |a| {
        if (std.mem.eql(u8, a, "-t")) {
            const v = args.next() orelse return std.debug.print(usage, .{});
            nt = try std.fmt.parseInt(usize, v, 10);
        } else {
            cpar = try std.fmt.parseFloat(f64, a);
        }
    }

    const gpa = std.heap.page_allocator;
    const t_start = nowNs();

    var tio = std.Io.Threaded.init(gpa, .{});
    const data = try std.Io.Dir.cwd().readFileAlloc(tio.io(), zpath, gpa, .limited(1 << 33));
    var zeros = try std.ArrayList(f64).initCapacity(gpa, 1 << 21);
    defer zeros.deinit(gpa);
    var lines = std.mem.tokenizeScalar(u8, data, '\n');
    while (lines.next()) |line| {
        const s = std.mem.trim(u8, line, " \t\r");
        if (s.len == 0) continue;
        try zeros.append(gpa, try std.fmt.parseFloat(f64, s));
    }
    gpa.free(data);
    const t_load = nowNs();

    const T = zeros.items[zeros.items.len - 1];
    const lambda = cpar / T;
    const l2 = lambda * lambda;
    const xf: f64 = @floatFromInt(x);
    const L = @log(xf);
    const sx = @sqrt(xf);
    const inv_s2l = 1.0 / (std.math.sqrt2 * lambda);
    const half_w = ERFC_CUT * std.math.sqrt2 * lambda;
    const lo: u64 = @intFromFloat(xf * @exp(-half_w) - 2.0);
    const hi: u64 = @intFromFloat(xf * @exp(half_w) + 2.0);

    std.debug.print("x = {d}  zeros = {d} (T = {d:.3})  lambda = {e:.4} (c = {d})  threads = {d}\n", .{ x, zeros.items.len, T, lambda, cpar, nt });

    // ---- zero sum: static zero-range fragments
    const zsums = try gpa.alloc(f64, nt);
    defer gpa.free(zsums);
    var zctx = ZeroCtx{ .zeros = zeros.items, .nt = nt, .sx = sx, .l2 = l2, .theta_l = L + l2 / 2.0, .L = L, .sums = zsums };
    {
        const threads = try gpa.alloc(std.Thread, nt);
        defer gpa.free(threads);
        for (threads, 0..) |*t, i| t.* = try std.Thread.spawn(.{}, zeroWorker, .{ &zctx, i });
        for (threads) |t| t.join();
    }
    var zk = Kahan{};
    for (zsums) |s| zk.add(s);
    const zerosum = zk.val();
    const t_zeros = nowNs();

    // ---- pole term (smoothed li(x)) and constants
    const li_x = realEi(L);
    const main_corr = l2 * xf * (L - 1.0) / (2.0 * L * L);
    const ln2 = @log(2.0);
    const pistar_smooth = li_x + main_corr + zerosum - ln2;
    std.debug.print("li(x) = {d:.6}  kernel-corr = {d:.6}  zerosum = {d:.6}  -ln2\n", .{ li_x, main_corr, zerosum });
    std.debug.print("pi*_smooth(zeros) = {d:.6}   ({d:.1}s load, {d:.1}s zerosum)\n", .{ pistar_smooth, @as(f64, @floatFromInt(t_load - t_start)) / 1e9, @as(f64, @floatFromInt(t_zeros - t_load)) / 1e9 });

    // ---- window: segmented, dispenser-parallel; (1/m)(chi - phi) corrections
    const wsq = iroot(hi, 2);
    std.debug.assert(wsq < lo); // base primes never lie inside the window
    const base = try rs.basePrimes(gpa, wsq);
    defer gpa.free(base);
    const wsums = try gpa.alloc(f64, nt);
    defer gpa.free(wsums);
    const wcounts = try gpa.alloc(u64, nt);
    defer gpa.free(wcounts);
    var wctx = WinCtx{ .lo = lo, .hi = hi, .x = x, .xf = xf, .inv_s2l = inv_s2l, .wsq = wsq, .base = base, .sums = wsums, .counts = wcounts, .gpa = gpa };
    {
        const threads = try gpa.alloc(std.Thread, nt);
        defer gpa.free(threads);
        for (threads, 0..) |*t, i| t.* = try std.Thread.spawn(.{}, winWorker, .{ &wctx, i });
        for (threads) |t| t.join();
    }
    var wk = Kahan{};
    var wprimes: u64 = 0;
    for (wsums) |s| wk.add(s);
    for (wcounts) |c| wprimes += c;
    // prime powers p^m (and base primes p) inside the window — serial, tiny
    for (base) |p| {
        var pk: u64 = 1;
        var m: f64 = 0;
        while (pk <= hi / p) {
            pk *= p;
            m += 1;
            if (pk <= lo) continue;
            const chi: f64 = if (pk <= x) 1.0 else 0.0;
            const u = std.math.log1p((@as(f64, @floatFromInt(pk)) - xf) / xf) * inv_s2l;
            wk.add((chi - 0.5 * erfc(u)) / m);
        }
    }
    const wcorr = wk.val();
    const pistar_sharp = pistar_smooth + wcorr;
    const t_win = nowNs();
    std.debug.print("window [{d}, {d}]: {d} ints, {d} primes  corr = {d:.6}   ({d:.1}s window)\n", .{ lo + 1, hi, hi - lo, wprimes, wcorr, @as(f64, @floatFromInt(t_win - t_zeros)) / 1e9 });
    std.debug.print("pi*_sharp = {d:.6}\n", .{pistar_sharp});

    // ---- Moebius unwind with exact tiny pi*(x^(1/m))
    var pk = Kahan{};
    pk.add(pistar_sharp);
    var m: u32 = 2;
    while (m <= 64) : (m += 1) {
        const mu = mobius(m);
        if (mu == 0) continue;
        const v = iroot(x, m);
        if (v < 2) break;
        const pse = try pistarExact(gpa, v);
        const term = @as(f64, @floatFromInt(mu)) * pse / @as(f64, @floatFromInt(m));
        pk.add(term);
        std.debug.print("  m={d:<2} x^(1/m)={d:<8} mu/m*pi* = {d:.6}\n", .{ m, v, term });
    }
    const pi_real = pk.val();
    const pi_round: u64 = @intFromFloat(@round(pi_real));
    const margin = @abs(pi_real - @round(pi_real));
    const secs = @as(f64, @floatFromInt(nowNs() - t_start)) / 1e9;

    std.debug.print("\npi(x) = {d:.6}  ->  {d}   (margin to half-integer: {d:.4})   {d:.1}s\n", .{ pi_real, pi_round, 0.5 - margin, secs });
    for (KNOWN) |kv| {
        if (kv.x == x) {
            const ok = kv.pi == pi_round;
            std.debug.print("published pi(x) = {d}   {s}\n", .{ kv.pi, if (ok) "MATCH" else "*** MISMATCH ***" });
            std.process.exit(if (ok) 0 else 1);
        }
    }
    std.debug.print("(no published reference in table)\n", .{});
}
