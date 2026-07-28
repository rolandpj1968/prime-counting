//! Rung 4+: an integer pi(x) from the zeros, via Riemann's pi* formula —
//! parallel, segmented window, f128 aggregates, and a pluggable KERNEL seam.
//!
//!   pi*(x) = sum_{p^m<=x} 1/m
//!          = li(x) - sum_rho li(x^rho) - ln 2 + integral_x^inf dt/(t(t^2-1)ln t)
//!
//! A kernel is the Mellin pair that smooths the sharp count: a prime-side
//! weight phi(n) (sharp step outside its window) and the matching zero-side
//! per-zero term. Comptime duck-typed interface:
//!
//!   init(x, T, c) K       — freeze all kernel parameters
//!   .lo, .hi: u64         — prime-side window support (phi = 1 below, 0 above)
//!   phi(n) f64            — smoothed weight, called only for n in (lo, hi]
//!   zeroTerm(g) f64       — full contribution of the conjugate zero pair at
//!                            gamma = g to the smoothed pi* (Kahan-added)
//!   cutoff(g) bool        — true when this zero and all above are negligible
//!   poleCorr() f64        — kernel correction to the li(x) pole term
//!
//! Kernels: gaussian (Galway/Platt, the proven baseline), with logan (FKBJ
//! band-limited: exact zero-sum truncation) and beurling-selberg
//! (majorant/minorant sandwich: certified bounds) planned — see
//! literature/buthe-*.pdf and docs/ANALYTIC.md.
//!
//! Parallel geometry mirrors the eventual distribution seam:
//!  - zero sum: static contiguous zero ranges, one Kahan per worker;
//!  - window: atomic dispenser over fixed-width segments, per-worker strips.
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

/// ln in f128 built from correctly-rounded arithmetic only: Zig's @log on
/// f128 is silently 53-bit (compiler_rt stub), and dEi/dL = x/L amplifies
/// ln(x) error by ~3e13 at x = 1e15. Exact power-of-two reduction + atanh
/// series (z <= 0.2 after reduction, ~26 terms to 1e-40).
const LN2_128: f128 = 0.69314718055994530941723212145817657;
fn ln128(v0: f128) f128 {
    var v = v0;
    var k: f128 = 0;
    while (v > 1.5) : (k += 1) v *= 0.5;
    while (v < 0.75) : (k -= 1) v *= 2.0;
    const z = (v - 1.0) / (v + 1.0);
    const z2 = z * z;
    var sum: f128 = 0;
    var term: f128 = z;
    var n: f128 = 1;
    while (n < 80) : (n += 2) {
        sum += term / n;
        if (term < 1e-40 and term > -1e-40) break;
        term *= z2;
    }
    return 2.0 * sum + k * LN2_128;
}

/// Ei(t) for real t > 0 via the convergent series gamma + ln t + sum t^k/(k k!)
/// (all terms positive -- no cancellation). f128 throughout: the f64
/// version's eps*e^L/L error (~0.03 at x = 1e15) was a dominant systematic.
fn realEi(t: f128) f128 {
    const euler: f128 = 0.57721566490153286060651209008240243;
    var s: f128 = euler + ln128(t);
    var p: f128 = 1;
    var i: f128 = 1;
    while (i < 500) : (i += 1) {
        p *= t / i;
        const term = p / i;
        s += term;
        if (term < 1e-36 * s) break;
    }
    return s;
}

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
    .{ .x = 100_000_000_000_000_000, .pi = 2_623_557_157_654_233 },
    .{ .x = 1_000_000_000_000_000_000, .pi = 24_739_954_287_740_860 },
};

// ---------------------------------------------------------------------------
// Kernels
// ---------------------------------------------------------------------------

const ERFC_CUT = 7.5; // erfc(7.5) ~ 4e-26: outside this, phi is exactly 0 or 1

/// Galway's Gaussian pair (Platt eq 3.3): phi_hat(s) = (x^s/s) e^(l^2 s^2/2),
/// phi(t) = (1/2) erfc(log(t/x)/(sqrt(2) l)), lambda = c/T. Per-zero term:
///   F(rho) = e^(l^2 rho^2/2) x^rho [ S/(rho L) - l^2/L^2 ],
/// S = asymptotic li series sum k!/(rho L)^k, K = 13 (|rho L| >= 195 at
/// x >= 1e6 => ~1e-25); -l^2/L^2 is the first-order kernel correction in its
/// cancelled form. Pole correction l^2 x(L-1)/(2L^2); trivial-zero integral
/// (~1/(2 x^2 L)) dropped.
const Gaussian = struct {
    xf: f64,
    sx: f64,
    L: f64,
    lambda: f64,
    l2: f64,
    inv_s2l: f64,
    theta_l: f64,
    lo: u64,
    hi: u64,

    fn init(x: u64, T: f64, c: f64) Gaussian {
        const xf: f64 = @floatFromInt(x);
        const L = @log(xf);
        const lambda = c / T;
        const l2 = lambda * lambda;
        const half_w = ERFC_CUT * std.math.sqrt2 * lambda;
        return .{
            .xf = xf,
            .sx = @sqrt(xf),
            .L = L,
            .lambda = lambda,
            .l2 = l2,
            .inv_s2l = 1.0 / (std.math.sqrt2 * lambda),
            .theta_l = L + l2 / 2.0,
            .lo = @intFromFloat(xf * @exp(-half_w) - 2.0),
            .hi = @intFromFloat(xf * @exp(half_w) + 2.0),
        };
    }

    fn phi(k: *const Gaussian, n: u64) f64 {
        // log(n/x) via log1p on the exact integer offset: no cancellation
        const u = std.math.log1p((@as(f64, @floatFromInt(n)) - k.xf) / k.xf) * k.inv_s2l;
        return 0.5 * erfc(u);
    }

    fn zeroTerm(k: *const Gaussian, g: f64) f64 {
        const damp = @exp(k.l2 * (0.25 - g * g) / 2.0);
        const rl = C{ .re = 0.5 * k.L, .im = g * k.L }; // rho L
        const irl = rl.inv();
        var S = C{ .re = 1, .im = 0 };
        var kk: f64 = 13;
        while (kk >= 1) : (kk -= 1) {
            S = irl.scale(kk).mul(S);
            S.re += 1;
        }
        var br = S.mul(irl); // S/(rho L)
        br.re -= k.l2 / (k.L * k.L); // first-order kernel correction
        const th = g * k.theta_l;
        const e = C{ .re = k.sx * damp * @cos(th), .im = k.sx * damp * @sin(th) };
        return -2.0 * e.mul(br).re;
    }

    fn cutoff(k: *const Gaussian, g: f64) bool {
        return @exp(k.l2 * (0.25 - g * g) / 2.0) < 1e-18;
    }

    fn poleCorr(k: *const Gaussian) f64 {
        return k.l2 * k.xf * (k.L - 1.0) / (2.0 * k.L * k.L);
    }
};

// ---------------------------------------------------------------------------
// Pipeline, generic over the kernel
// ---------------------------------------------------------------------------

fn run(comptime K: type, kern: *const K, x: u64, zeros: []const f64, nt: usize, gpa: std.mem.Allocator, t_start: u64, t_load: u64) !void {
    const ZCtx = struct {
        zeros: []const f64,
        nt: usize,
        kern: *const K,
        sums: []f64,
        fn work(ctx: *@This(), wid: usize) void {
            const n = ctx.zeros.len;
            const chunk = (n + ctx.nt - 1) / ctx.nt;
            const a = @min(wid * chunk, n);
            const b = @min(a + chunk, n);
            var k = Kahan{};
            for (ctx.zeros[a..b]) |g| {
                if (ctx.kern.cutoff(g)) break; // gammas ascend within a range
                k.add(ctx.kern.zeroTerm(g));
            }
            ctx.sums[wid] = k.val();
        }
    };
    const WCtx = struct {
        kern: *const K,
        x: u64,
        base: []const u64,
        next: std.atomic.Value(u64) = std.atomic.Value(u64).init(0),
        sums: []f64,
        counts: []u64,
        gpa: std.mem.Allocator,
        fn work(ctx: *@This(), wid: usize) void {
            const seg = ctx.gpa.alloc(bool, SEGW) catch unreachable;
            defer ctx.gpa.free(seg);
            const lo = ctx.kern.lo;
            const hi = ctx.kern.hi;
            var k = Kahan{};
            var cnt: u64 = 0;
            while (true) {
                const si = ctx.next.fetchAdd(1, .monotonic);
                const s_lo = lo + 1 + si * SEGW;
                if (s_lo > hi) break;
                const s_hi = @min(s_lo + SEGW - 1, hi);
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
                    k.add(chi - ctx.kern.phi(n));
                }
            }
            ctx.sums[wid] = k.val();
            ctx.counts[wid] = cnt;
        }
    };

    // ---- zero sum: static zero-range fragments
    const zsums = try gpa.alloc(f64, nt);
    defer gpa.free(zsums);
    var zctx = ZCtx{ .zeros = zeros, .nt = nt, .kern = kern, .sums = zsums };
    {
        const threads = try gpa.alloc(std.Thread, nt);
        defer gpa.free(threads);
        for (threads, 0..) |*t, i| t.* = try std.Thread.spawn(.{}, ZCtx.work, .{ &zctx, i });
        for (threads) |t| t.join();
    }
    var zk = Kahan{};
    for (zsums) |s| zk.add(s);
    const zerosum = zk.val();
    const t_zeros = nowNs();

    // ---- pole term (smoothed li(x)) and constants — f128 aggregates from
    // here down: at x = 1e15 the f64 ulp of the running ~3e13 totals is
    // already ~4e-3, and it reaches ~0.5 by 1e17.
    const li_x: f128 = realEi(ln128(@floatFromInt(x)));
    const pole_corr = kern.poleCorr();
    const ln2 = @log(2.0);
    const pistar_smooth: f128 = li_x + @as(f128, pole_corr) + @as(f128, zerosum) - @as(f128, ln2);
    std.debug.print("li(x) = {d:.6}  kernel-corr = {d:.6}  zerosum = {d:.6}  -ln2\n", .{ @as(f64, @floatCast(li_x)), pole_corr, zerosum });
    std.debug.print("pi*_smooth(zeros) = {d:.6}   ({d:.1}s load, {d:.1}s zerosum)\n", .{ @as(f64, @floatCast(pistar_smooth)), @as(f64, @floatFromInt(t_load - t_start)) / 1e9, @as(f64, @floatFromInt(t_zeros - t_load)) / 1e9 });

    // ---- window: segmented, dispenser-parallel; (1/m)(chi - phi) corrections
    const wsq = iroot(kern.hi, 2);
    std.debug.assert(wsq < kern.lo); // base primes never lie inside the window
    const base = try rs.basePrimes(gpa, wsq);
    defer gpa.free(base);
    const wsums = try gpa.alloc(f64, nt);
    defer gpa.free(wsums);
    const wcounts = try gpa.alloc(u64, nt);
    defer gpa.free(wcounts);
    var wctx = WCtx{ .kern = kern, .x = x, .base = base, .sums = wsums, .counts = wcounts, .gpa = gpa };
    {
        const threads = try gpa.alloc(std.Thread, nt);
        defer gpa.free(threads);
        for (threads, 0..) |*t, i| t.* = try std.Thread.spawn(.{}, WCtx.work, .{ &wctx, i });
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
        while (pk <= kern.hi / p) {
            pk *= p;
            m += 1;
            if (pk <= kern.lo) continue;
            const chi: f64 = if (pk <= x) 1.0 else 0.0;
            wk.add((chi - kern.phi(pk)) / m);
        }
    }
    const wcorr = wk.val();
    const pistar_sharp: f128 = pistar_smooth + @as(f128, wcorr);
    const t_win = nowNs();
    std.debug.print("window [{d}, {d}]: {d} ints, {d} primes  corr = {d:.6}   ({d:.1}s window)\n", .{ kern.lo + 1, kern.hi, kern.hi - kern.lo, wprimes, wcorr, @as(f64, @floatFromInt(t_win - t_zeros)) / 1e9 });
    std.debug.print("pi*_sharp = {d:.6}\n", .{@as(f64, @floatCast(pistar_sharp))});

    // ---- Moebius unwind with exact tiny pi*(x^(1/m)); the small terms sum
    // in f64 Kahan (magnitude ~1e5), the aggregate stays f128
    var pk = Kahan{};
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
    const pi_real: f128 = pistar_sharp + @as(f128, pk.val());
    const pi_round: u64 = @intFromFloat(@round(pi_real));
    const err: f64 = @floatCast(pi_real - @round(pi_real));
    const secs = @as(f64, @floatFromInt(nowNs() - t_start)) / 1e9;

    std.debug.print("\npi(x) = {d}  err = {d:.6}  (margin to half-integer: {d:.4})   {d:.1}s\n", .{ pi_round, err, 0.5 - @abs(err), secs });
    for (KNOWN) |kv| {
        if (kv.x == x) {
            const ok = kv.pi == pi_round;
            std.debug.print("published pi(x) = {d}   {s}\n", .{ kv.pi, if (ok) "MATCH" else "*** MISMATCH ***" });
            std.process.exit(if (ok) 0 else 1);
        }
    }
    std.debug.print("(no published reference in table)\n", .{});
}

pub fn main(init: std.process.Init) !void {
    var args = try init.minimal.args.iterateAllocator(init.gpa);
    _ = args.skip();
    const usage = "Usage: pistar <x> <zeros-file> [c] [-t nthreads] [--kernel gaussian]\n";
    const x_str = args.next() orelse return std.debug.print(usage, .{});
    const x = try std.fmt.parseInt(u64, x_str, 10);
    const zpath = args.next() orelse return std.debug.print(usage, .{});
    var cpar: f64 = 7.5;
    var nt: usize = 6;
    var kname: []const u8 = "gaussian";
    while (args.next()) |a| {
        if (std.mem.eql(u8, a, "-t")) {
            const v = args.next() orelse return std.debug.print(usage, .{});
            nt = try std.fmt.parseInt(usize, v, 10);
        } else if (std.mem.eql(u8, a, "--kernel")) {
            kname = args.next() orelse return std.debug.print(usage, .{});
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
    std.debug.print("x = {d}  zeros = {d} (T = {d:.3})  kernel = {s} (c = {d})  threads = {d}\n", .{ x, zeros.items.len, T, kname, cpar, nt });

    if (std.mem.eql(u8, kname, "gaussian")) {
        const kern = Gaussian.init(x, T, cpar);
        try run(Gaussian, &kern, x, zeros.items, nt, gpa, t_start, t_load);
    } else {
        std.debug.print("unknown kernel '{s}' (have: gaussian; logan and beurling-selberg planned)\n", .{kname});
        std.process.exit(2);
    }
}
