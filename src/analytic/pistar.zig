//! Rung 3: an integer pi(x) from the zeros, via Riemann's pi* formula.
//!
//!   pi*(x) = sum_{p^m<=x} 1/m
//!          = li(x) - sum_rho li(x^rho) - ln 2 + integral_x^inf dt/(t(t^2-1)ln t)
//!
//! smoothed with Galway's Gaussian pair exactly as rung 2. Each zero's term is
//!   F(rho) = int_{-inf}^{rho} x^z e^(l^2 z^2/2) / z dz   (horizontal path)
//!          = e^(l^2 rho^2/2) x^rho [ S/(rho L) - l^2/L^2 + O(small) ]
//! where S is the (convergent-in-practice) asymptotic li series
//!   S = sum_{k>=0} k!/(rho L)^k,  |rho L| >= 195 at x >= 1e6 => K=13 suffices,
//! and the -l^2/L^2 term is the first-order kernel correction
//! (int x^z(z-rho)/z dz = x^rho/L - rho li(x^rho) ~ -x^rho/(rho L^2), times
//! l^2 rho). Next order is O(1e-3) coherent at 1e12 -- budget-noted.
//! The pole term becomes li(x) + l^2 x(L-1)/(2L^2) + O(l^4); the trivial-zero
//! integral is ~1/(2 x^2 L) and is dropped.
//!
//! Sharp pi*(x) = smoothed zero sum + sum_{p^m in window} (1/m)(chi - phi),
//! window sieved exactly. Then Moebius-unwind with EXACT tiny pi*(x^(1/m)):
//!   pi(x) = sum_m mu(m)/m pi*(x^(1/m)),
//! round, report margin, and check against published pi(x).

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

const ERFC_CUT = 7.5;

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

const MU = [_]i8{ 0, 1, -1, -1, 0, -1, 1, -1, 0, 0, 1, -1, 0, -1, 1, 1, 0, -1, 0, -1, 0, 1, 1, -1, 0, 0, 1, 0, 0, -1, -1, -1, 0, 1, 1, 1, 0, -1, 1, 1, 0 };

const KNOWN = [_]struct { x: u64, pi: u64 }{
    .{ .x = 1_000_000, .pi = 78_498 },
    .{ .x = 10_000_000, .pi = 664_579 },
    .{ .x = 100_000_000, .pi = 5_761_455 },
    .{ .x = 1_000_000_000, .pi = 50_847_534 },
    .{ .x = 10_000_000_000, .pi = 455_052_511 },
    .{ .x = 100_000_000_000, .pi = 4_118_054_813 },
    .{ .x = 1_000_000_000_000, .pi = 37_607_912_018 },
    .{ .x = 10_000_000_000_000, .pi = 346_065_536_839 },
};

pub fn main(init: std.process.Init) !void {
    var args = try init.minimal.args.iterateAllocator(init.gpa);
    _ = args.skip();
    const usage = "Usage: pistar <x> <zeros-file> [c: lambda = c/T, default 7.5]\n";
    const x_str = args.next() orelse return std.debug.print(usage, .{});
    const x = try std.fmt.parseInt(u64, x_str, 10);
    const zpath = args.next() orelse return std.debug.print(usage, .{});
    const cpar: f64 = if (args.next()) |s| try std.fmt.parseFloat(f64, s) else 7.5;

    const gpa = std.heap.page_allocator;
    const t0 = nowNs();

    var tio = std.Io.Threaded.init(gpa, .{});
    const data = try std.Io.Dir.cwd().readFileAlloc(tio.io(), zpath, gpa, .limited(1 << 27));
    defer gpa.free(data);
    var zeros = try std.ArrayList(f64).initCapacity(gpa, 1 << 21);
    defer zeros.deinit(gpa);
    var lines = std.mem.tokenizeScalar(u8, data, '\n');
    while (lines.next()) |line| {
        const s = std.mem.trim(u8, line, " \t\r");
        if (s.len == 0) continue;
        try zeros.append(gpa, try std.fmt.parseFloat(f64, s));
    }

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

    std.debug.print("x = {d}  zeros = {d} (T = {d:.3})  lambda = {e:.4} (c = {d})\n", .{ x, zeros.items.len, T, lambda, cpar });

    // ---- zero sum: -2 Re F(rho) per positive-gamma zero
    var zk = Kahan{};
    for (zeros.items) |g| {
        const damp = @exp(l2 * (0.25 - g * g) / 2.0);
        if (damp < 1e-18) break;
        const rl = C{ .re = 0.5 * L, .im = g * L }; // rho L
        const irl = rl.inv();
        var S = C{ .re = 1, .im = 0 };
        var k: f64 = 13;
        while (k >= 1) : (k -= 1) {
            S = irl.scale(k).mul(S);
            S.re += 1;
        }
        var br = S.mul(irl); // S/(rho L)
        br.re -= l2 / (L * L); // first-order kernel correction
        const th = g * (L + l2 / 2.0);
        const e = C{ .re = sx * damp * @cos(th), .im = sx * damp * @sin(th) };
        zk.add(-2.0 * e.mul(br).re);
    }
    const zerosum = zk.val();

    // ---- pole term (smoothed li(x)) and constants
    const li_x = realEi(L);
    const main_corr = l2 * xf * (L - 1.0) / (2.0 * L * L);
    const ln2 = @log(2.0);
    const pistar_smooth = li_x + main_corr + zerosum - ln2;
    std.debug.print("li(x) = {d:.6}  kernel-corr = {d:.6}  zerosum = {d:.6}  -ln2\n", .{ li_x, main_corr, zerosum });
    std.debug.print("pi*_smooth(zeros) = {d:.6}\n", .{pistar_smooth});

    // ---- window: sieve [lo+1, hi] exactly, correct (1/m)(chi - phi)
    const wlen: usize = @intCast(hi - lo);
    const wsq = iroot(hi, 2);
    const base = try rs.basePrimes(gpa, wsq);
    defer gpa.free(base);
    const comp = try gpa.alloc(bool, wlen); // comp[i] : n = lo+1+i composite
    defer gpa.free(comp);
    @memset(comp, false);
    for (base) |p| {
        var q = (lo / p + 1) * p;
        while (q <= hi) : (q += p) {
            if (q > lo) comp[@intCast(q - lo - 1)] = true;
        }
    }
    var wk = Kahan{};
    var wprimes: u64 = 0;
    var i: usize = 0;
    while (i < wlen) : (i += 1) {
        if (comp[i]) continue;
        const n = lo + 1 + i;
        if (n <= wsq) continue; // base primes handled among prime powers below
        wprimes += 1;
        const chi: f64 = if (n <= x) 1.0 else 0.0;
        const u = std.math.log1p((@as(f64, @floatFromInt(n)) - xf) / xf) * inv_s2l;
        wk.add(chi - 0.5 * erfc(u));
    }
    // prime powers p^m (and base primes p) inside the window
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
    std.debug.print("window [{d}, {d}]: {d} ints, {d} primes  corr = {d:.6}\n", .{ lo + 1, hi, wlen, wprimes, wcorr });
    std.debug.print("pi*_sharp = {d:.6}\n", .{pistar_sharp});

    // ---- Moebius unwind with exact tiny pi*(x^(1/m))
    var pk = Kahan{};
    pk.add(pistar_sharp);
    var m: u32 = 2;
    while (m < MU.len) : (m += 1) {
        if (MU[m] == 0) continue;
        const v = iroot(x, m);
        if (v < 2) break;
        const pse = try pistarExact(gpa, v);
        const term = @as(f64, @floatFromInt(MU[m])) * pse / @as(f64, @floatFromInt(m));
        pk.add(term);
        std.debug.print("  m={d:<2} x^(1/m)={d:<8} mu/m*pi* = {d:.6}\n", .{ m, v, term });
    }
    const pi_real = pk.val();
    const pi_round: u64 = @intFromFloat(@round(pi_real));
    const margin = @abs(pi_real - @round(pi_real));
    const secs = @as(f64, @floatFromInt(nowNs() - t0)) / 1e9;

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
