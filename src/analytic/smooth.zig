//! Rung 2: Gaussian-smoothed explicit formula (Galway's kernel, per Platt eq 3.3),
//! applied to psi. Mellin pair: phi_hat(s) = (x^s/s) exp(lambda^2 s^2 / 2),
//! phi(t) = (1/2) erfc(log(t/x) / (sqrt(2) lambda)).
//!
//! Smoothed psi both ways:
//!   direct:  psi~(x) = sum_n Lambda(n) phi(n)            (sieve referee)
//!   zeros:   psi~(x) = x e^(l^2/2) - sum_rho phi_hat(rho) - ln 2pi - tail
//! Pair term for rho = 1/2 + i gamma:
//!   amp   = sqrt(x) exp(l^2 (1/4 - g^2) / 2)      <- Gaussian damping in gamma
//!   theta = g (L + l^2/2)                          <- phase shift from exp(l^2 rho^2/2)
//!   pair  = 2 amp (cos(theta)/2 + g sin(theta)) / (1/4 + g^2)
//! Trivial-zero tail: exp(2 l^2 k^2) factors differ from 1 by ~1e-10 at our
//! lambda, so the unsmoothed -ln 2pi - (1/2)ln(1-x^-2) is kept as-is.
//!
//! lambda = c/T annihilates the dropped tail (exp(-c^2/2) at gamma = T); the
//! price is a prime-side window of ~ x * 2*sqrt(2)*erfc_cut*lambda integers.

const std = @import("std");
// zig build-exe -O ReleaseFast -mcpu=native -lc --dep rs -Mroot=smooth.zig -Mrs=../rangesieve.zig -femit-bin=smooth
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

const ERFC_CUT = 7.5; // erfc(7.5) ~ 4e-26: outside this, phi is exactly 0 or 1

const Smoother = struct {
    xf: f64,
    inv_s2l: f64, // 1 / (sqrt(2) lambda)
    lo: u64, // below: phi = 1
    hi: u64, // above: phi = 0

    fn init(x: u64, lambda: f64) Smoother {
        const xf: f64 = @floatFromInt(x);
        const half_w = ERFC_CUT * std.math.sqrt2 * lambda;
        return .{
            .xf = xf,
            .inv_s2l = 1.0 / (std.math.sqrt2 * lambda),
            .lo = @intFromFloat(xf * @exp(-half_w) - 2.0),
            .hi = @intFromFloat(xf * @exp(half_w) + 2.0),
        };
    }
    fn phi(sm: *const Smoother, n: u64) f64 {
        if (n <= sm.lo) return 1.0;
        if (n > sm.hi) return 0.0;
        // log(n/x) via log1p on the exact integer offset: no cancellation
        const u = std.math.log1p((@as(f64, @floatFromInt(n)) - sm.xf) / sm.xf) * sm.inv_s2l;
        return 0.5 * erfc(u);
    }
};

fn psiSmoothDirect(gpa: std.mem.Allocator, sm: *const Smoother, win_count: *u64) !f64 {
    const primes = try rs.basePrimes(gpa, sm.hi);
    defer gpa.free(primes);
    var k = Kahan{};
    for (primes) |p| {
        const w = sm.phi(p);
        if (w != 0.0) k.add(@log(@as(f64, @floatFromInt(p))) * w);
        if (p > sm.lo and p <= sm.hi) win_count.* += 1;
    }
    for (primes) |p| {
        if (p * p > sm.hi) break;
        const lp = @log(@as(f64, @floatFromInt(p)));
        var pk = p * p;
        while (pk <= sm.hi) : (pk *= p) {
            const w = sm.phi(pk);
            if (w != 0.0) k.add(lp * w);
        }
    }
    return k.val();
}

pub fn main(init: std.process.Init) !void {
    var args = try init.minimal.args.iterateAllocator(init.gpa);
    _ = args.skip();
    const usage = "Usage: smooth <x> <zeros-file> [c: lambda = c/T, default 7.5]\n";
    const x_str = args.next() orelse return std.debug.print(usage, .{});
    const x = try std.fmt.parseInt(u64, x_str, 10);
    const zpath = args.next() orelse return std.debug.print(usage, .{});
    const c: f64 = if (args.next()) |s| try std.fmt.parseFloat(f64, s) else 7.5;

    const gpa = std.heap.page_allocator;

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
    const lambda = c / T;
    const l2 = lambda * lambda;
    const sm = Smoother.init(x, lambda);

    const xf: f64 = @floatFromInt(x);
    const L = @log(xf);
    const sx = @sqrt(xf);
    const ln2pi = @log(2.0 * std.math.pi);
    const tail = 0.5 * @log(1.0 - 1.0 / (xf * xf));

    var win_count: u64 = 0;
    const t0 = nowNs();
    const psi_ref = try psiSmoothDirect(gpa, &sm, &win_count);
    const t_sieve = @as(f64, @floatFromInt(nowNs() - t0)) / 1e9;

    std.debug.print("x = {d}  zeros = {d} (T = {d:.3})  lambda = {e:.4} (c = {d})\n", .{ x, zeros.items.len, T, lambda, c });
    std.debug.print("window [{d}, {d}]: {d} integers, {d} primes\n", .{ sm.lo + 1, sm.hi, sm.hi - sm.lo, win_count });
    std.debug.print("psi~_direct = {d:.9}   ({d:.1}s sieve)\n\n", .{ psi_ref, t_sieve });
    std.debug.print("{s:>9} {s:>12} {s:>22} {s:>14}\n", .{ "N", "T", "psi~_zeros", "err" });

    var k = Kahan{};
    var next_cp: usize = 100;
    const theta_l = L + l2 / 2.0;
    for (zeros.items, 0..) |g, i| {
        const amp = sx * @exp(l2 * (0.25 - g * g) / 2.0);
        const th = g * theta_l;
        k.add(2.0 * amp * (0.5 * @cos(th) + g * @sin(th)) / (0.25 + g * g));
        const n = i + 1;
        if (n == next_cp or n == zeros.items.len) {
            const psi_z = xf * @exp(l2 / 2.0) - k.val() - ln2pi - tail;
            std.debug.print("{d:>9} {d:>12.1} {d:>22.6} {d:>14.6}\n", .{ n, g, psi_z, psi_z - psi_ref });
            next_cp *= 10;
        }
    }
}
