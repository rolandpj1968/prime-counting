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

/// Büthe's Logan-function kernel (arXiv:1410.7008, Thms 3.1 + 4.1) — the FKBJ
/// band-limited family, third-generation version. Prime-side support is
/// EXACTLY [x e^-eps, x e^eps] (eps = c/T); the zero side uses every zero in
/// the table with weight l_c(eps*gamma), tail beyond T controlled by ~e^-c.
/// Window cost ~ c/T vs the Gaussian's ~ c^2/T — the structural win.
///
///   l_c(t)   = (c/sinh c) sin(sqrt(t^2-c^2))/sqrt(t^2-c^2)   (sinh form |t|<c)
///   eta_c(y) = (c/(2 sinh c)) I0(c sqrt(1-y^2)) on [-1,1]     (= l_c^hat)
///   mu(t<0)  = -int_{-inf}^t eta_{c,eps};  mu(t>0) = -mu(-t);  nu = int mu
///   Psi(rho) = lam^-1 (Ei1(rho L) + A (rho Ei2 - 2 rho^2 Ei3)(rho L)) l_c(eps g)
///   pi*_{c,eps} = li(x) + A x/L^2 - sum* Psi - ln 2 + trivial + Theta(35 eps)
///   pi*(x)      = pi*_{c,eps} - sum_window (1/m) M(p^m),
///   M(t) = lam^-1 [ mu(u) + (1/log t - 1/2)(mu(u) log(t/x) - nu(u)) ], u = log(t/x)
/// with lam = l_{c,eps}(i/2) = (c/sinh c) sinh(s)/s, s = sqrt(c^2 + eps^2/4),
/// and A = -l''_{c,eps}(0)/2 = eps^2 (coth c / c - 1/c^2)/2.
/// Seam mapping: phi(n) = chi(n<=x) + M(n), so (chi - phi) = -M as required.
const GRID = 1 << 16; // eta/mu/nu cumulative-integral grid over [-1, 1]

const LoganButhe = struct {
    x: u64,
    xf: f64,
    sx: f64,
    L: f64,
    c: f64,
    eps: f64,
    lam: f64,
    acorr: f64, // A_{c,eps}
    lo: u64,
    hi: u64,
    mu_tab: []f64, // mu(eps*y) on y-grid [-1,0], GRID+1 points
    nu_tab: []f64, // nu(eps*y) likewise

    fn bessI0(z: f64) f64 {
        var term: f64 = 1;
        var sum: f64 = 1;
        var n: f64 = 1;
        while (n < 200) : (n += 1) {
            term *= (z * z / 4.0) / (n * n);
            sum += term;
            if (term < 1e-17 * sum) break;
        }
        return sum;
    }

    /// l_c on the real line (even); sinh form below the crossover, sin above,
    /// series at the removable singularity |t| ~ c.
    fn ell(c: f64, t: f64) f64 {
        const norm = c / std.math.sinh(c);
        const w2 = t * t - c * c;
        const aw = @sqrt(@abs(w2));
        if (aw < 1e-3) return norm * (1.0 + w2 / 6.0 + w2 * w2 / 120.0);
        return if (w2 < 0) norm * std.math.sinh(aw) / aw else norm * @sin(aw) / aw;
    }

    fn init(gpa: std.mem.Allocator, x: u64, T: f64, c: f64) !LoganButhe {
        const xf: f64 = @floatFromInt(x);
        const L = @log(xf);
        const eps = c / T;
        const s = @sqrt(c * c + eps * eps / 4.0);
        const lam = (c / std.math.sinh(c)) * std.math.sinh(s) / s;
        const acorr = eps * eps * (std.math.cosh(c) / std.math.sinh(c) / c - 1.0 / (c * c)) / 2.0;
        // eta_c on [-1,1], cumulative trapezoid -> mu, then nu (in y-units;
        // the eps dilation scales mu by 1 (density 1/eps times dy=eps) and nu
        // by eps). mu(y) = -int_{-1}^{y} eta_c ; checked: mu(0-) -> -1/2.
        const h = 2.0 / @as(f64, GRID);
        const norm = c / (2.0 * std.math.sinh(c));
        const mu_tab = try gpa.alloc(f64, GRID / 2 + 1);
        const nu_tab = try gpa.alloc(f64, GRID / 2 + 1);
        var acc: f64 = 0; // integral of eta from -1
        var nacc: f64 = 0; // integral of mu (y-units)
        var prev_eta: f64 = 0; // eta(-1) = norm * I0(0)... = norm; fixed below
        var prev_mu: f64 = 0;
        var gi: usize = 0;
        while (gi <= GRID / 2) : (gi += 1) {
            const y = -1.0 + @as(f64, @floatFromInt(gi)) * h;
            const eta = norm * bessI0(c * @sqrt(@max(0.0, 1.0 - y * y)));
            if (gi > 0) acc += 0.5 * (prev_eta + eta) * h;
            const mu = -acc;
            if (gi > 0) nacc += 0.5 * (prev_mu + mu) * h;
            mu_tab[gi] = mu;
            nu_tab[gi] = nacc * eps; // nu carries one eps dilation factor
            prev_eta = eta;
            prev_mu = mu;
        }
        return .{
            .x = x,
            .xf = xf,
            .sx = @sqrt(xf),
            .L = L,
            .c = c,
            .eps = eps,
            .lam = lam,
            .acorr = acorr,
            .lo = @intFromFloat(xf * @exp(-eps) - 2.0),
            .hi = @intFromFloat(xf * @exp(eps) + 2.0),
            .mu_tab = mu_tab,
            .nu_tab = nu_tab,
        };
    }

    fn lookup(tab: []const f64, y_neg: f64) f64 {
        // y_neg in [-1, 0]; linear interp on the half-grid
        const f = (y_neg + 1.0) * @as(f64, GRID) / 2.0;
        const idx: usize = @intFromFloat(@min(f, @as(f64, GRID / 2) - 1e-9));
        const fr = f - @as(f64, @floatFromInt(idx));
        return tab[idx] * (1.0 - fr) + tab[idx + 1] * fr;
    }

    /// M_{x,c,eps}(n) per Thm 3.1; odd reflection handles n > x.
    fn bigM(k: *const LoganButhe, n: u64) f64 {
        const u = std.math.log1p((@as(f64, @floatFromInt(n)) - k.xf) / k.xf); // log(n/x)
        const y = u / k.eps; // in [-1, 1]
        var mu: f64 = undefined;
        var nu: f64 = undefined;
        // Branch on the integer cut, not sign(y): above 2^53, n in
        // (x - ulp/2, x] rounds to y == 0 and must still take the mu(0-)
        // branch to agree with chi's n <= x — else phi jumps by 1/lam there.
        if (n <= k.x) {
            mu = lookup(k.mu_tab, y);
            nu = lookup(k.nu_tab, y);
        } else {
            // mu odd => int_0^t mu = nu(-t) - nu(0) => nu is EVEN
            mu = -lookup(k.mu_tab, -y);
            nu = lookup(k.nu_tab, -y);
        }
        const logn = k.L + u;
        return (mu + (1.0 / logn - 0.5) * (mu * u - nu)) / k.lam;
    }

    fn phi(k: *const LoganButhe, n: u64) f64 {
        const chi: f64 = if (n <= k.x) 1.0 else 0.0;
        return chi + k.bigM(n);
    }

    /// -2 Re Psi(rho) for the conjugate pair at gamma = g.
    fn zeroTerm(k: *const LoganButhe, g: f64) f64 {
        const w = ell(k.c, k.eps * g);
        if (w == 0.0) return 0.0;
        const rl = C{ .re = 0.5 * k.L, .im = g * k.L }; // rho L
        const irl = rl.inv();
        // Ei_k(rl) ~ e^rl / rl^k * sum_j (k+j-1)!/(k-1)! / rl^j  — Horner
        var s1 = C{ .re = 1, .im = 0 };
        var s2 = C{ .re = 1, .im = 0 };
        var s3 = C{ .re = 1, .im = 0 };
        var j: f64 = 13;
        while (j >= 1) : (j -= 1) {
            s1 = irl.scale(j).mul(s1);
            s1.re += 1;
            s2 = irl.scale(j + 1.0).mul(s2);
            s2.re += 1;
            s3 = irl.scale(j + 2.0).mul(s3);
            s3.re += 1;
        }
        // Ei1 = e^rl/rl * s1; rho Ei2 = e^rl * rho/rl^2 * s2 = e^rl/(L^2 rho) s2...
        // assemble B = s1/rl + A (rho s2/rl^2 - 2 rho^2 s3/rl^3), all over e^rl
        const rho = C{ .re = 0.5, .im = g };
        const irl2 = irl.mul(irl);
        var b2 = rho.mul(irl2).mul(s2);
        const b3 = rho.mul(rho).mul(irl2).mul(irl).mul(s3);
        b2.re -= 2.0 * b3.re;
        b2.im -= 2.0 * b3.im;
        var B = s1.mul(irl);
        B.re += k.acorr * b2.re;
        B.im += k.acorr * b2.im;
        // e^rl = sqrt(x) e^(i g L)
        const e = C{ .re = k.sx * @cos(g * k.L), .im = k.sx * @sin(g * k.L) };
        return -2.0 * (w / k.lam) * e.mul(B).re;
    }

    fn cutoff(k: *const LoganButhe, g: f64) bool {
        _ = k;
        _ = g;
        return false; // every tabulated zero contributes; tail is ~e^-c by design
    }

    fn poleCorr(k: *const LoganButhe) f64 {
        return k.acorr * k.xf / (k.L * k.L);
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
    var cset = false;
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
            cset = true;
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
    } else if (std.mem.eql(u8, kname, "logan")) {
        // tail ~ sqrt(x) polylog e^-c: c must track (1/2) ln x (linear price
        // vs the Gaussian's quadratic — the structural win)
        const cl = if (cset) cpar else 0.5 * @log(@as(f64, @floatFromInt(x))) + 9.0;
        const kern = try LoganButhe.init(gpa, x, T, cl);
        std.debug.print("logan: c = {d:.2}  eps = {e:.4}  lam = {d:.6}\n", .{ cl, kern.eps, kern.lam });
        try run(LoganButhe, &kern, x, zeros.items, nt, gpa, t_start, t_load);
    } else {
        std.debug.print("unknown kernel '{s}' (have: gaussian, logan; beurling-selberg planned)\n", .{kname});
        std.process.exit(2);
    }
}
