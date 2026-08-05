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

const SEGW: u64 = 1 << 26; // window segment width (64M per strip: at 1e19 the base-prime
// scan per segment costs pi(sqrt hi) ~ 1.5e8 visits — fewer, fatter segments
// keep that overhead below the marking work itself

/// floor(x^(1/m)) exactly, by float guess + integer correction. The old
/// ipow saturated at 1<<63, "only compared against x" — true until x
/// reached 1e19 > 2^63, when the saturation sat BELOW x and the upward
/// correction spun forever (3.5h of v += 1 at the first 1e19 attempt,
/// diagnosed by fan). powLe compares in u128 with early exit: total by
/// construction, r stays <= x < 2^64 before every multiply.
fn iroot(x: u64, m: u32) u64 {
    if (m == 1) return x;
    var v: u64 = @intFromFloat(std.math.pow(f64, @floatFromInt(x), 1.0 / @as(f64, @floatFromInt(m))));
    while (v < std.math.maxInt(u64) and powLe(v + 1, m, x)) v += 1;
    while (v > 1 and !powLe(v, m, x)) v -= 1;
    return v;
}
/// v^m <= x, exact.
fn powLe(v: u64, m: u32, x: u64) bool {
    var r: u128 = 1;
    var i: u32 = 0;
    while (i < m) : (i += 1) {
        r *= v;
        if (r > x) return false;
    }
    return true;
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
    .{ .x = 10_000_000_000_000_000_000, .pi = 234_057_667_276_344_607 },
};

// ---------------------------------------------------------------------------
// Kernels
// ---------------------------------------------------------------------------

// ---- double-double phase: theta = gamma * thetaCoeff mod 2pi.
// gamma*L reaches ~4e9, so plain f64 leaves |dtheta| ~ gL*2^-53 ~ 4e-7 rad —
// a coherent worst-case bound of whole UNITS on the zero sum. TwoProd plus
// dd reduction brings |dtheta| to ~3e-16 rad absolute, collapsing that
// bound to ~1e-5. Needs hi/lo zeros (.bin table) to also kill the
// gamma-rounding half; with a text table the lo half is 0 and gamma's own
// f64 rounding remains the noise floor.
const TWO_PI_HI: f64 = 6.283185307179586; // fl(2pi)
const TWO_PI_LO: f64 = 2.4492935982947064e-16; // 2pi - TWO_PI_HI

fn twoProd(a: f64, b: f64) [2]f64 {
    const p = a * b;
    return .{ p, @mulAdd(f64, a, b, -p) };
}

fn ddPhase(g_hi: f64, g_lo: f64, c_hi: f64, c_lo: f64) f64 {
    const p = twoProd(g_hi, c_hi);
    const e = p[1] + (g_hi * c_lo + g_lo * c_hi);
    const n = @round(p[0] / (2.0 * std.math.pi));
    const q = twoProd(n, TWO_PI_HI);
    // p[0] - q[0] is exact (Sterbenz); the small pieces then rank-order
    return (p[0] - q[0]) - q[1] - n * TWO_PI_LO + e;
}

fn ddOfF128(v: f128) [2]f64 {
    const hi: f64 = @floatCast(v);
    return .{ hi, @floatCast(v - hi) };
}

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
    tc_hi: f64, // dd phase coefficient L + l^2/2
    tc_lo: f64,
    lo: u64,
    hi: u64,

    fn init(x: u64, T: f64, c: f64) Gaussian {
        const xf: f64 = @floatFromInt(x);
        const L = @log(xf);
        const lambda = c / T;
        const l2 = lambda * lambda;
        const half_w = ERFC_CUT * std.math.sqrt2 * lambda;
        const ldd = ddOfF128(ln128(@floatFromInt(x)));
        return .{
            .xf = xf,
            .sx = @sqrt(xf),
            .L = L,
            .lambda = lambda,
            .l2 = l2,
            .inv_s2l = 1.0 / (std.math.sqrt2 * lambda),
            .tc_hi = ldd[0],
            .tc_lo = ldd[1] + l2 / 2.0, // both ~1e-15-scale, safe to add
            .lo = @intFromFloat(xf * @exp(-half_w) - 2.0),
            .hi = @intFromFloat(xf * @exp(half_w) + 2.0),
        };
    }

    fn phi(k: *const Gaussian, n: u64) f64 {
        // log(n/x) via log1p on the exact integer offset: no cancellation
        const u = std.math.log1p((@as(f64, @floatFromInt(n)) - k.xf) / k.xf) * k.inv_s2l;
        return 0.5 * erfc(u);
    }

    fn zeroTerm(k: *const Gaussian, g: f64, th: f64) f64 {
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
        const e = C{ .re = k.sx * damp * @cos(th), .im = k.sx * damp * @sin(th) };
        return -2.0 * e.mul(br).re;
    }

    /// Oscillation amplitude bound of this pair's term — the lever arm a
    /// phase error acts on. Feeds the dual float-error books in run().
    fn ampEnv(k: *const Gaussian, g: f64) f64 {
        const damp = @exp(k.l2 * (0.25 - g * g) / 2.0);
        return 2.4 * k.sx * damp / (g * k.L);
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
/// GRID at 2^20 keeps the table h^2 error ~1e-11 per lookup so the certified
/// window-error component stays ~1e-2 even at 1.5e9 window terms.
const GRID = 1 << 22; // eta/mu/nu cumulative-integral grid over [-1, 1]:
// h^2 table error is the leading shared floor suspect AND the leading
// cert component; 2^22 puts R_window ~ 0.05 even at 1e19

const LoganButhe = struct {
    x: u64,
    xf: f64,
    sx: f64,
    L: f64,
    c: f64,
    eps: f64,
    lam: f64,
    acorr: f64, // A_{c,eps}
    tc_hi: f64, // dd phase coefficient L
    tc_lo: f64,
    lo: u64,
    hi: u64,
    envM: f64, // sup_{u>c} |l_c(u)| u
    tailI: f64, // int |l_c(u)|/u dens(u/eps) du — cert integrand, true weight
    l1tail: f64, // int_c^inf |l_c(u)|/u du, two-sided: should equal 4 e^-c
    mu_tab: []f64, // mu(eps*y) on y-grid [-1,0], GRID+1 points
    nu_tab: []f64, // nu(eps*y) likewise
    // measured-in-init bounds for the certified radius (y-units, x1.2 safety):
    etamax: f64, // max eta            (= max |mu'|, bounds nu'' for interp)
    detamax: f64, // max |eta'|         (= max |mu''|, bounds mu interp error)
    tvpp: f64, // total variation of eta' (~ int |eta''|, bounds trapezoid drift)

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
        var etamax: f64 = 0;
        var detamax: f64 = 0;
        var tvpp: f64 = 0;
        var prev_slope: f64 = 0;
        var gi: usize = 0;
        while (gi <= GRID / 2) : (gi += 1) {
            const y = -1.0 + @as(f64, @floatFromInt(gi)) * h;
            const eta = norm * bessI0(c * @sqrt(@max(0.0, 1.0 - y * y)));
            if (gi > 0) acc += 0.5 * (prev_eta + eta) * h;
            const mu = -acc;
            if (gi > 0) nacc += 0.5 * (prev_mu + mu) * h;
            mu_tab[gi] = mu;
            nu_tab[gi] = nacc * eps; // nu carries one eps dilation factor
            etamax = @max(etamax, eta);
            if (gi > 0) {
                const slope = (eta - prev_eta) / h;
                detamax = @max(detamax, @abs(slope));
                if (gi > 1) tvpp += @abs(slope - prev_slope);
                prev_slope = slope;
            }
            prev_eta = eta;
            prev_mu = mu;
        }
        const ldd = ddOfF128(ln128(@floatFromInt(x)));
        // ---- certified-tail integral from the TRUE weight (not its
        // envelope: the envelope discards the oscillation zeros and is
        // ~pi/2 loose). Same grid and cushions as Slepian's, so the two
        // certified tails are like-for-like and extremality is testable.
        var envM: f64 = 0;
        var envTail: f64 = 0; // sup over the last octaves: the far-field envelope
        var l1: f64 = 0;
        var tailI: f64 = 0;
        const du = std.math.pi / 256.0; // 8x finer: trapezoid on the concave |sin| arches undershot the exact L1 by 0.8% at pi/32
        const dens = struct {
            fn f(t: f64) f64 {
                return @log(t / (2.0 * std.math.pi)) / (2.0 * std.math.pi);
            }
        }.f;
        var prev_g: f64 = @abs(ell(c, c)) / c;
        var prev_d: f64 = prev_g * dens(c / eps);
        var u = c + du;
        while (u <= 50.0 * c) : (u += du) {
            const eh = @abs(ell(c, u));
            envM = @max(envM, eh * u);
            if (u >= 25.0 * c) envTail = @max(envTail, eh * u);
            const g = eh / u;
            const d = g * dens(u / eps);
            l1 += 0.5 * (prev_g + g) * du;
            tailI += 0.5 * (prev_d + d) * du;
            prev_g = g;
            prev_d = d;
        }
        // envM converges to |w(c)|c from below (the sup sits at the band
        // edge); the 1.2 cushion covers the last ~2% and the sampling.
        envM *= 1.2;
        // beyond 50c: |w(u)| <= envM/u, so int_{50c}^inf |w|/u du <= envM/(50c).
        // WITHOUT this the sampled integral sits 1.3% BELOW Logan's exact
        // optimum — impossible for an upper bound, and the tell that the
        // truncation, not the grid or the trig, was the deficit.
        // the far-field envelope, sampled over [25c, 50c] where |w|u is
        // already decreasing, bounds |w(u)|u for all u > 50c — 14x tighter
        // than the near-band peak envM the remainder used before
        const rem = 1.2 * envTail / (50.0 * c);
        l1 += rem;
        tailI += 2.0 * rem * dens(50.0 * c / eps);
        tailI *= 1.2;
        return .{
            .x = x,
            .xf = xf,
            .sx = @sqrt(xf),
            .L = L,
            .c = c,
            .eps = eps,
            .lam = lam,
            .acorr = acorr,
            .tc_hi = ldd[0],
            .tc_lo = ldd[1],
            .lo = @intFromFloat(xf * @exp(-eps) - 2.0),
            .hi = @intFromFloat(xf * @exp(eps) + 2.0),
            .envM = envM,
            .tailI = tailI,
            .l1tail = 2.0 * l1,
            .mu_tab = mu_tab,
            .nu_tab = nu_tab,
            .etamax = etamax * 1.2,
            .detamax = detamax * 1.2,
            .tvpp = tvpp * 1.2,
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
        // clamp: the lo/hi cushion admits a couple of integers just outside
        // the exact support, where M = 0 — and where an unclamped y would
        // hand lookup a NEGATIVE f (@intFromFloat to usize = UB, garbage
        // index; it bit exactly when a boundary-belt integer was prime)
        const y = @max(-1.0, @min(1.0, u / k.eps));
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
    fn zeroTerm(k: *const LoganButhe, g: f64, th: f64) f64 {
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
        // e^rl = sqrt(x) e^(i g L), phase dd-reduced by the caller
        const e = C{ .re = k.sx * @cos(th), .im = k.sx * @sin(th) };
        return -2.0 * (w / k.lam) * e.mul(B).re;
    }

    fn ampEnv(k: *const LoganButhe, g: f64) f64 {
        return 2.4 * k.sx * @abs(ell(k.c, k.eps * g)) / (k.lam * g * k.L);
    }

    fn cutoff(k: *const LoganButhe, g: f64) bool {
        _ = k;
        _ = g;
        return false; // every tabulated zero contributes; tail is ~e^-c by design
    }

    fn poleCorr(k: *const LoganButhe) f64 {
        return k.acorr * k.xf / (k.L * k.L);
    }

    // ---- certified-radius components (analytic terms only; f64 rounding
    // and the Thm 4.1 Theta-constants are EXCLUDED and said so in output —
    // a bracket must state what it covers).

    /// Bound on the dropped zeros gamma > T. Per conjugate pair:
    /// |term| <= 2 sqrt(x) * 1.3/(lam*g*L) * env(eps*g), where beyond the
    /// support crossover |l_c(t)| <= (c/sinh c)*min(1, 1/sqrt(t^2-c^2)) —
    /// the sinh-normalisation IS the e^-c tail. Zero counting via the
    /// Rosser bound |N(t) - M(t)| <= 0.137 ln t + 0.443 lnln t + 4.35.
    fn certTail(k: *const LoganButhe, T: f64) f64 {
        const cf = 2.0 * k.sx * 1.3 / (k.lam * k.L);
        const qT = 0.137 * @log(T) + 0.443 * @log(@log(T)) + 4.35;
        return cf * k.tailI + 3.0 * qT * cf * (k.envM / k.c) / T;
    }

    /// Coarse bound on the K=13 Ei asymptotic truncation: first omitted
    /// term at the lowest zero (largest relative error), sector-cushioned
    /// by 4|rho|, times the zero count. ~1e-10 territory; kept honest.
    fn certSeries(k: *const LoganButhe, g1: f64, nzeros: f64) f64 {
        const r1 = @sqrt(0.25 + g1 * g1) * k.L;
        var fac: f64 = 1; // 16!/2! = product 3..16
        var j: f64 = 3;
        while (j <= 16) : (j += 1) fac *= j;
        const rel = fac / std.math.pow(f64, r1, 14.0) * 4.0 * g1;
        return nzeros * 2.0 * k.sx / (k.lam * g1 * k.L) * rel;
    }

    /// Buthe Thm 4.1's Theta(35 eps) explicit-formula remainder, taken in
    /// its decomposed form from Lemmas 4.3 (3.2 A/log x, the pole/li
    /// comparison), 4.4 (440 eps^4/c^2, prime powers via Brun-Titchmarsh)
    /// and 4.5 (34.9 eps, the w_infinity comparison — itself six explicit
    /// pieces). Hypotheses: x > 30000, c >= 1, 0 < eps <= 0.01; asserted.
    fn certTheta(k: *const LoganButhe) f64 {
        std.debug.assert(k.x > 30000 and k.c >= 1.0 and k.eps <= 0.01);
        const e2 = k.eps * k.eps;
        return 34.9 * k.eps + 3.2 * k.acorr / k.L + 440.0 * e2 * e2 / (k.c * k.c);
    }

    /// Per-window-term bound on the mu/nu table error (trapezoid drift
    /// h^2/12 * int|eta''| plus linear-interp h^2/8 * max|mu''|); the nu
    /// contribution carries the eps dilation and is negligible but kept.
    fn certWindowPerTerm(k: *const LoganButhe) f64 {
        const h = 2.0 / @as(f64, GRID);
        const dmu = h * h * (k.tvpp / 12.0 + k.detamax / 8.0);
        const dnu = k.eps * h * h * (2.0 * k.etamax / 12.0 + k.etamax / 8.0);
        return (dmu * (1.0 + 0.5 * k.eps) + 0.5 * dnu) / k.lam;
    }
};

/// Slepian/prolate kernel — the L2-extremal experiment. The bump is the
/// zeroth prolate spheroidal wave function psi_0 on [-1,1] with bandwidth
/// parameter c (Bouwkamp Legendre eigen-solve of
/// ((1-y^2)u')' + (chi - c^2 y^2)u = 0, ground state, even chain).
/// psi_0 maximizes the in-band fraction lambda_0 of its transform's L2
/// mass; leakage 1-lambda_0 ~ e^{-2c} vs Logan's L1-extremal e^{-c} tail.
/// Logan is the prover's kernel (worst-case signs on dropped zeros);
/// this is the empiricist's (RMS signs). Same window support, same
/// mu/nu/M dressing — Buthe's A-correction is the bump's second moment
/// (verified numerically to 2e-16 against his closed form for Logan),
/// so the Ei machinery is moment-generic and psi_0 inherits it.
///
/// Zero-side weight via the self-transform property (no Bessel sums):
///   int psi_0(y) cos(c u y) dy = mu_0 psi_0(u)  =>  eta_hat(t) = psi_0(t/c)/psi_0(0)
/// tabulated as ln psi_0 on [0,1] (positive ground state), exp on lookup.
const Slepian = struct {
    const KMAXE = 300; // even-chain length cap (k up to 598)

    x: u64,
    xf: f64,
    sx: f64,
    L: f64,
    c: f64,
    eps: f64,
    lam: f64,
    acorr: f64, // eps^2 m2(psi_0)/2
    tc_hi: f64,
    tc_lo: f64,
    lo: u64,
    hi: u64,
    mu_tab: []f64,
    nu_tab: []f64,
    lnw_tab: []f64, // ln psi_0 on u-grid [0,1], GRID/2+1 points
    etamax: f64,
    detamax: f64,
    tvpp: f64,
    chi0: f64, // prolate ODE eigenvalue (diagnostic)
    envM: f64, // sup_{u>c} |eta_hat(u)| u — the L1 envelope constant
    l1tail: f64, // int_c^inf |eta_hat(u)|/u du (vs Logan's proven 4 e^-c)
    tailI: f64, // int_c^inf |eta_hat(u)|/u density(u/eps) du — cert integrand

    fn init(gpa: std.mem.Allocator, x: u64, T: f64, c: f64) !Slepian {
        const xf: f64 = @floatFromInt(x);
        const L = @log(xf);
        const eps = c / T;

        // ---- Bouwkamp: symmetric tridiagonal over even normalized-Legendre
        // coefficients; ground state by Sturm bisection + inverse iteration
        const ne: usize = @min(KMAXE, @as(usize, @intFromFloat(2.0 * c + 40.0)));
        var diag: [KMAXE]f64 = undefined;
        var off: [KMAXE]f64 = undefined; // off[i] couples entries i, i+1
        for (0..ne) |i| {
            const k: f64 = @floatFromInt(2 * i);
            const bkk = if (i == 0) 1.0 / 3.0 else (2.0 * k * (k + 1.0) - 1.0) / ((2.0 * k + 3.0) * (2.0 * k - 1.0));
            diag[i] = k * (k + 1.0) + c * c * bkk;
            if (i + 1 < ne) {
                off[i] = c * c * (k + 2.0) * (k + 1.0) / ((2.0 * k + 3.0) * @sqrt((2.0 * k + 1.0) * (2.0 * k + 5.0)));
            }
        }
        var blo: f64 = 0.0;
        var bhi: f64 = diag[0];
        for (0..ne) |i| bhi = @max(bhi, diag[i] + 2.0 * (if (i + 1 < ne) off[i] else 0.0));
        var it: usize = 0;
        while (it < 200) : (it += 1) { // Sturm count of eigenvalues < mid
            const mid = 0.5 * (blo + bhi);
            var cnt: usize = 0;
            var d = diag[0] - mid;
            if (d < 0) cnt += 1;
            for (1..ne) |i| {
                const dd = if (d != 0) d else 1e-300;
                d = diag[i] - mid - off[i - 1] * off[i - 1] / dd;
                if (d < 0) cnt += 1;
            }
            if (cnt >= 1) bhi = mid else blo = mid;
        }
        const chi0 = 0.5 * (blo + bhi);
        var a: [KMAXE]f64 = undefined;
        var v: [KMAXE]f64 = undefined;
        for (0..ne) |i| a[i] = 1.0 / @sqrt(@as(f64, @floatFromInt(ne)));
        var round: usize = 0;
        while (round < 4) : (round += 1) { // (A - chi I) w = a, Thomas
            const sh = chi0 - 1e-9;
            var cp: [KMAXE]f64 = undefined;
            var dp: [KMAXE]f64 = undefined;
            const b0 = diag[0] - sh;
            cp[0] = off[0] / b0;
            dp[0] = a[0] / b0;
            for (1..ne) |i| {
                const m = (diag[i] - sh) - off[i - 1] * cp[i - 1];
                cp[i] = if (i + 1 < ne) off[i] / m else 0.0;
                dp[i] = (a[i] - off[i - 1] * dp[i - 1]) / m;
            }
            v[ne - 1] = dp[ne - 1];
            var ii = ne - 1;
            while (ii > 0) : (ii -= 1) v[ii - 1] = dp[ii - 1] - cp[ii - 1] * v[ii];
            var nrm: f64 = 0;
            for (0..ne) |i| nrm += v[i] * v[i];
            nrm = @sqrt(nrm);
            for (0..ne) |i| a[i] = v[i] / nrm;
        }
        if (a[0] < 0) for (0..ne) |i| {
            a[i] = -a[i];
        };
        const n0 = a[0] * @sqrt(2.0); // int psi_0 (only P_0 integrates)

        // ---- tables: eta = psi_0/n0 on y in [-1,0] (even mirror), same
        // trapezoid machinery as Logan, plus ln psi_0 on u = -y in [0,1],
        // second moment m2 and lam = int eta cosh(eps y / 2) * 2
        const h = 2.0 / @as(f64, GRID);
        const mu_tab = try gpa.alloc(f64, GRID / 2 + 1);
        const nu_tab = try gpa.alloc(f64, GRID / 2 + 1);
        const lnw_tab = try gpa.alloc(f64, GRID / 2 + 1);
        var acc: f64 = 0;
        var nacc: f64 = 0;
        var m2acc: f64 = 0;
        var lamacc: f64 = 0;
        var prev_eta: f64 = 0;
        var prev_mu: f64 = 0;
        var prev_y2e: f64 = 0;
        var prev_che: f64 = 0;
        var etamax: f64 = 0;
        var detamax: f64 = 0;
        var tvpp: f64 = 0;
        var prev_slope: f64 = 0;
        var gi: usize = 0;
        while (gi <= GRID / 2) : (gi += 1) {
            const y = -1.0 + @as(f64, @floatFromInt(gi)) * h;
            // psi_0(y) = sum a_i sqrt(2i + 1/2) P_{2i}(y), upward recurrence
            var p0: f64 = 1.0;
            var p1: f64 = y;
            var psi: f64 = a[0] * @sqrt(0.5);
            var k: usize = 1;
            while (k < 2 * ne) : (k += 1) {
                const kf: f64 = @floatFromInt(k);
                const p2 = ((2.0 * kf + 1.0) * y * p1 - kf * p0) / (kf + 1.0);
                if ((k + 1) % 2 == 0 and (k + 1) / 2 < ne)
                    psi += a[(k + 1) / 2] * @sqrt(@as(f64, @floatFromInt(k + 1)) + 0.5) * p2;
                p0 = p1;
                p1 = p2;
            }
            const eta = @max(psi, 1e-300) / n0;
            lnw_tab[GRID / 2 - gi] = @log(@max(psi, 1e-300));
            const y2e = y * y * eta;
            const che = eta * std.math.cosh(eps * y / 2.0);
            if (gi > 0) {
                acc += 0.5 * (prev_eta + eta) * h;
                m2acc += 0.5 * (prev_y2e + y2e) * h;
                lamacc += 0.5 * (prev_che + che) * h;
            }
            const mu = -acc;
            if (gi > 0) nacc += 0.5 * (prev_mu + mu) * h;
            mu_tab[gi] = mu;
            nu_tab[gi] = nacc * eps;
            etamax = @max(etamax, eta);
            if (gi > 0) {
                const slope = (eta - prev_eta) / h;
                detamax = @max(detamax, @abs(slope));
                if (gi > 1) tvpp += @abs(slope - prev_slope);
                prev_slope = slope;
            }
            prev_eta = eta;
            prev_mu = mu;
            prev_y2e = y2e;
            prev_che = che;
        }
        const m2 = 2.0 * m2acc;
        const lam = 2.0 * lamacc;
        const ldd = ddOfF128(ln128(@floatFromInt(x)));
        // ---- out-of-band envelope, the L1 way. Beyond the band the
        // self-transform identity no longer holds, so use the exact
        // Legendre-Bessel continuation. Sampled on a pi/32 grid (the
        // oscillation period is 2pi) out to 50c: a SAMPLED sup, not a
        // proven one — flagged in the output. This replaces a Parseval +
        // Cauchy-Schwarz bound whose 1-lambda_0 ~ e^{-2c} input is below
        // f64's ability to measure at all, and which therefore carried a
        // 1000x ignorance cushion (= 32x after the square root).
        const kmax: usize = 2 * (ne - 1);
        const jb = try gpa.alloc(f64, kmax + 1);
        defer gpa.free(jb);
        var envM: f64 = 0;
        var envTail: f64 = 0; // sup over the last octaves: the far-field envelope
        var l1: f64 = 0;
        var tailI: f64 = 0; // int |eta^(u)|/u * density(u/eps) du — the cert integrand
        const du = std.math.pi / 256.0; // 8x finer: trapezoid on the concave |sin| arches undershot the exact L1 by 0.8% at pi/32
        const dens = struct {
            fn f(t: f64) f64 {
                return @log(t / (2.0 * std.math.pi)) / (2.0 * std.math.pi);
            }
        }.f;
        var prev_g: f64 = @abs(etaHat(a[0..ne], n0, c, jb)) / c;
        var prev_d: f64 = prev_g * dens(c / eps);
        var u = c + du;
        while (u <= 50.0 * c) : (u += du) {
            const eh = @abs(etaHat(a[0..ne], n0, u, jb));
            envM = @max(envM, eh * u);
            if (u >= 25.0 * c) envTail = @max(envTail, eh * u);
            const g = eh / u;
            const d = g * dens(u / eps);
            l1 += 0.5 * (prev_g + g) * du;
            tailI += 0.5 * (prev_d + d) * du;
            prev_g = g;
            prev_d = d;
        }
        envM *= 1.2;
        // see LoganButhe.init
        // the far-field envelope, sampled over [25c, 50c] where |w|u is
        // already decreasing, bounds |w(u)|u for all u > 50c — 14x tighter
        // than the near-band peak envM the remainder used before
        const rem = 1.2 * envTail / (50.0 * c);
        l1 += rem;
        tailI += 2.0 * rem * dens(50.0 * c / eps);
        tailI *= 1.2;
        return .{
            .x = x,
            .xf = xf,
            .sx = @sqrt(xf),
            .L = L,
            .c = c,
            .eps = eps,
            .lam = lam,
            .acorr = eps * eps * m2 / 2.0,
            .tc_hi = ldd[0],
            .tc_lo = ldd[1],
            .lo = @intFromFloat(xf * @exp(-eps) - 2.0),
            .hi = @intFromFloat(xf * @exp(eps) + 2.0),
            .mu_tab = mu_tab,
            .nu_tab = nu_tab,
            .lnw_tab = lnw_tab,
            .etamax = etamax * 1.2,
            .detamax = detamax * 1.2,
            .tvpp = tvpp * 1.2,
            .chi0 = chi0,
            .envM = envM,
            .l1tail = 2.0 * l1, // both signs of u (Logan's constant is two-sided)
            .tailI = tailI,
        };
    }

    /// j_k(t), k = 0..kmax, by downward (Miller) recurrence. Normalized on
    /// whichever of the analytic j_0, j_1 is larger — they are pi/2 out of
    /// phase, so this never divides by a near-zero.
    fn sphBessel(kmax: usize, t: f64, out: []f64) void {
        @memset(out, 0);
        var m: usize = kmax + 25 + @as(usize, @intFromFloat(t));
        var jp: f64 = 0;
        var j: f64 = 1e-300;
        while (m > 0) : (m -= 1) {
            const kf: f64 = @floatFromInt(m);
            const jm = (2.0 * kf + 1.0) / t * j - jp;
            jp = j;
            j = jm;
            if (m - 1 <= kmax) out[m - 1] = j;
            if (@abs(j) > 1e200) {
                jp *= 1e-200;
                j *= 1e-200;
                for (out) |*v| v.* *= 1e-200;
            }
        }
        const a0 = @sin(t) / t;
        const a1 = @sin(t) / (t * t) - @cos(t) / t;
        const scale = if (@abs(a0) >= @abs(a1)) a0 / out[0] else a1 / out[1];
        for (out) |*v| v.* *= scale;
    }

    /// eta_hat(u) for ANY real u, from the exact continuation
    ///   eta_hat(u) = (1/N0) sum_i a_i sqrt(2i+1/2) 2 (-1)^i j_{2i}(u),
    /// since int_{-1}^{1} P_k(y) cos(uy) dy = 2 (-1)^{k/2} j_k(u) for even k.
    /// Validated against direct quadrature in prolate-proto.py (1e-16).
    fn etaHat(a: []const f64, n0: f64, u: f64, jb: []f64) f64 {
        sphBessel(jb.len - 1, u, jb);
        var s: f64 = 0;
        for (a, 0..) |ai, i| {
            const k = 2 * i;
            const sgn: f64 = if (i % 2 == 0) 1.0 else -1.0;
            s += ai * @sqrt(@as(f64, @floatFromInt(k)) + 0.5) * 2.0 * sgn * jb[k];
        }
        return s / n0;
    }

    /// eta_hat(t) = psi_0(t/c)/psi_0(0) by self-transform; ln-interp.
    fn weight(k: *const Slepian, t: f64) f64 {
        const u = @max(0.0, @min(1.0, t / k.c));
        const f = u * @as(f64, GRID) / 2.0;
        const idx: usize = @intFromFloat(@min(f, @as(f64, GRID / 2) - 1e-9));
        const fr = f - @as(f64, @floatFromInt(idx));
        const lnw = k.lnw_tab[idx] * (1.0 - fr) + k.lnw_tab[idx + 1] * fr;
        return @exp(lnw - k.lnw_tab[0]);
    }

    fn bigM(k: *const Slepian, n: u64) f64 {
        const u = std.math.log1p((@as(f64, @floatFromInt(n)) - k.xf) / k.xf);
        const y = @max(-1.0, @min(1.0, u / k.eps)); // clamp: see LoganButhe
        var mu: f64 = undefined;
        var nu: f64 = undefined;
        if (n <= k.x) { // integer cut, not sign(y): see LoganButhe
            mu = LoganButhe.lookup(k.mu_tab, y);
            nu = LoganButhe.lookup(k.nu_tab, y);
        } else {
            mu = -LoganButhe.lookup(k.mu_tab, -y);
            nu = LoganButhe.lookup(k.nu_tab, -y);
        }
        const logn = k.L + u;
        return (mu + (1.0 / logn - 0.5) * (mu * u - nu)) / k.lam;
    }

    fn phi(k: *const Slepian, n: u64) f64 {
        const chi: f64 = if (n <= k.x) 1.0 else 0.0;
        return chi + k.bigM(n);
    }

    fn zeroTerm(k: *const Slepian, g: f64, th: f64) f64 {
        const w = k.weight(k.eps * g);
        if (w == 0.0) return 0.0;
        const rl = C{ .re = 0.5 * k.L, .im = g * k.L };
        const irl = rl.inv();
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
        const rho = C{ .re = 0.5, .im = g };
        const irl2 = irl.mul(irl);
        var b2 = rho.mul(irl2).mul(s2);
        const b3 = rho.mul(rho).mul(irl2).mul(irl).mul(s3);
        b2.re -= 2.0 * b3.re;
        b2.im -= 2.0 * b3.im;
        var B = s1.mul(irl);
        B.re += k.acorr * b2.re;
        B.im += k.acorr * b2.im;
        const e = C{ .re = k.sx * @cos(th), .im = k.sx * @sin(th) };
        return -2.0 * (w / k.lam) * e.mul(B).re;
    }

    fn ampEnv(k: *const Slepian, g: f64) f64 {
        return 2.4 * k.sx * k.weight(k.eps * g) / (k.lam * g * k.L);
    }

    fn cutoff(k: *const Slepian, g: f64) bool {
        _ = k;
        _ = g;
        return false;
    }

    fn poleCorr(k: *const Slepian) f64 {
        return k.acorr * k.xf / (k.L * k.L);
    }

    /// Dropped zeros gamma > T, bounded the L1 way — structurally identical
    /// to Logan's, with |eta_hat(u)| <= envM/u replacing his
    /// (c/sinh c)/sqrt(u^2-c^2). Like-for-like, so the two certified tails
    /// differ only through envM vs c/sinh c, and Logan's extremality
    /// (he minimizes int_{|t|>c}|f(t)/t| dt) says he must win: if this
    /// ever prints lower than Logan's, the implementation is wrong.
    fn certTail(k: *const Slepian, T: f64) f64 {
        // substitute u = eps t: sum_{gamma>T} amp |w| <= cf * int_c^inf
        // |eta^(u)|/u dens(u/eps) du, precomputed as tailI in init from the
        // TRUE weight (not its envelope — that alone was 4.5x loose).
        const cf = 2.0 * k.sx * 1.3 / (k.lam * k.L);
        const qT = 0.137 * @log(T) + 0.443 * @log(@log(T)) + 4.35;
        return cf * k.tailI + 3.0 * qT * cf * (k.envM / k.c) / T;
    }

    fn certSeries(k: *const Slepian, g1: f64, nzeros: f64) f64 {
        const r1 = @sqrt(0.25 + g1 * g1) * k.L;
        var fac: f64 = 1;
        var j: f64 = 3;
        while (j <= 16) : (j += 1) fac *= j;
        const rel = fac / std.math.pow(f64, r1, 14.0) * 4.0 * g1;
        return nzeros * 2.0 * k.sx / (k.lam * g1 * k.L) * rel;
    }

    /// Same FORM as Buthe Thm 4.1's Theta(35 eps), but his constants are
    /// derived from Logan-specific bounds (Lemma 3.3 bounds
    /// |phi_inf,c,eps - f1| using ell_c). The Weil-Barner argument carries
    /// over to any even bump, so the O(eps) shape is right and the
    /// constants are plausibly close — but they are NOT proven for psi_0.
    /// Flagged as borrowed in the output; do not cite as certified.
    const THETA_BORROWED = true;
    fn certTheta(k: *const Slepian) f64 {
        std.debug.assert(k.x > 30000 and k.c >= 1.0 and k.eps <= 0.01);
        const e2 = k.eps * k.eps;
        return 34.9 * k.eps + 3.2 * k.acorr / k.L + 440.0 * e2 * e2 / (k.c * k.c);
    }

    fn certWindowPerTerm(k: *const Slepian) f64 {
        const h = 2.0 / @as(f64, GRID);
        const dmu = h * h * (k.tvpp / 12.0 + k.detamax / 8.0);
        const dnu = k.eps * h * h * (2.0 * k.etamax / 12.0 + k.etamax / 8.0);
        return (dmu * (1.0 + 0.5 * k.eps) + 0.5 * dnu) / k.lam;
    }
};

// ---------------------------------------------------------------------------
// Pipeline, generic over the kernel
// ---------------------------------------------------------------------------

fn run(comptime K: type, kern: *const K, x: u64, zeros: []const f64, zlo: []const f64, nt: usize, gpa: std.mem.Allocator, t_start: u64, t_load: u64) !void {
    const ZCtx = struct {
        zeros: []const f64,
        zlo: []const f64,
        nt: usize,
        kern: *const K,
        sums: []f64,
        // dual float-error books: worst-case (coherent) vs stochastic-model
        // (RMS, mean-zero independent). The ratio is the measured "sqrt(N)
        // of sign-mixing" this run is implicitly relying on.
        amp1: []f64, // sum |amp|
        amp2: []f64, // sum amp^2
        abst: []f64, // sum |term|
        fn work(ctx: *@This(), wid: usize) void {
            const n = ctx.zeros.len;
            const chunk = (n + ctx.nt - 1) / ctx.nt;
            const a = @min(wid * chunk, n);
            const b = @min(a + chunk, n);
            var k = Kahan{};
            var a1: f64 = 0;
            var a2: f64 = 0;
            var at: f64 = 0;
            for (ctx.zeros[a..b], ctx.zlo[a..b]) |g, gl| {
                if (ctx.kern.cutoff(g)) break; // gammas ascend within a range
                const th = ddPhase(g, gl, ctx.kern.tc_hi, ctx.kern.tc_lo);
                const t = ctx.kern.zeroTerm(g, th);
                k.add(t);
                const amp = ctx.kern.ampEnv(g);
                a1 += amp;
                a2 += amp * amp;
                at += @abs(t);
            }
            ctx.sums[wid] = k.val();
            ctx.amp1[wid] = a1;
            ctx.amp2[wid] = a2;
            ctx.abst[wid] = at;
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
    const zamp1 = try gpa.alloc(f64, nt);
    defer gpa.free(zamp1);
    const zamp2 = try gpa.alloc(f64, nt);
    defer gpa.free(zamp2);
    const zabst = try gpa.alloc(f64, nt);
    defer gpa.free(zabst);
    var zctx = ZCtx{ .zeros = zeros, .zlo = zlo, .nt = nt, .kern = kern, .sums = zsums, .amp1 = zamp1, .amp2 = zamp2, .abst = zabst };
    {
        const threads = try gpa.alloc(std.Thread, nt);
        defer gpa.free(threads);
        for (threads, 0..) |*t, i| t.* = try std.Thread.spawn(.{}, ZCtx.work, .{ &zctx, i });
        for (threads) |t| t.join();
    }
    var zk = Kahan{};
    var s_amp1: f64 = 0;
    var s_amp2: f64 = 0;
    var s_abst: f64 = 0;
    for (zsums) |s| zk.add(s);
    for (zamp1) |s| s_amp1 += s;
    for (zamp2) |s| s_amp2 += s;
    for (zabst) |s| s_abst += s;
    const zerosum = zk.val();
    // dual books: per-term error = amp * dtheta + |term| * ops*u, plus
    // Neumaier 2u*sum|t|. Worst case adds moduli; the stochastic model
    // (mean-zero independent, uniform-ish) adds variances. dtheta ~ 6u
    // covers ddPhase's final adds + one ulp each for sin/cos (with a
    // text table gamma's own rounding is NOT included — .bin is the
    // instrument here). MODEL, not certificate, and labeled so.
    const uu = 1.1102230246251565e-16;
    const dtheta = 6.0 * uu;
    const ops = 42.0 * uu;
    const f_worst = dtheta * s_amp1 + ops * s_abst + 2.0 * uu * s_abst;
    const f_rms = @sqrt(dtheta * dtheta * s_amp2 / 3.0 + ops * ops * s_abst * s_abst / @max(1.0, @as(f64, @floatFromInt(zeros.len))));
    std.debug.print("float(zerosum): worst-case {e:.2} | stochastic-model rms {e:.2} | mixing-leverage {d:.0}x\n", .{ f_worst, f_rms, f_worst / @max(f_rms, 1e-300) });
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
    var npow: u64 = 0;
    for (base) |p| {
        var pk: u64 = 1;
        var m: f64 = 0;
        while (pk <= kern.hi / p) {
            pk *= p;
            m += 1;
            if (pk <= kern.lo) continue;
            const chi: f64 = if (pk <= x) 1.0 else 0.0;
            wk.add((chi - kern.phi(pk)) / m);
            npow += 1;
        }
    }
    const wcorr = wk.val();
    const pistar_sharp: f128 = pistar_smooth + @as(f128, wcorr);
    const t_win = nowNs();
    // certified-radius report (kernels that can price their own errors)
    if (comptime @hasDecl(K, "certTail")) {
        const T = zeros[zeros.len - 1];
        const r_tail = kern.certTail(T);
        const r_series = kern.certSeries(zeros[0], @floatFromInt(zeros.len));
        const r_win = kern.certWindowPerTerm() * @as(f64, @floatFromInt(wprimes + npow));
        const r_theta = kern.certTheta();
        const borrowed = comptime @hasDecl(K, "THETA_BORROWED");
        std.debug.print("cert: R_tail {e:.2}  R_series {e:.2}  R_window {e:.2}  R_theta {e:.2}{s}  => R_analytic {e:.2}\n", .{ r_tail, r_series, r_win, r_theta, if (borrowed) " (BORROWED)" else "", r_tail + r_series + r_win + r_theta });
        // the kernel-independent systematic measured in the 2026-08-04
        // sweeps — real, reproducible, mechanism unknown, NOT in the radius
        const floor_est = 4e-4 * @sqrt(@as(f64, @floatFromInt(x)) / 1e15) * @exp(-0.8 * (kern.c - 15.0));
        std.debug.print("cert: shared-floor estimate {e:.2} (EMPIRICAL, unpriced — mechanism on worklist)\n", .{floor_est});
        if (borrowed) std.debug.print("cert: NOTE R_theta constants are Buthe's, proven for the Logan bump only\n", .{});
        std.debug.print("cert: EXCLUDED: f64 rounding (see float() line; dd/ball rung to fold in)\n", .{});
    }
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
    var zlos = try std.ArrayList(f64).initCapacity(gpa, 1 << 21);
    defer zlos.deinit(gpa);
    if (data.len >= 16 and std.mem.eql(u8, data[0..8], "ZROSDD01")) {
        // lmfdb2bin.py hi/lo pairs (little-endian, matches x86); lo feeds
        // the dd phase so gamma's own rounding stops being the noise floor
        const n = std.mem.bytesToValue(u64, data[8..16]);
        std.debug.assert(data.len >= 16 + 16 * n);
        try zeros.ensureTotalCapacityPrecise(gpa, n);
        try zlos.ensureTotalCapacityPrecise(gpa, n);
        var i: usize = 0;
        while (i < n) : (i += 1) {
            zeros.appendAssumeCapacity(std.mem.bytesToValue(f64, data[16 + 16 * i ..][0..8]));
            zlos.appendAssumeCapacity(std.mem.bytesToValue(f64, data[24 + 16 * i ..][0..8]));
        }
    } else {
        var lines = std.mem.tokenizeScalar(u8, data, '\n');
        while (lines.next()) |line| {
            const s = std.mem.trim(u8, line, " \t\r");
            if (s.len == 0) continue;
            try zeros.append(gpa, try std.fmt.parseFloat(f64, s));
        }
        try zlos.appendNTimes(gpa, 0.0, zeros.items.len);
    }
    gpa.free(data);
    const t_load = nowNs();

    const T = zeros.items[zeros.items.len - 1];
    std.debug.print("x = {d}  zeros = {d} (T = {d:.3})  kernel = {s} (c = {d})  threads = {d}\n", .{ x, zeros.items.len, T, kname, cpar, nt });

    if (std.mem.eql(u8, kname, "gaussian")) {
        const kern = Gaussian.init(x, T, cpar);
        try run(Gaussian, &kern, x, zeros.items, zlos.items, nt, gpa, t_start, t_load);
    } else if (std.mem.eql(u8, kname, "logan")) {
        // empirical tuning rule (2026-08-04 c-sweep, envelope < 0.05):
        // tail e^{-1.25c} governs to ~1e16, the shared floor (0.625 ln x
        // slope) above; old 0.5 ln x + 9 was ~2x conservative
        const lx = @log(@as(f64, @floatFromInt(x)));
        const cl = if (cset) cpar else @max(8.0, @max(0.4 * lx - 2.9, 0.625 * lx - 11.3));
        const kern = try LoganButhe.init(gpa, x, T, cl);
        const l1_opt = 2.0 * @log((1.0 + @exp(-cl)) / (1.0 - @exp(-cl)));
        std.debug.print("logan: c = {d:.2}  eps = {e:.4}  lam = {d:.6}  L1tail = {e:.3} vs proven optimum {e:.3} ({d:.3}x)\n", .{ cl, kern.eps, kern.lam, kern.l1tail, l1_opt, kern.l1tail / l1_opt });
        try run(LoganButhe, &kern, x, zeros.items, zlos.items, nt, gpa, t_start, t_load);
    } else if (std.mem.eql(u8, kname, "slepian")) {
        // e^{-2.3c} tail never binds above 1e14 — the shared floor sets c
        const lx = @log(@as(f64, @floatFromInt(x)));
        const cl = if (cset) cpar else @max(8.0, 0.625 * lx - 11.3);
        const kern = try Slepian.init(gpa, x, T, cl);
        std.debug.print("slepian: c = {d:.2}  eps = {e:.4}  lam = {d:.6}  chi0 = {d:.4}  m2*2/eps^2 = {e:.4}\n", .{ cl, kern.eps, kern.lam, kern.chi0, kern.acorr * 2.0 / (kern.eps * kern.eps) });
        // like-for-like against Logan's proven optimum for the same functional
        const l1_logan = 2.0 * @log((1.0 + @exp(-cl)) / (1.0 - @exp(-cl)));
        std.debug.print("slepian: envM = {e:.3} (sampled sup |eta^(u)|u; logan c/sinh c = {e:.3})  L1tail = {e:.3} vs logan {e:.3} ({d:.2}x)\n", .{ kern.envM, cl / std.math.sinh(cl), kern.l1tail, l1_logan, kern.l1tail / l1_logan });
        try run(Slepian, &kern, x, zeros.items, zlos.items, nt, gpa, t_start, t_load);
    } else {
        std.debug.print("unknown kernel '{s}' (have: gaussian, logan, slepian)\n", .{kname});
        std.process.exit(2);
    }
}
