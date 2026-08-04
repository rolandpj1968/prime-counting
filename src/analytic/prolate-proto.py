#!/usr/bin/env python3
"""Prototype: prolate psi_0 via Bouwkamp Legendre eigen-solve, pure python.

Validates the formulas before the Zig port:
 1. ground-state solve of the prolate ODE ((1-y^2)u')' + (chi - c^2 y^2)u = 0
    in the normalized-Legendre basis (even chain, symmetric tridiagonal)
 2. ODE residual check on a grid
 3. self-transform check: int psi0(y) cos(c x y) dy = mu0 psi0(x), ratio
    constancy over x in [0,1]; lambda0 = c mu0^2 / (2 pi)
 4. m2 identity: int y^2 eta_logan dy =? coth(c)/c - 1/c^2  (validates that
    Buthe's A-correction is the second moment => moment-generic dressing)
 5. eta-hat via Legendre/spherical-Bessel vs direct quadrature
"""
import math

def solve_prolate(c, kmax=200):
    # even chain k = 0,2,...,kmax; normalized Legendre Pbar_k = sqrt(k+1/2) P_k
    ks = list(range(0, kmax + 1, 2))
    diag = []
    off = []  # off[i] couples ks[i] and ks[i+1]
    for i, k in enumerate(ks):
        bkk = (2.0 * k * (k + 1) - 1.0) / ((2 * k + 3) * (2 * k - 1)) if k > 0 else 1.0 / 3.0
        diag.append(k * (k + 1) + c * c * bkk)
        if i + 1 < len(ks):
            bk2 = (k + 2.0) * (k + 1.0) / ((2 * k + 3) * math.sqrt((2 * k + 1) * (2 * k + 5)))
            off.append(c * c * bk2)
    n = len(ks)

    def sturm(x):
        # number of eigenvalues < x
        cnt = 0
        d = diag[0] - x
        if d < 0: cnt += 1
        for i in range(1, n):
            d = diag[i] - x - (off[i-1] * off[i-1]) / (d if d != 0 else 1e-300)
            if d < 0: cnt += 1
        return cnt

    lo, hi = 0.0, max(diag) + 2 * max(off)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if sturm(mid) >= 1: hi = mid
        else: lo = mid
    chi = 0.5 * (lo + hi)

    # inverse iteration for the ground state
    v = [1.0 / math.sqrt(n)] * n
    for _ in range(4):
        # solve (A - (chi - 1e-10) I) w = v, Thomas
        sh = chi - 1e-9
        cp = [0.0] * n; dp = [0.0] * n
        b0 = diag[0] - sh
        cp[0] = off[0] / b0; dp[0] = v[0] / b0
        for i in range(1, n):
            m = (diag[i] - sh) - off[i-1] * cp[i-1]
            cp[i] = (off[i] / m) if i + 1 < n else 0.0
            dp[i] = (v[i] - off[i-1] * dp[i-1]) / m
        w = [0.0] * n
        w[n-1] = dp[n-1]
        for i in range(n - 2, -1, -1):
            w[i] = dp[i] - cp[i] * w[i+1]
        nrm = math.sqrt(sum(x * x for x in w))
        v = [x / nrm for x in w]
    if v[0] < 0: v = [-x for x in v]
    return chi, ks, v

def legendre_eval(ks, a, y):
    # sum a_i Pbar_{ks[i]}(y) by upward recurrence over ALL degrees
    kmax = ks[-1]
    p0, p1 = 1.0, y
    s = 0.0
    ai = {k: x for k, x in zip(ks, a)}
    for k in range(0, kmax + 1):
        pk = p0 if k == 0 else p1
        if k in ai:
            s += ai[k] * math.sqrt(k + 0.5) * pk
        if k >= 1:
            p2 = ((2 * k + 1) * y * p1 - k * p0) / (k + 1)
            p0, p1 = p1, p2
        else:
            p1 = y
    return s

def sph_jn(nmax, t):
    # downward (Miller) recurrence, normalized by j0
    if t == 0.0:
        return [1.0 if k == 0 else 0.0 for k in range(nmax + 1)]
    m = nmax + 20 + int(t)
    jp, j = 0.0, 1e-30
    out = [0.0] * (nmax + 1)
    for k in range(m, 0, -1):
        jm = (2 * k + 1) / t * j - jp
        jp, j = j, jm
        if k - 1 <= nmax:
            out[k - 1] = j
        if abs(j) > 1e250:
            jp *= 1e-250; j *= 1e-250
            for i in range(len(out)):
                out[i] *= 1e-250
    scale = (math.sin(t) / t) / out[0]
    return [x * scale for x in out]

C = 28.57  # the c at x = 1e17
chi, ks, a = solve_prolate(C)
print(f"c = {C}: chi0 = {chi:.6f}, {len(ks)} even Legendre coeffs, a0 = {a[0]:.6f}")
print(f"tail coeffs: a[k=100] = {a[50]:.3e}  a[k=160] = {a[80]:.3e}")

# 2. ODE residual at sample points (finite differences)
h = 1e-4
worst = 0.0
for y in [0.0, 0.3, 0.6, 0.9, 0.99]:
    u = lambda t: legendre_eval(ks, a, t)
    upp = (u(y + h) - 2 * u(y) + u(y - h)) / (h * h)
    up = (u(y + h) - u(y - h)) / (2 * h)
    res = (1 - y * y) * upp - 2 * y * up + (chi - C * C * y * y) * u(y)
    worst = max(worst, abs(res) / (abs(u(0.0)) * C * C))
print(f"ODE relative residual (fd h={h}): worst {worst:.2e}")

# 3. self-transform: F(x) = int psi0(y) cos(Cxy) dy vs psi0(x)
def quad(f, a_, b_, n=4000):
    hh = (b_ - a_) / n
    s = 0.5 * (f(a_) + f(b_))
    for i in range(1, n):
        s += f(a_ + i * hh)
    return s * hh
ratios = []
for x in [0.1, 0.35, 0.6, 0.85, 1.0]:
    F = quad(lambda y: legendre_eval(ks, a, y) * math.cos(C * x * y), -1, 1)
    ratios.append(F / legendre_eval(ks, a, x))
mu0 = sum(ratios) / len(ratios)
spread = max(abs(r - mu0) for r in ratios) / abs(mu0)
lam0 = C * mu0 * mu0 / (2 * math.pi)
print(f"self-transform: mu0 = {mu0:.6e}, ratio spread {spread:.2e}")
print(f"lambda0 = {lam0:.15f}, leakage 1-lambda0 = {1 - lam0:.3e}  (e^-2c = {math.exp(-2 * C):.3e})")

# 4. m2 identity for the LOGAN bump (validates moment-generic A-correction)
def bessi0(z):
    t, s, n = 1.0, 1.0, 1
    while n < 300:
        t *= (z * z / 4) / (n * n)
        s += t
        if t < 1e-18 * s: break
        n += 1
    return s
norm = C / (2 * math.sinh(C))
m2_logan = quad(lambda y: y * y * norm * bessi0(C * math.sqrt(max(0.0, 1 - y * y))), -1, 1, 20000)
m2_closed = math.cosh(C) / math.sinh(C) / C - 1.0 / (C * C)
print(f"m2(logan): numeric {m2_logan:.12f} vs coth(c)/c - 1/c^2 = {m2_closed:.12f}  diff {abs(m2_logan-m2_closed):.2e}")

# 5. eta-hat: Legendre-Bessel vs direct quadrature; also m2 and lambda for prolate
N0 = a[0] * math.sqrt(2.0)  # int psi0 = a0 * sqrt(2)
def etahat_lb(t):
    jn = sph_jn(ks[-1], t)
    s = 0.0
    for i, k in enumerate(ks):
        s += a[i] * math.sqrt(k + 0.5) * (2.0 * (-1) ** (k // 2)) * jn[k]
    return s / N0
for t in [0.0, 5.0, 14.0, 25.0, C]:
    direct = quad(lambda y: legendre_eval(ks, a, y) * math.cos(t * y), -1, 1) / N0
    lb = etahat_lb(t)
    print(f"eta-hat({t:５.2f}): legendre-bessel {lb:+.9e}  quadrature {direct:+.9e}  diff {abs(lb-direct):.1e}")
m2_pro = quad(lambda y: y * y * legendre_eval(ks, a, y), -1, 1, 20000) / N0
print(f"m2(prolate) = {m2_pro:.9f}  (logan's: {m2_closed:.9f})")
eps = 2.822e-7
lam_pro = quad(lambda y: legendre_eval(ks, a, y) * math.exp(eps * y / 2), -1, 1, 20000) / N0
print(f"lambda(prolate, eps={eps}) = {lam_pro:.9f}")
# compare zero-side weights at matched band: logan ell vs prolate eta-hat
def ell(t):
    w2 = t * t - C * C
    aw = math.sqrt(abs(w2))
    nrm = C / math.sinh(C)
    if aw < 1e-6: return nrm
    return nrm * (math.sinh(aw) / aw if w2 < 0 else math.sin(aw) / aw)
print("\nweight comparison (zero-side damping):")
for t in [0.0, 0.25 * C, 0.5 * C, 0.75 * C, 0.9 * C, C]:
    print(f"  t/c={t / C:4.2f}: logan {ell(t):+.6e}   prolate {etahat_lb(t):+.6e}")
