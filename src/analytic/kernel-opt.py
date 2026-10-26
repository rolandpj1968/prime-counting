#!/usr/bin/env python3
"""Kernel optimisation for pistar: optimise the functional the code actually pays.

Every pistar kernel is an even bump eta on [-1, 1] with int eta = 1, dilated
to the window |log(n/x)| <= eps, eps = c/T. Zero gamma is weighted by
eta_hat(eps*gamma), eta_hat(u) = int eta(y) cos(u y) dy, so the zeros dropped
above T are exactly u > c. Per conjugate pair a dropped zero contributes about
2 sqrt(x) eta_hat(u) / (gamma L), and zeros arrive with density
w(gamma) = log(gamma/2pi)/2pi. So, with u = eps*gamma:

  certified tail   R1 = (sqrt(x) / (pi L)) * J1,   J1 = int_c^inf |eta_hat(u)| log(u/(2 pi eps)) / u du
  RMS (random phases)
                   R2 = sqrt(x eps / (pi L^2) * J2), J2 = int_c^inf eta_hat(u)^2 log(u/(2 pi eps)) / u^2 du

Logan minimises int_c^inf |eta_hat|/u (J1 without the density weight); the
prolate psi_0 maximises in-band L2 energy relative to total energy. Neither is
the optimum of J1 or J2. Here:

  - kernels are eta = sum_k a_k Pbar_{2k}(y), Pbar_n = sqrt(n + 1/2) P_n, so
    eta_hat(u) = sum_k a_k sqrt(2k + 1/2) 2 (-1)^k j_{2k}(u)  (exact), and
    int eta = 1  <=>  a_0 = 1/sqrt(2)
  - J2-optimal: a quadratic form under one linear constraint -> linear solve
  - J1-optimal: weighted L1 -> linear program (scipy HiGHS)
  - Logan and prolate evaluated in the same machinery (Logan in closed form,
    prolate by the Bouwkamp eigen-solve, as in prolate-proto.py)

Usage: kernel-opt.py [x] [T]   (defaults: the 1e15 / 251M-zero runs, T = 1.01e8)
"""
import math
import sys

import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import linprog
from scipy.special import i0, spherical_jn

X = float(sys.argv[1]) if len(sys.argv) > 1 else 1e15
T = float(sys.argv[2]) if len(sys.argv) > 2 else 1.01e8
L = math.log(X)


# ---------------------------------------------------------------- basis
def basis_hat(K, u):
    """G[i, k] = transform of Pbar_{2k} at u[i]."""
    G = np.empty((len(u), K))
    for k in range(K):
        G[:, k] = math.sqrt(2 * k + 0.5) * 2.0 * (-1) ** k * spherical_jn(2 * k, u)
    return G


def basis_val(K, y):
    """V[i, k] = Pbar_{2k}(y[i])."""
    V = np.empty((len(y), K))
    P = np.polynomial.legendre
    for k in range(K):
        e = np.zeros(2 * k + 1)
        e[2 * k] = 1.0
        V[:, k] = math.sqrt(2 * k + 0.5) * P.legval(y, e)
    return V


# ---------------------------------------------------------------- kernels
def logan_hat(c, u):
    w2 = u * u - c * c
    out = np.empty_like(u)
    neg = w2 < 0
    aw = np.sqrt(np.abs(w2))
    small = aw < 1e-6
    out[~small & neg] = np.sinh(aw[~small & neg]) / aw[~small & neg]
    out[~small & ~neg] = np.sin(aw[~small & ~neg]) / aw[~small & ~neg]
    out[small] = 1.0
    return out * c / math.sinh(c)


def logan_val(c, y):
    return c / (2 * math.sinh(c)) * i0(c * np.sqrt(np.maximum(0.0, 1 - y * y)))


def prolate_coeffs(c, K):
    """Ground state of the prolate ODE on the even chain, scaled so int = 1."""
    ks = np.arange(0, 2 * K, 2, dtype=float)
    bkk = np.where(ks > 0, (2 * ks * (ks + 1) - 1) / ((2 * ks + 3) * (2 * ks - 1)), 1 / 3)
    diag = ks * (ks + 1) + c * c * bkk
    kk = ks[:-1]
    off = c * c * (kk + 2) * (kk + 1) / ((2 * kk + 3) * np.sqrt((2 * kk + 1) * (2 * kk + 5)))
    _, v = eigh_tridiagonal(diag, off, select="i", select_range=(0, 0))
    a = v[:, 0]
    return a / (a[0] * math.sqrt(2))  # int Pbar_0 = sqrt(2)


# ---------------------------------------------------------------- functionals
class Grid:
    """u-grid on [c, U] fine enough for the 2pi oscillation, plus an
    asymptotic tail beyond U: eta_hat ~ 2 eta(1) sin(u)/u."""

    def __init__(self, c, K, du=0.02, umax_mult=60.0):
        self.c = c
        self.eps = c / T
        self.U = umax_mult * c
        self.u = np.arange(c, self.U, du)
        self.du = du
        self.w = np.log(self.u / (2 * math.pi * self.eps))  # = log(gamma / 2pi)
        self.G = basis_hat(K, self.u)

    def J(self, hat, eta1):
        u, w, du = self.u, self.w, self.du
        wU = math.log(self.U / (2 * math.pi * self.eps))
        j1 = np.trapezoid(np.abs(hat) * w / u, dx=du) + (4 / math.pi) * abs(eta1) * wU / self.U
        j2 = np.trapezoid(hat * hat * w / (u * u), dx=du) + 2 * eta1 ** 2 * wU / (3 * self.U ** 3)
        return j1, j2


def R1(j1):
    return math.sqrt(X) / (math.pi * L) * j1


def R2(j2, eps):
    return math.sqrt(X * eps / (math.pi * L * L) * j2)


def j2_opt(gr, K, ridge):
    """min a^T Q a + ridge * |eta|^2 (Parseval: |eta|_2^2 = |a|^2) s.t. a_0 = 1/sqrt 2."""
    G, u, w, du = gr.G, gr.u, gr.w, gr.du
    W = (w / (u * u)) * du
    Q = G.T @ (G * W[:, None]) + ridge * np.eye(K)
    e0 = np.zeros(K)
    e0[0] = 1.0
    z = np.linalg.solve(Q, e0)
    return z / (z[0] * math.sqrt(2))


def j1_opt(c, K, du=0.025, umax_mult=40.0):
    """Weighted L1 by LP on the FULL u-grid (sparse, HiGHS), with the tail
    beyond U priced through eta(1). A subsampled grid lets the LP park zeros
    on the samples and oscillate in between, and a short U lets it hide mass
    beyond U: at U = 20c the 'gain' over Logan read 15-20%, at U = 40c 3-5%.
    Returns coefficients; evaluate them on a finer grid than this one."""
    from scipy import sparse
    gr = Grid(c, K, du=du, umax_mult=umax_mult)
    G = sparse.csr_matrix(gr.G)
    q = (gr.w / gr.u) * du
    n = len(q)
    wU = math.log(gr.U / (2 * math.pi * gr.eps))
    tailc = (4 / math.pi) * wU / gr.U
    one = np.array([[math.sqrt(2 * k + 0.5) for k in range(K)]])  # eta(1) = one . a
    I = sparse.identity(n, format="csr")
    z = sparse.csr_matrix((n, 1))
    A = sparse.vstack([
        sparse.hstack([G, -I, z]), sparse.hstack([-G, -I, z]),
        sparse.hstack([sparse.csr_matrix(one), sparse.csr_matrix((1, n)), sparse.csr_matrix([[-1.0]])]),
        sparse.hstack([sparse.csr_matrix(-one), sparse.csr_matrix((1, n)), sparse.csr_matrix([[-1.0]])]),
    ]).tocsr()
    cvec = np.concatenate([np.zeros(K), q, [tailc]])
    Aeq = sparse.csr_matrix(([1.0], ([0], [0])), shape=(1, K + n + 1))
    res = linprog(cvec, A_ub=A, b_ub=np.zeros(2 * n + 2), A_eq=Aeq, b_eq=[1 / math.sqrt(2)],
                  bounds=[(None, None)] * K + [(0, None)] * (n + 1), method="highs")
    if not res.success:
        raise RuntimeError(res.message)
    return res.x[:K]


def report(name, gr, hat_on_grid, eta1, eta_samples):
    j1, j2 = gr.J(hat_on_grid, eta1)
    return dict(name=name, j1=j1, j2=j2, r1=R1(j1), r2=R2(j2, gr.eps),
                etamax=float(np.max(np.abs(eta_samples))), etaneg=float(np.min(eta_samples)), eta1=eta1)


def main():
    ys = np.linspace(-1, 1, 2001)
    print(f"x = {X:.3g}  T = {T:.3g}  L = {L:.3f}   (R1 = certified tail, R2 = RMS model; window ints = 2 x c / T)")
    print(f"{'c':>4} {'kernel':<9} {'J1':>10} {'J2':>10} {'R1 cert':>10} {'R2 rms':>10} {'max|eta|':>9} {'min eta':>9} {'eta(1)':>10}")
    rows = []
    # J2-opt is solved in double precision: its J2 floors near 1e-16 from
    # c ~ 20, so rows past that measure the solver, not the kernel. J1-opt is
    # run to c = 16 (the LP's 1e-7 tolerances meet J1 ~ 1e-6 beyond that).
    for c in [6, 8, 10, 12, 14, 16, 18, 20]:
        K = int(c + 30)
        gr = Grid(c, K, du=0.005)
        out = []
        # Logan
        out.append(report("logan", gr, logan_hat(c, gr.u), float(logan_val(c, np.array([1.0]))[0]), logan_val(c, ys)))
        # prolate
        V = basis_val(K, ys)
        kers = [("prolate", prolate_coeffs(c, K)), ("J2-opt", j2_opt(gr, K, ridge=1e-14))]
        if c <= 16:
            kers.append(("J1-opt", j1_opt(c, K)))
        for name, a in kers:
            eta = V @ a
            eta1 = float(sum(a[k] * math.sqrt(2 * k + 0.5) for k in range(K)))
            out.append(report(name, gr, gr.G @ a, eta1, eta))
        for r in out:
            print(f"{c:>4} {r['name']:<9} {r['j1']:>10.3e} {r['j2']:>10.3e} {r['r1']:>10.3e} {r['r2']:>10.3e} {r['etamax']:>9.3g} {r['etaneg']:>9.2g} {r['eta1']:>10.2e}")
            rows.append((c, r))
        print()
    # c needed for targets
    print("c needed (log-interpolated) for each target:")
    names = ["logan", "prolate", "J2-opt", "J1-opt"]
    for key, tgts in [("r2", [0.05, 0.01, 1e-3]), ("r1", [0.5, 0.1])]:
        for tgt in tgts:
            cs = []
            for nm in names:
                pts = [(c, r[key]) for c, r in rows if r["name"] == nm]
                cneed = None
                for (c0, v0), (c1, v1) in zip(pts, pts[1:]):
                    if v0 >= tgt > v1 and v1 > 0:
                        f = (math.log(v0) - math.log(tgt)) / (math.log(v0) - math.log(v1))
                        cneed = c0 + f * (c1 - c0)
                        break
                cs.append(cneed)
            print(f"  {key} <= {tgt:<6g}: " + "  ".join(f"{nm} {('%.2f' % cv) if cv else '  -  '}" for nm, cv in zip(names, cs)))


if __name__ == "__main__":
    main()
