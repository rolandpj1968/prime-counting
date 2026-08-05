#!/usr/bin/env python3
"""Independent referee for Zig's f64 transcendentals.

Reads "name arg_bits result_bits" triples from trigemit.zig, reconstructs the
exact f64 values from their bit patterns, recomputes each function in Python's
decimal module at 80 significant digits (pure software, no libm involvement —
unlike numpy/R, which would re-call the very library under test), and reports
the error in ulps of the true value.

Correctly-rounded means <= 0.5 ulp. glibc claims <= 1 ulp for these and is
believed sub-ulp in practice; anything materially above 1 ulp is a finding,
and anything above ~2^29 ulp would indicate a silently-downcast stub of the
kind we found in @log(f128) (compiler_rt logq).
"""
import struct
import sys
from decimal import Decimal, getcontext, localcontext

getcontext().prec = 80

PI = Decimal(
    "3.14159265358979323846264338327950288419716939937510"
    "58209749445923078164062862089986280348253421170679"
)


def f64(bits):
    return struct.unpack("<d", struct.pack("<Q", bits))[0]


def ulp(x):
    """Width of one f64 ulp at |x| (x a Decimal, exact)."""
    if x == 0:
        return Decimal(2) ** -1074
    e = x.copy_abs().adjusted()  # decimal exponent
    b = Decimal(2) ** (int(Decimal(x.copy_abs()).ln() / Decimal(2).ln()))
    # refine to the true binary exponent
    while b > x.copy_abs():
        b /= 2
    while b * 2 <= x.copy_abs():
        b *= 2
    _ = e
    return b * Decimal(2) ** -52


def dsin(x):
    with localcontext() as ctx:
        ctx.prec = 100
        # reduce into [-pi, pi] (args already are, but be safe)
        n = int(x / (2 * PI))
        x = x - 2 * PI * n
        term = x
        s = x
        k = 1
        while abs(term) > Decimal(10) ** -95:
            term = -term * x * x / ((2 * k) * (2 * k + 1))
            s += term
            k += 1
        return +s


def dcos(x):
    with localcontext() as ctx:
        ctx.prec = 100
        n = int(x / (2 * PI))
        x = x - 2 * PI * n
        term = Decimal(1)
        s = Decimal(1)
        k = 1
        while abs(term) > Decimal(10) ** -95:
            term = -term * x * x / ((2 * k - 1) * (2 * k))
            s += term
            k += 1
        return +s


def dsinh(x):
    with localcontext() as ctx:
        ctx.prec = 100
        e = x.exp()
        return +((e - 1 / e) / 2)


def dcosh(x):
    with localcontext() as ctx:
        ctx.prec = 100
        e = x.exp()
        return +((e + 1 / e) / 2)


def dexp(x):
    with localcontext() as ctx:
        ctx.prec = 100
        return +x.exp()


def dlog(x):
    with localcontext() as ctx:
        ctx.prec = 100
        return +x.ln()


def dlog1p(x):
    with localcontext() as ctx:
        ctx.prec = 100
        return +(1 + x).ln()


def dsqrt(x):
    with localcontext() as ctx:
        ctx.prec = 100
        return +x.sqrt()


def derfc(x):
    """erfc via the Maclaurin series for erf; |x| <= 7.5 needs care for
    large |x| where erf -> +-1 and cancellation eats the series, so switch
    to the continued-fraction-free asymptotic-free route: use the series
    for erf at 200 digits, which holds the cancellation."""
    with localcontext() as ctx:
        ctx.prec = 220
        t = x
        s = x
        k = 1
        x2 = x * x
        while True:
            t = -t * x2 / k
            term = t / (2 * k + 1)
            s += term
            if abs(term) < Decimal(10) ** -215:
                break
            k += 1
            if k > 5000:
                break
        two_over_sqrtpi = 2 / PI.sqrt()
        erf = two_over_sqrtpi * s
        return +(1 - erf)


FUNCS = {
    "sin": dsin,
    "cos": dcos,
    "sin_big": dsin,
    "sinh": dsinh,
    "cosh": dcosh,
    "exp": dexp,
    "log": dlog,
    "log1p": dlog1p,
    "sqrt": dsqrt,
    "erfc": derfc,
}

worst = {}
for line in sys.stdin:
    parts = line.split()
    if len(parts) != 3:
        continue
    name, ab, rb = parts
    if name not in FUNCS:
        continue
    a = Decimal(f64(int(ab, 16)))
    got = Decimal(f64(int(rb, 16)))
    exact = FUNCS[name](a)
    if exact == 0:
        continue
    err = abs(got - exact) / ulp(exact)
    cur = worst.get(name)
    if cur is None or err > cur[0]:
        worst[name] = (err, float(a))

print(f"{'function':<10} {'max err (ulp)':>16}  {'at argument':>18}   verdict")
for name in sorted(worst):
    e, a = worst[name]
    ef = float(e)
    verdict = (
        "correctly rounded" if ef <= 0.5
        else "sub-ulp (normal libm)" if ef <= 1.0
        else "LOOSE" if ef < 100 else "*** SUSPECT STUB ***"
    )
    print(f"{name:<10} {ef:>16.4f}  {a:>18.10g}   {verdict}")
