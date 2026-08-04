#!/usr/bin/env python3
"""Decode LMFDB bulk zeta-zeros files (Platt's format) to a binary hi/lo table.

Same source format as lmfdb2txt.py (see its docstring: 104-bit deltas in
units of 2^-101, exact integer accumulation mandatory). Output is the
double-double table pistar consumes directly:

  file  = [8-byte magic "ZROSDD01"] [u64 count] pair^count
  pair  = [f64 hi][f64 lo]   little-endian

where hi is the correctly-rounded f64 of the exact zero ordinate gamma and
lo is the correctly-rounded f64 of the exact remainder gamma - hi. The pair
represents gamma to ~2^-79 absolute at gamma ~ 1e8 (|lo| <= ulp(hi)/2, and
lo itself carries 53 bits) — far below the 2^-101 source precision and far
below need: the dd phase gamma*L needs delta(gamma) ~ 1e-11 to make the
coherent worst-case float bound negligible; 2^-79 leaves ~1e15 headroom.

Exactness: gamma = t0 + G/2^101 with t0 an f64 (dyadic rational), so
gamma = N/2^101 for an exact integer N. Python int/int division is
correctly rounded, giving hi; N - hi*2^101 is computed exactly via
as_integer_ratio (hi is dyadic, its denominator divides 2^101), giving lo.

Usage: lmfdb2bin.py out.bin zeros_14.dat zeros_5000.dat ...
       (pass .dat files in ascending <t> order; monotonicity and N(t)
       continuity are asserted, trailer bytes skipped as in lmfdb2txt.py)
"""
import struct
import sys

D = 1 << 101

out_path = sys.argv[1]
out = open(out_path, 'wb')
out.write(b"ZROSDD01")
out.write(struct.pack('<Q', 0))  # count, patched at the end

count = 0
last = 0.0
ncont = None  # N(t) continuity across blocks and files
buf = bytearray()
for fn in sys.argv[2:]:
    d = open(fn, 'rb').read()
    nb, = struct.unpack_from('<Q', d, 0)
    off = 8
    for _ in range(nb):
        t0, t1 = struct.unpack_from('<dd', d, off)
        n0, n1 = struct.unpack_from('<QQ', d, off + 16)
        assert ncont is None or n0 == ncont, (fn, n0, ncont)
        ncont = n1
        off += 32
        p, q = t0.as_integer_ratio()
        base = p * (D // q)  # exact t0 * 2^101 (q divides 2^101)
        G = 0
        for _ in range(n1 - n0):
            zlo, zmid, zhi = struct.unpack_from('<QIB', d, off)
            off += 13
            G += zlo | (zmid << 64) | (zhi << 96)
            N = base + G
            hi = N / D  # correctly-rounded f64 of the exact gamma
            assert hi > last, (fn, hi, last)
            last = hi
            hp, hq = hi.as_integer_ratio()
            lo = (N - hp * (D // hq)) / D  # exact remainder, rounded once
            buf += struct.pack('<dd', hi, lo)
            count += 1
        if len(buf) >= 1 << 22:
            out.write(buf)
            buf.clear()
    if off != len(d):
        print(f"{fn}: NOTE skipping {len(d) - off} trailer bytes", file=sys.stderr)
    print(f"{fn}: done, N={ncont}, last={last:.6f}", file=sys.stderr)
out.write(buf)
out.seek(8)
out.write(struct.pack('<Q', count))
out.close()
print(f"{out_path}: {count} zeros, {16 * count + 16} bytes", file=sys.stderr)
