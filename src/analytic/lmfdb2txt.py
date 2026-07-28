#!/usr/bin/env python3
"""Decode LMFDB bulk zeta-zeros files (Platt's format) to a plain-text table.

Source: https://beta.lmfdb.org/riemann-zeta-zeros/ (send cookie "human=1" to
pass the JS gate; be polite — space out requests). Files data/zeros_<t>.dat,
named by starting height. Format (from README.dvi, validated against Odlyzko):

  file  = [u64 nblocks] block*
  block = [f64 t0][f64 t1][u64 N(t0)][u64 N(t1)] zero^(N(t1)-N(t0))
  zero  = [u64 lo][u32 mid][u8 hi]   (13 bytes, little-endian)

Each 104-bit Z = lo | mid<<64 | hi<<96 is the DELTA to the previous zero (the
first is relative to t0), in units of 2^-101 — absolute precision ~2.5e-31.
Deltas must be accumulated in exact integer arithmetic (a float cumulative sum
random-walks to ~1e-6 over millions of zeros); we round once per zero on
output. Output text (15 decimals) is exact to f64 ulp, i.e. ~gamma*2^-52 —
full 2^-101 precision needs a hi/lo pair format (the double-double rung).

Validation history (2026-07-28): first 4520 zeros vs Odlyzko zeros1 worst
diff 2.9e-9; first 2,001,052 vs Odlyzko zeros6 worst diff 3.0e-9 (both within
Odlyzko's own +-4e-9).

Usage: lmfdb2txt.py zeros_14.dat zeros_5000.dat ... > prefix.txt
       (pass files in ascending <t> order; monotonicity is asserted)
"""
import struct
import sys

last = 0.0
ncont = None  # N(t) continuity across blocks and files
for fn in sys.argv[1:]:
    d = open(fn, 'rb').read()
    nb, = struct.unpack_from('<Q', d, 0)
    off = 8
    for _ in range(nb):
        t0, t1 = struct.unpack_from('<dd', d, off)
        n0, n1 = struct.unpack_from('<QQ', d, off + 16)
        assert ncont is None or n0 == ncont, (fn, n0, ncont)
        ncont = n1
        off += 32
        G = 0
        for _ in range(n1 - n0):
            lo, mid, hi = struct.unpack_from('<QIB', d, off)
            off += 13
            G += lo | (mid << 64) | (hi << 96)
            g = t0 + G / 2**101
            assert g > last, (fn, g, last)
            last = g
            sys.stdout.write(f"{g:.15f}\n")
    # Some files carry trailing bytes beyond the advertised block count
    # (md5-canonical; decodes as a non-monotone pseudo-block — not zero
    # data). Honor the count, note and skip the trailer.
    if off != len(d):
        print(f"{fn}: NOTE skipping {len(d) - off} trailer bytes", file=sys.stderr)
    print(f"{fn}: done, N={ncont}, last={last:.6f}", file=sys.stderr)
