#!/bin/bash
# zig build-exe -O ReleaseFast -mcpu=native -lc --dep rs -Mroot=smooth.zig -Mrs=../rangesieve.zig -femit-bin=smooth
# Usage: smooth <x> <zeros-file> [c]  with lambda = c/T (default 7.5)

time ./smooth 1000000    ../../zeros/zeros1.txt
echo
time ./smooth 100000000  ../../zeros/zeros6.txt
echo
time ./smooth 1000000000 ../../zeros/zeros6.txt
