#!/bin/bash
# zig build-exe -O ReleaseFast -mcpu=native -lc --dep rs -Mroot=pistar.zig -Mrs=../rangesieve.zig -femit-bin=pistar
# Usage: pistar <x> <zeros-file> [c]  with lambda = c/T (default 7.5)

for x in 1e6 1e8 1e9 1e10 1e11 1e12 1e13; do
  time ./pistar ${x/e*/}$(printf '0%.0s' $(seq 1 ${x#*e})) ../../zeros/zeros6.txt | tail -2
  echo
done
