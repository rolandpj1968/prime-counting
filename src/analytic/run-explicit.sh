#!/bin/bash
# zig build-exe -O ReleaseFast -mcpu=native --dep rs -Mroot=explicit.zig -Mrs=../rangesieve.zig -femit-bin=explicit
# Zeros tables (Odlyzko, +-4e-9): ../../zeros/zeros1.txt (first 1e5), zeros6.txt (first 2e6)

time ./explicit 1000000    ../../zeros/zeros6.txt
echo
time ./explicit 100000000  ../../zeros/zeros6.txt
echo
time ./explicit 1000000000 ../../zeros/zeros6.txt
