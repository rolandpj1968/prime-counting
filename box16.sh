#!/bin/bash
# 16-core Graviton characterization queue: scaling knee, alpha portrait, 1e20 rung.
# Invocation, from the repo dir: git pull && bash box16.sh   (inside tmux)
set -e
zig build-exe -O ReleaseFast -mcpu=native src/pi.zig -femit-bin=./pi

# 1. Core-scaling sweep at 1e19 — where is Graviton's bandwidth knee?
for t in 1 2 4 8 16; do
  pins=$(seq -s, 0 $((t - 1)))
  echo "=== t=$t ===" | tee -a ~/scale19.log
  ./pi 1e19 --pin-list "$pins" --check 2>&1 | tail -3 | tee -a ~/scale19.log
done

# 2. Machine portrait: fit alpha(x) for THIS box at 16 threads
echo "=== calibrate 16t ===" | tee ~/calib16.log
./pi --calibrate --budget 900 -t 16 2>&1 | tee -a ~/calib16.log

# 3. First ARM parallel ladder rung (laptop 6-core reference: 739 s)
echo "=== 1e20 16t ===" | tee ~/pi20_arm.log
./pi 1e20 --pin-list "$(seq -s, 0 15)" --check -v 2>&1 | tail -6 | tee -a ~/pi20_arm.log

echo "ALL DONE"
