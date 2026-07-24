#!/bin/bash
# 2xl (8-core) bandwidth-partition probe: scaling legs at 1e19 with the SAME
# default alpha fit as the 4xl sweep — matched mistuning, apples-to-apples.
# Partitioned slice predicts knee ~3-4 cores and t=8 well slower than the
# 4xl's 158.8 s; unpartitioned+quiet predicts t=8 matching or beating it.
# (t=1 skipped: 960.8/962.6 s established on two instances.)
# Invocation, from the repo dir: git pull && bash box8.sh   (inside tmux)
set -e
zig build-exe -O ReleaseFast -mcpu=native src/pi.zig -femit-bin=./pi

for t in 2 4 8; do
  pins=$(seq -s, 0 $((t - 1)))
  echo "=== t=$t ===" | tee -a ~/scale19_2xl.log
  ./pi 1e19 --pin-list "$pins" --check 2>&1 | tail -3 | tee -a ~/scale19_2xl.log
done

echo "ALL DONE"
