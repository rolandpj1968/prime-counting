#!/bin/bash
# 1-core box follow-up queue (after the TLB probe):
#   1. serial alpha portrait -> ~/calib1.log   (pairs with 16-core calib for 1/sqrt(t) check)
#   2. first real x>2^64 ARM rung, serial 1e20 -> ~/pi20_serial_arm.log (~90 min)
# NOTE: pulls + rebuilds first. Run inside tmux: bash box1.sh
set -e
cd ~/prime-counting
git pull
zig build-exe -O ReleaseFast -mcpu=native src/pi.zig -femit-bin=./pi

echo "=== calibrate t=1 ===" | tee ~/calib1.log
./pi --calibrate -t 1 --budget 900 2>&1 | tee -a ~/calib1.log

echo "=== 1e20 serial, segw 524160 ===" | tee ~/pi20_serial_arm.log
./pi 1e20 -t 1 --segw 524160 -v --check 2>&1 | tee -a ~/pi20_serial_arm.log

echo "ALL DONE"
