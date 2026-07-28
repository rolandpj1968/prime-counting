# zig build-exe  -O ReleaseFast --dep rs -Mroot=play.zig -Mrs=../rangesieve.zig -femit-bin=play

# ./play 1000 10
# ./play 1000000 100
# ./play 1000000000 1000

# time ./play 8 2
# time ./play 64 4
# time ./play 512 8
# time ./play 4096 16
# time ./play 32768 32
# time ./play 262144 64
# time ./play 2097152 128
# time ./play 16777216 256
# time ./play 134217728 512
# time ./play 1073741824 1024
# time ./play 8589934592 2048
# time ./play 68719476736 4096
# time ./play 549755813888 8192

echo
echo

# time ./play 4 2
# time ./play 16 4
# time ./play 64 8
# time ./play 256 16
# time ./play 1024 32
# time ./play 4096 64
# time ./play 16384 128
# time ./play 65536 256
# time ./play 262144 512
# time ./play 1048576 1024
# time ./play 4194304 2048
# time ./play 16777216 4096
# time ./play 67108864 8192
# time ./play 268435456 16384
# time ./play 1073741824 32768
# time ./play 4294967296 65536
# time ./play 17179869184 131072
# time ./play 68719476736 262144
# time ./play 274877906944 524288
# time ./play 1099511627776 1048576
# time ./play 4398046511104 2097152
# time ./play 17592186044416 4194304


echo "10^10"
time ./play 10000000000 100000

echo
echo
echo "10^12"
time ./play 1000000000000 1000000

echo
echo
echo "10^14"
time ./play 100000000000000 10000000
