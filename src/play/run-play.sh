# zig build-exe  -O ReleaseFast --dep rs -Mroot=play.zig -Mrs=../rangesieve.zig -femit-bin=play

# ./play 1000 10
# ./play 1000000 100
# ./play 1000000000 1000

./play 8 2
./play 64 4
./play 512 8
./play 4096 16
./play 32768 32
./play 262144 64
./play 2097152 128
./play 16777216 256
./play 134217728 512
./play 1073741824 1024
./play 8589934592 2048
./play 68719476736 4096
./play 549755813888 8192

echo
echo

./play 4 2
./play 16 4
./play 64 8
./play 256 16
./play 1024 32
./play 4096 64
./play 16384 128
./play 65536 256
./play 262144 512
./play 1048576 1024
./play 4194304 2048
./play 16777216 4096
./play 67108864 8192
./play 268435456 16384
./play 1073741824 32768
./play 4294967296 65536
./play 17179869184 131072
./play 68719476736 262144
./play 274877906944 524288
