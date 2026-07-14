//! Backing store: bit-packed, parameterised over the word type.
//! `BitPacked(u64)` is the old default; u8/u16/u32 vary the RMW/popcount
//! granularity (density is always 8 flags/byte — the word doesn't change that).

const std = @import("std");

pub fn BitPacked(comptime Word: type) type {
    return struct {
        const wbits: u64 = @bitSizeOf(Word);
        const Shift = std.math.Log2Int(Word);

        pub const name = std.fmt.comptimePrint("[]u{d}", .{wbits});
        pub const flags_per_byte: u64 = 8;

        pub const State = struct { w: []Word };

        fn numWords(n: u64) usize {
            return @intCast((n + wbits) / wbits); // ceil((n+1)/wbits)
        }

        pub fn footprintBytes(n: u64) usize {
            return numWords(n) * @sizeOf(Word);
        }

        pub fn init(gpa: std.mem.Allocator, n: u64) !State {
            return .{ .w = try gpa.alloc(Word, numWords(n)) };
        }

        pub fn deinit(st: *State, gpa: std.mem.Allocator) void {
            gpa.free(st.w);
        }

        pub fn clearAll(st: *State, n: u64) void {
            _ = n;
            @memset(st.w, 0);
        }

        pub fn get(st: *const State, i: u64) bool {
            return ((st.w[@intCast(i / wbits)] >> @as(Shift, @intCast(i % wbits))) & 1) != 0;
        }

        pub fn set(st: *State, i: u64) void {
            st.w[@intCast(i / wbits)] |= @as(Word, 1) << @as(Shift, @intCast(i % wbits));
        }

        pub fn countComposites(st: *const State, n: u64) u64 {
            _ = n;
            var c: u64 = 0;
            for (st.w) |word| c += @popCount(word);
            return c;
        }
    };
}

test "BitPacked round-trip across word widths" {
    const a = std.testing.allocator;
    inline for (.{ u8, u16, u32, u64 }) |Word| {
        const S = BitPacked(Word);
        var st = try S.init(a, 1000);
        defer S.deinit(&st, a);
        S.clearAll(&st, 1000);
        S.set(&st, 500);
        S.set(&st, 999);
        try std.testing.expect(S.get(&st, 500));
        try std.testing.expect(!S.get(&st, 501));
        try std.testing.expectEqual(@as(u64, 2), S.countComposites(&st, 1000));
    }
}
