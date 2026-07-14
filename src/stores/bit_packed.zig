//! Backing store: hand-rolled []u64, 1 bit per integer, raw word access.
//! bit i lives at word[i >> 6], position i & 63.

const std = @import("std");

pub const name = "hand-rolled []u64";

pub const State = struct { w: []u64 };

fn numWords(n: u64) usize {
    return @intCast((n + 64) / 64); // ceil((n+1)/64)
}

pub fn footprintBytes(n: u64) usize {
    return numWords(n) * 8;
}

pub fn init(gpa: std.mem.Allocator, n: u64) !State {
    return .{ .w = try gpa.alloc(u64, numWords(n)) };
}

pub fn deinit(st: *State, gpa: std.mem.Allocator) void {
    gpa.free(st.w);
}

pub fn clearAll(st: *State, n: u64) void {
    _ = n;
    @memset(st.w, 0);
}

pub fn get(st: *const State, i: u64) bool {
    return ((st.w[@intCast(i >> 6)] >> @as(u6, @intCast(i & 63))) & 1) != 0;
}

pub fn set(st: *State, i: u64) void {
    st.w[@intCast(i >> 6)] |= @as(u64, 1) << @as(u6, @intCast(i & 63));
}

pub fn countComposites(st: *const State, n: u64) u64 {
    _ = n; // padding bits past n are always 0 (never struck), so popcount is exact
    var c: u64 = 0;
    for (st.w) |word| c += @popCount(word);
    return c;
}
