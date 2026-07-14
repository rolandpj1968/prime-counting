//! Backing store: std.DynamicBitSet — the idiomatic library option. Same 1
//! bit/integer as bit_packed, but through the std abstraction.

const std = @import("std");

pub const name = "std.DynamicBitSet";

pub const flags_per_byte: u64 = 8; // store-size knob converts bytes -> flags

pub const State = struct { bits: std.DynamicBitSet };

pub fn footprintBytes(n: u64) usize {
    const words = (n + 1 + 63) / 64;
    return @intCast(words * 8);
}

pub fn init(gpa: std.mem.Allocator, n: u64) !State {
    return .{ .bits = try std.DynamicBitSet.initEmpty(gpa, @intCast(n + 1)) };
}

pub fn deinit(st: *State, gpa: std.mem.Allocator) void {
    _ = gpa; // managed bitset stores its own allocator
    st.bits.deinit();
}

pub fn clearAll(st: *State, n: u64) void {
    st.bits.setRangeValue(.{ .start = 0, .end = @intCast(n + 1) }, false);
}

pub fn get(st: *const State, i: u64) bool {
    return st.bits.isSet(@intCast(i));
}

pub fn set(st: *State, i: u64) void {
    st.bits.set(@intCast(i));
}

pub fn countComposites(st: *const State, n: u64) u64 {
    _ = n;
    return @intCast(st.bits.count());
}
