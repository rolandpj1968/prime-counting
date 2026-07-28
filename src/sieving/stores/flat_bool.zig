//! Backing store: one byte per integer ([]bool). The dumb baseline — the
//! correctness oracle and performance floor.

const std = @import("std");

pub const name = "flat []bool";

pub const flags_per_byte: u64 = 1; // store-size knob converts bytes -> flags

pub const State = struct { b: []bool };

pub fn footprintBytes(n: u64) usize {
    return @intCast(n + 1);
}

pub fn init(gpa: std.mem.Allocator, n: u64) !State {
    return .{ .b = try gpa.alloc(bool, @intCast(n + 1)) };
}

pub fn deinit(st: *State, gpa: std.mem.Allocator) void {
    gpa.free(st.b);
}

pub fn clearAll(st: *State, n: u64) void {
    _ = n;
    @memset(st.b, false);
}

pub fn get(st: *const State, i: u64) bool {
    return st.b[@intCast(i)];
}

pub fn set(st: *State, i: u64) void {
    st.b[@intCast(i)] = true;
}

pub fn countComposites(st: *const State, n: u64) u64 {
    var c: u64 = 0;
    var i: u64 = 0;
    while (i <= n) : (i += 1) {
        if (st.b[@intCast(i)]) c += 1;
    }
    return c;
}
