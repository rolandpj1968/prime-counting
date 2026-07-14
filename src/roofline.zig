//! Roofline listener — a deliberately-paced fan experiment. Alternates a
//! compute-bound (L1-resident) and a memory-bound (DRAM) sieve phase, each held
//! long enough for the fan to respond, announcing each transition. Run in the
//! FOREGROUND and listen: compute-bound = core burning watts at full IPC = fan
//! rises; memory-bound = core stalled on DRAM = fan eases. Correlate ear to
//! regime — the thing the fast background sweep didn't let you do.
//!
//!   zig build-exe -O ReleaseFast -mcpu=native src/roofline.zig -femit-bin=./sieve-roofline
//!   ./sieve-roofline        # Ctrl-C to stop
//!
//! Its own root/main, so it doesn't touch the main sieve build.

const std = @import("std");
const common = @import("common.zig");
const wheel = @import("wheel.zig");
const sieve = @import("sieve.zig");
const store_bit_packed = @import("stores/bit_packed.zig");

const W = wheel.Wheel(&[_]u64{}); // all-wheel: strike-heavy → biggest power swing
const N: u64 = 1_000_000_000;
const HOLD_NS: u64 = 15 * std.time.ns_per_s;

fn phase(comptime S: type, gpa: std.mem.Allocator, banner: []const u8, hint: []const u8) !void {
    std.debug.print("\n============================================================\n", .{});
    std.debug.print(">>> {s}\n>>> {s}\n", .{ banner, hint });
    std.debug.print("============================================================\n", .{});
    var st = try S.init(gpa, N);
    defer S.deinit(&st, gpa);
    const t_start = common.nowNs();
    var runs: u64 = 0;
    var best: u64 = std.math.maxInt(u64);
    while (common.nowNs() - t_start < HOLD_NS) {
        const t0 = common.nowNs();
        S.sieve(&st, N);
        best = @min(best, common.nowNs() - t0);
        runs += 1;
    }
    const rate = @as(f64, @floatFromInt(N)) / (@as(f64, @floatFromInt(best)) / 1e9) / 1e6;
    std.debug.print("    held ~15s: {d} runs, {d:.0} M ints/s\n", .{ runs, rate });
}

pub fn main() !void {
    const gpa = std.heap.page_allocator;
    const Compute = sieve.Sieve(W, store_bit_packed, 8 * 1024); // L1-resident, compute-bound
    const Memory = sieve.Sieve(W, store_bit_packed, 128 * 1024 * 1024); // DRAM, memory-bound

    std.debug.print("Roofline listener. Ctrl-C to stop. Listen to the fan.\n", .{});
    while (true) {
        try phase(Compute, gpa, "COMPUTE-BOUND  (8 KiB, L1-resident)", "fan should RISE  (core at full IPC, burning watts)");
        try phase(Memory, gpa, "MEMORY-BOUND   (128 MiB, DRAM)", "fan should EASE  (core stalled waiting on RAM)");
    }
}
