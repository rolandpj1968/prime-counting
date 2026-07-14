//! The benchmark harness and the sieve-implementation contract.
//!
//! A sieve *implementation* is any namespace (typically its own file under
//! src/impls/) exposing these decls — this is Zig's "interface": duck-typed at
//! comptime, zero runtime cost.
//!
//!   pub const name: []const u8;                       // human label
//!   pub const State: type;                            // impl-owned buffers
//!   pub fn footprintBytes(n: u64) usize;              // for the memory report
//!   pub fn init(gpa, n) !State;                       // allocate — untimed
//!   pub fn sieve(*State, n) void;                     // THE timed operation
//!   pub fn count(*const State, n) u64;                // separate count pass
//!   pub fn deinit(*State, gpa) void;
//!
//! Knob-tuned variants (segment size, wheel, ...) are `fn(comptime knobs) type`
//! families *inside* an impl's file — a different axis from this interface.

const std = @import("std");
const common = @import("common.zig");

/// Turn a missing/mis-named decl into an early, readable compile error instead
/// of a cryptic failure deep inside `run`. (Zig doesn't check the interface
/// until instantiation — this is the idiomatic workaround.)
pub fn assertSieveImpl(comptime Impl: type) void {
    const required = [_][]const u8{ "name", "State", "footprintBytes", "init", "sieve", "count", "deinit" };
    inline for (required) |decl| {
        if (!@hasDecl(Impl, decl)) {
            @compileError(@typeName(Impl) ++ ": missing required decl `" ++ decl ++ "` (see sieve contract in bench.zig)");
        }
    }
}

pub fn run(comptime Impl: type, gpa: std.mem.Allocator, n: u64, repeats: usize) !void {
    comptime assertSieveImpl(Impl);

    const mib = @as(f64, @floatFromInt(Impl.footprintBytes(n))) / (1024.0 * 1024.0);
    std.debug.print("\n=== {s} ===\n", .{Impl.name});
    std.debug.print("N = {d}   footprint = {d:.1} MiB   repeats = {d}\n", .{ n, mib, repeats });

    var st = try Impl.init(gpa, n);
    defer Impl.deinit(&st, gpa);

    var best_ns: u64 = std.math.maxInt(u64);
    std.debug.print("sieve runs (ms):", .{});
    var r: usize = 0;
    while (r < repeats) : (r += 1) {
        const t0 = common.nowNs();
        Impl.sieve(&st, n);
        const dt = common.nowNs() - t0;
        best_ns = @min(best_ns, dt);
        std.debug.print(" {d:.1}", .{common.ms(dt)});
    }
    const rate = @as(f64, @floatFromInt(n)) / (@as(f64, @floatFromInt(best_ns)) / 1e9) / 1e6;
    std.debug.print("  -> best {d:.1} ms  ({d:.0} M ints/s)\n", .{ common.ms(best_ns), rate });

    const tc0 = common.nowNs();
    const pi = Impl.count(&st, n);
    std.debug.print("count: {d:.1} ms (separate)\n", .{common.ms(common.nowNs() - tc0)});

    if (common.expectedPi(n)) |e| {
        std.debug.print("pi({d}) = {d}  (expected {d}) {s}\n", .{ n, pi, e, if (pi == e) "OK" else "FAIL" });
    } else {
        std.debug.print("pi({d}) = {d}\n", .{ n, pi });
    }
}
