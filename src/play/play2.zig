//! RPJ playground

const std = @import("std");
const assert = std.debug.assert;
// zig run --dep rs -Mroot=play2.zig -Mrs=../rangesieve.zig
const rs = @import("rs");

fn comptime_prime(comptime a: u64) u64 {
    return switch (a) {
        // 0 => 1, // for p_(b-1) - not needed
        1 => 2,
        2 => 3,
        3 => 5,
        4 => 7,
        5 => 11,
        6 => 13,
        7 => 17,
        8 => 19,
        else => unreachable,
    };
}

fn phi_comptime_a(x: u64, comptime a: u64) u64 {
    if (a == 0) return x;

    return phi_comptime_a(x, a - 1) - phi_comptime_a(x / comptime_prime(a), a - 1);
}

const TBL_8B_X_LIMIT = 256;
const TBL_8B_A_LIMIT = 8;

pub const Counters = struct { n: u64 = 0, x0: u64 = 0, a0: u64 = 0, a_min: u64 = 0, phi1: u64 = 0, x8b: u64 = 0, x8b_: u64 = 0, x8b_a: u64 = 0 };

fn prime(a: u64, primes: []const u64) u64 {
    return primes[a - 1];
}

fn phi(x: u64, a: u64, primes: []const u64, c: *Counters) u64 {
    c.n += 1;

    if (x == 0) {
        c.x0 += 1;
        return 0;
    }
    if (a == 0) {
        c.a0 += 1;
        return x;
    }

    const p_a = prime(a, primes);

    return phi(x, a - 1, primes, c) - phi(x / p_a, a - 1, primes, c);
}

fn phi_x_gt_0(x: u64, a: u64, primes: []const u64, c: *Counters) u64 {
    assert(0 < x);

    c.n += 1;

    if (x == 0) {
        c.x0 += 1;
        return 0;
    }
    if (a == 0) {
        c.a0 += 1;
        return x;
    }

    const p_a = prime(a, primes);

    var phi_val = phi_x_gt_0(x, a - 1, primes, c);
    if (p_a <= x) {
        phi_val -= phi_x_gt_0(x / p_a, a - 1, primes, c);
    }
    return phi_val;
}

fn phi_x0(x: u64, a: u64, primes: []const u64, c: *Counters) u64 {
    if (x == 0) {
        c.x0 += 1;
        return 0;
    }

    return phi_x_gt_0(x, a, primes, c);
}

fn phi_a_ge_a_min(x: u64, a: u64, comptime a_min: u64, primes: []const u64, c: *Counters) u64 {
    assert(a_min <= a);

    c.n += 1;

    if (x == 0) {
        c.x0 += 1;
        return 0;
    }

    if (a == a_min) {
        c.a_min += 1;
        return phi_comptime_a(x, a_min);
    }

    const p_a = prime(a, primes);

    return phi_a_ge_a_min(x, a - 1, a_min, primes, c) - phi_a_ge_a_min(x / p_a, a - 1, a_min, primes, c);
}

fn phi_a_min(x: u64, a: u64, comptime a_min: u64, primes: []const u64, c: *Counters) u64 {
    if (a < a_min) {
        c.a0 += 1; // TODO
        inline for (0..a_min) |comptime_a| {
            if (a == comptime_a)
                return phi_comptime_a(x, comptime_a);
        }
        unreachable;
    }

    return phi_a_ge_a_min(x, a, a_min, primes, c);
}

fn phi_simple(x: u64, a: u64, primes: []const u64, pis: []const u64, tbl8Bit: *const [TBL_8B_X_LIMIT][TBL_8B_A_LIMIT]u8, c: *Counters) u64 {
    c.n += 1;

    if (x == 0) {
        c.x0 += 1;
        return 0;
    }
    if (a == 0) {
        c.a0 += 1;
        return x;
    }

    const p_a = primes[a - 1];

    if (x <= p_a) {
        c.phi1 += 1;
        return 1;
    }

    if (x < TBL_8B_X_LIMIT) {
        if (a < TBL_8B_A_LIMIT) {
            c.x8b += 1;
            //return tbl8Bit[x][a];
        } else {
            c.x8b_ += 1;
            c.x8b_a += a;
        }
    }

    return phi_simple(x, a - 1, primes, pis, tbl8Bit, c) - phi_simple(x / p_a, a - 1, primes, pis, tbl8Bit, c);
}

fn pisToY(gpa: std.mem.Allocator, y: u64, primes: []u64) ![]u64 {
    const pis = try gpa.alloc(u64, y + 1);
    var a: u64 = 0;
    var i: u64 = 0;
    var pi: u64 = 0;
    while (i <= y) : (i += 1) {
        if (a < primes.len and primes[a] == i) {
            pi += 1;
            a += 1;
        }
        pis[i] = pi;
    }
    return pis;
}

fn genTbl8Bit(primes: []u64) [TBL_8B_X_LIMIT][TBL_8B_A_LIMIT]u8 {
    var tbl: [TBL_8B_X_LIMIT][TBL_8B_A_LIMIT]u8 = undefined;

    for (0..TBL_8B_X_LIMIT) |x| {
        for (0..TBL_8B_A_LIMIT) |a| {
            if (x == 0) {
                tbl[x][a] = 0;
            } else if (a == 0) {
                tbl[x][a] = @intCast(x);
            } else {
                const p_a: u8 = @intCast(primes[a - 1]);
                tbl[x][a] = tbl[x][a - 1] - tbl[x / p_a][a - 1];
            }
        }
    }

    return tbl;
}

pub fn main(init: std.process.Init) !void {
    var args = try init.minimal.args.iterateAllocator(init.gpa);
    _ = args.skip();
    const x_str = args.next() orelse {
        std.debug.print("Usage: zig run play.zig -- <x> <y>\n", .{});
        return;
    };
    const x = std.fmt.parseInt(u64, x_str, 10) catch |err| {
        std.debug.print("Error: '{s}' is not a valid u64 number (Error: {any})\n", .{ x_str, err });
        return;
    };
    const y_str = args.next() orelse {
        std.debug.print("Usage: zig run play.zig -- <x> <y>\n", .{});
        return;
    };
    const y = std.fmt.parseInt(u64, y_str, 10) catch |err| {
        std.debug.print("Error: '{s}' is not a valid u64 number (Error: {any})\n", .{ y_str, err });
        return;
    };

    const gpa = std.heap.page_allocator;

    const primes = try rs.basePrimes(gpa, y); // y+1?
    defer gpa.free(primes);

    const pis = try pisToY(gpa, y, primes);
    defer gpa.free(pis);

    // const tbl8Bit = genTbl8Bit(primes);

    const a = primes.len;

    var c = Counters{};

    // const phi_x_y = phi_simple(x, a, primes, pis, &tbl8Bit, &c);
    // const phi_x_y = phi(x, a, primes, &c);
    // const phi_x_y = phi_x0(x, a, primes, &c);
    const phi_x_y = phi_a_min(x, a, 16, primes, &c);

    const n_f: f64 = @floatFromInt(c.n);
    const x_f: f64 = @floatFromInt(x);
    const ln_n: f64 = @log(n_f);
    const ln_x = @log(x_f);

    std.debug.print("x: {d:>12} | y: {d:>6} ------ nodes/x: {d:>7.4} ln: {d:5.3} ------- phi(x,y): {d:>12} | {any}\n", .{ x, y, n_f / x_f, ln_n / ln_x, phi_x_y, c });
}
