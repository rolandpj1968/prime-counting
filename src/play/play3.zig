//! RPJ playground

const std = @import("std");
const assert = std.debug.assert;
// zig run --dep rs -Mroot=play2.zig -Mrs=../rangesieve.zig
const rs = @import("rs");

pub const Counters = struct { n: u64 = 0, x0: u64 = 0, a0: u64 = 0, phi1: u64 = 0, pi: u64 = 0, pi_: u64 = 0 };

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

    const p_a = primes[a];

    return phi(x, a - 1, primes, c) - phi(x / p_a, a - 1, primes, c);
}

fn phi1(x: u64, a: u64, primes: []const u64, c: *Counters) u64 {
    c.n += 1;

    if (x == 0) {
        c.x0 += 1;
        return 0;
    }
    if (a == 0) {
        c.a0 += 1;
        return x;
    }

    const p_a = primes[a];

    if (x <= p_a) {
        c.phi1 += 1;
        return 1;
    }

    return phi1(x, a - 1, primes, c) - phi1(x / p_a, a - 1, primes, c);
}

/// floor(sqrt(n)), exact for u64.
pub fn isqrt(n: u64) u64 {
    if (n < 2) return n;
    var x: u64 = @intFromFloat(@sqrt(@as(f64, @floatFromInt(n))));
    while (x * x > n) x -= 1;
    while ((x + 1) * (x + 1) <= n) x += 1;
    return x;
}

fn phi1pi(x: u64, a: u64, primes: []const u64, pis: []const u64, c: *Counters) u64 {
    c.n += 1;

    if (x == 0) {
        c.x0 += 1;
        return 0;
    }
    if (a == 0) {
        c.a0 += 1;
        return x;
    }

    const p_a = primes[a];

    if (x <= p_a) {
        c.phi1 += 1;
        return 1;
    }

    // Could be tighter, but avoids infinite recursion
    if (x < p_a * p_a) {
        const pi_x = pi_by_phi1p1(x, primes, pis, c);
        return 1 + pi_x - a;
    }

    return phi1pi(x, a - 1, primes, pis, c) - phi1pi(x / p_a, a - 1, primes, pis, c);
}

fn pi_by_phi1p1(x: u64, primes: []const u64, pis: []const u64, c: *Counters) u64 {
    if (x < pis.len) {
        c.pi += 1;
        return pis[x];
    }

    c.pi_ += 1;
    const sqrt_x = isqrt(x);
    const pi_sqrt_x = pis[sqrt_x];
    return phi1pi(x, pi_sqrt_x, primes, pis, c) - 1 + pi_sqrt_x;
}

fn phil(x: u64, a: u64, primes: []const u64, c: *Counters) u64 {
    c.n += 1;

    if (x == 0) {
        c.x0 += 1;
        return 0;
    }
    if (a == 0) {
        c.a0 += 1;
        return x;
    }

    var phi_val = x;

    for (1..(a + 1)) |b| {
        const p_b = primes[b];
        phi_val -= phi(x / p_b, b - 1, primes, c);
    }

    return phi_val;
}

fn phil1(x: u64, a: u64, primes: []const u64, c: *Counters) u64 {
    c.n += 1;

    if (x == 0) {
        c.x0 += 1;
        return 0;
    }
    if (a == 0) {
        c.a0 += 1;
        return x;
    }

    var phi_val = x;

    for (1..(a + 1)) |b| {
        const p_b = primes[b];
        phi_val -= phi(x / p_b, b - 1, primes, c);
    }

    return phi_val;
}

fn pisToY(gpa: std.mem.Allocator, y: u64, primes: []u64) ![]u64 {
    const pis = try gpa.alloc(u64, y + 1);
    var a: u64 = 1;
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

    const primes_tmp = try rs.basePrimes(gpa, y); // y+1?
    defer gpa.free(primes_tmp);

    // prepend 1 for natural access
    const primes = try gpa.alloc(u64, primes_tmp.len + 1);
    defer gpa.free(primes);

    primes[0] = 1;
    @memcpy(primes[1..], primes_tmp);

    const a = primes.len - 1;

    var c = Counters{};

    // const phi_x_y = phi(x, a, primes, &c);
    // const phi_x_y = phi1(x, a, primes, &c);

    // const phi_x_y = phil(x, a, primes, &c);

    const pis = try pisToY(gpa, y, primes);
    defer gpa.free(pis);

    const phi_x_y = phi1pi(x, a, primes, pis, &c);

    const n_f: f64 = @floatFromInt(c.n);
    const x_f: f64 = @floatFromInt(x);
    const ln_n: f64 = @log(n_f);
    const ln_x = @log(x_f);

    std.debug.print("x: {d:>12} | y: {d:>6} ------ nodes/x: {d:>7.4} ln: {d:5.3} ------- phi(x,y): {d:>12} | {any}\n", .{ x, y, n_f / x_f, ln_n / ln_x, phi_x_y, c });
}
