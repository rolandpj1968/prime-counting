//! RPJ playground

const std = @import("std");
const assert = std.debug.assert;
// zig run --dep rs -Mroot=play.zig -Mrs=../rangesieve.zig
const rs = @import("rs");

// pub fn main() !void {
//     std.debug.print("Hello world!\n", .{});
// }

pub const Counters = struct { n: u64 = 0, it: u64 = 0, x0: u64 = 0, ta: u64 = 0, one: u64 = 0, ones: u64 = 0, pi: u64 = 0, hipi: u64 = 0, ch: u64 = 0, lt_pxp: u64 = 0, l: u64 = 0, b: u64 = 0, la: u64 = 0, ba: u64 = 0 };

fn phi_simple(x: u64, a: u64, primes: []const u64, pis: []const u64, c: *Counters) u64 {
    c.n += 1;

    if (x == 0) {
        c.x0 += 1;
        return 0;
    }
    if (a == 0) {
        //c.a0 += 1;
        return x;
    }

    const p_a = primes[a - 1];

    return phi_simple(x, a - 1, primes, pis, c) - phi_simple(x / p_a, a - 1, primes, pis, c);
}

fn phi_simple_ltpa(x: u64, a: u64, primes: []const u64, pis: []const u64, c: *Counters) u64 {
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
        c.one += 1;
        return 1;
    }

    return phi_simple_ltpa(x, a - 1, primes, pis, c) - phi_simple_ltpa(x / p_a, a - 1, primes, pis, c);
}

fn phi_simple_pi(x: u64, a: u64, primes: []const u64, pis: []const u64, c: *Counters) u64 {
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
        c.pi += 1;
        c.one += 1;
        return 1;
    }

    if (x < p_a * p_a and x < pis.len) {
        const pi_x = pis[x];
        return 1 + pi_x - a;
    }

    return phi_simple_pi(x, a - 1, primes, pis, c) - phi_simple_pi(x / p_a, a - 1, primes, pis, c);
}

fn phi_linear(x: u64, a: u64, primes: []const u64, pis: []const u64, c: *Counters) u64 {
    c.n += 1;

    if (x == 0) {
        c.x0 += 1;
        return 0;
    }
    if (a == 0) {
        c.a0 += 1;
        return x;
    }

    var phi = x;

    var b: u64 = 0;
    while (b < a) : (b += 1) {
        c.it += 1;
        phi -= phi_linear(x / primes[b], b, primes, pis, c);
    }

    return phi;
}

fn phi_linear_ltpa(x: u64, a: u64, primes: []const u64, pis: []const u64, c: *Counters) u64 {
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
        c.one += 1;
        return 1;
    }

    var phi = x;

    var b: u64 = 0;
    while (b < a) : (b += 1) {
        c.it += 1;
        phi -= phi_linear_ltpa(x / primes[b], b, primes, pis, c);
    }

    return phi;
}

fn phi_linear_chop(x: u64, a: u64, primes: []const u64, pis: []const u64, c: *Counters) u64 {
    c.n += 1;

    if (x == 0) {
        c.x0 += 1;
        return 0;
    }
    if (a == 0) {
        c.a0 += 1;
        return x;
    }

    var phi = x;

    var b: u64 = 0;
    while (b < a) : (b += 1) {
        const p_b = primes[b];

        if (x < p_b) {
            c.ch += 1;
            break;
        }

        if (x < p_b * p_b) c.lt_pxp += 1;

        c.it += 1;
        phi -= phi_linear_chop(x / p_b, b, primes, pis, c);
    }

    return phi;
}

fn phi_linear_chop_ltpa(x: u64, a: u64, primes: []const u64, pis: []const u64, c: *Counters) u64 {
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
        c.one += 1;
        return 1;
    }

    var phi = x;

    var b: u64 = 0;
    while (b < a) : (b += 1) {
        const p_b = primes[b];

        // Never triggers, cos we've already bailed above on more general condition x <= p_a
        if (x < p_b) {
            c.ch += 1;
            break;
        }

        if (x < p_b * p_b) c.lt_pxp += 1;

        c.it += 1;
        phi -= phi_linear_chop_ltpa(x / p_b, b, primes, pis, c);
    }

    return phi;
}

fn phi_linear_chop2(x: u64, a: u64, primes: []const u64, pis: []const u64, c: *Counters) u64 {
    c.n += 1;

    if (x == 0) {
        c.x0 += 1;
        return 0;
    }
    if (a == 0) {
        c.a0 += 1;
        return x;
    }

    var phi = x;

    var b: u64 = 0;
    while (b < a) : (b += 1) {
        const p_b = primes[b];

        if (x < p_b) {
            c.ch += 1;
            break;
        }

        // <= p_b*2 is probably fine too - this is a very weak condition (even x <= p_b * p_bm1 is fine I think)
        if (x < p_b * 2) {
            // x/p_b = 1, so phi(x/p_b, b) = 1
            phi -= 1;
            continue;
        }

        c.it += 1;
        phi -= phi_linear_chop2(x / p_b, b, primes, pis, c);
    }

    return phi;
}

fn phi_linear_chop2_ltpa(x: u64, a: u64, primes: []const u64, pis: []const u64, c: *Counters) u64 {
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
        c.one += 1;
        return 1;
    }

    var phi = x;

    var b: u64 = 0;
    while (b < a) : (b += 1) {
        const p_b = primes[b];

        // Never triggers, cos we've already bailed above on more general condition x <= p_a
        if (x < p_b) {
            c.ch += 1;
            break;
        }

        // <= p_b*2 is probably fine too - this is a very weak condition (even x <= p_b * p_bm1 is fine I think)
        if (x < p_b * 2) {
            // x/p_b = 1, so phi(x/p_b, b) = 1
            phi -= 1;
            continue;
        }

        c.it += 1;
        phi -= phi_linear_chop2_ltpa(x / p_b, b, primes, pis, c);
    }

    return phi;
}

fn phi_linear_chop2_pi(x: u64, a: u64, primes: []const u64, pis: []const u64, c: *Counters) u64 {
    c.n += 1;

    if (x == 0) {
        c.x0 += 1;
        return 0;
    }
    if (a == 0) {
        c.a0 += 1;
        return x;
    }
    // Hrmm, this should be primes[a-1] but it works?
    const p_a = primes[a];
    if (x < pis.len and x < p_a * p_a) {
        const pi_x = pis[x];
        return 1 + (if (pi_x > a) pi_x - a else 0);
    }

    var phi = x;

    var b: u64 = 0;
    while (b < a) : (b += 1) {
        const p_b = primes[b];

        if (x < p_b) {
            c.ch += 1;
            break;
        }

        if (x < p_b * 2) {
            // x/p_b = 1, so phi(x/p_b, b) = 1
            phi -= 1;
            continue;
        }

        c.it += 1;
        phi -= phi_linear_chop2_pi(x / p_b, b, primes, pis, c);
    }

    return phi;
}

const MAX_TINY_A = 4;
const MAX_TINY_A_2 = 4;

fn explicit_prime(comptime a: u64) u64 {
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

fn phi_explicit_a(x: u64, comptime a: u64) u64 {
    if (a == 0) return x;

    return phi_explicit_a(x, a - 1) - phi_explicit_a(x / explicit_prime(a), a - 1);
}

fn phi_tiny_a(x: u64, a: u64) u64 {
    inline for (0..(MAX_TINY_A + 1)) |explicit_a| {
        if (a == explicit_a)
            return phi_explicit_a(x, explicit_a);
    }

    unreachable;
}

fn phi_linear_all(x: u64, a: u64, primes: []const u64, pis: []const u64, c: *Counters) u64 {
    c.n += 1;

    if (x == 0) {
        std.debug.print("                                                              x == 0 | a: {:>6}\n", .{a});
        c.x0 += 1;
        return 0;
    }

    if (a <= MAX_TINY_A) {
        c.ta += 1;
        return phi_tiny_a(x, a);
    }

    const p_a = primes[a - 1];

    if (x <= p_a) {
        std.debug.print("                                                              ONE x: {:>12} | a: {:>6}\n", .{ x, a });
        c.one += 1;
        return 1;
    }

    if (x <= p_a * p_a) {
        if (x < pis.len) {
            c.pi += 1;
            //std.debug.print("                                                              LO PI: x: {d:>12} | a: {d:>6} | p_a: {d:>6} | p_(a_1): {d:>6} | 10000/p_a(+1): {d:>6} | phi(x,a): {d:>6}\n", .{ x, a, p_a, primes[a], 10000 / primes[a], phi_simple(x, a, primes, pis, c) });
            const pi_x = pis[x];
            return 1 + pi_x - a;
        } else {
            c.hipi += 1;
            //std.debug.print("                                                              HI PI: x: {d:>12} | a: {d:>6} | p_a: {d:>6} | p_(a_1): {d:>6} | 10000/p_a(+1): {d:>6} | phi(x,a): {d:>6}\n", .{ x, a, p_a, primes[a], 10000 / primes[a], phi_simple(x, a, primes, pis, c) });
        }
    }

    if (x <= primes.len) {
        c.l += 1;
        c.la += a;
    } else {
        c.b += 1;
        c.ba += a;
    }

    var phi = x;

    inline for (1..(MAX_TINY_A_2 + 1)) |explicit_a| {
        phi -= phi_explicit_a(x / explicit_prime(explicit_a), explicit_a - 1);
    }

    var b: u64 = MAX_TINY_A_2 + 1;
    var p_bm1: u64 = explicit_prime(MAX_TINY_A_2);

    while (b <= a) : (b += 1) {
        const p_b = primes[b - 1];

        //if (x <= p_bm1 * p_b) {
        if (x / p_b <= p_bm1) {
            c.ones += 1;
            // x/p_b <= p_(b-1), so phi(x/p_b, b-1) = 1
            //    ... and same holds for the rest of the indices up to 'a'
            return phi - (a - (b - 1));
        }

        c.it += 1;
        phi -= phi_linear_all(x / p_b, b - 1, primes, pis, c);

        p_bm1 = p_b;
    }

    return phi;
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

pub fn main(init: std.process.Init) !void {
    var args = try init.minimal.args.iterateAllocator(init.gpa);
    //var args = init.args.iterate();
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

    const a = primes.len;

    // std.debug.print("\n", .{});
    // std.debug.print("x: {d:>6} / y: {d:>2} / a: {d:>2} / primes: {any}\n", .{ x, y, a, primes });

    var c = Counters{};

    //const phi_fn = phi_simple;
    //const phi_fn = phi_simple_pi;
    //const phi_fn = phi_simple_ltpa;
    //const phi_fn = phi_linear;
    //const phi_fn = phi_linear_ltpa;
    //const phi_fn = phi_linear_chop;
    //const phi_fn = phi_linear_chop_ltpa;
    //const phi_fn = phi_linear_chop2;
    //const phi_fn = phi_linear_chop2_ltpa;
    //const phi_fn = phi_linear_chop2_pi;
    const phi_fn = phi_linear_all;

    const phi_x_y = phi_fn(x, a, primes, pis, &c);

    const n_f: f64 = @floatFromInt(c.n);
    const x_f: f64 = @floatFromInt(x);
    const ln_n: f64 = @log(n_f);
    const ln_x = @log(x_f);

    std.debug.print("x: {d:>12} | y: {d:>6} ------ nodes/x: {d:>7.4} ln: {d:5.3} ------- phi(x,y): {d:>12} | {any}\n", .{ x, y, n_f / x_f, ln_n / ln_x, phi_x_y, c });

    const la_f: f64 = @floatFromInt(c.la);
    const l_f: f64 = @floatFromInt(c.l);
    std.debug.print("                                                    ------------------> la/l: {d:>5.3}\n", .{la_f / l_f});

    const ba_f: f64 = @floatFromInt(c.ba);
    const b_f: f64 = @floatFromInt(c.b);
    std.debug.print("                                                    ------------------> ba/b: {d:>5.3}\n", .{ba_f / b_f});

    // std.debug.print("\n", .{});
    // std.debug.print("phi({d:>6}, {d:>2}) = {d:>6} / [{any}]\n", .{ x, y, phi_x_y, c });

    //std.debug.print("phi({}, {}) = {}\n", .{ x, y, phi_x_y });

    if (x == y * y) {
        std.debug.print("\n", .{});
        std.debug.print("pi({} = {})\n", .{ x, phi_x_y + primes.len - 1 });
    }
}
