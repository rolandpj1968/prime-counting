//! RPJ playground

const std = @import("std");
const assert = std.debug.assert;
// zig run --dep rs -Mroot=play4.zig -Mrs=../rangesieve.zig
const rs = @import("rs");

pub const Counters = struct { n: u64 = 0, x0: u64 = 0, a0: u64 = 0, a1: u64 = 0, phi1: u64 = 0, pi: u64 = 0, pi_: u64 = 0, cb: u64 = 0, cb_: u64 = 0, memo: u64 = 0 };

pub const NodeKey = struct { x: u64, a: u64 };

pub const NodeCounts = std.AutoHashMap(NodeKey, u64);
pub const NodeMemo = std.AutoHashMap(NodeKey, u64);

/// floor(sqrt(n)), exact for u64.
pub fn isqrt(n: u64) u64 {
    if (n < 2) return n;
    var x: u64 = @intFromFloat(@sqrt(@as(f64, @floatFromInt(n))));
    while (x * x > n) x -= 1;
    while ((x + 1) * (x + 1) <= n) x += 1;
    return x;
}

fn icbrt(x: u64) u64 {
    if (x == 0) return 0;
    var r: u64 = @intFromFloat(std.math.pow(f64, @floatFromInt(x), 1.0 / 3.0));
    if (r == 0) r = 1;
    while (r * r * r > x) r -= 1;
    while ((r + 1) * (r + 1) * (r + 1) <= x) r += 1;
    return r;
}

pub const PhiContext = struct { primes: []const u64, pis: []const u64, memo: *NodeMemo, c: *Counters, nc: *NodeCounts, opt_phi1: bool, opt_pi: bool, opt_cb: bool, opt_memo: bool };

fn phi(x: u64, a: u64, ctx: *const PhiContext) !u64 {
    ctx.c.n += 1;

    const xa = NodeKey{ .x = x, .a = a };
    const result = try ctx.nc.getOrPut(xa);
    if (!result.found_existing) {
        result.value_ptr.* = 0;
    }
    result.value_ptr.* += 1;

    if (x == 0) {
        ctx.c.x0 += 1;
        return 0;
    }
    if (a == 0) {
        ctx.c.a0 += 1;
        return x;
    }
    // if (a == 1) {
    //     ctx.c.a1 += 1;
    //     return x - x / 2;
    // }

    if (ctx.opt_memo) {
        if (ctx.memo.get(xa)) |phi_val| {
            ctx.c.memo += 1;
            return phi_val;
        }
    }

    var b_max = a;

    var phi_val = x;

    if (ctx.opt_phi1) {
        const p_a = ctx.primes[a];

        if (x <= p_a) {
            ctx.c.phi1 += 1;
            return 1;
        }

        if (ctx.opt_pi) {
            if (x < p_a * p_a) {
                assert(x > p_a);
                const sqrt_x = isqrt(x);
                const pi_sqrt_x = ctx.pis[sqrt_x];
                assert(pi_sqrt_x < a);
                if (x < ctx.pis.len) {
                    ctx.c.pi += 1;
                    // pi's are just special-case phi memo-isation
                    return ctx.pis[x] - a + 1;
                } else {
                    ctx.c.pi_ += 1;
                    return try phi(x, pi_sqrt_x, ctx) - (a - pi_sqrt_x);
                }
            }

            if (ctx.opt_cb) {
                if (x < p_a * p_a * p_a) {
                    ctx.c.cb += 1;
                    assert(x >= p_a * p_a);
                    const cbrt_x = icbrt(x);
                    const pi_cbrt_x = ctx.pis[cbrt_x];
                    assert(pi_cbrt_x < a);

                    const n_lim = a - (pi_cbrt_x + 1);
                    phi_val += (n_lim * (n_lim + 1)) / 2;

                    // const p_pi_cbrt_x = ctx.primes[pi_cbrt_x];
                    // var x_o_p_bm1: u64 = 0;
                    // for ((pi_cbrt_x + 1)..(a + 1)) |b| {
                    //     const p_b = ctx.primes[b];
                    //     const x_o_p_b = x / p_b;
                    //     if (x_o_p_b / p_pi_cbrt_x == x_o_p_bm1 / p_pi_cbrt_x) {
                    //         ctx.c.cb_ += 1;
                    //     }
                    //     phi_val -= try phi(x / p_b, pi_cbrt_x, ctx);
                    //     x_o_p_bm1 = x_o_p_b;
                    // }

                    // Inverted loop, collecting common grandchildren together.
                    for ((pi_cbrt_x + 1)..(a + 1)) |b| {
                        const p_b = ctx.primes[b];
                        phi_val -= x / p_b;
                    }
                    for (1..(pi_cbrt_x + 1)) |d| {
                        const p_d = ctx.primes[d];

                        var x_o_p_bm1: u64 = 0;
                        var phi_val_last: u64 = 0;
                        for ((pi_cbrt_x + 1)..(a + 1)) |b| {
                            const p_b = ctx.primes[b];
                            const x_o_p_b = x / p_b;

                            if (x_o_p_bm1 == 0 or x_o_p_bm1 / p_d != x_o_p_b / p_d) {
                                phi_val_last = try phi(x_o_p_b / p_d, d - 1, ctx);
                            }
                            phi_val += phi_val_last;

                            x_o_p_bm1 = x_o_p_b;
                        }
                    }

                    b_max = pi_cbrt_x;
                }
            }
        }
    }

    // assert(b_max >= 1);
    // phi_val -= x / 2;

    for (1..(b_max + 1)) |b| {
        const p_b = ctx.primes[b];
        phi_val -= try phi(x / p_b, b - 1, ctx);
    }

    if (ctx.opt_memo) {
        try ctx.memo.put(xa, phi_val);
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

const NodeCountPair = struct {
    xa: NodeKey,
    count: u64,
};

fn lessThanByNodeCountPair(_: void, lhs: NodeCountPair, rhs: NodeCountPair) bool {
    if (lhs.xa.x != rhs.xa.x) {
        return lhs.xa.x < rhs.xa.x;
    }
    return lhs.xa.a < rhs.xa.a;
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

    var memo = std.AutoHashMap(NodeKey, u64).init(gpa);
    defer memo.deinit();

    const a = primes.len - 1;

    var c = Counters{};
    var nc = std.AutoHashMap(NodeKey, u64).init(gpa);
    defer nc.deinit();

    const pis = try pisToY(gpa, y, primes);
    defer gpa.free(pis);

    const ctx: PhiContext = .{ .primes = primes, .pis = pis, .memo = &memo, .c = &c, .nc = &nc, .opt_phi1 = false, .opt_pi = false, .opt_cb = false, .opt_memo = true };
    const phi_x_y = try phi(x, a, &ctx);

    const n_f: f64 = @floatFromInt(c.n);
    const x_f: f64 = @floatFromInt(x);
    const ln_n: f64 = @log(n_f);
    const ln_x = @log(x_f);

    std.debug.print("x: {d:>12} | y: {d:>6} ------ nodes/x: {d:>7.4} ln: {d:5.3} ------- phi(x,y): {d:>12} | {any} | nc: {d:>12} | mc: {d:>12}\n", .{ x, y, n_f / x_f, ln_n / ln_x, phi_x_y, c, nc.count(), memo.count() });

    // var ncs = std.ArrayList(NodeCountPair).empty;
    // defer ncs.deinit(gpa);

    // var it1 = nc.iterator();
    // while (it1.next()) |kv| {
    //     const xa = kv.key_ptr.*;
    //     const count = kv.value_ptr.*;
    //     try ncs.append(
    //         gpa,
    //         .{
    //             .xa = xa,
    //             .count = count,
    //         },
    //     );
    // }

    // std.sort.block(NodeCountPair, ncs.items, {}, lessThanByNodeCountPair);

    // for (ncs.items) |xac| {
    //     const xa = xac.xa;
    //     const count = xac.count;
    //     std.debug.print("    x: {d:>12} | a: {d:>6}    --x-- {d:>6}\n", .{ xa.x, xa.a, count });
    // }
    // std.debug.print("\n\n", .{});
    // std.debug.print("Nodes: \n", .{});
}
