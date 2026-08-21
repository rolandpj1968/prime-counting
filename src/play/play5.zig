//! RPJ playground
// zig run --dep rs -Mroot=play5.zig -Mrs=../rangesieve.zig -- 100

const std = @import("std");
const assert = std.debug.assert;
const rs = @import("rs");

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

const CHECK_FIND_PRED = true;

// Check predicate is all-true then all-false over the slice
fn checkFindFirstTruePred(
    comptime T: type,
    slice: []const T,
    ctx: anytype,
    pred: fn (@TypeOf(ctx), T) bool,
) void {
    var found_true = false;
    var i: usize = 0;
    while (i < slice.len) {
        if (pred(ctx, slice[i])) {
            found_true = true;
        } else {
            assert(found_true == false);
        }
        i += 1;
    }
}

// Assume predicate is all-false then all-true over the slice
fn findFirstTrueIndex(
    comptime T: type,
    slice: []const T,
    ctx: anytype,
    pred: fn (@TypeOf(ctx), T) bool,
) usize {
    if (CHECK_FIND_PRED) {
        checkFindFirstTruePred(T, slice, ctx, pred);
    }

    var i: usize = slice.len;
    while (i > 0) {
        i -= 1;
        if (!pred(ctx, slice[i])) {
            return i + 1;
        }
    }
    return 0;
}

const IsUnitCtx = struct {
    x: u64,
    b: u64,
    primes: []const u64,
    primes2: []const u64,
};

fn isUnit(ctx: IsUnitCtx, p_a: u64) bool {
    const x = ctx.x;
    const b = ctx.b;
    const p_b = ctx.primes[b];
    const p_bm1 = ctx.primes[b - 1];

    const x_o_p_a_o_p_b = x / p_a / p_b;
    assert(x_o_p_a_o_p_b == x / p_b / p_a);

    return x_o_p_a_o_p_b <= p_bm1;
}

fn findFirstNonUnitIndex(x: u64, b: u64, a_limit: u64, primes: []const u64, primes2: []const u64) usize {
    const index = findFirstTrueIndex(u64, primes2[b + 1 .. a_limit], IsUnitCtx{ .x = x, .b = b, .primes = primes, .primes2 = primes2 }, isUnit);

    return b + 1 + index;
}

fn to_f(i: u64) f64 {
    return @floatFromInt(i);
}

fn pc(v: u64, n: u64) f64 {
    return to_f(v) / to_f(n) * 100.0;
}

pub fn main(init: std.process.Init) !void {
    var args = try init.minimal.args.iterateAllocator(init.gpa);
    _ = args.skip();
    const x_str = args.next() orelse {
        std.debug.print("Usage: zig run play5.zig -- <x>\n", .{});
        return;
    };
    const x = std.fmt.parseInt(u64, x_str, 10) catch |err| {
        std.debug.print("Error: '{s}' is not a valid u64 number (Error: {any})\n", .{ x_str, err });
        return;
    };

    const y = isqrt(x);
    const z = icbrt(x);

    const gpa = std.heap.page_allocator;

    const primes2 = try rs.basePrimes(gpa, y);
    defer gpa.free(primes2);

    const primes = try rs.basePrimes(gpa, z);
    defer gpa.free(primes);

    std.debug.print("x: {d:>16} | y = x^1/2: {d:>8} | pi(y): {d:>8} | z = x^1/3: {d:>8} | pi(z): {d:>8}\n\n", .{ x, y, primes2.len, z, primes.len });

    var i: usize = 0;
    while (i < 5 and i < primes2.len) {
        const p_a = primes2[primes2.len - 1 - i];
        std.debug.print("    x / {d:>8} : {d:>8}\n", .{ p_a, x / p_a });

        var j: usize = 0;
        while (j < 1 and j < primes.len) {
            const p_b = primes[primes.len - 1 - j];
            std.debug.print("    x / {d:>8} / {d:>8} : {d:>8}\n", .{ p_a, p_b, x / p_a / p_b });

            j += 1;
        }

        i += 1;
    }

    {
        i = 0;
        const p_b = primes[primes.len - 1];
        while (i < primes2.len and x / primes2[primes2.len - 1 - i] / p_b < p_b) {
            i += 1;
        }
        std.debug.print("\n\n{d:>8} out of {d:>8} p <= x^1/2 have 2nd order phi = 1 --- {d:>6.2}%\n", .{ i, primes2.len, pc(i, primes2.len) });

        const lim = i + 5;
        while (i < lim and i < primes2.len) {
            const p_a = primes2[primes2.len - 1 - i];
            std.debug.print("    x / {d:>8} / {d:>8} : {d:>8}\n", .{ p_a, p_b, x / p_a / p_b });

            i += 1;
        }
    }

    std.debug.print("\n", .{});

    var total_nodes: u64 = 0;
    var total_ones: u64 = 0;
    var total_pis: u64 = 0;
    var total_pi_dups: u64 = 0;
    var total_pi_dups2: u64 = 0;
    var total_rest: u64 = 0;

    var b: usize = primes.len - 1;
    var loop_no: usize = 0;

    while (true) {
        const ones_limit2 = findFirstNonUnitIndex(x, b, primes2.len, primes, primes2);

        const p_b = primes[b];
        const p_bm1 = primes[b - 1];
        const total = primes2.len - 1 - b;

        var a = primes2.len - 1;

        while (b < a) {
            const p_a = primes2[a];

            const x_o_p_a_o_p_b = x / p_a / p_b;
            assert(x_o_p_a_o_p_b == x / p_b / p_a);
            if (x_o_p_a_o_p_b > p_bm1) {
                break;
            }

            a -= 1;
        }
        const ones_limit = a;

        var last_x_o_p_a_o_p_b: u64 = 0;
        var pi_dups_count: u64 = 0;
        var pi_dups2_count: u64 = 0;

        while (b < a) {
            const p_a = primes2[a];

            const x_o_p_a_o_p_b = x / p_a / p_b;
            assert(x_o_p_a_o_p_b == x / p_b / p_a);
            if (x_o_p_a_o_p_b >= p_bm1 * p_bm1) {
                break;
            }

            if (last_x_o_p_a_o_p_b != 0 and last_x_o_p_a_o_p_b == x_o_p_a_o_p_b) {
                pi_dups_count += 1;
            }
            if (last_x_o_p_a_o_p_b != 0 and (last_x_o_p_a_o_p_b == x_o_p_a_o_p_b or (x_o_p_a_o_p_b % 2 == 0 and last_x_o_p_a_o_p_b == x_o_p_a_o_p_b - 1))) {
                pi_dups2_count += 1;
            }

            last_x_o_p_a_o_p_b = x_o_p_a_o_p_b;

            a -= 1;
        }
        const pis_limit = a;

        const ones_count = primes2.len - 1 - ones_limit;
        const pis_count = ones_limit - pis_limit;
        const rest_count = pis_limit - b;

        if (loop_no < 5 or b < 5) {
            std.debug.print("        b: {d:>8} | p_b: {d:>8} | total p_a's {d:>8} | 1's: {d:>8} - {d:>6.2}% | pi's: {d:>8} - {d:>6.2}% | rest: {d:>8} - {d:>6.2}%\n", .{ b, p_b, total, ones_count, pc(ones_count, total), pis_count, pc(pis_count, total), rest_count, pc(rest_count, total) });
            std.debug.print("                 ones_limit: {d:>8} ones_limit2 {d:>8}\n", .{ ones_limit, ones_limit2 });
        }

        total_nodes += total;
        total_ones += ones_count;
        total_pis += pis_count;
        total_pi_dups += pi_dups_count;
        total_pi_dups2 += pi_dups2_count;
        total_rest += rest_count;

        if (b == 1) {
            break;
        }
        b -= 1;
        loop_no += 1;
    }

    std.debug.print("\n", .{});
    std.debug.print("x: {d:>16} | pi(x^1/2): {d:>8} | pi(x^1/3): {d:>8} | nodes: {d:>16} | ones: {d:>12} - {d:>6.2}% | pis: {d:>12} - {d:>6.2}% | rest: {d:>12} - {d:>6.2}%\n", .{ x, primes2.len, primes.len, total_nodes, total_ones, pc(total_ones, total_nodes), total_pis, pc(total_pis, total_nodes), total_rest, pc(total_rest, total_nodes) });
    std.debug.print("    pis: {d:>12} | pi-dups: {d:>12} - {d:>6.2}%| pi-dups2: {d:>12} - {d:>6.2}%\n", .{ total_pis, total_pi_dups, pc(total_pi_dups, total_pis), total_pi_dups2, pc(total_pi_dups2, total_pis) });
}
