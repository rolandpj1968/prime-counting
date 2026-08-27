//! RPJ playground
// zig run --dep rs -Mroot=pi-phi-batch.zig -Mrs=../rangesieve.zig -- 100

const std = @import("std");
const assert = std.debug.assert;
const rs = @import("rs");

fn to_f(i: u64) f64 {
    return @floatFromInt(i);
}

fn pc(v: u64, n: u64) f64 {
    return to_f(v) / to_f(n) * 100.0;
}

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

// Generate primes with "natural" phi indices - primes[0] == 1, primes[1] == 2, ...
fn genPrimesTo(gpa: std.mem.Allocator, x: u64) ![]u64 {
    // primes0[0] == 2
    const primes0 = try rs.basePrimes(gpa, x);
    defer gpa.free(primes0);

    // prepend 1 for natural phi param access
    const primes = try gpa.alloc(u64, primes0.len + 1);
    primes[0] = 1;
    @memcpy(primes[1..], primes0);

    // primes[1] == 2
    return primes;
}

fn genPisTo(gpa: std.mem.Allocator, y: u64, primes: []u64) ![]u64 {
    const pis = try gpa.alloc(u64, y + 1);

    var a: u64 = 1;
    var pi: u64 = 0;

    var i: u64 = 0;
    while (i <= y) : (i += 1) {
        if (a < primes.len and primes[a] == i) {
            pi += 1;
            a += 1;
        }
        pis[i] = pi;
    }
    return pis;
}

pub const PhiCtx = struct {
    x: u64,
    sqrt_x: u64,
    pi_sqrt_x: u64,
    cbrt_x: u64,
    pi_cbrt_x: u64,
    primes: []const u64,
    pis: []const u64,
};

pub const PhiStats = struct {
    phi_simples: u64 = 0,
    phi_simple_x0s: u64 = 0,
    phi_simple_a0s: u64 = 0,
};

fn phiSimple(ctx: PhiCtx, x: u64, a: u64, stats: *PhiStats) u64 {
    stats.phi_simples += 1;

    if (x == 0) {
        stats.phi_simple_x0s += 1;
        return 0;
    }

    if (a == 0) {
        stats.phi_simple_a0s += 1;
        return x;
    }

    const p_a = ctx.primes[a];

    return phiSimple(ctx, x, a - 1, stats) - phiSimple(ctx, x / p_a, a - 1, stats);
}

// degenerate case - phi(x,a) == 1
fn isUnitPhi(x: u64, a: u64, primes: []const u64) bool {
    return x > 0 and x <= primes[a];
}

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
fn findFirstTrueIndexBackwards(
    comptime T: type,
    slice: []const T,
    ctx: anytype,
    pred: fn (@TypeOf(ctx), T) bool,
) usize {
    var i: usize = slice.len;
    while (i > 0) {
        i -= 1;
        if (!pred(ctx, slice[i])) {
            return i + 1;
        }
    }
    return 0;
}

// Assume predicate is all-false then all-true over the slice
fn findFirstTrueIndexForwards(
    comptime T: type,
    slice: []const T,
    ctx: anytype,
    pred: fn (@TypeOf(ctx), T) bool,
) usize {
    var i: usize = 0;
    while (i < slice.len) {
        if (pred(ctx, slice[i])) {
            return i;
        }
        i += 1;
    }
    return slice.len;
}

// Assume predicate is all-false then all-true over the slice
fn findFirstTrueIndexBinary(
    comptime T: type,
    slice: []const T,
    ctx: anytype,
    pred: fn (@TypeOf(ctx), T) bool,
) usize {
    var lo: usize = 0;
    var hi: usize = slice.len;

    while (lo != hi) {
        assert(lo < hi);

        const mid = (lo + hi) / 2;

        assert(mid < slice.len);
        assert(lo <= mid);
        assert(mid < hi);

        if (pred(ctx, slice[mid])) {
            hi = mid;
        } else {
            lo = mid + 1;
        }
    }

    assert(lo == hi);
    return lo;
}

const CHECK_FIND_PRED = false;
const CHECK_FIND_INDEX = false;

fn findFirstTrueIndex(
    comptime T: type,
    slice: []const T,
    ctx: anytype,
    pred: fn (@TypeOf(ctx), T) bool,
) usize {
    if (CHECK_FIND_PRED) {
        checkFindFirstTruePred(T, slice, ctx, pred);
    }

    const index = findFirstTrueIndexBinary(T, slice, ctx, pred);

    if (CHECK_FIND_INDEX) {
        const index1 = findFirstTrueIndexBackwards(T, slice, ctx, pred);
        const index2 = findFirstTrueIndexForwards(T, slice, ctx, pred);

        assert(index == index1);
        assert(index == index2);
    }

    return index;
}

const IsUnitCtx = struct {
    x: u64,
    b: u64,
    primes: []const u64,
    primes2: []const u64,
};

fn isUnit(ctx: IsUnitCtx, p_a: u64) bool {
    const x = ctx.x;
    assert(x > 0);
    const b = ctx.b;
    assert(b > 0);
    const p_b = ctx.primes[b];

    const x_o_p_a_o_p_b = x / p_a / p_b;
    assert(x_o_p_a_o_p_b == x / p_b / p_a);
    assert(x_o_p_a_o_p_b == x / (p_b * p_a));

    return isUnitPhi(x_o_p_a_o_p_b, b - 1, ctx.primes);
}

fn findFirstUnitIndex(x: u64, b: u64, a_limit: u64, primes: []const u64, primes2: []const u64) usize {
    const index = findFirstTrueIndex(u64, primes2[b + 1 .. a_limit], IsUnitCtx{ .x = x, .b = b, .primes = primes, .primes2 = primes2 }, isUnit);

    return b + 1 + index;
}

const IsPiCtx = struct {
    x: u64,
    b: u64,
    primes: []const u64,
    primes2: []const u64,
};

fn isPi(ctx: IsPiCtx, p_a: u64) bool {
    const x = ctx.x;
    const b = ctx.b;
    const p_b = ctx.primes[b];
    const p_bm1 = ctx.primes[b - 1];

    const x_o_p_a_o_p_b = x / p_a / p_b;
    assert(x_o_p_a_o_p_b == x / p_b / p_a);
    assert(x_o_p_a_o_p_b == x / (p_b * p_a));

    return x_o_p_a_o_p_b < p_bm1 * p_bm1;
}

fn findFirstPiIndex(x: u64, b: u64, a_limit: u64, primes: []const u64, primes2: []const u64) usize {
    const index = findFirstTrueIndex(u64, primes2[b + 1 .. a_limit], IsPiCtx{ .x = x, .b = b, .primes = primes, .primes2 = primes2 }, isPi);

    return b + 1 + index;
}

// return sum(a_min <= a < a_lim){ phi(x/p_bs/p_a, b)}
//fn phi_batch(x: u64, p_bs: u64, b: u64, primes: []const u64, i_min: u64, i_lim: u64) u64 {}

fn phiBatch1(ctx: PhiCtx, b: u64, a_min: u64, a_lim: u64, stats: *PhiStats) u64 {
    var phi: u64 = 0;

    var a = a_min;
    while (a < a_lim) : (a += 1) {
        const p_a = ctx.primes[a];

        phi += phiSimple(ctx, ctx.x / p_a, b, stats);
    }

    return phi;
}

// return phi(x, pi(sqrt(x)))
fn phiBatchTop(ctx: PhiCtx, stats: *PhiStats) u64 {
    var phi = ctx.x;

    // Gourdan
    const top_phi_ones = (ctx.pi_sqrt_x - ctx.pi_cbrt_x - 1) * (ctx.pi_sqrt_x - ctx.pi_cbrt_x) / 2;
    phi += top_phi_ones;

    // First batch - still need to invert one level
    phi -= phiBatch1(ctx, ctx.pi_cbrt_x, ctx.pi_cbrt_x + 1, ctx.pi_sqrt_x + 1, stats);

    var a = ctx.pi_cbrt_x;
    while (a > 0) : (a -= 1) {
        const p_a = ctx.primes[a];

        phi -= phiSimple(ctx, ctx.x / p_a, a - 1, stats);
    }

    const CHECK_PHI = true;
    if (CHECK_PHI) {
        var dummy = PhiStats{};
        assert(phi == phiSimple(ctx, ctx.x, ctx.pi_sqrt_x, &dummy));
    }

    return phi;
}

fn piBatch(gpa: std.mem.Allocator, x: u64) !u64 {
    const sqrt_x = isqrt(x);
    const cbrt_x = icbrt(x);

    const primes = try genPrimesTo(gpa, sqrt_x);
    defer gpa.free(primes);

    const pis = try genPisTo(gpa, sqrt_x, primes);
    defer gpa.free(pis);

    const pi_sqrt_x = primes.len - 1;
    assert(pis[sqrt_x] == pi_sqrt_x);

    const primes3 = try genPrimesTo(gpa, cbrt_x);
    defer gpa.free(primes3);

    const pi_cbrt_x = primes3.len - 1;
    assert(pis[cbrt_x] == pi_cbrt_x);

    std.debug.print("x: {d:>16} | x^1/2: {d:>8} | pi(x^1/2): {d:>8} | x^1/3: {d:>8} | pi(x^1/3): {d:>8}\n\n", .{ x, sqrt_x, pi_sqrt_x, cbrt_x, pi_cbrt_x });

    std.debug.print("\n", .{});

    const ctx = PhiCtx{
        .x = x,
        .sqrt_x = sqrt_x,
        .pi_sqrt_x = pi_sqrt_x,
        .cbrt_x = cbrt_x,
        .pi_cbrt_x = pi_cbrt_x,
        .primes = primes,
        .pis = pis,
    };

    var stats = PhiStats{};

    const phi = phiBatchTop(ctx, &stats);

    std.debug.print("\n                          {any}\n\n", .{stats});

    return phi + pi_sqrt_x - 1;
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

    const gpa = std.heap.page_allocator;

    const pi_x = try piBatch(gpa, x);

    std.debug.print("\npi({d}) = {d}", .{ x, pi_x });

    const CHECK_PI = true;
    if (CHECK_PI) {
        const primes_to_x = try genPrimesTo(gpa, x);
        defer gpa.free(primes_to_x);
        std.debug.print(" [expecting {d}]", .{primes_to_x.len - 1});
    }
    std.debug.print("\n\n", .{});
}
