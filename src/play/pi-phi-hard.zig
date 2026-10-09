//! RPJ playground
// zig build-exe -O ReleaseFast --dep rs -Mroot=pi-phi-hard.zig -Mrs=../rangesieve.zig -femit-bin=pi-phi-hard

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
    phi_batches: u64 = 0,
    phi_batch_1s: u64 = 0,
    phi_batch_others: u64 = 0,
    phi_batch_uniques: u64 = 0,
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

// ---------------------------------------------------------------------------
// Hard-leaf census.
//
// The top-level batch expands phi(x, pi(sqrt x)) into level-2 leaves
//   phi(u, c),  u = x / (p_a * p_b),  c = b - 1,  b <= pi(cbrt x) < a <= pi(sqrt x)
// (b > pi(cbrt x) pairs are all phi == 1 and counted in closed form by
// phiBatchTop). Each leaf is classified by u against p = p_b, DR-style:
//
//   C0     c == 0                       phi = u
//   ONE    u <  p                       phi = 1
//   EASY_T p <= u < p^2, u <= sqrt x    phi = pi(u) - c + 1, pi from table
//   EASY_X p <= u < p^2, u >  sqrt x    phi = pi(u) - c + 1, pi NOT in table
//   H2     p^2 <= u < p^3               phi = pi(u) - c + 1 + P2(u, c): pi values only
//   H3     u >= p^3                     needs genuine phi recursion
//
// Classes are evaluated as the playground does today (phiSimple with the
// pi(sqrt u) collapse) so "work" = phiSimple calls actually spent per class,
// and the sum is checked against pi(x).
//
// For the analytic candidates (EASY_X, H2, H3) we also price an HKM-style
// per-query correction interval [u, u + S(u)], S(u) = 0.18 sqrt(u) ln(u) log2(u)
// (the constant measured on main, docs/COMBINATORIAL.md R1), and measure how
// much consecutive leaves at the same level b share: sum S vs union length.
// ---------------------------------------------------------------------------

const Class = enum { C0, ONE, EASY_T, EASY_X, H2, H3 };
const NCLASS = @typeInfo(Class).@"enum".fields.len;

const ClassStats = struct {
    leaves: u64 = 0,
    distinct: u64 = 0,
    work: u64 = 0,
    sum_sqrt_u: f64 = 0,
    max_u: u64 = 0,
    sum_s: f64 = 0,
    union_s: f64 = 0,
    pi_queries: u64 = 0, // H2: pi(u/p_i) terms in P2(u, c)
};

const NBAND = 64;
const Band = struct {
    nb: u64 = 0,
    leaves: [NCLASS]u64 = [_]u64{0} ** NCLASS,
    work: [NCLASS]u64 = [_]u64{0} ** NCLASS,
};

const Census = struct {
    cls: [NCLASS]ClassStats = [_]ClassStats{.{}} ** NCLASS,
    band: [NBAND]Band = [_]Band{.{}} ** NBAND,
};

fn corrWidth(u: u64) f64 {
    const uf = to_f(u);
    if (u < 3) return 0;
    return 0.18 * @sqrt(uf) * @log(uf) * std.math.log2(uf);
}

fn classify(ctx: PhiCtx, u: u64, b: u64) Class {
    const p = ctx.primes[b];
    if (b == 1) return .C0;
    if (u < p) return .ONE;
    if (u / p < p) return if (u <= ctx.sqrt_x) .EASY_T else .EASY_X;
    if (u / p / p < p) return .H2;
    return .H3;
}

fn evalLeaf(ctx: PhiCtx, u: u64, c: u64, cls: Class, stats: *PhiStats) u64 {
    switch (cls) {
        .C0 => return u,
        .ONE => return 1,
        .EASY_T => return ctx.pis[u] + 1 - c,
        else => {
            const k = ctx.pis[isqrt(u)];
            if (k < c) return phiSimple(ctx, u, k, stats) + k - c;
            return phiSimple(ctx, u, c, stats);
        },
    }
}

const Open = struct { lo: u64 = 0, hi: f64 = 0, active: bool = false };

// all level-2 leaves at level b, a in (b, pi(sqrt x)]
fn censusBatch(ctx: PhiCtx, b: u64, census: *Census, stats: *PhiStats) u64 {
    const c = b - 1;
    const p_b = ctx.primes[b];
    const band = &census.band[std.math.log2_int(u64, p_b)];
    band.nb += 1;

    var open = [_]Open{.{}} ** NCLASS;
    var sum: u64 = 0;
    var last_u: u64 = 0;
    var last_phi: u64 = 0;
    var a = b + 1;
    while (a <= ctx.pi_sqrt_x) : (a += 1) {
        const u = ctx.x / ctx.primes[a] / p_b;
        const cls = classify(ctx, u, b);
        const ci = @intFromEnum(cls);
        const cs = &census.cls[ci];
        cs.leaves += 1;
        band.leaves[ci] += 1;
        if (u == last_u) {
            sum += last_phi;
            continue;
        }
        cs.distinct += 1;
        cs.max_u = @max(cs.max_u, u);
        cs.sum_sqrt_u += @sqrt(to_f(u));

        const w0 = stats.phi_simples;
        const phi = evalLeaf(ctx, u, c, cls, stats);
        cs.work += stats.phi_simples - w0;
        band.work[ci] += stats.phi_simples - w0;
        if (cls == .H2) cs.pi_queries += ctx.pis[isqrt(u)] - c;

        if (cls == .EASY_X or cls == .H2 or cls == .H3) {
            // u descends with a: merge [u, u+S] into the open run if it reaches it
            const s = corrWidth(u);
            cs.sum_s += s;
            const o = &open[ci];
            if (o.active and to_f(u) + s >= to_f(o.lo)) {
                o.lo = u;
            } else {
                if (o.active) cs.union_s += o.hi - to_f(o.lo);
                o.* = .{ .lo = u, .hi = to_f(u) + s, .active = true };
            }
        }
        last_u = u;
        last_phi = phi;
        sum += phi;
    }
    for (&open, 0..) |*o, ci| {
        if (o.active) census.cls[ci].union_s += o.hi - to_f(o.lo);
    }
    return sum;
}

fn piCensus(gpa: std.mem.Allocator, x: u64) !u64 {
    const sqrt_x = isqrt(x);
    const cbrt_x = icbrt(x);
    const primes = try genPrimesTo(gpa, sqrt_x);
    defer gpa.free(primes);
    const pis = try genPisTo(gpa, sqrt_x, primes);
    defer gpa.free(pis);
    const pi_sqrt_x = primes.len - 1;
    const pi_cbrt_x = pis[cbrt_x];

    const ctx = PhiCtx{ .x = x, .sqrt_x = sqrt_x, .pi_sqrt_x = pi_sqrt_x, .cbrt_x = cbrt_x, .pi_cbrt_x = pi_cbrt_x, .primes = primes, .pis = pis };

    var stats = PhiStats{};
    const census = try gpa.create(Census);
    defer gpa.destroy(census);
    census.* = .{};

    var phi = x + (pi_sqrt_x - pi_cbrt_x - 1) * (pi_sqrt_x - pi_cbrt_x) / 2;
    var base: u64 = 0;
    var a: u64 = 1;
    while (a <= pi_sqrt_x) : (a += 1) base += x / primes[a];
    var leaf_sum: u64 = 0;
    var b: u64 = 1;
    while (b <= pi_cbrt_x) : (b += 1) leaf_sum += censusBatch(ctx, b, census, &stats);
    phi -= base - leaf_sum;

    const xf = to_f(x);
    const x23 = std.math.pow(f64, xf, 2.0 / 3.0);
    std.debug.print("x = {d}  sqrt x = {d}  x^2/3 = {e:.3}  pi(sqrt x) = {d}  pi(cbrt x) = {d}\n\n", .{ x, sqrt_x, x23, pi_sqrt_x, pi_cbrt_x });

    var tot_leaves: u64 = 0;
    for (census.cls) |cs| tot_leaves += cs.leaves;
    std.debug.print("{s:<7} {s:>12} {s:>7} {s:>12} {s:>12} {s:>10} {s:>10} {s:>10} {s:>10} {s:>10}\n", .{ "class", "leaves", "%", "distinct", "work", "max u", "sum sqrtu", "sum S", "union S", "share" });
    inline for (@typeInfo(Class).@"enum".fields, 0..) |f, ci| {
        const cs = census.cls[ci];
        const share = if (cs.union_s > 0) cs.sum_s / cs.union_s else 0;
        std.debug.print("{s:<7} {d:>12} {d:>6.2}% {d:>12} {d:>12} {e:>10.2} {e:>10.2} {e:>10.2} {e:>10.2} {d:>9.2}x\n", .{ f.name, cs.leaves, pc(cs.leaves, tot_leaves), cs.distinct, cs.work, to_f(cs.max_u), cs.sum_sqrt_u, cs.sum_s, cs.union_s, share });
    }
    std.debug.print("H2 P2 pi-queries: {d}   total work: {d}   (work / x^2/3 = {d:.3})\n\n", .{ census.cls[@intFromEnum(Class.H2)].pi_queries, stats.phi_simples, to_f(stats.phi_simples) / x23 });

    std.debug.print("{s:>6} {s:>6} | {s:>11} {s:>12} | {s:>11} {s:>12} | {s:>11} {s:>12} | {s:>11}\n", .{ "p_b~2^", "#b", "H3 leaves", "H3 work", "H2 leaves", "H2 work", "EX leaves", "EX work", "ET leaves" });
    for (census.band, 0..) |bd, k| {
        if (bd.nb == 0) continue;
        const H3 = @intFromEnum(Class.H3);
        const H2 = @intFromEnum(Class.H2);
        const EX = @intFromEnum(Class.EASY_X);
        const ET = @intFromEnum(Class.EASY_T);
        std.debug.print("{d:>6} {d:>6} | {d:>11} {d:>12} | {d:>11} {d:>12} | {d:>11} {d:>12} | {d:>11}\n", .{ k, bd.nb, bd.leaves[H3], bd.work[H3], bd.leaves[H2], bd.work[H2], bd.leaves[EX], bd.work[EX], bd.leaves[ET] });
    }
    return phi + pi_sqrt_x - 1;
}

const KNOWN = [_]u64{ 0, 4, 25, 168, 1229, 9592, 78498, 664579, 5761455, 50847534, 455052511, 4118054813, 37607912018, 346065536839 };

pub fn main(init: std.process.Init) !void {
    var args = try init.minimal.args.iterateAllocator(init.gpa);
    _ = args.skip();
    const x_str = args.next() orelse return std.debug.print("Usage: pi-phi-hard <x>\n", .{});
    const x = try std.fmt.parseInt(u64, x_str, 10);
    const pi_x = try piCensus(std.heap.page_allocator, x);
    var check: []const u8 = "";
    var p10: u64 = 1;
    for (KNOWN) |k| {
        if (p10 == x) check = if (k == pi_x) "  MATCH" else "  *** MISMATCH ***";
        p10 *%= 10;
    }
    std.debug.print("\npi({d}) = {d}{s}\n", .{ x, pi_x, check });
}
