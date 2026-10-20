//! RPJ playground: phi via level-synchronous residual pools.
// zig build-exe -O ReleaseFast --dep rs -Mroot=pi-phi-pool.zig -Mrs=../rangesieve.zig -femit-bin=pi-phi-pool
//
// Instead of recursing into each set's head (pi-phi-sets), record the head's
// residual leaves phi(u, c), u = floor(x/m), in a pool for level c, and let
// leaves from DIFFERENT denominators m that share the same u combine: the pool
// at level c holds one signed coefficient k(u) per distinct u. One step down:
//
//   k(u, c-1) += k(u, c),   k(floor(u/p_c), c-1) -= k(u, c)
//
// (floor(floor(x/m)/p) = floor(x/(m p)), so every u stays a floor quotient of
// x and indexes densely: u <= sqrt x -> lo[u], else hi[x/u].) Before stepping,
// each pool entry is resolved if it can be, once, with its combined
// coefficient, with p = p_{c+1}:
//   u <  p                 phi = 1
//   u <  p^2, u <= L       phi = pi(u) - c + 1
//   u <  p^3, u <= L       phi = pi(u) - c + 1 + sum_{c<i<=pi(sqrt u)} (pi(u/p_i) - i + 1)
//   c == 0                 phi = u
// The top level is the same set peeling as pi-phi-sets, with heads fed into
// the pool at their level instead of recursed.

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

const Ctx = struct {
    x: u64,
    sqrt_x: u64,
    tab: u64,
    primes: []const u64,
    pis: []const u64,
};

const Stats = struct {
    inserted: u64 = 0, // head elements fed into pools by the top level
    live: u64 = 0, // sum over levels of distinct live u
    steps: u64 = 0, // entries stepped down a level
    cancelled: u64 = 0, // entries whose combined coefficient hit 0
    r_one: u64 = 0,
    r_pi: u64 = 0,
    r_p2: u64 = 0,
    r_c0: u64 = 0,
    p2_lookups: u64 = 0,
    top_pi: u64 = 0,
    top_p2: u64 = 0,
    max_live: u64 = 0,
    // representation census (REP_STATS), summed over levels
    rep_live: u64 = 0, // nonzero entries at level start
    rep_adj_runs: u64 = 0, // maximal runs of consecutive floor-quotient indices
    rep_const_runs: u64 = 0, // ... further split where the count changes
    rep_unit: u64 = 0, // |k| == 1
    rep_small: u64 = 0, // 2 <= |k| <= 8
    rep_max_k: u64 = 0,
    rep_big_u: u64 = 0, // entries with u > sqrt x (the "hi" half)
    // node census (NODE_STATS)
    top_nodes: u64 = 0, // distinct (u, c) resolved inside top-level sets
    uniq_u: u64 = 0, // distinct u over all levels (<= 2 sqrt x + 1)
    band_arr: [64]u64 = [_]u64{0} ** 64, // pool arrivals at level c, by log2 p_c
    band_live: [64]u64 = [_]u64{0} ** 64,
    band_step: [64]u64 = [_]u64{0} ** 64,
    runs: u64 = 0, // per u, maximal runs of consecutive levels with the same nonzero count
};

// Node census: distinct phi nodes, distinct u, and where arrivals collide.
const NODE_STATS = false; // instrumentation only
var seen: []bool = &.{};
var prev_k: []i64 = &.{};
var prev_c: []u64 = &.{};

fn markSeen(ctx: Ctx, u: u64, st: *Stats) void {
    const i = qIndex(ctx, u);
    if (!seen[i]) {
        seen[i] = true;
        st.uniq_u += 1;
    }
}

fn band(ctx: Ctx, c: u64) usize {
    return if (c == 0) 0 else std.math.log2_int(u64, ctx.primes[c]);
}

// Representation census: is the per-level count vector compressible?
const REP_STATS = false; // sorts every level: instrumentation only

// position of u in the increasing floor-quotient order
fn qIndex(ctx: Ctx, u: u64) u64 {
    return if (u <= ctx.sqrt_x) u else 2 * ctx.sqrt_x + 1 - ctx.x / u;
}

fn repCensus(gpa: std.mem.Allocator, ctx: Ctx, pool: *Pool, st: *Stats) !void {
    var es: std.ArrayList([2]i64) = .empty;
    defer es.deinit(gpa);
    for (pool.touched.items) |u| {
        const k = pool.slot(ctx, u).*;
        if (k != 0) try es.append(gpa, .{ @intCast(qIndex(ctx, u)), k });
        if (k != 0 and u > ctx.sqrt_x) st.rep_big_u += 1;
    }
    std.mem.sort([2]i64, es.items, {}, struct {
        fn lt(_: void, a: [2]i64, b: [2]i64) bool {
            return a[0] < b[0];
        }
    }.lt);
    for (es.items, 0..) |e, i| {
        const ak = @abs(e[1]);
        if (ak == 1) st.rep_unit += 1 else if (ak <= 8) st.rep_small += 1;
        st.rep_max_k = @max(st.rep_max_k, ak);
        const adj = i > 0 and es.items[i - 1][0] + 1 == e[0];
        if (!adj) st.rep_adj_runs += 1;
        if (!adj or es.items[i - 1][1] != e[1]) st.rep_const_runs += 1;
    }
    st.rep_live += es.items.len;
}

// Dense pool over the floor-quotient set, with a touched list.
const Pool = struct {
    lo: []i64,
    hi: []i64,
    touched: std.ArrayList(u64), // u values with a nonzero-at-some-point slot

    fn init(gpa: std.mem.Allocator, sqrt_x: u64) !Pool {
        const lo = try gpa.alloc(i64, sqrt_x + 1);
        @memset(lo, 0);
        const hi = try gpa.alloc(i64, sqrt_x + 2);
        @memset(hi, 0);
        return .{ .lo = lo, .hi = hi, .touched = .empty };
    }

    fn slot(pool: *Pool, ctx: Ctx, u: u64) *i64 {
        return if (u <= ctx.sqrt_x) &pool.lo[u] else &pool.hi[ctx.x / u];
    }

    fn add(pool: *Pool, gpa: std.mem.Allocator, ctx: Ctx, u: u64, k: i64) !void {
        const s = pool.slot(ctx, u);
        if (s.* == 0) try pool.touched.append(gpa, u);
        s.* += k;
    }
};

fn firstBelow(ctx: Ctx, P: u64, lo: u64, hi: u64, bound: u64) u64 {
    var l = lo;
    var h = hi;
    while (l < h) {
        const m = l + (h - l) / 2;
        if (ctx.x / P / ctx.primes[m] < bound) h = m else l = m + 1;
    }
    return l;
}

fn phiP2(ctx: Ctx, u: u64, c: u64, st: *Stats) u64 {
    const k = ctx.pis[isqrt(u)];
    var s: u64 = ctx.pis[u] + 1 - c;
    var i = c + 1;
    while (i <= k) : (i += 1) s += ctx.pis[u / ctx.primes[i]] + 1 - i;
    st.p2_lookups += k - c;
    return s;
}

// Resolve phi(u, c) in closed form if a case applies.
fn resolve(ctx: Ctx, u: u64, c: u64, st: *Stats) ?u64 {
    if (c == 0) {
        st.r_c0 += 1;
        return u;
    }
    const p = ctx.primes[c + 1];
    if (u < p) {
        st.r_one += 1;
        return if (u == 0) 0 else 1;
    }
    if (u > ctx.tab) return null;
    if (u < p *| p) {
        st.r_pi += 1;
        return ctx.pis[u] + 1 - c;
    }
    if (u < p *| p *| p) {
        st.r_p2 += 1;
        return phiP2(ctx, u, c, st);
    }
    return null;
}

// Top-level set S(p_b, c = b-1, (b, A]): peel tails, feed head into pool.
fn topSet(gpa: std.mem.Allocator, ctx: Ctx, b: u64, A: u64, pool: *Pool, st: *Stats) !i64 {
    const c = b - 1;
    const P = ctx.primes[b];
    const a0 = b + 1;
    const a1 = A + 1;
    if (a0 >= a1) return 0;
    if (c == 0) {
        var s: i64 = 0;
        var a = a0;
        while (a < a1) : (a += 1) s += @intCast(ctx.x / P / ctx.primes[a]);
        st.r_c0 += a1 - a0;
        return s;
    }
    const p = ctx.primes[c + 1];
    const a_one = firstBelow(ctx, P, a0, a1, p);
    var sum: i64 = @intCast(a1 - a_one);
    const a_root = firstBelow(ctx, P, a0, a_one, ctx.tab + 1);
    const a_pi = @max(a_root, firstBelow(ctx, P, a0, a_one, p *| p));
    const a_p2 = @max(a_root, firstBelow(ctx, P, a0, a_pi, p *| p *| p));
    var a = a_pi;
    var last_u: u64 = 0;
    while (a < a_one) : (a += 1) {
        const u = ctx.x / P / ctx.primes[a];
        if (NODE_STATS and u != last_u) {
            st.top_nodes += 1;
            markSeen(ctx, u, st);
            last_u = u;
        }
        sum += @intCast(ctx.pis[u] + 1 - c);
    }
    st.top_pi += a_one - a_pi;
    a = a_p2;
    last_u = 0;
    while (a < a_pi) : (a += 1) {
        const u = ctx.x / P / ctx.primes[a];
        if (NODE_STATS and u != last_u) {
            st.top_nodes += 1;
            markSeen(ctx, u, st);
            last_u = u;
        }
        sum += @intCast(phiP2(ctx, u, c, st));
    }
    st.top_p2 += a_pi - a_p2;
    a = a0;
    while (a < a_p2) : (a += 1) try pool.add(gpa, ctx, ctx.x / P / ctx.primes[a], 1);
    st.inserted += a_p2 - a0;
    if (NODE_STATS) st.band_arr[band(ctx, c)] += a_p2 - a0;
    return sum;
}

fn piPool(gpa: std.mem.Allocator, x: u64, tab: u64, st: *Stats) !u64 {
    const sqrt_x = isqrt(x);
    const primes = try genPrimesTo(gpa, sqrt_x);
    defer gpa.free(primes);
    const pis = try genPisTo(gpa, tab, primes);
    defer gpa.free(pis);
    const A = primes.len - 1;
    const ext = try gpa.alloc(u64, primes.len + 1);
    defer gpa.free(ext);
    @memcpy(ext[0..primes.len], primes);
    ext[primes.len] = std.math.maxInt(u64);
    const ctx = Ctx{ .x = x, .sqrt_x = sqrt_x, .tab = tab, .primes = ext, .pis = pis };

    var cur = try Pool.init(gpa, sqrt_x);
    var nxt = try Pool.init(gpa, sqrt_x);
    if (NODE_STATS) {
        seen = try gpa.alloc(bool, 2 * sqrt_x + 2);
        @memset(seen, false);
        prev_k = try gpa.alloc(i64, 2 * sqrt_x + 2);
        @memset(prev_k, 0);
        prev_c = try gpa.alloc(u64, 2 * sqrt_x + 2);
        @memset(prev_c, std.math.maxInt(u64));
    }

    var total: i64 = @intCast(x);
    var a: u64 = 1;
    while (a <= A) : (a += 1) total -= @intCast(x / primes[a]);

    // top-level sets at c = b - 1 >= pi(cbrt x) are all-ones: count directly
    const cbrt_x = icbrt(x);
    var pc3: u64 = 0;
    while (pc3 + 1 <= A and primes[pc3 + 1] <= cbrt_x) pc3 += 1;
    var b: u64 = A - 1;
    while (b > pc3) : (b -= 1) total += @intCast(A - b); // phi(u, b-1) = 1 for all a > b
    var c: u64 = pc3; // highest level that can hold a pool entry: b = pc3 + 1

    while (true) {
        // top-level set b = c + 1 feeds level c
        if (c + 1 <= pc3) total += try topSet(gpa, ctx, c + 1, A, &cur, st);

        if (REP_STATS) try repCensus(gpa, ctx, &cur, st);
        st.live += cur.touched.items.len;
        st.max_live = @max(st.max_live, cur.touched.items.len);
        if (NODE_STATS) {
            st.band_live[band(ctx, c)] += cur.touched.items.len;
            for (cur.touched.items) |u| markSeen(ctx, u, st);
        }
        for (cur.touched.items) |u| {
            const s = cur.slot(ctx, u);
            const k = s.*;
            s.* = 0;
            if (k == 0) {
                st.cancelled += 1;
                continue;
            }
            if (NODE_STATS) {
                const i = qIndex(ctx, u);
                if (!(prev_k[i] == k and prev_c[i] == c + 1)) st.runs += 1;
                prev_k[i] = k;
                prev_c[i] = c;
            }
            if (resolve(ctx, u, c, st)) |phi| {
                total += k * @as(i64, @intCast(phi));
                continue;
            }
            // c >= 1 here (c == 0 always resolves)
            st.steps += 1;
            try nxt.add(gpa, ctx, u, k);
            const v = u / ctx.primes[c];
            if (v > 0) try nxt.add(gpa, ctx, v, -k);
            if (NODE_STATS) {
                st.band_step[band(ctx, c)] += 1;
                st.band_arr[band(ctx, c - 1)] += if (v > 0) 2 else 1;
            }
        }
        cur.touched.clearRetainingCapacity();
        if (c == 0) break;
        std.mem.swap(Pool, &cur, &nxt);
        c -= 1;
    }
    return @as(u64, @intCast(total)) + A - 1;
}

const KNOWN = [_]u64{ 0, 4, 25, 168, 1229, 9592, 78498, 664579, 5761455, 50847534, 455052511, 4118054813, 37607912018, 346065536839 };

pub fn main(init: std.process.Init) !void {
    var args = try init.minimal.args.iterateAllocator(init.gpa);
    _ = args.skip();
    const x_str = args.next() orelse return std.debug.print("Usage: pi-phi-pool <x> [table exponent, default 0.5]\n", .{});
    const x = try std.fmt.parseInt(u64, x_str, 10);
    const e: f64 = if (args.next()) |s| try std.fmt.parseFloat(f64, s) else 0.5;
    const sqrt_x = isqrt(x);
    var tab: u64 = @intFromFloat(@floor(std.math.pow(f64, to_f(x), e)));
    tab = @max(1, @min(tab, sqrt_x));
    var st = Stats{};
    const pi_x = try piPool(std.heap.page_allocator, x, tab, &st);

    var check: []const u8 = "";
    var p10: u64 = 1;
    for (KNOWN) |k| {
        if (p10 == x) check = if (k == pi_x) "  MATCH" else "  *** MISMATCH ***";
        p10 *%= 10;
    }
    const x23 = std.math.pow(f64, to_f(x), 2.0 / 3.0);
    std.debug.print("pi({d}) = {d}{s}\n", .{ x, pi_x, check });
    std.debug.print("L = {d}  inserted {d}  live {d} ({d:.3} x^2/3, max/level {d})  steps {d}  cancelled {d}\n", .{ tab, st.inserted, st.live, to_f(st.live) / x23, st.max_live, st.steps, st.cancelled });
    if (REP_STATS) std.debug.print("rep: live {d}  adjacent-runs {d} ({d:.2}/run)  constant-count runs {d} ({d:.2}/run)  |k|=1 {d:.1}%  2..8 {d:.1}%  max |k| {d}  u > sqrt x {d:.1}%\n", .{ st.rep_live, st.rep_adj_runs, to_f(st.rep_live) / to_f(st.rep_adj_runs), st.rep_const_runs, to_f(st.rep_live) / to_f(st.rep_const_runs), pc(st.rep_unit, st.rep_live), pc(st.rep_small, st.rep_live), st.rep_max_k, pc(st.rep_big_u, st.rep_live) });
    if (NODE_STATS) {
        const sq = to_f(sqrt_x);
        const nodes = st.live + st.top_nodes;
        std.debug.print("nodes: distinct (u,c) {d} = pool {d} + top-resolved {d}  ({d:.3} sqrt x, {d:.4} x^2/3)  stepped (interior) {d} ({d:.3} sqrt x)  distinct u {d} ({d:.3} sqrt x of max 2)\n", .{ nodes, st.live, st.top_nodes, to_f(nodes) / sq, to_f(nodes) / x23, st.steps, to_f(st.steps) / sq, st.uniq_u, to_f(st.uniq_u) / sq });
        std.debug.print("  coefficient runs over levels: {d} ({d:.3} sqrt x) vs live nodes {d}: {d:.1} levels per run\n", .{ st.runs, to_f(st.runs) / sq, st.live, to_f(st.live) / @max(1.0, to_f(st.runs)) });
        std.debug.print("  log2 p_c |   arrivals       live   arr/live    stepped\n", .{});
        for (0..64) |k| {
            if (st.band_arr[k] == 0 and st.band_live[k] == 0) continue;
            std.debug.print("  {d:>8} | {d:>10} {d:>10} {d:>9.2}x {d:>10}\n", .{ k, st.band_arr[k], st.band_live[k], to_f(st.band_arr[k]) / @max(1.0, to_f(st.band_live[k])), st.band_step[k] });
        }
    }
    std.debug.print("resolved: one {d}  pi {d}  p2 {d}  c0 {d}   p2 lookups {d} ({d:.3} x^2/3)   top: pi {d} p2 {d}\n", .{ st.r_one, st.r_pi, st.r_p2, st.r_c0, st.p2_lookups, to_f(st.p2_lookups) / x23, st.top_pi, st.top_p2 });
}
