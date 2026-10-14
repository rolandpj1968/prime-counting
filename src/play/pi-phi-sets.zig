//! RPJ playground: phi as recursion over SETS of leaves, not single leaves.
// zig build-exe -O ReleaseFast --dep rs -Mroot=pi-phi-sets.zig -Mrs=../rangesieve.zig -femit-bin=pi-phi-sets
//
// A set S(P, c, [a0, a1)) stands for  sum_{a0 <= a < a1} phi(u_a, c),
// u_a = floor(x / (P p_a)), strictly non-increasing in a. One step:
//
//   - zero tail  (u_a == 0)             binary search, dropped
//   - one tail   (1 <= u_a < p_{c+1})   binary search, counted in bulk
//   - pi run     (u_a < p_c^2, u_a <= sqrt x)
//                                       phi = pi(u) - c + 1, one lookup per element
//   - the rest   (head of the range)    descend b -> b-1 on the whole range:
//         S(P, c, R) = S(P, c-1, R) - S(P p_c, c-1, R)
//   - c == 0                            sum of u_a, one division per element
//
// So all per-element work is either a pi lookup or a c == 0 division; every
// other element is resolved in bulk by a binary search on its set. Counted:
// sets created, binary searches, element touches by kind. Top level is the
// same expansion as phiBatchTop:
//   phi(x, A) = x - sum_a x/p_a + sum_{b < A} S(p_b, b-1, (b, A]).

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
    tab: u64, // pi table bound L: pi(n) known for n <= L only
    primes: []const u64,
    pis: []const u64,
};

const Stats = struct {
    sets: u64 = 0,
    searches: u64 = 0,
    ones: u64 = 0, // elements resolved in bulk as phi == 1
    zeros: u64 = 0,
    pi_touches: u64 = 0,
    c0_touches: u64 = 0,
    c0_runs: u64 = 0, // distinct-quotient runs actually visited when GROUPED
    pi_runs: u64 = 0,
    p2_elems: u64 = 0, // elements resolved by the P2 formula
    p2_runs: u64 = 0,
    p2_lookups: u64 = 0, // pi(u / p_i) lookups spent on them
    head_elems: u64 = 0, // elements that had to descend
    max_depth: u64 = 0,
    sets_by_depth: [64]u64 = [_]u64{0} ** 64,
    touches_by_depth: [64]u64 = [_]u64{0} ** 64,
};

// first a in [lo, hi) with u_a < bound (u_a non-increasing), else hi
fn firstBelow(ctx: Ctx, P: u64, lo: u64, hi: u64, bound: u64, st: *Stats) u64 {
    st.searches += 1;
    var l = lo;
    var h = hi;
    while (l < h) {
        const m = l + (h - l) / 2;
        if (ctx.x / P / ctx.primes[m] < bound) h = m else l = m + 1;
    }
    return l;
}

// Grouped per-element work. Over a range of primes, q = floor(N / p_a) is
// constant on runs, and the run holding p_a ends at a' = pi(N / q) + 1
// (p <= N/q  <=>  N/p >= q). Every p_a <= sqrt x, so N/q is clamped there
// and the pi table answers it. One O(1) step per distinct quotient instead
// of one division per element.
const GROUPED = false; // fewer ops, but slower wall-clock: pis[] random access (8 MB at 1e12)

fn runEnd(ctx: Ctx, N: u64, q: u64, a1: u64) u64 {
    const lim = N / q;
    if (lim >= ctx.tab) return a1; // only valid for tab == sqrt x (all p_a <= sqrt x)
    return @min(a1, ctx.pis[lim] + 1);
}

// sum_{a0 <= a < a1} floor(N / p_a)
fn sumQuot(ctx: Ctx, N: u64, a0: u64, a1: u64, st: *Stats) u64 {
    var s: u64 = 0;
    var a = a0;
    while (a < a1) {
        const q = N / ctx.primes[a];
        const e = runEnd(ctx, N, q, a1);
        s += q * (e - a);
        st.c0_runs += 1;
        a = e;
    }
    return s;
}

// sum_{a0 <= a < a1} (pi(floor(N / p_a)) + 1 - c); every quotient <= sqrt x
fn sumPiQuot(ctx: Ctx, N: u64, c: u64, a0: u64, a1: u64, st: *Stats) u64 {
    var s: u64 = 0;
    var a = a0;
    while (a < a1) {
        const q = N / ctx.primes[a];
        const e = runEnd(ctx, N, q, a1);
        s += (ctx.pis[q] + 1 - c) * (e - a);
        st.pi_runs += 1;
        a = e;
    }
    return s;
}

// Third case: p^2 <= u < p^3 (p = p_{c+1}) and u <= sqrt x, so every pi it
// needs is in the table:
//   phi(u, c) = pi(u) - c + 1 + sum_{c < i <= pi(sqrt u)} (pi(u/p_i) - i + 1)
const P2_CASE = true;

fn phiP2(ctx: Ctx, u: u64, c: u64, st: *Stats) u64 {
    const k = ctx.pis[isqrt(u)];
    var s: u64 = ctx.pis[u] + 1 - c;
    var i = c + 1;
    while (i <= k) : (i += 1) s += ctx.pis[u / ctx.primes[i]] + 1 - i;
    st.p2_lookups += k - c;
    return s;
}

fn sumP2(ctx: Ctx, N: u64, c: u64, a0: u64, a1: u64, st: *Stats) u64 {
    var s: u64 = 0;
    var a = a0;
    while (a < a1) {
        const q = N / ctx.primes[a];
        const e = if (GROUPED) runEnd(ctx, N, q, a1) else a + 1;
        s += phiP2(ctx, q, c, st) * (e - a);
        st.p2_runs += 1;
        a = e;
    }
    st.p2_elems += a1 - a0;
    return s;
}

// depth = number of primes in P beyond the top-level p_b
fn evalSet(ctx: Ctx, P: u64, c: u64, a0: u64, a1_in: u64, depth: u64, st: *Stats) u64 {
    if (a0 >= a1_in) return 0;
    st.sets += 1;
    st.sets_by_depth[depth] += 1;
    st.max_depth = @max(st.max_depth, depth);

    const a1 = firstBelow(ctx, P, a0, a1_in, 1, st);
    st.zeros += a1_in - a1;
    if (a0 >= a1) return 0;

    if (c == 0) {
        st.c0_touches += a1 - a0;
        if (GROUPED) {
            const s = sumQuot(ctx, ctx.x / P, a0, a1, st);
            st.touches_by_depth[depth] += a1 - a0;
            return s;
        }
        var s: u64 = 0;
        var a = a0;
        while (a < a1) : (a += 1) s += ctx.x / P / ctx.primes[a];
        st.touches_by_depth[depth] += a1 - a0;
        return s;
    }

    const p_c = ctx.primes[c];
    const p = ctx.primes[c + 1];
    const a_one = firstBelow(ctx, P, a0, a1, p, st);
    var sum: u64 = a1 - a_one;
    st.ones += a1 - a_one;

    // everything below is table-bound: u <= L
    const a_root = firstBelow(ctx, P, a0, a_one, ctx.tab + 1, st);
    // [a_pi, a_one): u < p^2  =>  phi = pi(u) - c + 1
    const a_pi = @max(a_root, firstBelow(ctx, P, a0, a_one, p *| p, st));
    // [a_p2, a_pi): p^2 <= u < p^3  =>  P2 formula
    const a_p2 = if (P2_CASE) @max(a_root, firstBelow(ctx, P, a0, a_pi, p *| p *| p, st)) else a_pi;

    if (GROUPED) {
        sum += sumPiQuot(ctx, ctx.x / P, c, a_pi, a_one, st);
    } else {
        var a = a_pi;
        while (a < a_one) : (a += 1) sum += ctx.pis[ctx.x / P / ctx.primes[a]] + 1 - c;
    }
    st.pi_touches += a_one - a_pi;
    st.touches_by_depth[depth] += a_one - a_pi;
    if (a_p2 < a_pi) sum += sumP2(ctx, ctx.x / P, c, a_p2, a_pi, st);

    // head [a0, a_p2): descend one prime
    if (a0 < a_p2) {
        st.head_elems += a_p2 - a0;
        sum += evalSet(ctx, P, c - 1, a0, a_p2, depth, st);
        if (P <= ctx.x / p_c) sum -= evalSet(ctx, P * p_c, c - 1, a0, a_p2, depth + 1, st);
    }
    return sum;
}

fn piSets(gpa: std.mem.Allocator, x: u64, tab: u64, st: *Stats) !u64 {
    const sqrt_x = isqrt(x);
    const primes = try genPrimesTo(gpa, sqrt_x);
    defer gpa.free(primes);
    const pis = try genPisTo(gpa, tab, primes);
    defer gpa.free(pis);
    const A = primes.len - 1;
    // primes[A + 1] is read as p_{c+1} when c == A - 1; extend by one sentinel
    const ext = try gpa.alloc(u64, primes.len + 1);
    defer gpa.free(ext);
    @memcpy(ext[0..primes.len], primes);
    ext[primes.len] = std.math.maxInt(u64);
    const ctx = Ctx{ .x = x, .sqrt_x = sqrt_x, .tab = tab, .primes = ext, .pis = pis };

    var phi: u64 = x;
    var a: u64 = 1;
    while (a <= A) : (a += 1) phi -= x / primes[a];
    var b: u64 = 1;
    while (b < A) : (b += 1) phi += evalSet(ctx, primes[b], b - 1, b + 1, A + 1, 0, st);
    return phi + A - 1;
}

const KNOWN = [_]u64{ 0, 4, 25, 168, 1229, 9592, 78498, 664579, 5761455, 50847534, 455052511, 4118054813, 37607912018, 346065536839 };

pub fn main(init: std.process.Init) !void {
    var args = try init.minimal.args.iterateAllocator(init.gpa);
    _ = args.skip();
    const x_str = args.next() orelse return std.debug.print("Usage: pi-phi-sets <x>\n", .{});
    const x = try std.fmt.parseInt(u64, x_str, 10);
    // optional pi-table exponent e: L = x^e, clamped to [1, sqrt x]
    const e: f64 = if (args.next()) |s| try std.fmt.parseFloat(f64, s) else 0.5;
    const sqrt_x = isqrt(x);
    var tab: u64 = @intFromFloat(@floor(std.math.pow(f64, to_f(x), e)));
    tab = @max(1, @min(tab, sqrt_x));
    if (GROUPED and tab != sqrt_x) return std.debug.print("GROUPED needs the full sqrt x table\n", .{});
    var st = Stats{};
    const pi_x = try piSets(std.heap.page_allocator, x, tab, &st);

    var check: []const u8 = "";
    var p10: u64 = 1;
    for (KNOWN) |k| {
        if (p10 == x) check = if (k == pi_x) "  MATCH" else "  *** MISMATCH ***";
        p10 *%= 10;
    }
    const x23 = std.math.pow(f64, to_f(x), 2.0 / 3.0);
    const touches = st.pi_touches + st.c0_touches;
    std.debug.print("pi({d}) = {d}{s}\n", .{ x, pi_x, check });
    std.debug.print("pi table L = {d} (x^{d:.3}), prime list to sqrt x = {d}\n", .{ tab, e, sqrt_x });
    std.debug.print("sets {d} ({d:.3} x^2/3)  searches {d}  touches {d} ({d:.3} x^2/3: pi {d}, c0 {d})  bulk ones {d}  zeros {d}  max depth {d}\n", .{ st.sets, to_f(st.sets) / x23, st.searches, touches, to_f(touches) / x23, st.pi_touches, st.c0_touches, st.ones, st.zeros, st.max_depth });
    if (GROUPED) std.debug.print("grouped: runs {d} ({d:.3} x^2/3: pi {d}, c0 {d}) vs {d} element touches = {d:.1}x fewer\n", .{ st.pi_runs + st.c0_runs, to_f(st.pi_runs + st.c0_runs) / x23, st.pi_runs, st.c0_runs, touches, to_f(touches) / to_f(st.pi_runs + st.c0_runs) });
    std.debug.print("P2 case: {d} elements in {d} runs, {d} pi lookups ({d:.3} x^2/3)   head elements descended: {d}\n", .{ st.p2_elems, st.p2_runs, st.p2_lookups, to_f(st.p2_lookups) / x23, st.head_elems });
    std.debug.print("depth:  sets / touches\n", .{});
    var d: usize = 0;
    while (d <= st.max_depth) : (d += 1) std.debug.print("  {d:>2}: {d:>14} / {d:>14}\n", .{ d, st.sets_by_depth[d], st.touches_by_depth[d] });
}
