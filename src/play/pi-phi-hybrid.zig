//! RPJ playground: hybrid of 3-item sets and per-u count pools.
// zig build-exe -O ReleaseFast --dep rs -Mroot=pi-phi-hybrid.zig -Mrs=../rangesieve.zig -femit-bin=pi-phi-hybrid
//
// Residuals live in two forms, both advanced level by level (c descending):
//
//   set   (P, [a0, a1), sign): sum_a sign * phi(floor(x/(P p_a)), c),
//         u non-increasing in a, so every case below is a binary-searched tail
//   pool  k(u) per u <= T: one signed count per floor quotient, merged across
//         every set/denominator that produced it
//
// At level c (p = p_{c+1}), each set:
//   - drops its zero tail, counts its one tail (u < p) in bulk
//   - sends its tail with u <= T into the pool, one insert per element  <- the combine point
//   - resolves T < u <= L elements with u < p^2 (pi) or u < p^3 (P2) one by one
//   - descends the head as two sets: (P, c-1) and (P p_c, c-1, -sign)
// Then each pool entry is resolved once (one / pi / P2 / c0) or stepped:
//   k(u, c-1) += k,  k(u/p_c, c-1) -= k     (u/p_c <= T: the pool is closed downward)
//
// Usage: pi-phi-hybrid <x> [t = combine exponent, T = x^t] [e = pi-table exponent, L = x^e]

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
    tab: u64, // L
    thr: u64, // T (<= sqrt x)
    primes: []const u64,
    pis: []const u64,
};

const Set = struct { P: u64, a0: u64, a1: u64, sign: i64 };

const Stats = struct {
    sets: u64 = 0,
    set_pi: u64 = 0, // elementwise pi resolutions inside sets
    set_p2: u64 = 0,
    set_c0: u64 = 0, // elementwise c == 0 divisions inside sets
    inserted: u64 = 0, // set elements sent to the pool
    merged: u64 = 0, // ... that landed on an occupied slot
    live: u64 = 0, // pool entries processed, summed over levels
    steps: u64 = 0,
    p2_lookups: u64 = 0,
    max_sets_level: u64 = 0,
    ins_range: u64 = 0, // sum over inserted tails of min(N, u_first - u_last + 1): cost if inserted per u
    ins_tails: u64 = 0,
    dense_tails: u64 = 0, // tails with N > range (more elements than u-values)
    dense_elems: u64 = 0,
    dense_range: u64 = 0,
    // rough-walk census (ROUGH_STATS): per inserted tail (u_last, u_first]
    rw_tails: u64 = 0,
    rw_elems: u64 = 0, // E: elements (current insert cost)
    rw_range: u64 = 0, // R: u-values in (u_last, u_first]
    rw_rough: u64 = 0, // W: c-rough n in (u_last, u_first] (rough-walk lookups)
    rw_union: u64 = 0, // c-rough n in the union of a level's intervals
    rw_pred: f64 = 0, // Mertens: sum range * prod_{p <= p_c} (1 - 1/p)
    rw_levels: u64 = 0, // levels with any pooled tail
    rw_tau_w: u64 = 0, // sum_c tails_c * W-union_c: plain S_c(x/n) query over all tails
    rw_max_tau: u64 = 0,
};

// Rough-walk census: what would walking c-rough n (instead of elements) cost?
const ROUGH_STATS = false; // O(T) per level: instrumentation only

const Pool = struct {
    k: []i64,
    touched: std.ArrayList(u64) = .empty,

    fn add(pool: *Pool, gpa: std.mem.Allocator, u: u64, v: i64) !bool {
        const was = pool.k[u] != 0;
        if (!was) try pool.touched.append(gpa, u);
        pool.k[u] += v;
        return was;
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

// phi(u, c) in closed form if a case applies (pool entries)
fn resolve(ctx: Ctx, u: u64, c: u64, st: *Stats) ?u64 {
    if (c == 0) return u;
    const p = ctx.primes[c + 1];
    if (u < p) return if (u == 0) 0 else 1;
    if (u > ctx.tab) return null;
    if (u < p *| p) return ctx.pis[u] + 1 - c;
    if (u < p *| p *| p) return phiP2(ctx, u, c, st);
    return null;
}

// one set at level c; children for level c-1 appended to `next`
fn stepSet(gpa: std.mem.Allocator, ctx: Ctx, s: Set, c: u64, pool: *Pool, next: *std.ArrayList(Set), ivals: *std.ArrayList([2]u64), st: *Stats) !i64 {
    st.sets += 1;
    const P = s.P;
    const a1 = firstBelow(ctx, P, s.a0, s.a1, 1);
    if (s.a0 >= a1) return 0;
    var sum: i64 = 0;
    if (c == 0) {
        var a = s.a0;
        while (a < a1) : (a += 1) sum += @intCast(ctx.x / P / ctx.primes[a]);
        st.set_c0 += a1 - s.a0;
        return s.sign * sum;
    }
    const p = ctx.primes[c + 1];
    const a_one = firstBelow(ctx, P, s.a0, a1, p);
    sum += @intCast(a1 - a_one);

    // combine point: u <= T goes to the pool
    const a_t = firstBelow(ctx, P, s.a0, a_one, ctx.thr + 1);
    var a = a_t;
    while (a < a_one) : (a += 1) {
        if (try pool.add(gpa, ctx.x / P / ctx.primes[a], s.sign)) st.merged += 1;
    }
    st.inserted += a_one - a_t;
    if (a_t < a_one) {
        const n = a_one - a_t;
        const rng = ctx.x / P / ctx.primes[a_t] - ctx.x / P / ctx.primes[a_one - 1] + 1;
        st.ins_range += @min(n, rng);
        if (ROUGH_STATS) try ivals.append(gpa, .{ ctx.x / P / ctx.primes[a_one - 1], ctx.x / P / ctx.primes[a_t] });
        st.rw_elems += n;
        st.ins_tails += 1;
        if (n > rng) {
            st.dense_tails += 1;
            st.dense_elems += n;
            st.dense_range += rng;
        }
    }

    // T < u <= L: pi / P2 elementwise
    const a_root = firstBelow(ctx, P, s.a0, a_t, ctx.tab + 1);
    const a_pi = @max(a_root, firstBelow(ctx, P, s.a0, a_t, p *| p));
    const a_p2 = @max(a_root, firstBelow(ctx, P, s.a0, a_pi, p *| p *| p));
    a = a_pi;
    while (a < a_t) : (a += 1) sum += @intCast(ctx.pis[ctx.x / P / ctx.primes[a]] + 1 - c);
    st.set_pi += a_t - a_pi;
    a = a_p2;
    while (a < a_pi) : (a += 1) sum += @intCast(phiP2(ctx, ctx.x / P / ctx.primes[a], c, st));
    st.set_p2 += a_pi - a_p2;

    // head descends as two sets
    if (s.a0 < a_p2) {
        try next.append(gpa, .{ .P = P, .a0 = s.a0, .a1 = a_p2, .sign = s.sign });
        const p_c = ctx.primes[c];
        if (P <= ctx.x / p_c) try next.append(gpa, .{ .P = P * p_c, .a0 = s.a0, .a1 = a_p2, .sign = -s.sign });
    }
    return s.sign * sum;
}

// For level c: rough[v] = #{1 <= n <= v : n has no prime factor <= p_c}.
fn roughLevel(gpa: std.mem.Allocator, ctx: Ctx, c: u64, lpf: []const u32, rough: []u32, mert: f64, ivals: *std.ArrayList([2]u64), st: *Stats) !void {
    _ = gpa;
    const p_c = ctx.primes[c];
    var cnt: u32 = 0;
    rough[0] = 0;
    var v: u64 = 1;
    while (v <= ctx.thr) : (v += 1) {
        if (v == 1 or lpf[v] > p_c) cnt += 1;
        rough[v] = cnt;
    }
    const iv = ivals.items;
    const union0 = st.rw_union;
    defer {
        st.rw_levels += 1;
        st.rw_tau_w += iv.len * (st.rw_union - union0);
        st.rw_max_tau = @max(st.rw_max_tau, iv.len);
    }
    for (iv) |e| {
        st.rw_tails += 1;
        st.rw_range += e[1] - e[0];
        st.rw_rough += rough[e[1]] - rough[e[0]];
        st.rw_pred += to_f(e[1] - e[0]) * mert;
    }
    std.mem.sort([2]u64, iv, {}, struct {
        fn lt(_: void, a: [2]u64, b: [2]u64) bool {
            return a[0] < b[0];
        }
    }.lt);
    var lo = iv[0][0];
    var hi = iv[0][1];
    for (iv[1..]) |e| {
        if (e[0] <= hi) {
            hi = @max(hi, e[1]);
        } else {
            st.rw_union += rough[hi] - rough[lo];
            lo = e[0];
            hi = e[1];
        }
    }
    st.rw_union += rough[hi] - rough[lo];
}

fn piHybrid(gpa: std.mem.Allocator, x: u64, tab: u64, thr: u64, st: *Stats) !u64 {
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
    const ctx = Ctx{ .x = x, .sqrt_x = sqrt_x, .tab = tab, .thr = thr, .primes = ext, .pis = pis };

    var cur = Pool{ .k = try gpa.alloc(i64, thr + 1) };
    var nxt = Pool{ .k = try gpa.alloc(i64, thr + 1) };
    @memset(cur.k, 0);
    @memset(nxt.k, 0);
    var sets: std.ArrayList(Set) = .empty;
    var next_sets: std.ArrayList(Set) = .empty;
    var ivals: std.ArrayList([2]u64) = .empty;

    // rough-walk census support: least prime factor to T, per-level prefix
    // counts of c-rough n, Mertens product per level
    const lpf = try gpa.alloc(u32, thr + 1);
    const rough = try gpa.alloc(u32, thr + 1);
    const mertens = try gpa.alloc(f64, A + 2);
    if (ROUGH_STATS) {
        @memset(lpf, 0);
        var i: u64 = 2;
        while (i <= thr) : (i += 1) {
            if (lpf[i] != 0) continue;
            var j = i;
            while (j <= thr) : (j += i) {
                if (lpf[j] == 0) lpf[j] = @intCast(i);
            }
        }
        mertens[0] = 1;
        var cc: u64 = 1;
        while (cc <= A) : (cc += 1) mertens[cc] = mertens[cc - 1] * (1.0 - 1.0 / to_f(primes[cc]));
    }

    var total: i64 = @intCast(x);
    var a: u64 = 1;
    while (a <= A) : (a += 1) total -= @intCast(x / primes[a]);

    // top-level sets with p_b > cbrt x are all-ones
    const cbrt_x = icbrt(x);
    var pc3: u64 = 0;
    while (pc3 + 1 <= A and primes[pc3 + 1] <= cbrt_x) pc3 += 1;
    var b: u64 = A - 1;
    while (b > pc3) : (b -= 1) total += @intCast(A - b);

    var c: u64 = pc3;
    while (true) {
        if (c + 1 <= pc3 and c + 2 <= A + 1) try sets.append(gpa, .{ .P = primes[c + 1], .a0 = c + 2, .a1 = A + 1, .sign = 1 });
        st.max_sets_level = @max(st.max_sets_level, sets.items.len);
        for (sets.items) |s| total += try stepSet(gpa, ctx, s, c, &cur, &next_sets, &ivals, st);
        if (ROUGH_STATS and ivals.items.len > 0) try roughLevel(gpa, ctx, c, lpf, rough, mertens[c], &ivals, st);
        ivals.clearRetainingCapacity();
        sets.clearRetainingCapacity();
        std.mem.swap(std.ArrayList(Set), &sets, &next_sets);

        st.live += cur.touched.items.len;
        for (cur.touched.items) |u| {
            const k = cur.k[u];
            cur.k[u] = 0;
            if (k == 0) continue;
            if (resolve(ctx, u, c, st)) |phi| {
                total += k * @as(i64, @intCast(phi));
                continue;
            }
            st.steps += 1;
            _ = try nxt.add(gpa, u, k);
            const v = u / ctx.primes[c];
            if (v > 0) _ = try nxt.add(gpa, v, -k);
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
    const x_str = args.next() orelse return std.debug.print("Usage: pi-phi-hybrid <x> [t] [e]\n", .{});
    const x = try std.fmt.parseInt(u64, x_str, 10);
    const t: f64 = if (args.next()) |s| try std.fmt.parseFloat(f64, s) else 0.5;
    const e: f64 = if (args.next()) |s| try std.fmt.parseFloat(f64, s) else 0.5;
    const sqrt_x = isqrt(x);
    const xf = to_f(x);
    const thr: u64 = @max(1, @min(sqrt_x, @as(u64, @intFromFloat(@floor(std.math.pow(f64, xf, t))))));
    const tab: u64 = @max(1, @min(sqrt_x, @as(u64, @intFromFloat(@floor(std.math.pow(f64, xf, e))))));
    var st = Stats{};
    const pi_x = try piHybrid(std.heap.page_allocator, x, tab, thr, &st);

    var check: []const u8 = "";
    var p10: u64 = 1;
    for (KNOWN) |k| {
        if (p10 == x) check = if (k == pi_x) "  MATCH" else "  *** MISMATCH ***";
        p10 *%= 10;
    }
    const x23 = std.math.pow(f64, xf, 2.0 / 3.0);
    const elem = st.set_pi + st.set_p2 + st.set_c0;
    std.debug.print("pi({d}) = {d}{s}\n", .{ x, pi_x, check });
    if (ROUGH_STATS) std.debug.print("rough walk: {d} tails  E elements {d}  R u-range {d}  W rough n {d} (Mertens {d:.0}, ratio {d:.3})  W-union {d}  |  E/W {d:.2}x  R/W {d:.2}x  W/W-union {d:.1}x\n", .{ st.rw_tails, st.rw_elems, st.rw_range, st.rw_rough, st.rw_pred, to_f(st.rw_rough) / @max(1.0, st.rw_pred), st.rw_union, to_f(st.rw_elems) / @max(1.0, to_f(st.rw_rough)), to_f(st.rw_range) / @max(1.0, to_f(st.rw_rough)), to_f(st.rw_rough) / @max(1.0, to_f(st.rw_union)) });
    if (ROUGH_STATS) std.debug.print("plain query: {d} levels, max tails/level {d}, sum tails*W-union {d}  vs E {d} ({d:.3}x)  vs W {d} ({d:.3}x)\n", .{ st.rw_levels, st.rw_max_tau, st.rw_tau_w, st.rw_elems, to_f(st.rw_tau_w) / @max(1.0, to_f(st.rw_elems)), st.rw_rough, to_f(st.rw_tau_w) / @max(1.0, to_f(st.rw_rough)) });
    std.debug.print("compactness: {d} tails, {d} elements over {d} u-slots ({d:.1} elems/slot); dense tails (N > range) {d}: {d} elems over {d} slots ({d:.1}x); overlap: {d} slots -> {d} live ({d:.1}x)\n", .{ st.ins_tails, st.inserted, st.ins_range, to_f(st.inserted) / @max(1.0, to_f(st.ins_range)), st.dense_tails, st.dense_elems, st.dense_range, to_f(st.dense_elems) / @max(1.0, to_f(st.dense_range)), st.ins_range, st.live, to_f(st.ins_range) / @max(1.0, to_f(st.live)) });
    std.debug.print("T = {d}  L = {d}  sets {d} (max/level {d})  set elementwise {d} (pi {d} p2 {d} c0 {d})  inserted {d} merged {d} ({d:.1}%)  pool live {d} steps {d}  p2 lookups {d}  [x^2/3 = {e:.2}]\n", .{ thr, tab, st.sets, st.max_sets_level, elem, st.set_pi, st.set_p2, st.set_c0, st.inserted, st.merged, pc(st.merged, @max(st.inserted, 1)), st.live, st.steps, st.p2_lookups, x23 });
}
