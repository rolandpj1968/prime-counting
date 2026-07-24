//! gourdon.zig — Xavier Gourdon's 2001 reorganisation of Deléglise–Rivat.
//!
//!   π(x) = A − B + ω + φ₀ + Σ          (Theorem 1 / eq (5) of gourdon.ps)
//!
//! Term status (top-down build; lmo.zig referees the total):
//!   ω+B — FUSED on ONE segmented O(1)-kill counter folded to √z. ω's leaves query
//!         φ(v,bi) at stage bi<π(x*) DURING the fold; after the full fold the counter
//!         is φ(·,π(√z)), so B's π(x/p)=φ(x/p,π(√z))+π(√z)−1 is read straight off it
//!         (Legendre; z<y² ⇒ P₂ term vanishes). Gourdon folds π(√z) primes vs lmo's
//!         π(y) — fewer kills AND fewer prefix queries → beats lmo (op-counts confirm).
//!         ω uses Gourdon's C/D split (x/pm<p² ⇒ closed form 1+max(0,π(v)−bi), π via
//!         the π oracle) and DR/LMO dense/sparse leaf enumeration (p>√y ⇒ m prime).
//!   A/Σ — query π only at points ≤ √x (proven); monotone cursor / binary search on
//!         the prime list — no table, no sweep.
//!   φ₀  — closed form φ(x/n,1) = ⌈(x/n)/2⌉ (k=1).
//!   Σ   — seven closed forms (a=π(y), b=π(x^1/3), c=π(√(x/y)), d=π(x*)).
//!
//! Memory O(√x)+O(y)+O(SEG), no O(z) structure — ~12 MB at 10^15. One O(z)-time pass.
//! x* = max(x^1/4, x/y²); primes 0-indexed; k=1. (computeB is kept only as B's ref.)

const std = @import("std");
const common = @import("common.zig");
const lmo = @import("lmo.zig");

var g_bmarks: u64 = 0; // INST-only: B sieve marks (composite bit-clears in answerPi)

// ---------------------------------------------------------------- integer roots
fn isqrt(n: u64) u64 {
    return common.isqrt(n);
}
// Generic-width roots of x (result always fits u64). Copied from lmo.zig.
fn icbrtG(comptime T: type, x: T) u64 {
    if (x == 0) return 0;
    var r: u64 = @intFromFloat(std.math.pow(f64, @floatFromInt(x), 1.0 / 3.0));
    if (r == 0) r = 1;
    while (@as(T, r) * @as(T, r) * @as(T, r) > x) r -= 1;
    while (@as(T, r + 1) * @as(T, r + 1) * @as(T, r + 1) <= x) r += 1;
    return r;
}
fn isqrtG(comptime T: type, n: T) u64 {
    if (n < 2) return @intCast(n);
    var x: u64 = @intFromFloat(@sqrt(@as(f64, @floatFromInt(n))));
    while (@as(T, x) * @as(T, x) > n) x -= 1;
    while (@as(T, x + 1) * @as(T, x + 1) <= n) x += 1;
    return x;
}
/// ⌊x/d⌋ with x wide (X) and d a u64 that fits — result fits u64 (it's ≤ z or √x).
/// ⌊x / (a·b)⌋ with the product formed in X, so it is safe when a·b overflows u64
/// (reachable past ~10²³). b = 0 ⇒ no bound; a·b > x ⇒ 0.
inline fn mBound(comptime X: type, x: X, a: u64, b: u64) u64 {
    if (b == 0) return std.math.maxInt(u64);
    const d: X = @as(X, a) * @as(X, b);
    return if (d > x) 0 else @intCast(x / d);
}

inline fn xdiv(comptime X: type, x: X, d: u64) u64 {
    return @intCast(x / @as(X, d));
}
/// Pin the calling thread to one logical CPU (best-effort). Copied from lmo.zig.
fn pinToCpu(cpu: u32) void {
    var set: std.os.linux.cpu_set_t = @splat(0);
    set[cpu / @bitSizeOf(usize)] |= @as(usize, 1) << @intCast(cpu % @bitSizeOf(usize));
    std.os.linux.sched_setaffinity(0, &set) catch {};
}

// --------------------------------------------------------------- small sieve
/// Primes ≤ nmax by a SEGMENTED Eratosthenes: base primes ≤ √nmax mark each cache-
/// sized segment in turn, so the transient sieve buffer is O(SEG), not O(nmax/8).
/// (The returned prime list is still O(nmax/ln) — that's the persistent cost.)
fn sievePrimes(comptime P: type, gpa: std.mem.Allocator, nmax: u64) ![]P {
    if (nmax < 2) return gpa.alloc(P, 0);
    const rt: usize = @intCast(isqrt(nmax)); // base primes needed up to √nmax
    // base primes ≤ rt (simple; rt = x^(1/4) here, tiny)
    const bcomp = try gpa.alloc(bool, rt + 1);
    defer gpa.free(bcomp);
    @memset(bcomp, false);
    var i: usize = 2;
    while (i * i <= rt) : (i += 1) {
        if (!bcomp[i]) {
            var j = i * i;
            while (j <= rt) : (j += i) bcomp[j] = true;
        }
    }
    var nb: usize = 0;
    i = 2;
    while (i <= rt) : (i += 1) {
        if (!bcomp[i]) nb += 1;
    }
    const base = try gpa.alloc(u32, nb);
    defer gpa.free(base);
    {
        var k: usize = 0;
        i = 2;
        while (i <= rt) : (i += 1) {
            if (!bcomp[i]) {
                base[k] = @intCast(i);
                k += 1;
            }
        }
    }
    // upper bound on π(nmax) — Dusart: π(n) ≤ (n/ln n)(1 + 1.2762/ln n). Using ⌊ln⌋
    // (< ln, via log2·ln2) only loosens it, so this is a safe cap with ~3% slop.
    const ln: u64 = @max(1, (@as(u64, log2Floor(@intCast(nmax))) * 693) / 1000);
    const t1: u64 = nmax / ln;
    const cap: usize = @max(@as(usize, 16), @as(usize, @intCast(t1 + (t1 * 1276) / (ln * 1000) + 64)));
    var primes = try gpa.alloc(P, cap); // base primes stay u32 (≤ √nmax = x^(1/4))
    errdefer gpa.free(primes);
    var cnt: usize = 0;

    const SEG: usize = 1 << 18;
    const seg = try gpa.alloc(bool, SEG);
    defer gpa.free(seg);
    var lo: u64 = 0;
    while (lo <= nmax) : (lo += SEG) {
        const hi = @min(lo + SEG, nmax + 1);
        const len: usize = @intCast(hi - lo);
        @memset(seg[0..len], false);
        for (base) |p32| {
            const p: u64 = p32;
            if (p * p >= hi) break;
            var j: u64 = @max(p * p, ((lo + p - 1) / p) * p);
            while (j < hi) : (j += p) seg[@intCast(j - lo)] = true;
        }
        if (lo == 0) {
            seg[0] = true; // 0
            if (len > 1) seg[1] = true; // 1
        }
        var k: usize = 0;
        while (k < len) : (k += 1) {
            if (!seg[k]) {
                primes[cnt] = @intCast(lo + @as(u64, k));
                cnt += 1;
            }
        }
    }
    return gpa.realloc(primes, cnt);
}

/// Persistent tables: the fused μ+δ leaf table only to y (ω leaves, φ₀ — never read
/// past y), an explicit
/// prime list only to plist_max (the largest prime any loop enumerates by value), and
/// the O(1) π oracle over [0, √x] as a coprime-30 bitset. Storing the √x primes as a
/// u64 list instead costs 8·π(√x) — 3.5 GB at x=10²⁰, ~10× the oracle.
/// δ saturation point. √y is the largest p the leaf test δ[m] > p ever compares
/// against, so saturating is exact iff √y < DELTA_MAX, i.e. y < 4.295e9. With
/// y = α·x^(1/3) and α from chooseAlpha that holds to x ≈ 10^25 and fails at 10^26
/// (y = 8.7e9, √y = 93529). Sieve.init enforces it rather than trusting the bound.
/// Fused μ+δ leaf table, one u16 per m ∈ [0..y]:
///   0        ⇒ μ(m) = 0 (some p² | m) — rejected everywhere, so μ's zero bit is free
///   bit 15   ⇒ μ(m) = +1 (even # prime factors); clear ⇒ −1. (+1 sets the bit so
///              m = 1 encodes as 0x8000 and cannot collide with the 0 sentinel.)
///   bits 0-14⇒ min(π(spf(m)), LEAF_SAT) — the INDEX of the smallest prime factor,
///              not its value. The only consumer is the leaf test spf(m) > p with
///              p = primes[bi], and spf > p ⟺ π(spf) > bi+1, so the compare is
///              against the loop index and needs no prime value at all. Indices
///              compress ~ln p better than values: 15 bits reaches p_32766 =
///              386,051, i.e. √y < 386k ⇒ x ≈ 10^29 — past the runtime horizon —
///              where the old u16 VALUE encoding died at √y > 65535 (x ≈ 10^25).
/// One table and one stream where μ (i8) + δ (u16) were two; 2 bytes/m instead of 3.
const LEAF_SAT: u16 = 0x7FFF;
const LEAF_PLUS: u16 = 0x8000;
const LeafIdxTooNarrow = error.YExceedsLeafIndexRange;

fn Sieve(comptime P: type) type {
    return struct {
        const Self = @This();
        leaf: []u16, // fused μ + spf-index table, [0..y] — see LEAF_SAT above
        primes: []P, // ascending, ≤ plist_max (enumeration only)
        pio: PiOracle, // π(v) and descending prime walk over [0, √x]

    fn init(gpa: std.mem.Allocator, sqx: u64, y: u64, plist_max: u64, nthreads: usize, pins: ?[]const u32) !Self {
        const primes = try sievePrimes(P, gpa, plist_max);
        // Capped at plist_max (≈ y): every surviving count() consumer is ≤ y —
        // A's v < y pairs, Σ₄/₅/₆ (p > x* ≥ x/y² bounds each), the scalars, the
        // leaf-table guard, and the window base-prime count x^(1/4). π above y
        // comes from buildBoundaryPi + sieved windows. This is where the √x/30
        // resident term (3.5 GB at 10^22) actually disappears.
        _ = sqx;
        const pio = try buildPiOraclePar(gpa, plist_max, nthreads, pins);
        const Ny: usize = @intCast(y + 1);

        // Guard: saturation is exact iff no compared index reaches LEAF_SAT, i.e.
        // π(√y) < LEAF_SAT. An error, not a comment — the u16-value guard's lesson.
        if (pio.count(isqrt(y)) >= LEAF_SAT) return LeafIdxTooNarrow;

        const leaf = try gpa.alloc(u16, Ny);
        errdefer gpa.free(leaf);
        @memset(leaf, LEAF_PLUS); // μ=+1, spf-index unset (0)

        // Pass 1, primes ascending: flip the sign bit per prime factor; first flip
        // also stamps the spf index (first-write-wins ⇒ smallest, since q ascends).
        for (primes, 0..) |q32, qi| {
            const q: usize = @intCast(q32);
            if (q >= Ny) break;
            const qidx1: u16 = @intCast(@min(qi + 1, LEAF_SAT)); // 1-based π(q), saturated
            var j: usize = q;
            while (j < Ny) : (j += q) {
                var v = leaf[j] ^ LEAF_PLUS;
                if (v & LEAF_SAT == 0) v |= qidx1;
                leaf[j] = v;
            }
        }
        // Pass 2, AFTER all sign/idx writes so nothing revives a zeroed entry:
        // q² | m ⇒ μ(m) = 0 ⇒ the 0 sentinel.
        for (primes) |q32| {
            const q: usize = @intCast(q32);
            const qq = q * q;
            if (qq >= Ny) break;
            var j: usize = qq;
            while (j < Ny) : (j += qq) leaf[j] = 0;
        }
        if (Ny > 0) leaf[0] = 0;
        // leaf[1] stays LEAF_PLUS: μ(1) = +1 (φ₀ needs it), idx 0 ⇒ never a leaf.
        return .{ .leaf = leaf, .primes = primes, .pio = pio };
    }
    fn deinit(self: *Self, gpa: std.mem.Allocator) void {
        gpa.free(self.leaf);
        gpa.free(self.primes);
        gpa.free(self.pio.bits);
        gpa.free(self.pio.pref);
    }
    };
}

// ============================================================================
// COPIED from lmo.zig: mod-30 wheel + 3-level O(1)-kill counter (Counter3P).
// ============================================================================
inline fn log2Floor(n: usize) u6 {
    return @intCast(63 - @clz(@as(u64, @max(n, 1))));
}
const W30 = [8]u8{ 1, 7, 11, 13, 17, 19, 23, 29 };
const W30GAP = [8]u8{ 6, 4, 2, 4, 2, 4, 6, 2 };
const COP30: [30]bool = blk: {
    var c: [30]bool = undefined;
    for (0..30) |r| c[r] = (r % 2 != 0 and r % 3 != 0 and r % 5 != 0);
    break :blk c;
};
const W30IDX: [30]u8 = blk: { // residue r (coprime to 30) → its index in W30
    var t: [30]u8 = @splat(0);
    for (W30, 0..) |res, i| t[res] = @intCast(i);
    break :blk t;
};
const MASK30: [15]u64 = blk: {
    @setEvalBranchQuota(20000);
    var m: [15]u64 = @splat(0);
    for (0..15) |w| {
        for (0..64) |i| {
            const n = w * 64 + i;
            const r = n % 30;
            var coprime = true;
            for ([_]u64{ 2, 3, 5 }) |q| {
                if (r % q == 0) coprime = false;
            }
            if (coprime) m[w] |= @as(u64, 1) << @as(u6, @intCast(i));
        }
    }
    break :blk m;
};

inline fn phiSmall(v: u64, b: usize) i64 {
    return switch (b) {
        0 => @intCast(v),
        1 => @intCast(v - v / 2),
        2 => @intCast(v - v / 2 - v / 3 + v / 6),
        else => unreachable,
    };
}

const Counter3P = struct {
    bits: []u64,
    cnt1: []u32,
    cnt2: []u32,
    s1: u6,
    s2: u6,
    nwords: usize,
    total: i64,

    fn init(gpa: std.mem.Allocator, seg: usize) !Counter3P {
        const nwords = (seg + 63) / 64;
        const t: u6 = @intCast(@max(1, (log2Floor(nwords) + 2) / 3));
        const nblocks = (nwords >> t) + 1;
        const nsuper = (nblocks >> t) + 1;
        return .{
            .bits = try gpa.alloc(u64, nwords),
            .cnt1 = try gpa.alloc(u32, nblocks),
            .cnt2 = try gpa.alloc(u32, nsuper),
            .s1 = t,
            .s2 = t,
            .nwords = nwords,
            .total = 0,
        };
    }
    fn deinit(self: *Counter3P, gpa: std.mem.Allocator) void {
        gpa.free(self.bits);
        gpa.free(self.cnt1);
        gpa.free(self.cnt2);
    }
    fn reset(self: *Counter3P, len: usize) void {
        const nw = (len + 63) / 64;
        for (0..nw) |w| self.bits[w] = MASK30[w % 15];
        if (len % 64 != 0) self.bits[nw - 1] &= (@as(u64, 1) << @as(u6, @intCast(len % 64))) - 1;
        @memset(self.bits[nw..self.nwords], 0);
        @memset(self.cnt1, 0);
        @memset(self.cnt2, 0);
        var tot: i64 = 0;
        for (self.bits, 0..) |word, w| {
            const c: u32 = @popCount(word);
            self.cnt1[w >> self.s1] += c;
            self.cnt2[(w >> self.s1) >> self.s2] += c;
            tot += c;
        }
        self.total = tot;
    }
    inline fn kill(self: *Counter3P, i: usize) void {
        const w = i >> 6;
        const b = @as(u64, 1) << @as(u6, @intCast(i & 63));
        const alive: u32 = @intFromBool(self.bits[w] & b != 0);
        self.bits[w] &= ~b;
        const blk = w >> self.s1;
        self.cnt1[blk] -= alive;
        self.cnt2[blk >> self.s2] -= alive;
        self.total -= @as(i64, alive);
    }
    /// Word-parallel strike: clear all bits of `t` from word `w` in one op (delta via popcount).
    inline fn killWord(self: *Counter3P, w: usize, t: u64) void {
        const old = self.bits[w];
        self.bits[w] = old & ~t;
        const d: u32 = @popCount(old & t);
        const blk = w >> self.s1;
        self.cnt1[blk] -= d;
        self.cnt2[blk >> self.s2] -= d;
        self.total -= @as(i64, d);
    }
    fn prefix(self: *const Counter3P, i: usize) i64 {
        const w = i >> 6;
        const blk = w >> self.s1;
        const sblk = blk >> self.s2;
        var s: i64 = 0;
        for (self.cnt2[0..sblk]) |c| s += c;
        for (self.cnt1[sblk << self.s2 .. blk]) |c| s += c;
        for (self.bits[blk << self.s1 .. w]) |word| s += @popCount(word);
        const r: u6 = @intCast(i & 63);
        const mask: u64 = if (r == 63) ~@as(u64, 0) else (@as(u64, 1) << (r + 1)) - 1;
        s += @popCount(self.bits[w] & mask);
        return s;
    }
};

/// 2-level variant (one count array ≈ √nwords blocks): O(1) kill with ONE decrement
/// (vs 3P's two), O(√nwords) prefix (vs 3P's ~nwords^(1/3)). Same MASK30 reset +
/// branchless kill as Counter3P — the ONLY difference is the count level, to A/B the
/// kill-heavy vs query-heavy tradeoff. gourdon's ω is kill-dominated (C/D offload).
const Counter2P = struct {
    bits: []u64,
    cnt: []u32,
    s1: u6,
    nwords: usize,
    total: i64,

    fn init(gpa: std.mem.Allocator, seg: usize) !Counter2P {
        const nwords = (seg + 63) / 64;
        const s1: u6 = @intCast((log2Floor(nwords) + 1) / 2); // block ≈ √nwords
        const nblocks = (nwords >> s1) + 1;
        return .{ .bits = try gpa.alloc(u64, nwords), .cnt = try gpa.alloc(u32, nblocks), .s1 = s1, .nwords = nwords, .total = 0 };
    }
    fn deinit(self: *Counter2P, gpa: std.mem.Allocator) void {
        gpa.free(self.bits);
        gpa.free(self.cnt);
    }
    fn reset(self: *Counter2P, len: usize) void {
        const nw = (len + 63) / 64;
        for (0..nw) |w| self.bits[w] = MASK30[w % 15];
        if (len % 64 != 0) self.bits[nw - 1] &= (@as(u64, 1) << @as(u6, @intCast(len % 64))) - 1;
        @memset(self.bits[nw..self.nwords], 0);
        @memset(self.cnt, 0);
        var tot: i64 = 0;
        for (self.bits, 0..) |word, w| {
            const c: u32 = @popCount(word);
            self.cnt[w >> self.s1] += c;
            tot += c;
        }
        self.total = tot;
    }
    inline fn kill(self: *Counter2P, i: usize) void {
        const w = i >> 6;
        const b = @as(u64, 1) << @as(u6, @intCast(i & 63));
        const alive: u32 = @intFromBool(self.bits[w] & b != 0);
        self.bits[w] &= ~b;
        self.cnt[w >> self.s1] -= alive;
        self.total -= @as(i64, alive);
    }
    /// Word-parallel strike: clear all bits of `t` from word `w` in one op (delta via popcount).
    inline fn killWord(self: *Counter2P, w: usize, t: u64) void {
        const old = self.bits[w];
        self.bits[w] = old & ~t;
        const d: u32 = @popCount(old & t);
        self.cnt[w >> self.s1] -= d;
        self.total -= @as(i64, d);
    }
    /// Counted kills maintain cnt/total on every strike so the counter can be
    /// queried at any point mid-fold. The fold stages ABOVE π(x*) are never queried
    /// mid-fold (leaves only exist for bi < nax), so they can clear bits and skip
    /// the bookkeeping entirely, with one rebuild() before the counter is next read.
    inline fn strike(self: *Counter2P, i: usize) void {
        self.bits[i >> 6] &= ~(@as(u64, 1) << @as(u6, @intCast(i & 63)));
    }
    inline fn strikeWord(self: *Counter2P, w: usize, t: u64) void {
        self.bits[w] &= ~t;
    }
    /// Recompute cnt/total from bits — one popcount pass, same shape as reset()'s.
    fn rebuild(self: *Counter2P) void {
        @memset(self.cnt, 0);
        var tot: i64 = 0;
        for (self.bits, 0..) |word, w| {
            const c: u32 = @popCount(word);
            self.cnt[w >> self.s1] += c;
            tot += c;
        }
        self.total = tot;
    }
    fn prefix(self: *const Counter2P, i: usize) i64 {
        const w = i >> 6;
        const blk = w >> self.s1;
        var s: i64 = 0;
        for (self.cnt[0..blk]) |c| s += c;
        for (self.bits[blk << self.s1 .. w]) |word| s += @popCount(word);
        const r: u6 = @intCast(i & 63);
        const mask: u64 = if (r == 63) ~@as(u64, 0) else (@as(u64, 1) << (r + 1)) - 1;
        s += @popCount(self.bits[w] & mask);
        return s;
    }
};

/// ω+B counter: MEASURED Counter2P > Counter3P for gourdon (~8.5% faster @10^14).
/// gourdon's ω is kill-dominated (C/D offloads most queries), so 3P's extra per-kill
/// decrement costs more than its cheaper query saves — the opposite of lmo's S2, where
/// the query-heavy pattern makes 3P win. Swap to Counter3P to re-measure.
const Ctr = Counter2P;

// ------------------------------------------------------------- segmented π
/// Answer π at every point in `pts` (order preserved in `out`), no O(z) table: sort
/// the queries, sweep [1,z] in cache-sized segments carrying a running π, read each
/// query off as the walk passes it. base = primes ≤ √z. O(z) time, O(SEG+n) memory.
///
/// mod-30 WHEEL, BIT-PACKED: one bit per integer, each word initialised to MASK30
/// (2,3,5 pre-struck). Composites of p≥7 are cleared by wheel-stepping (bit ops, no
/// div); primes are counted 64-at-a-time by popcount. Survivors in [0,v] = {1} ∪
/// {coprime-30 primes ≤ v}, so π(v) = 2 + popcount([0,v]) for v ≥ 7 (the +2 = 3 for
/// {2,3,5} − 1 for the spurious "1"). All B query points are ≥ √x ≥ 7.
fn answerPi(comptime INST: bool, gpa: std.mem.Allocator, base: []const u32, pts: []const u64, out: []u64, z: u64) !void {
    const n = pts.len;
    const ord = try gpa.alloc(usize, n);
    defer gpa.free(ord);
    for (ord, 0..) |*o, i| o.* = i;
    const Ctx = struct {
        p: []const u64,
        fn lt(self: @This(), a: usize, b: usize) bool {
            return self.p[a] < self.p[b];
        }
    };
    std.mem.sort(usize, ord, Ctx{ .p = pts }, Ctx.lt);

    const SEG: usize = 960 * 256; // multiple of lcm(30,64)=960 ⇒ MASK30[w%15] tiles
    const NW = SEG / 64;
    const bits = try gpa.alloc(u64, NW);
    defer gpa.free(bits);

    var qi: usize = 0;
    // small v (< 7): π(0..6) = 0,0,1,2,2,3,3 — the wheel formula (2+…) only holds v≥7.
    while (qi < n and pts[ord[qi]] < 7) : (qi += 1) {
        const v = pts[ord[qi]];
        out[ord[qi]] = if (v < 2) 0 else if (v < 3) 1 else if (v < 5) 2 else 3;
    }
    var gc: u64 = 0; // survivor popcount in [0, lo)
    var lo: u64 = 0; // multiple of SEG (hence of 960)
    while (lo <= z) : (lo += SEG) {
        const hi = @min(lo + SEG, z + 1);
        const len: usize = @intCast(hi - lo);
        const nw: usize = (len + 63) / 64;
        for (0..nw) |w| bits[w] = MASK30[w % 15];
        if (len % 64 != 0) bits[nw - 1] &= (@as(u64, 1) << @as(u6, @intCast(len % 64))) - 1;
        for (base) |p32| {
            const p: u64 = p32;
            if (p < 7) continue; // 2,3,5 already struck by MASK30
            if (p * p >= hi) break;
            const L = @max(p * p, lo);
            var k = (L + p - 1) / p;
            while (!COP30[@intCast(k % 30)]) k += 1; // first coprime-30 multiplier
            var w = p * k;
            var widx: u8 = W30IDX[@intCast(k % 30)];
            while (w < hi) {
                if (INST) g_bmarks += 1;
                const i: usize = @intCast(w - lo);
                bits[i >> 6] &= ~(@as(u64, 1) << @as(u6, @intCast(i & 63)));
                w += p * W30GAP[widx];
                widx = (widx + 1) & 7;
            }
        }
        var rc: u64 = 0; // survivor popcount in [lo, lo+w·64)
        var w: usize = 0;
        while (w < nw) : (w += 1) {
            const word = bits[w];
            const wbase = lo + @as(u64, w) * 64;
            while (qi < n and pts[ord[qi]] < wbase + 64 and pts[ord[qi]] < hi) : (qi += 1) {
                const b: u6 = @intCast(pts[ord[qi]] - wbase);
                const mask: u64 = if (b == 63) ~@as(u64, 0) else (@as(u64, 1) << (b + 1)) - 1;
                out[ord[qi]] = 2 + gc + rc + @popCount(word & mask);
            }
            rc += @popCount(word);
        }
        gc += rc;
    }
    while (qi < n) : (qi += 1) out[ord[qi]] = 2 + gc; // defensive
}

// ------------------------------------------------------------------- π oracle
const W30CNT: [30]u8 = blk: { // r → #{s ∈ [1,r] : gcd(s,30)=1}
    var t: [30]u8 = @splat(0);
    var c: u8 = 0;
    for (0..30) |r| {
        if (r > 0 and COP30[r]) c += 1;
        t[r] = c;
    }
    break :blk t;
};

/// π(v) and descending prime enumeration over [0, n], as a coprime-30 primality
/// bitset plus block prefix counts. One bit per coprime-30 residue ⇒ n/30 bytes,
/// against 8·π(n) for an explicit u64 prime list — 10.4× less at n=10¹⁰ — and π(v)
/// costs a prefix load plus ≤8 popcounts on one cache line of bits, against ~29
/// dependent misses for a binary search striding a 3.5 GB array.
///
/// Word w covers exactly the integers [240w, 240w+240): 64 bits ÷ 8 bits-per-30.
/// Bit j ↔ value 30·(j/8) + W30[j%8], so bit 0 ↔ 1 (cleared: 1 is not prime) and
/// 2, 3, 5 have no bit at all — count() adds them back by inspection.
const PiOracle = struct {
    bits: []u64, // 1 = prime (coprime-30 candidates only)
    pref: []u64, // pref[b] = # set bits strictly below word b·BLKW
    n: u64,
    nbits: usize,

    const BLKW: usize = 8; // words per prefix block = one 64 B cache line

    inline fn kidx(v: u64) usize { // bit index of coprime-30 v
        return @intCast((v / 30) * 8 + W30IDX[@intCast(v % 30)]);
    }
    inline fn kpos(v: u64) usize { // # coprime-30 integers in [1, v] = one past v's bit
        return @intCast((v / 30) * 8 + W30CNT[@intCast(v % 30)]);
    }
    inline fn val(j: usize) u64 { // inverse of kidx
        return @as(u64, j / 8) * 30 + W30[j % 8];
    }

    fn countBelow(self: *const PiOracle, k0: usize) u64 {
        const k = @min(k0, self.nbits);
        if (k == 0) return 0;
        const w = k >> 6;
        var c: u64 = self.pref[w / BLKW];
        for (self.bits[(w / BLKW) * BLKW .. w]) |word| c += @popCount(word);
        const r: u32 = @intCast(k & 63);
        if (r != 0) c += @popCount(self.bits[w] & ((@as(u64, 1) << @intCast(r)) - 1));
        return c;
    }

    /// π(v), exact for every v ≤ n.
    fn count(self: *const PiOracle, v: u64) u64 {
        if (v < 2) return 0;
        var c: u64 = 1; // 2
        if (v >= 3) c += 1;
        if (v >= 5) c += 1;
        return c + self.countBelow(kpos(v));
    }

    /// Largest prime ≤ v (v ≤ n). Returns 0 when v < 2.
    fn prevPrime(self: *const PiOracle, v: u64) u64 {
        if (v < 7) return if (v >= 5) 5 else if (v >= 3) 3 else if (v >= 2) 2 else 0;
        const k = @min(kpos(v), self.nbits);
        if (k == 0) return 5;
        var w = (k - 1) >> 6;
        const r: u32 = @intCast(k & 63);
        var m = self.bits[w];
        if (r != 0) m &= (@as(u64, 1) << @intCast(r)) - 1;
        while (true) {
            if (m != 0) return val((w << 6) + 63 - @clz(m));
            if (w == 0) return 5;
            w -= 1;
            m = self.bits[w];
        }
    }
};

const ORACLE_SEGW: usize = 1 << 15; // 32768 words = 256 KB of bits = 7.86M integers

/// The ω sweep's segment width (bits = 32 KB = L1d) — also the granularity of the
/// bpi boundary array, so window bases index it directly.
const SWEEP_SEGW: usize = 273 * 960;

/// First index with list[i] ≥ val (the prime list is ascending).
fn lowerBound(comptime P: type, list: []const P, val: u64) usize {
    var lo: usize = 0;
    var hi: usize = list.len;
    while (lo < hi) {
        const mid = (lo + hi) / 2;
        if (@as(u64, @intCast(list[mid])) < val) lo = mid + 1 else hi = mid;
    }
    return lo;
}

/// Strike one oracle segment (words [w0, w0+ORACLE_SEGW)). A word covers exactly
/// 240 integers, so segments touch disjoint words — the parallel build needs no
/// synchronisation beyond handing out segment indices.
fn oracleStrikeSeg(bits: []u64, base: []const u32, n: u64, nw: usize, w0: usize) void {
    const w1 = @min(w0 + ORACLE_SEGW, nw);
    const lo: u64 = @as(u64, w0) * 240;
    const hi: u64 = @min(@as(u64, w1) * 240, n + 1);
    for (base) |p32| {
        const p: u64 = p32;
        if (p < 7) continue;
        if (p * p >= hi) break;
        var m: u64 = @max(p, (lo + p - 1) / p);
        while (!COP30[@intCast(m % 30)]) m += 1;
        var wp: u8 = W30IDX[@intCast(m % 30)];
        var v: u64 = p * m;
        while (v < hi) {
            const k = PiOracle.kidx(v);
            bits[k >> 6] &= ~(@as(u64, 1) << @intCast(k & 63));
            v += p * W30GAP[wp];
            wp = (wp + 1) & 7;
        }
    }
}

/// Streamed boundary-π over [0, √x]: bpi[k] = π(k·SWEEP_SEGW − 1) plus π(√x)
/// itself, built by sieving windows and retaining ONLY the per-segment counts —
/// O(√x/segw) memory where the resident oracle kept √x/30 bytes of bits. This is
/// what lets the oracle be capped at plist_max (≈ y): every π value above y now
/// comes from bpi + a freshly sieved window.
// ---------------------------------------------- verbose progress
// Periodic intra-phase progress lines for multi-hour runs. Enabled only by
// Config.verbose (set per call in piGourdonV); every tick is a load-and-branch
// when disabled, and at most one clock read + atomic add per work unit when on.
// Work units are coarse (segments/windows/chunks), so the cost is noise.
var g_progress: bool = false;
const PROG_GATE_NS: u64 = 60 * 1_000_000_000; // min interval between lines

const Prog = struct {
    name: []const u8 = "",
    total: u64 = 1,
    done: std.atomic.Value(u64) = std.atomic.Value(u64).init(0),
    t0: u64 = 0,
    last: std.atomic.Value(u64) = std.atomic.Value(u64).init(0),

    fn begin(self: *Prog, name: []const u8, total: u64) void {
        if (!g_progress) return;
        const now = common.nowNs();
        self.name = name;
        self.total = @max(total, 1);
        self.done.store(0, .monotonic);
        self.t0 = now;
        self.last.store(now, .monotonic);
    }
    /// `pos` is the calling worker's current position in value space (sweep u,
    /// window base, prime value, …) — a representative coordinate, not a frontier:
    /// workers hold different dispatch units, so successive lines can wobble.
    fn tick(self: *Prog, pos: u64) void {
        if (!g_progress) return;
        const d = self.done.fetchAdd(1, .monotonic) + 1;
        const now = common.nowNs();
        // anchor the rate clock at the first completed unit, not at begin(): some
        // phases (A/Σ v-windows) only start after another phase drains per-worker.
        // Re-arm the gate too, so a late-starting phase doesn't print a zero-rate line.
        if (d == 1) {
            self.t0 = now;
            self.last.store(now, .monotonic);
            return;
        }
        const prev = self.last.load(.monotonic);
        if (now -% prev < PROG_GATE_NS) return;
        // one printer per gate window; losers just keep working
        if (self.last.cmpxchgStrong(prev, now, .monotonic, .monotonic) != null) return;
        const dd = @min(d, self.total);
        const el = @as(f64, @floatFromInt(now - self.t0)) / 1e9;
        const pct = 100.0 * @as(f64, @floatFromInt(dd)) / @as(f64, @floatFromInt(self.total));
        const eta = el * @as(f64, @floatFromInt(self.total - dd)) / @as(f64, @floatFromInt(dd));
        std.debug.print("    [{s:>6}] {d:5.1}%  {d}/{d}  @ {d}  elapsed {d:.0} s  eta {d:.0} s\n", .{ self.name, pct, dd, self.total, pos, el, eta });
    }
};
var prog_bpi: Prog = .{};
var prog_ap: Prog = .{};
var prog_aw: Prog = .{};
var prog_om: Prog = .{};

fn buildBoundaryPi(gpa: std.mem.Allocator, sqx: u64, segw: usize, nthreads: usize, pins: ?[]const u32) !struct { bpi: []u64, total: u64 } {
    const nbpi = @as(usize, @intCast(sqx / segw)) + 2;
    const bpi = try gpa.alloc(u64, nbpi);
    errdefer gpa.free(bpi);
    @memset(bpi, 0);
    const base = try sievePrimes(u32, gpa, isqrt(sqx));
    defer gpa.free(base);

    const AW: u64 = 16 * @as(u64, segw);
    const nwin: usize = @intCast(sqx / AW + 1);
    const nww = (16 * segw) / 240;
    const nsegs_total = nbpi - 1; // count slots: segment k = [k·segw, (k+1)·segw)

    const Ctx = struct {
        disp: std.atomic.Value(usize),
        bpi: []u64,
        base: []const u32,
        sqx: u64,
        nwin: usize,
        nsegs: usize,
        segw: usize,
    };
    var ctx = Ctx{
        .disp = std.atomic.Value(usize).init(0),
        .bpi = bpi,
        .base = base,
        .sqx = sqx,
        .nwin = nwin,
        .nsegs = nsegs_total,
        .segw = segw,
    };
    const Worker = struct {
        fn run(cx: *Ctx, bits: []u64, cpu: ?u32) void {
            if (cpu) |c| pinToCpu(c);
            const segww = cx.segw / 240; // words per segment
            while (true) {
                const wi = cx.disp.fetchAdd(1, .monotonic);
                if (wi >= cx.nwin) break;
                const wlo = @as(u64, wi) * (16 * @as(u64, cx.segw));
                const whi = @min(wlo + 16 * @as(u64, cx.segw), cx.sqx + 1);
                bwinFill(u32, bits, cx.base, cx.base.len, wlo, whi);
                if (whi == cx.sqx + 1) {
                    // zero every bit at values > √x so counts stay exact
                    const kend = PiOracle.kpos(cx.sqx) - @as(usize, @intCast(wlo / 30)) * 8;
                    const wend = kend >> 6;
                    const rend: u32 = @intCast(kend & 63);
                    if (rend != 0 and wend < bits.len)
                        bits[wend] &= (@as(u64, 1) << @intCast(rend)) - 1;
                    const zfrom = if (rend != 0) wend + 1 else wend;
                    if (zfrom < bits.len) @memset(bits[zfrom..], 0);
                }
                var si: usize = 0;
                while (si < 16) : (si += 1) {
                    const seg = wi * 16 + si;
                    if (seg >= cx.nsegs) break;
                    var c: u64 = 0;
                    for (bits[si * segww ..][0..segww]) |word| c += @popCount(word);
                    cx.bpi[seg] = c;
                }
                prog_bpi.tick(wlo);
            }
        }
    };
    prog_bpi.begin("bpi", nwin);
    const wbits = try gpa.alloc(u64, nthreads * nww);
    defer gpa.free(wbits);
    const threads = try gpa.alloc(std.Thread, nthreads);
    defer gpa.free(threads);
    var spawned: usize = 0;
    for (1..nthreads) |i| {
        const cpu: ?u32 = if (pins) |pp| pp[i] else null;
        threads[i] = std.Thread.spawn(.{}, Worker.run, .{ &ctx, wbits[i * nww ..][0..nww], cpu }) catch break;
        spawned = i;
    }
    Worker.run(&ctx, wbits[0..nww], if (pins) |pp| pp[0] else null);
    var j: usize = 1;
    while (j <= spawned) : (j += 1) threads[j].join();

    // exclusive scan: slot k becomes π(k·segw − 1) = 3 (for 2,3,5) + coprimes below
    var acc: u64 = 0;
    for (0..nbpi) |k| {
        const c = bpi[k];
        bpi[k] = if (k == 0) 0 else acc + 3;
        acc += c;
    }
    return .{ .bpi = bpi, .total = acc + 3 };
}

/// Segmented wheel sieve into the oracle bitset. Marks only coprime-30 multiples
/// (3.75× fewer stores than a byte sieve) and writes no prime list, so this both
/// costs less time than sievePrimes(n) and leaves n/30 bytes behind instead of 8·π(n).
fn buildPiOracle(gpa: std.mem.Allocator, n: u64) !PiOracle {
    return buildPiOraclePar(gpa, n, 1, null);
}

/// Parallel oracle build: the strike phase fans segments out via an atomic
/// dispenser (disjoint words, no locks), and the prefix pass becomes per-block
/// popcounts in parallel followed by one serial exclusive scan. The build is the
/// largest SERIAL phase at scale — Amdahl's first bite on a many-core box.
fn buildPiOraclePar(gpa: std.mem.Allocator, n: u64, nthreads: usize, pins: ?[]const u32) !PiOracle {
    const nbits: usize = if (n < 1) 0 else PiOracle.kpos(n);
    const nw = (nbits + 63) / 64;
    const nwp = ((nw + PiOracle.BLKW - 1) / PiOracle.BLKW) * PiOracle.BLKW;
    const bits = try gpa.alloc(u64, @max(nwp, PiOracle.BLKW));
    errdefer gpa.free(bits);
    @memset(bits, ~@as(u64, 0));
    @memset(bits[nw..], 0);
    if (nw > 0 and nbits % 64 != 0) bits[nw - 1] &= (@as(u64, 1) << @intCast(nbits % 64)) - 1;
    if (nbits > 0) bits[0] &= ~@as(u64, 1); // 1 is not prime

    const base = try sievePrimes(u32, gpa, isqrt(n));
    defer gpa.free(base);

    const nseg = (nw + ORACLE_SEGW - 1) / ORACLE_SEGW;
    if (nthreads <= 1 or nseg < 2 * nthreads) {
        var w0: usize = 0;
        while (w0 < nw) : (w0 += ORACLE_SEGW) oracleStrikeSeg(bits, base, n, nw, w0);
    } else {
        const Ctx = struct {
            disp: std.atomic.Value(usize),
            bits: []u64,
            base: []const u32,
            n: u64,
            nw: usize,
            nseg: usize,
        };
        var ctx = Ctx{
            .disp = std.atomic.Value(usize).init(0),
            .bits = bits,
            .base = base,
            .n = n,
            .nw = nw,
            .nseg = nseg,
        };
        const Worker = struct {
            fn run(cx: *Ctx, cpu: ?u32) void {
                if (cpu) |c| pinToCpu(c);
                while (true) {
                    const t = cx.disp.fetchAdd(1, .monotonic);
                    if (t >= cx.nseg) break;
                    oracleStrikeSeg(cx.bits, cx.base, cx.n, cx.nw, t * ORACLE_SEGW);
                }
            }
        };
        const threads = try gpa.alloc(std.Thread, nthreads);
        defer gpa.free(threads);
        var spawned: usize = 0;
        for (1..nthreads) |i| {
            const cpu: ?u32 = if (pins) |pp| pp[i] else null;
            threads[i] = std.Thread.spawn(.{}, Worker.run, .{ &ctx, cpu }) catch break;
            spawned = i;
        }
        Worker.run(&ctx, if (pins) |pp| pp[0] else null);
        var j: usize = 1;
        while (j <= spawned) : (j += 1) threads[j].join();
    }

    // pref[b] = # set bits strictly below block b. Phase 1 (parallel-friendly):
    // per-block popcounts. Phase 2: one serial exclusive scan in place.
    const nblk = bits.len / PiOracle.BLKW;
    const pref = try gpa.alloc(u64, nblk + 1);
    errdefer gpa.free(pref);
    if (nthreads <= 1 or nblk < 1 << 16) {
        for (0..nblk) |b| {
            var c: u64 = 0;
            for (bits[b * PiOracle.BLKW ..][0..PiOracle.BLKW]) |word| c += @popCount(word);
            pref[b] = c;
        }
    } else {
        const PCtx = struct { bits: []const u64, pref: []u64, b0: usize, b1: usize };
        const PW = struct {
            fn run(cx: PCtx) void {
                for (cx.b0..cx.b1) |b| {
                    var c: u64 = 0;
                    for (cx.bits[b * PiOracle.BLKW ..][0..PiOracle.BLKW]) |word| c += @popCount(word);
                    cx.pref[b] = c;
                }
            }
        };
        const chunk = (nblk + nthreads - 1) / nthreads;
        const threads = try gpa.alloc(std.Thread, nthreads);
        defer gpa.free(threads);
        var spawned: usize = 0;
        for (1..nthreads) |i| {
            const b0 = @min(i * chunk, nblk);
            const b1 = @min(b0 + chunk, nblk);
            threads[i] = std.Thread.spawn(.{}, PW.run, .{PCtx{ .bits = bits, .pref = pref, .b0 = b0, .b1 = b1 }}) catch break;
            spawned = i;
        }
        PW.run(.{ .bits = bits, .pref = pref, .b0 = 0, .b1 = @min(chunk, nblk) });
        var j: usize = 1;
        while (j <= spawned) : (j += 1) threads[j].join();
    }
    var acc: u64 = 0;
    for (0..nblk) |b| {
        const c = pref[b];
        pref[b] = acc;
        acc += c;
    }
    pref[nblk] = acc;
    return .{ .bits = bits, .pref = pref, .n = n, .nbits = nbits };
}

// ------------------------------------------------- windowed π (the segmented oracle)
// Increment toward retiring the resident [0, √x] bitset (35 GB at 10^24): π(v)
// answered from a per-thread WINDOW covering exactly the sweep segment [lo, hi).
// This works because the ω sweep already delivers every C-leaf query inside the
// current segment (the m-/q-walk guard is `if (v >= hi) break; if (v >= lo)`), so
// the sweep itself is the window pass — no deferral, no sorting. All C-leaves have
// v ≤ √x, i.e. they live in the first √x/z of the z-line, so the window is built
// only for those early segments (~0.2% of segments at 10^20).
//
// A window is segw = 262080 integers = segw/240 = 1092 words (240 ints per word),
// plus a per-word u32 prefix so a query is O(1): 13 KB per thread, L1-resident.
// The global offset comes from bpi[k] = π(k·segw − 1), a boundary array derived
// (for now) from the resident oracle; the streaming build replaces that source
// when the oracle is capped at y in the final increment.

/// Per-thread window over one sweep segment. count(vv) is exact for vv in
/// [lo, hi) with vv ≤ √x.
const PiWin = struct {
    bits: []u64, // segw/240 words, coprime-30, 1 = prime
    pref: []u32, // pref[w] = # set bits in words [0, w)
    lo: u64, // segment base (240-aligned; segw multiple)
    base: u64, // π(lo − 1), full (includes 2, 3, 5 once lo > 0)

    fn count(self: *const PiWin, vv: u64) u64 {
        var c: u64 = self.base;
        if (self.lo == 0) { // 2, 3, 5 have no wheel bit; oracle count() fixes up too
            if (vv >= 2) c += 1;
            if (vv >= 3) c += 1;
            if (vv >= 5) c += 1;
        }
        const k = PiOracle.kpos(vv) - @as(usize, @intCast(self.lo / 30)) * 8;
        const w = k >> 6;
        c += self.pref[w];
        const r: u32 = @intCast(k & 63);
        if (r != 0) c += @popCount(self.bits[w] & ((@as(u64, 1) << @intCast(r)) - 1));
        return c;
    }
};

/// Build the window for segment [lo, hi): MASK30 fill (phase from the GLOBAL word
/// index), cursor-driven strikes for base primes 7.. (indices 3..nwb — the cursors
/// live in the caller and carry across consecutive segments, so the per-window
/// fast-forward cost is paid once per block, not once per segment), then the
/// prefix. Bits beyond min(hi, √x) may be wrong (their striking primes exceed the
/// base set) — harmless, queries never reach them.
fn pwinBuild(comptime P: type, win: *PiWin, primes: []const P, nwb: usize, pw_next: []u64, pw_wpos: []u8, lo: u64, hi: u64) void {
    const nwin = win.bits.len;
    const gw0: usize = @intCast(lo / 240);
    // Wheel-COMPRESSED bits: every position is a coprime-30 integer, so the blank
    // slate is all-ones (NOT the counter's dense MASK30 pattern — that mistake
    // undercounts π everywhere and the suite differential exists to catch it).
    @memset(win.bits, ~@as(u64, 0));
    if (lo == 0) win.bits[0] &= ~@as(u64, 1); // 1 is not prime
    win.lo = lo;

    for (3..nwb) |bi| {
        const p: u64 = @intCast(primes[bi]);
        var j: u64 = pw_next[bi];
        var wp: u8 = pw_wpos[bi];
        while (j < hi) {
            const k = PiOracle.kidx(j) - gw0 * 64;
            win.bits[k >> 6] &= ~(@as(u64, 1) << @intCast(k & 63));
            j += p * W30GAP[wp];
            wp = (wp + 1) & 7;
        }
        pw_next[bi] = j;
        pw_wpos[bi] = wp;
    }

    var acc: u32 = 0;
    for (0..nwin) |w| {
        win.pref[w] = acc;
        acc += @popCount(win.bits[w]);
    }
    win.pref[nwin] = acc;
}

/// Fast-forward the window strike cursors to block_lo — the window twin of the
/// fold's next[]/wpos[] init.
fn pwinInitCursors(comptime P: type, primes: []const P, nwb: usize, pw_next: []u64, pw_wpos: []u8, block_lo: u64) void {
    for (3..nwb) |bi| {
        const p: u64 = @intCast(primes[bi]);
        var m0: u64 = if (block_lo == 0) p else @max(p, (block_lo + p - 1) / p);
        while (!COP30[@intCast(m0 % 30)]) m0 += 1;
        pw_next[bi] = p * m0;
        pw_wpos[bi] = W30IDX[@intCast(m0 % 30)];
    }
}

/// Descending chunked primality window for B's prevPrime walk over (y, √x].
/// Unlike PiWin (which rides the ascending sweep and carries strike cursors),
/// refills jump DOWNWARD, so each fill computes fresh per-prime offsets — paid
/// once per chunk, and B's cursor is monotone within a block so each block fills
/// its own p-subrange exactly once. lo > hi marks "invalid" (set at block start).
const BWin = struct {
    bits: []u64,
    lo: u64 = 1,
    hi: u64 = 0,
};

fn bwinFill(comptime P: type, bits: []u64, primes: []const P, nwb: usize, lo: u64, hi: u64) void {
    @memset(bits, ~@as(u64, 0));
    if (lo == 0) bits[0] &= ~@as(u64, 1); // 1 is not prime
    const gw0: usize = @intCast(lo / 240);
    for (3..nwb) |bi| {
        const p: u64 = @intCast(primes[bi]);
        if (p * p >= hi) break;
        var m0: u64 = @max(p, (lo + p - 1) / p);
        while (!COP30[@intCast(m0 % 30)]) m0 += 1;
        var wp: u8 = W30IDX[@intCast(m0 % 30)];
        var j: u64 = p * m0;
        while (j < hi) {
            const k = PiOracle.kidx(j) - gw0 * 64;
            bits[k >> 6] &= ~(@as(u64, 1) << @intCast(k & 63));
            j += p * W30GAP[wp];
            wp = (wp + 1) & 7;
        }
    }
}

/// Largest prime ≤ v0, walking the window and refilling downward as needed.
/// Exact while every composite ≤ v0 has a factor ≤ primes[nwb-1] (v0 ≤ x^(1/2)
/// with nwb = π(x^(1/4))). Once the window floor reaches `floor_hint` the search
/// stops (the caller's `pB > y` guard makes anything lower irrelevant) and 0 is
/// returned.
fn bwinPrev(comptime P: type, bw: *BWin, primes: []const P, nwb: usize, v0: u64, floor_hint: u64) u64 {
    if (v0 < 7) return if (v0 >= 5) 5 else if (v0 >= 3) 3 else if (v0 >= 2) 2 else 0;
    const span: u64 = @as(u64, @intCast(bw.bits.len)) * 240;
    var v = v0;
    while (true) {
        if (v < bw.lo or v >= bw.hi) {
            const whi = (v / 240) * 240 + 240;
            const wlo = whi -| span;
            bwinFill(P, bw.bits, primes, nwb, wlo, whi);
            bw.lo = wlo;
            bw.hi = whi;
        }
        const gk0 = @as(usize, @intCast(bw.lo / 30)) * 8; // global bit index of window start
        const k = PiOracle.kpos(v) - gk0; // one past v's bit, window-local
        if (k > 0) {
            var w = (k - 1) >> 6;
            const r: u32 = @intCast(k & 63);
            var m = bw.bits[w];
            if (r != 0) m &= (@as(u64, 1) << @intCast(r)) - 1;
            while (true) {
                if (m != 0) return PiOracle.val((w << 6) + 63 - @clz(m) + gk0);
                if (w == 0) break;
                w -= 1;
                m = bw.bits[w];
            }
        }
        if (bw.lo <= floor_hint) return 0; // nothing relevant remains below
        v = bw.lo - 1; // contiguous descent into the next lower chunk
    }
}

const Terms = struct {
    A: i128,
    sig4: i128,
    sig5: i128,
    sig6: i128,
    a: i128,
    b: i128,
    c: i128,
    d: i128,
    P: i128,
};

/// A, Σ₄/₅/₆ and the scalars query the π oracle (all query points ≤ √x — proven).
/// Only B has v = x/p up to z > √x, so B alone takes a segmented sweep over its
/// ~π(√x) points. π-free terms (ω, φ₀, Σ₀–₃) are handled in piGourdonV.
fn computeTerms(comptime X: type, s: anytype, x: X, y: u64, sqx: u64, x13: u64, sqz: u64, xstar: u64) Terms {
    const pr = s.primes;
    const pio = &s.pio;
    const a: i128 = @intCast(pio.count(y));
    const b: i128 = @intCast(pio.count(x13));
    const c: i128 = @intCast(pio.count(sqz));
    const d: i128 = @intCast(pio.count(xstar));
    const Pi: i128 = @intCast(pio.count(sqx));

    // A = Σ_{x*<p≤x^1/3} Σ_{p<q≤√(x/p)} χ(x/pq) π(x/pq),  χ = 2 if x/pq<y else 1.
    // (q ≤ √(x/x*) = y, so the explicit list suffices for the inner enumeration.)
    var A: i128 = 0;
    for (pr, 0..) |p32, pidx| {
        const p: u64 = @intCast(p32);
        if (p <= xstar) continue;
        if (p > x13) break;
        const qhi = isqrtG(X, x / @as(X, p));
        var qi = pidx + 1;
        while (qi < pr.len) : (qi += 1) {
            const qq: u64 = @intCast(pr[qi]);
            if (qq > qhi) break;
            const v = xdiv(X, x, p * qq); // < √x
            const chi: i128 = if (v < y) 2 else 1;
            A += chi * @as(i128, @intCast(pio.count(v)));
        }
    }
    // Σ₄ = a·Σ_{x*<p≤√(x/y)} π(x/py)   (v ≤ x^{5/12} < √x)
    var sig4: i128 = 0;
    for (pr) |p32| {
        const p: u64 = @intCast(p32);
        if (p <= xstar) continue;
        if (p > sqz) break;
        sig4 += @intCast(pio.count(xdiv(X, x, p * y)));
    }
    sig4 *= a;
    // Σ₅ = Σ_{√(x/y)<p≤x^1/3} π(x/p²)   (v < y ≤ √x)
    var sig5: i128 = 0;
    for (pr) |p32| {
        const p: u64 = @intCast(p32);
        if (p <= sqz) continue;
        if (p > x13) break;
        sig5 += @intCast(pio.count(xdiv(X, x, p * p)));
    }
    // Σ₆ = −Σ_{x*<p≤x^1/3} π(√(x/p))²   (v ≤ x^{3/8} < √x)
    var sig6: i128 = 0;
    for (pr) |p32| {
        const p: u64 = @intCast(p32);
        if (p <= xstar) continue;
        if (p > x13) break;
        const t: i128 = @intCast(pio.count(isqrt(xdiv(X, x, p))));
        sig6 += t * t;
    }
    return .{ .A = A, .sig4 = sig4, .sig5 = sig5, .sig6 = -sig6, .a = a, .b = b, .c = c, .d = d, .P = Pi };
}

const Partial = struct { A: i128 = 0, sig4: i128 = 0, sig5: i128 = 0, sig6: i128 = 0 };

/// Parallel A/Σ (Model-A phase 1), in two sub-phases per worker:
///   1. p-chunks from an atomic dispenser: Σ₄/Σ₅/Σ₆ (every query ≤ y, since
///      p > x* ≥ x/y² bounds them — proven per-term) and A's v < y pairs (χ = 2).
///   2. v-windows over [y, √x] from a second dispenser: A's v ≥ y pairs (χ = 1),
///      π(v) answered from a freshly sieved wheel window + per-word prefix, based
///      at bpi[window/segw]. Threads own disjoint windows; for each window the
///      (p, q) pairs mapping into it are enumerated off the prime list (q < y
///      always, so the list suffices).
/// The resident oracle serves only v ≤ y queries and the scalars — everything the
/// y-capped oracle of the final increment can still answer. π(√x) arrives as a
/// parameter (bpi's builder knows it).
fn computeTermsPar(comptime X: type, gpa: std.mem.Allocator, s: anytype, x: X, y: u64, sqx: u64, x13: u64, sqz: u64, xstar: u64, nthreads: usize, pins: ?[]const u32, bpi: []const u64, pi_sqx: i128, segw: usize) !Terms {
    const pio = &s.pio;
    const a: i128 = @intCast(pio.count(y));
    const b: i128 = @intCast(pio.count(x13));
    const c: i128 = @intCast(pio.count(sqz));
    const d: i128 = @intCast(pio.count(xstar));
    const Pi: i128 = pi_sqx;

    const i_lo: usize = @intCast(pio.count(xstar)); // first prime index with p > x*
    const i_hi: usize = @intCast(pio.count(x13)); // one past last prime ≤ x^1/3
    const n_primes = if (i_hi > i_lo) i_hi - i_lo else 0;
    const nchunks = @max(@as(usize, 1), @min(n_primes, nthreads * 16)); // over-partition: small p ⇒ more work
    const csize = (n_primes + nchunks - 1) / nchunks;

    const AWIN: u64 = 16 * @as(u64, segw); // 16 sweep segments: L2-scale windows

    const Ctx = struct {
        disp: std.atomic.Value(usize),
        disp2: std.atomic.Value(usize),
        s: @TypeOf(s),
        x: X,
        y: u64,
        x13: u64,
        sqz: u64,
        xstar: u64,
        sqx: u64,
        i_lo: usize,
        i_hi: usize,
        csize: usize,
        nchunks: usize,
        w0: u64,
        nwb: usize,
        bpi: []const u64,
        segw: usize,
    };
    var ctx = Ctx{
        .disp = std.atomic.Value(usize).init(0),
        .disp2 = std.atomic.Value(usize).init(0),
        .s = s,
        .x = x,
        .y = y,
        .x13 = x13,
        .sqz = sqz,
        .xstar = xstar,
        .sqx = sqx,
        .i_lo = i_lo,
        .i_hi = i_hi,
        .csize = csize,
        .nchunks = nchunks,
        .w0 = (y / AWIN) * AWIN,
        .nwb = @max(3, @as(usize, @intCast(pio.count(isqrt(sqx))))),
        .bpi = bpi,
        .segw = segw,
    };
    const partials = try gpa.alloc(Partial, nthreads);
    defer gpa.free(partials);
    @memset(partials, Partial{});

    const Worker = struct {
        fn run(cx: *Ctx, out: *Partial, wbits: []u64, wpref: []u32, cpu: ?u32) void {
            if (cpu) |cp| pinToCpu(cp);
            const prr = cx.s.primes;
            const PP = std.meta.Child(@TypeOf(prr));
            const po = &cx.s.pio;
            var A: i128 = 0;
            var s4: i128 = 0;
            var s5: i128 = 0;
            var s6: i128 = 0;
            while (true) {
                const ch = cx.disp.fetchAdd(1, .monotonic);
                if (ch >= cx.nchunks) break;
                const lo = cx.i_lo + ch * cx.csize;
                const hi = @min(cx.i_lo + (ch + 1) * cx.csize, cx.i_hi);
                var idx = lo;
                while (idx < hi) : (idx += 1) {
                    const p: u64 = @intCast(prr[idx]);
                    if (p > cx.xstar and p <= cx.x13) {
                        // A, v < y pairs only (χ = 2): q > x/(p·y) ⟺ v < y, up to
                        // division slop rechecked per q. v ≥ y pairs → window phase.
                        const qhi = isqrtG(X, cx.x / @as(X, p));
                        var qi = lowerBound(PP, prr, mBound(X, cx.x, p, cx.y) + 1);
                        if (qi <= idx) qi = idx + 1;
                        while (qi < prr.len) : (qi += 1) {
                            const qq: u64 = @intCast(prr[qi]);
                            if (qq > qhi) break;
                            const v = xdiv(X, cx.x, p * qq);
                            if (v >= cx.y) continue;
                            A += 2 * @as(i128, @intCast(po.count(v)));
                        }
                        const t: i128 = @intCast(po.count(isqrt(xdiv(X, cx.x, p)))); // Σ₆
                        s6 += t * t;
                    }
                    if (p > cx.xstar and p <= cx.sqz) s4 += @intCast(po.count(xdiv(X, cx.x, p * cx.y))); // Σ₄
                    if (p > cx.sqz and p <= cx.x13) s5 += @intCast(po.count(xdiv(X, cx.x, p * p))); // Σ₅
                }
                prog_ap.tick(if (hi > lo) @intCast(prr[hi - 1]) else 0);
            }
            // Phase 2 — A's v ∈ [y, √x] pairs (χ = 1) over disjoint v-windows.
            const awin: u64 = @intCast(wbits.len * 240);
            while (true) {
                const wi = cx.disp2.fetchAdd(1, .monotonic);
                const wlo = cx.w0 + @as(u64, wi) * awin;
                if (wlo > cx.sqx) break;
                const whi = @min(wlo + awin, cx.sqx + 1);
                bwinFill(PP, wbits, prr, cx.nwb, wlo, whi);
                var acc: u32 = 0;
                for (0..wbits.len) |w2| {
                    wpref[w2] = acc;
                    acc += @popCount(wbits[w2]);
                }
                wpref[wbits.len] = acc;
                const win = PiWin{ .bits = wbits, .pref = wpref, .lo = wlo, .base = cx.bpi[@intCast(wlo / cx.segw)] };
                const vmin = @max(wlo, cx.y);
                var pidx = cx.i_lo;
                while (pidx < cx.i_hi) : (pidx += 1) {
                    const p: u64 = @intCast(prr[pidx]);
                    // beyond this p even q = p + 1 lands below the window
                    if (@as(u128, p) * p * vmin > @as(u128, cx.x)) break;
                    const qmax = @min(isqrtG(X, cx.x / @as(X, p)), mBound(X, cx.x, p, vmin));
                    var qi = lowerBound(PP, prr, mBound(X, cx.x, p, whi) + 1);
                    if (qi <= pidx) qi = pidx + 1;
                    while (qi < prr.len) : (qi += 1) {
                        const qq: u64 = @intCast(prr[qi]);
                        if (qq > qmax) break;
                        const v = xdiv(X, cx.x, p * qq);
                        if (v < vmin or v >= whi) continue;
                        A += @intCast(win.count(v));
                    }
                }
                prog_aw.tick(wlo);
            }
            out.* = .{ .A = A, .sig4 = s4, .sig5 = s5, .sig6 = s6 };
        }
    };

    prog_ap.begin("A/Σ p", nchunks);
    prog_aw.begin("A/Σ v", @intCast((sqx - ctx.w0) / AWIN + 1));
    const nww = (16 * segw) / 240;
    const wbits_flat = try gpa.alloc(u64, nthreads * nww);
    defer gpa.free(wbits_flat);
    const wpref_flat = try gpa.alloc(u32, nthreads * (nww + 1));
    defer gpa.free(wpref_flat);

    const threads = try gpa.alloc(std.Thread, nthreads);
    defer gpa.free(threads);
    var spawned: usize = 0;
    for (1..nthreads) |i| {
        const cpu: ?u32 = if (pins) |pp| pp[i] else null;
        threads[i] = std.Thread.spawn(.{}, Worker.run, .{ &ctx, &partials[i], wbits_flat[i * nww ..][0..nww], wpref_flat[i * (nww + 1) ..][0 .. nww + 1], cpu }) catch break;
        spawned = i;
    }
    Worker.run(&ctx, &partials[0], wbits_flat[0..nww], wpref_flat[0 .. nww + 1], if (pins) |pp| pp[0] else null);
    var j: usize = 1;
    while (j <= spawned) : (j += 1) threads[j].join();

    var A: i128 = 0;
    var sig4: i128 = 0;
    var sig5: i128 = 0;
    var sig6: i128 = 0;
    for (partials) |pt| {
        A += pt.A;
        sig4 += pt.sig4;
        sig5 += pt.sig5;
        sig6 += pt.sig6;
    }
    sig4 *= a;
    return .{ .A = A, .sig4 = sig4, .sig5 = sig5, .sig6 = -sig6, .a = a, .b = b, .c = c, .d = d, .P = Pi };
}

/// B = Σ_{y<p≤√x} π(x/p). The only term needing π(v) for v > √x (v ∈ [√x, z]).
/// Answered by a running-π segmented sweep over just [√x, z], starting at π(√x).
fn computeB(comptime INST: bool, gpa: std.mem.Allocator, s: anytype, x: u64, y: u64, z: u64, sqx: u64) !i128 {
    const pr = s.primes;
    var nb: usize = 0;
    for (pr) |p32| {
        const p: u64 = @intCast(p32);
        if (p <= y) continue;
        if (p > sqx) break;
        nb += 1;
    }
    const bpts = try gpa.alloc(u64, nb);
    defer gpa.free(bpts);
    const bans = try gpa.alloc(u64, nb);
    defer gpa.free(bans);
    var k: usize = 0;
    for (pr) |p32| {
        const p: u64 = @intCast(p32);
        if (p <= y) continue;
        if (p > sqx) break;
        bpts[k] = x / p;
        k += 1;
    }
    try answerPi(INST, gpa, pr, bpts, bans, z);
    var B: i128 = 0;
    for (bans) |pv| B += @intCast(pv);
    return B;
}

// ------------------------------------------------------------------------ ω+B
/// FUSED ω and B on ONE segmented counter. The counter folds every prime ≤ √z (not
/// just ≤ x*): ω leaves are queried at their stage bi < π(x*) DURING the fold, and
/// after the full fold the counter holds φ(·,π(√z)), so B's π(x/p) = φ(x/p,π(√z)) +
/// π(√z) − 1 is read straight off it (Legendre; valid since z < y² ⇒ no product of
/// two primes > √z is ≤ z). Gourdon folds π(√z) primes vs lmo's π(y) — ~6× fewer at
/// 10^12 — so one pass serves both, beating lmo's separate S2-fold + fused-P₂.
/// C/D leaf evaluator: closed form (C, is_d=false) when x/pm < p², else the counter
/// (D, is_d=true — the leaf needs a cross-block φ correction). π(v) comes from the
/// oracle; the old y-sized pi_tab direct-index table was measured at 0.7% (inside
/// noise) against 196 MB at 10^20, so it is gone.
inline fn evalPhi(comptime INST: bool, st: *Stats, vv: u64, b: usize, pp: u64, pw: *const PiWin, pr_bi: i64, ct: *const Ctr, l: u64, sx: u64) struct { phi: i64, is_d: bool } {
    var piv: i64 = -1;
    if (pp * pp > vv) {
        // vv ∈ [lo, hi) by the walk guards, and vv ≤ √x ⇒ the window was built
        if (vv <= sx) piv = @intCast(pw.count(vv));
    }
    if (piv >= 0) {
        if (INST) st.n_easy += 1;
        return .{ .phi = 1 + @max(0, piv - @as(i64, @intCast(b))), .is_d = false };
    }
    if (INST) st.n_hard += 1;
    return .{ .phi = pr_bi + ct.prefix(@intCast(vv - l)), .is_d = true };
}

const OmegaB = struct { omega: i128, b: i128 };

/// ω+B instrumentation counters (INST-gated). Per-thread, summed after join.
const Stats = struct {
    n_seg: u64 = 0,
    n_mwalk: u64 = 0,
    n_small: u64 = 0,
    n_easy: u64 = 0,
    n_hard: u64 = 0,
    n_kill: u64 = 0,
    n_bq: u64 = 0,
    fn add(self: *Stats, o: Stats) void {
        self.n_seg += o.n_seg;
        self.n_mwalk += o.n_mwalk;
        self.n_small += o.n_small;
        self.n_easy += o.n_easy;
        self.n_hard += o.n_hard;
        self.n_kill += o.n_kill;
        self.n_bq += o.n_bq;
    }
};

/// Per-thread scratch: the counter + all cursors.
/// One pending strike for a bucketed prime. Copied into the bucket by VALUE rather
/// than referenced by prime index: a bucket of indices would touch next[]/wpos[]/
/// primes[] at scattered offsets and move MORE bytes than the sequential scan it
/// replaces (~2.2 MB vs 1.7 MB per segment at 10^20). Contiguous 16-byte entries
/// bring that to ~140 KB. p ≤ √z < 2^32 for any x this code can address.
/// One pending strike, packed into a u64: [p : 43][wp : 3][offset-in-slot-segment
/// : 18]. Offset relative to the SLOT's segment base (< segw = 2^18 always, and
/// √z < 2^32 fits 43 bits) — scale-safe with no guard.
///
/// Storage is contiguous per-slot VECTORS, not an arena of chains. The arena
/// (8 B/prime exact, chain-linked) measured time-neutral at 10^18-19 — where it
/// fit in L3 — and cost +33% at 10^22 (22671 s vs 17072 s pre-refactor, perf:
/// one indexed u32 chain-load = 37.6% of ALL cycles): ~3 MB/thread of chains ×6
/// threads crossed the 16 MB L3 exactly between 10^21 (9.6 MB, fast) and 10^22
/// (18 MB, cliff), turning ~10^12 refile events into dependent DRAM round-trips.
/// Vectors drain SEQUENTIALLY — the prefetcher's favorite shape — at the price
/// of doubling slack (tens of MB at 10^22; noise against the hours reclaimed).
/// Every performance verdict carries an "...at the scales measured" rider; this
/// one is stamped for 10^22.
const Ent = u64;

inline fn entPack(off: u64, wp: u8, p: u64) Ent {
    return off | (@as(u64, wp) << 21) | (p << 24);
}

const Bucket = struct {
    items: []Ent = &.{},
    len: usize = 0,
    inline fn push(self: *Bucket, gpa: std.mem.Allocator, e: Ent) !void {
        if (self.len == self.items.len) {
            self.items = try gpa.realloc(self.items, if (self.items.len == 0) 64 else self.items.len * 2);
        }
        self.items[self.len] = e;
        self.len += 1;
    }
};

const Scratch = struct {
    ctr: Ctr,
    cur: []u64,
    next: []u64,
    wpos: []u8,
    phi_run: []i64,
    buck: []Bucket, // ring of per-slot vectors, indexed by (segment index & ring_mask)
    pwin: PiWin, // per-thread π window over the current sweep segment
    pw_next: []u64, // window strike cursors (base primes ≤ x^(1/4))
    pw_wpos: []u8,
    bwin: BWin, // descending prime window for B's cursor
    fn init(gpa: std.mem.Allocator, nax: usize, naz: usize, segw: usize, nring: usize, nwb: usize) !Scratch {
        const nwin = segw / 240;
        return .{
            .ctr = try Ctr.init(gpa, segw),
            .cur = try gpa.alloc(u64, nax),
            .next = try gpa.alloc(u64, naz),
            .wpos = try gpa.alloc(u8, naz),
            .phi_run = try gpa.alloc(i64, nax),
            .buck = blk: {
                const b = try gpa.alloc(Bucket, nring);
                @memset(b, Bucket{});
                break :blk b;
            },
            .pwin = .{ .bits = try gpa.alloc(u64, nwin), .pref = try gpa.alloc(u32, nwin + 1), .lo = 0, .base = 0 },
            .pw_next = try gpa.alloc(u64, nwb),
            .pw_wpos = try gpa.alloc(u8, nwb),
            .bwin = .{ .bits = try gpa.alloc(u64, nwin) },
        };
    }
    fn deinit(self: *Scratch, gpa: std.mem.Allocator) void {
        self.ctr.deinit(gpa);
        gpa.free(self.cur);
        gpa.free(self.next);
        gpa.free(self.wpos);
        gpa.free(self.phi_run);
        for (self.buck) |*b| if (b.items.len > 0) gpa.free(b.items);
        gpa.free(self.buck);
        gpa.free(self.pwin.bits);
        gpa.free(self.pwin.pref);
        gpa.free(self.pw_next);
        gpa.free(self.pw_wpos);
        gpa.free(self.bwin.bits);
    }
};

/// Shared read-only context + per-block output arrays + block dispenser.
fn BlkCtx(comptime X: type, comptime P: type) type {
    return struct {
        leaf: []const u16,
        primes: []const P,
        pio: *const PiOracle,
        x: X,
        y: u64,
        sqx: u64,
        sqrt_y: u64,
        xstar: u64,
        total: u64, // z+1
        nax: usize,
        naz: usize,
        nay: usize,
        segw: usize,
        nseg: usize,
        nb: usize,
        naz_i: i64,
        blk_total: []i64,
        blk_mu: []i64,
        blk_total_full: []i64,
        blk_bcount: []u64,
        blk_omega: []i128,
        blk_b: []i128,
        bpi: []const u64, // π(k·segw − 1) at segment boundaries below √x (bpi[0] = 0)
        bpi_owned: bool, // tests derive their own; piGourdonV's is shared
        nwb: usize, // # base primes for window strikes (π(x^(1/4)); indices 3.. used)
        nbuck: usize, // first prime index handled by the bucket ring (naz ⇒ disabled)
        nring: usize, // ring slots, power of two ≥ max segments a prime can skip
        ring_mask: usize,
        nsmall: usize, // # of small primes (bi = 3 .. 3+nsmall) folded word-parallel
        small_tmpl: []const u64, // flat p-word templates per small prime
        small_off: []const usize, // small_off[k] = start of prime k's templates
        disp: std.atomic.Value(usize),
    };
}

/// Small primes p (excluding 2,3,5) get a word-parallel strike when their coprime-30
/// multiples are ≥ ~1 per 64-bit word, i.e. p ≲ 64·8/30 ≈ 17. A couple above break-even
/// are harmless. Template p, phase r has bit j set iff j ≡ r (mod p), j∈[0,63].
const SMALL_STRIKE_MAX: u64 = 30;

/// Build flat coprime-agnostic multiple-of-p templates for primes[3 .. 3+nsmall].
fn buildSmallTemplates(gpa: std.mem.Allocator, primes: anytype) !struct { tmpl: []u64, off: []usize, n: usize } {
    var n: usize = 0;
    var total: usize = 0;
    while (3 + n < primes.len and @as(u64, primes[3 + n]) <= SMALL_STRIKE_MAX) : (n += 1) {
        total += @intCast(primes[3 + n]);
    }
    const off = try gpa.alloc(usize, n + 1);
    const tmpl = try gpa.alloc(u64, total);
    var cur: usize = 0;
    for (0..n) |k| {
        off[k] = cur;
        const p: usize = @intCast(primes[3 + k]);
        for (0..p) |r| { // phase r: bits at j ≡ r (mod p)
            var word: u64 = 0;
            var j: usize = r;
            while (j < 64) : (j += p) word |= @as(u64, 1) << @as(u6, @intCast(j));
            tmpl[cur + r] = word;
        }
        cur += p;
    }
    off[n] = cur;
    return .{ .tmpl = tmpl, .off = off, .n = n };
}

inline fn nfold_c_init(nax: usize, naz: usize) usize {
    return @min(nax, naz);
}

/// Fold prime index bi out of the current segment. `counted` selects whether cnt/
/// total are maintained on every strike. That bookkeeping exists so the counter can
/// be READ mid-fold, which only the stages below π(x*) ever do — above it, clearing
/// the bit is the whole job and one rebuild() at the end restores the counts.
inline fn foldPrime(comptime INST: bool, comptime counted: bool, ctx: anytype, ctr: *Ctr, st: *Stats, bi: usize, p: u64, lo: u64, hi: u64, len: usize, next: []u64, wpos: []u8) void {
    const k = bi - 3;
    if (k < ctx.nsmall) {
        // word-parallel strike: AND out all multiples of p across the segment's
        // words (2,3,5-multiples already 0 ⇒ popcount delta counts only live kills)
        const tp = ctx.small_tmpl[ctx.small_off[k]..ctx.small_off[k + 1]];
        const pp: u64 = @intCast(tp.len); // = p
        const nwlen = (len + 63) >> 6;
        var r: usize = @intCast((pp - (lo % pp)) % pp); // phase of word 0 (pos 0 = int lo)
        const dr: usize = @intCast(pp - (64 % pp)); // r -= 64 mod p per word
        const ppz: usize = @intCast(pp);
        for (0..nwlen) |w| {
            if (counted) ctr.killWord(w, tp[r]) else ctr.strikeWord(w, tp[r]);
            r += dr;
            if (r >= ppz) r -= ppz;
        }
        if (INST) st.n_kill += nwlen;
    } else {
        var j: u64 = next[bi];
        var wp: u8 = wpos[bi];
        while (j < hi) {
            if (INST) st.n_kill += 1;
            if (counted) ctr.kill(@intCast(j - lo)) else ctr.strike(@intCast(j - lo));
            j += p * W30GAP[wp];
            wp = (wp + 1) & 7;
        }
        next[bi] = j;
        wpos[bi] = wp;
    }
}

/// Sweep block t (segw-aligned) with fast-forwarded cursors + LOCAL phi_run; write the
/// reduction data into ctx.blk_*[t] (disjoint per t ⇒ concurrency-safe across threads).
fn runOneBlock(comptime INST: bool, comptime X: type, comptime P: type, gpa: std.mem.Allocator, ctx: *BlkCtx(X, P), sc: *Scratch, st: *Stats, t: usize) !void {
    const primes = ctx.primes;
    const x = ctx.x;
    const y = ctx.y;
    const sqx = ctx.sqx;
    const sqrt_y = ctx.sqrt_y;
    const nax = ctx.nax;
    const naz = ctx.naz;
    const segw = ctx.segw;
    const nseg = ctx.nseg;
    const nb = ctx.nb;
    const total = ctx.total;
    const naz_i = ctx.naz_i;
    const cur = sc.cur;
    const next = sc.next;
    const wpos = sc.wpos;
    const phi_run = sc.phi_run;

    const block_lo = @min((t * nseg / nb) * segw, total);
    const block_hi = @min(((t + 1) * nseg / nb) * segw, total);

    // ── fast-forward all cursors to block_lo ──
    @memset(phi_run, 0);
    for (0..nax) |bi| {
        const p: u64 = @intCast(primes[bi]);
        const bnd: u64 = if (block_lo == 0) y else @min(y, xdiv(X, x, p * block_lo));
        if (p <= sqrt_y) {
            cur[bi] = bnd;
        } else {
            const pl: usize = @intCast(ctx.pio.count(bnd));
            cur[bi] = if (pl > 0) @as(u64, pl - 1) else 0;
        }
    }
    for (3..naz) |bi| {
        const p: u64 = @intCast(primes[bi]);
        var m0: u64 = if (block_lo == 0) 1 else (block_lo + p - 1) / p;
        if (m0 == 0) m0 = 1;
        while (!COP30[@intCast(m0 % 30)]) m0 += 1;
        next[bi] = p * m0;
        wpos[bi] = W30IDX[@intCast(m0 % 30)];
    }
    // Window strike cursors for the C-leaf π windows (only blocks touching [0, √x]
    // ever build one — all C-leaves have v ≤ √x).
    if (block_lo <= sqx) pwinInitCursors(P, primes, ctx.nwb, sc.pw_next, sc.pw_wpos, block_lo);
    sc.bwin.lo = 1; // invalidate: this block's cursor starts a fresh descent
    sc.bwin.hi = 0;

    // Seed the bucket ring: each sparse prime is filed under the segment its first
    // multiple lands in, so a segment later touches only the primes that hit it.
    const nbuck = ctx.nbuck;
    const rmask = ctx.ring_mask;
    for (sc.buck) |*b| b.len = 0;
    for (nbuck..naz) |bi| {
        if (next[bi] >= block_hi) continue;
        const si = (next[bi] - block_lo) / segw;
        const base = block_lo + @as(u64, @intCast(si)) * segw;
        try sc.buck[@as(usize, @intCast(si)) & rmask].push(gpa, entPack(next[bi] - base, wpos[bi], @intCast(primes[bi])));
    }

    // B's descending cursor over the primes of (y, √x] — the only loop that needs
    // prime VALUES above the explicit list, so it walks the oracle bitset instead.
    var pB: u64 = if (block_lo == 0) sqx else @min(sqx, xdiv(X, x, block_lo));
    pB = bwinPrev(P, &sc.bwin, primes, ctx.nwb, pB, y);

    const omega_mu = ctx.blk_mu[t * nax ..][0..nax];
    @memset(omega_mu, 0);
    var omega: i128 = 0;
    var b_sum: i128 = 0;
    var phi_run_full: i64 = 0;
    var b_count: u64 = 0;

    // x* ≤ √z always (y ≤ √x and x ≤ y³), so the counted stages are a prefix of the
    // folded ones; the @min is defensive.
    const nfold_c = @min(nax, naz);

    var lo: u64 = block_lo;
    while (lo < block_hi) : (lo += segw) {
        const hi = @min(lo + @as(u64, segw), block_hi);
        const len: usize = @intCast(hi - lo);
        sc.ctr.reset(len);
        if (INST) st.n_seg += 1;
        if (lo <= sqx) { // π window for this segment's C-leaves
            sc.pwin.base = ctx.bpi[@intCast(lo / segw)];
            pwinBuild(P, &sc.pwin, primes, ctx.nwb, sc.pw_next, sc.pw_wpos, lo, hi);
        }
        for (0..nfold_c) |bi| {
            const p: u64 = @intCast(primes[bi]);
            {
                if (p > 2 and p <= sqrt_y) {
                    // DENSE m-walk. Both segment bounds are algebraic, not arithmetic:
                    //   v = x/(mp) < hi  ⟺  m > x/(p·hi),   v ≥ lo  ⟺  m ≤ x/(p·lo)
                    // so the descent range costs two divisions per (p, segment) rather
                    // than one per m. That matters because ~63% of the m are rejected
                    // (μ=0 or spf ≤ p) and used to pay a full x/(mp) before being
                    // thrown away — a libcall on the u128 path. A rejected m is now a
                    // load and a compare, and only a surviving leaf divides.
                    var m: u64 = @min(cur[bi], mBound(X, x, p, lo));
                    const mlo: u64 = @max(y / p, mBound(X, x, p, hi));
                    while (m > mlo) {
                        if (INST) st.n_mwalk += 1;
                        // one load answers μ≠0 AND spf>p: (leaf & SAT) > bi+1
                        const lv = ctx.leaf[@intCast(m)];
                        if (@as(usize, lv & LEAF_SAT) > bi + 1) {
                            const v: u64 = xdiv(X, x, m * p);
                            const sign: i64 = if (lv & LEAF_PLUS != 0) -1 else 1; // -μ
                            if (bi <= 2) {
                                if (INST) st.n_small += 1;
                                omega += @as(i128, sign) * @as(i128, phiSmall(v, bi));
                            } else {
                                const r = evalPhi(INST, st, v, bi, p, &sc.pwin, phi_run[bi], &sc.ctr, lo, sqx);
                                omega += @as(i128, sign) * @as(i128, r.phi);
                                if (r.is_d) omega_mu[bi] += sign;
                            }
                        }
                        m -= 1;
                    }
                    cur[bi] = m;
                } else if (p > 2) {
                    // SPARSE q-walk (p>√y): every prime q∈(p,y] is a valid leaf (sign +1)
                    var qc: usize = @intCast(cur[bi]);
                    while (qc > bi) {
                        if (INST) st.n_mwalk += 1;
                        const q: u64 = @intCast(primes[qc]);
                        const v: u64 = xdiv(X, x, p * q);
                        if (v >= hi) break;
                        if (v >= lo) {
                            const r = evalPhi(INST, st, v, bi, p, &sc.pwin, phi_run[bi], &sc.ctr, lo, sqx);
                            omega += @as(i128, r.phi);
                            if (r.is_d) omega_mu[bi] += 1;
                        }
                        qc -= 1;
                    }
                    cur[bi] = qc;
                }
            }
            // phi_run accumulates HERE — after this stage's leaves (which never
            // touch the counter, so total equals the old pre-leaf snapshot) and
            // before this stage's fold. The seg_cnt snapshot array is gone.
            if (bi >= 3) phi_run[bi] += sc.ctr.total;
            if (bi >= 3) foldPrime(INST, true, ctx, &sc.ctr, st, bi, p, lo, hi, len, next, wpos);
        }
        // Stages above π(x*) carry no leaves, so nothing reads the counter until the
        // B queries below: strike bits only, then restore cnt/total in one pass.
        for (nfold_c..nbuck) |bi| {
            const p: u64 = @intCast(primes[bi]);
            if (bi >= 3) foldPrime(INST, false, ctx, &sc.ctr, st, bi, p, lo, hi, len, next, wpos);
        }
        if (nbuck < naz) {
            // Drain the primes filed under this segment, refiling each under the
            // segment of its next multiple. Entries are contiguous, so this streams.
            const si: usize = @intCast((lo - block_lo) / segw);
            const b = &sc.buck[si & rmask];
            const n = b.len;
            b.len = 0; // refiling targets a strictly later slot, never this one
            for (b.items[0..n]) |e| {
                var j: u64 = lo + (e & 0x1FFFFF);
                var wp: u8 = @intCast((e >> 21) & 7);
                const pe: u64 = e >> 24;
                while (j < hi) {
                    if (INST) st.n_kill += 1;
                    sc.ctr.strike(@intCast(j - lo));
                    j += pe * W30GAP[wp];
                    wp = (wp + 1) & 7;
                }
                if (j < block_hi) {
                    const sj: usize = @intCast((j - block_lo) / segw);
                    // audit flag 10a: a same-slot refile after realloc would alias
                    // the slice being drained; the ring-span invariant forbids it.
                    std.debug.assert((sj & rmask) != (si & rmask));
                    const sbase = block_lo + @as(u64, @intCast(sj)) * segw;
                    try sc.buck[sj & rmask].push(gpa, entPack(j - sbase, wp, pe));
                }
            }
        }
        if (naz > nfold_c) sc.ctr.rebuild();
        // B queries off the fully-folded counter (φ(·,π√z))
        while (pB > y and xdiv(X, x, pB) < hi) {
            const v: u64 = xdiv(X, x, pB);
            if (v >= lo) {
                if (INST) st.n_bq += 1;
                b_sum += @as(i128, phi_run_full + sc.ctr.prefix(@intCast(v - lo)) + naz_i - 1);
                b_count += 1;
            }
            pB = bwinPrev(P, &sc.bwin, primes, ctx.nwb, pB - 1, y);
        }
        phi_run_full += sc.ctr.total;
        prog_om.tick(lo);
    }
    for (0..nax) |bi| ctx.blk_total[t * nax + bi] = phi_run[bi];
    ctx.blk_total_full[t] = phi_run_full;
    ctx.blk_bcount[t] = b_count;
    ctx.blk_omega[t] = omega;
    ctx.blk_b[t] = b_sum;
}

/// Serial prefix-sum stitch: ω D-leaf correction (O(nb·π(x*))) + B correction.
fn reduceOmB(comptime X: type, comptime P: type, gpa: std.mem.Allocator, ctx: *BlkCtx(X, P)) !OmegaB {
    const nax = ctx.nax;
    const nb = ctx.nb;
    const prefix = try gpa.alloc(i64, nax);
    defer gpa.free(prefix);
    @memset(prefix, 0);
    var omega_total: i128 = 0;
    var b_total: i128 = 0;
    var prefix_full: i64 = 0;
    for (0..nb) |t| {
        omega_total += ctx.blk_omega[t];
        b_total += ctx.blk_b[t];
        var corr: i128 = 0;
        for (3..nax) |bi| corr += @as(i128, prefix[bi]) * @as(i128, ctx.blk_mu[t * nax + bi]);
        omega_total += corr;
        b_total += @as(i128, prefix_full) * @as(i128, @intCast(ctx.blk_bcount[t]));
        for (3..nax) |bi| prefix[bi] += ctx.blk_total[t * nax + bi];
        prefix_full += ctx.blk_total_full[t];
    }
    return .{ .omega = omega_total, .b = b_total };
}

/// Build the BlkCtx (derived sizes + output arrays). null ⇒ no fold primes (tiny x).
fn initBlkCtx(comptime X: type, comptime P: type, gpa: std.mem.Allocator, s: anytype, x: X, y: u64, z: u64, xstar: u64, nb: usize, bpi_in: ?[]const u64, segw_in: usize) !?BlkCtx(X, P) {
    const primes = s.primes;
    const sqx = isqrtG(X, x);
    const sqz = isqrt(z);
    const naz: usize = @intCast(s.pio.count(sqz));
    if (naz == 0) return null;
    const nax: usize = @intCast(s.pio.count(xstar));
    const segw: usize = segw_in;

    // Bucket ring: primes sparse enough to miss most segments. A prime contributes
    // segw·(8/30)/p strikes per segment, so bucketing pays once p > segw·8/30. Only
    // primes above x* qualify (the counted prefix keeps its cursor scan), and the
    // ring must span the furthest a prime can jump: coprime-30 gaps reach 6, so
    // ⌈6p/segw⌉ + 2 slots, rounded up to a power of two for masking.
    const bucket_min: u64 = @max(segw * 8 / 30, 1);
    var nbuck: usize = naz;
    {
        var i = nfold_c_init(nax, naz);
        while (i < naz and @as(u64, primes[i]) <= bucket_min) i += 1;
        if (i < naz) nbuck = i;
    }
    var nring: usize = 8;
    if (nbuck < naz) {
        const span = (6 * @as(u64, primes[naz - 1])) / segw + 4;
        while (@as(u64, nring) < span) nring *= 2;
    }
    const total = z + 1;
    // Window support: base primes for the strikes (p ≤ x^(1/4) ≤ y, so the capped
    // oracle answers) and π at segment boundaries below √x — passed in by
    // piGourdonV (streamed); the null path (tests, whose oracle spans √x) derives
    // it locally and owns the allocation.
    const nwb: usize = @max(3, @as(usize, @intCast(s.pio.count(isqrt(sqx)))));
    const bpi_owned = bpi_in == null;
    const bpi: []const u64 = bpi_in orelse blk: {
        const nbpi = @as(usize, @intCast(sqx / segw)) + 2;
        const b2 = try gpa.alloc(u64, nbpi);
        b2[0] = 0;
        for (1..nbpi) |k| b2[k] = s.pio.count(@as(u64, k) * segw - 1);
        break :blk b2;
    };

    const tm = try buildSmallTemplates(gpa, primes);
    return .{
        .leaf = s.leaf,
        .primes = primes,
        .pio = &s.pio,
        .x = x,
        .y = y,
        .sqx = sqx,
        .sqrt_y = isqrt(y),
        .xstar = xstar,
        .total = total,
        .nax = nax,
        .naz = naz,
        .nay = @intCast(s.pio.count(y)),
        .segw = segw,
        .nseg = (total + segw - 1) / segw,
        .nb = nb,
        .naz_i = @intCast(naz),
        .blk_total = try gpa.alloc(i64, nb * nax),
        .blk_mu = try gpa.alloc(i64, nb * nax),
        .blk_total_full = try gpa.alloc(i64, nb),
        .blk_bcount = try gpa.alloc(u64, nb),
        .blk_omega = try gpa.alloc(i128, nb),
        .blk_b = try gpa.alloc(i128, nb),
        .bpi = bpi,
        .bpi_owned = bpi_owned,
        .nwb = nwb,
        .nbuck = nbuck,
        .nring = nring,
        .ring_mask = nring - 1,
        .nsmall = tm.n,
        .small_tmpl = tm.tmpl,
        .small_off = tm.off,
        .disp = std.atomic.Value(usize).init(0),
    };
}
fn freeBlkCtx(comptime X: type, comptime P: type, gpa: std.mem.Allocator, ctx: *BlkCtx(X, P)) void {
    gpa.free(ctx.blk_total);
    gpa.free(ctx.blk_mu);
    gpa.free(ctx.blk_total_full);
    gpa.free(ctx.blk_bcount);
    gpa.free(ctx.blk_omega);
    gpa.free(ctx.blk_b);
    gpa.free(ctx.small_tmpl);
    gpa.free(ctx.small_off);
    if (ctx.bpi_owned) gpa.free(@constCast(ctx.bpi));
}

fn omegaCounter(comptime INST: bool, comptime X: type, gpa: std.mem.Allocator, s: anytype, x: X, y: u64, z: u64, xstar: u64) !OmegaB {
    return omegaBlocked(INST, X, gpa, s, x, y, z, xstar, 1, null, SWEEP_SEGW);
}

/// Serial ω+B over nb blocks (nb=1 = monolithic sweep). Threaded: omegaBlockedPar.
fn omegaBlocked(comptime INST: bool, comptime X: type, gpa: std.mem.Allocator, s: anytype, x: X, y: u64, z: u64, xstar: u64, nb: usize, bpi_in: ?[]const u64, segw: usize) !OmegaB {
    const P = std.meta.Child(@TypeOf(s.primes));
    var ctx = (try initBlkCtx(X, P, gpa, s, x, y, z, xstar, nb, bpi_in, segw)) orelse return .{ .omega = 0, .b = 0 };
    defer freeBlkCtx(X, P, gpa, &ctx);
    prog_om.begin("ω+B", ctx.nseg);
    var sc = try Scratch.init(gpa, ctx.nax, ctx.naz, ctx.segw, ctx.nring, ctx.nwb);
    defer sc.deinit(gpa);
    var st = Stats{};
    for (0..nb) |t| try runOneBlock(INST, X, P, gpa, &ctx, &sc, &st, t);
    const r = try reduceOmB(X, P, gpa, &ctx);
    if (INST) std.debug.print("  ωB-stats: nax={d} naz={d} nb={d} segs={d} mwalk={d} leaves(small/C/D)={d}/{d}/{d} kills={d} Bqueries={d}\n", .{ ctx.nax, ctx.naz, nb, st.n_seg, st.n_mwalk, st.n_small, st.n_easy, st.n_hard, st.n_kill, st.n_bq });
    return r;
}

/// Threaded ω+B (Model-A phase 2): nthreads workers pull blocks from an atomic
/// dispenser, each with its own Scratch + Stats; block outputs are disjoint so no
/// locks. Serial reduceOmB stitches after join. nb = nthreads·k_over (over-partition).
fn omegaBlockedPar(comptime INST: bool, comptime X: type, gpa: std.mem.Allocator, s: anytype, x: X, y: u64, z: u64, xstar: u64, nb: usize, nthreads: usize, pins: ?[]const u32, bpi_in: ?[]const u64, segw: usize) !OmegaB {
    const P = std.meta.Child(@TypeOf(s.primes));
    var ctx = (try initBlkCtx(X, P, gpa, s, x, y, z, xstar, nb, bpi_in, segw)) orelse return .{ .omega = 0, .b = 0 };
    defer freeBlkCtx(X, P, gpa, &ctx);
    prog_om.begin("ω+B", ctx.nseg);

    const scratches = try gpa.alloc(Scratch, nthreads);
    defer gpa.free(scratches);
    var ninit: usize = 0;
    defer for (0..ninit) |k| scratches[k].deinit(gpa);
    for (0..nthreads) |i| {
        scratches[i] = try Scratch.init(gpa, ctx.nax, ctx.naz, ctx.segw, ctx.nring, ctx.nwb);
        ninit = i + 1;
    }
    const stats = try gpa.alloc(Stats, nthreads);
    defer gpa.free(stats);
    @memset(stats, Stats{});

    const Worker = struct {
        fn run(g: std.mem.Allocator, cx: *BlkCtx(X, P), sc: *Scratch, stt: *Stats, cpu: ?u32, errp: *?anyerror) void {
            if (cpu) |c| pinToCpu(c);
            while (true) {
                const t = cx.disp.fetchAdd(1, .monotonic);
                if (t >= cx.nb) break;
                runOneBlock(INST, X, P, g, cx, sc, stt, t) catch |e| {
                    errp.* = e; // bucket growth OOM: record and stop this worker
                    return;
                };
            }
        }
    };
    const threads = try gpa.alloc(std.Thread, nthreads);
    defer gpa.free(threads);
    const werr = try gpa.alloc(?anyerror, nthreads);
    defer gpa.free(werr);
    @memset(werr, null);
    var spawned: usize = 0;
    for (1..nthreads) |i| {
        const cpu: ?u32 = if (pins) |pp| pp[i] else null;
        threads[i] = std.Thread.spawn(.{}, Worker.run, .{ gpa, &ctx, &scratches[i], &stats[i], cpu, &werr[i] }) catch break;
        spawned = i;
    }
    Worker.run(gpa, &ctx, &scratches[0], &stats[0], if (pins) |pp| pp[0] else null, &werr[0]);
    var j: usize = 1;
    while (j <= spawned) : (j += 1) threads[j].join();
    for (werr) |e| if (e) |ee| return ee;

    if (INST) {
        var st = Stats{};
        for (stats) |ss| st.add(ss);
        std.debug.print("  ωB-par-stats: nb={d} nthreads={d} segs={d} kills={d} Bqueries={d}\n", .{ nb, nthreads, st.n_seg, st.n_kill, st.n_bq });
    }
    return reduceOmB(X, P, gpa, &ctx);
}

// ---------------------------------------------- ω reference (for the diff-check)
fn phiRec(primes: []const u32, pi: []const u32, u: u64, b: usize) i64 {
    if (u == 0) return 0;
    if (b == 0) return @intCast(u);
    const pb: u64 = primes[b - 1];
    if (pb > u) return 1;
    if (pb * pb > u) return @as(i64, @intCast(pi[@intCast(u)])) - @as(i64, @intCast(b)) + 1;
    return phiRec(primes, pi, u, b - 1) - phiRec(primes, pi, u / pb, b - 1);
}
fn omegaNaive(primes: []const u32, leaf: []const u16, pi: []const u32, x: u64, y: u64, xstar: u64) i128 {
    var omega: i128 = 0;
    for (primes, 0..) |p32, pidx| {
        const p: u64 = p32;
        if (p <= 2) continue;
        if (p > xstar) break;
        var m: u64 = y / p + 1;
        while (m <= y) : (m += 1) {
            const lv = leaf[@intCast(m)];
            if (@as(usize, lv & LEAF_SAT) <= pidx + 1) continue; // μ=0 or spf ≤ p
            const u = x / (p * m);
            const mun: i128 = if (lv & LEAF_PLUS != 0) 1 else -1;
            omega += mun * @as(i128, phiRec(primes, pi, u, pidx));
        }
    }
    return -omega;
}

/// Full π-table to n (used only by the small-x ω diff-check).
fn buildPi(gpa: std.mem.Allocator, n: u64) ![]u32 {
    const N: usize = @intCast(n + 1);
    const comp = try gpa.alloc(bool, N);
    defer gpa.free(comp);
    @memset(comp, false);
    var i: usize = 2;
    while (i * i < N) : (i += 1) {
        if (!comp[i]) {
            var j: usize = i * i;
            while (j < N) : (j += i) comp[j] = true;
        }
    }
    const pi = try gpa.alloc(u32, N);
    var c: u32 = 0;
    pi[0] = 0;
    if (N > 1) pi[1] = 0;
    i = 2;
    while (i < N) : (i += 1) {
        if (!comp[i]) c += 1;
        pi[i] = c;
    }
    return pi;
}

// --------------------------------------------------------------------- driver
pub const GResult = struct {
    pi: i128,
    A: i128,
    B: i128,
    omega: i128,
    phi0: i128,
    sigma: i128,
    y: u64,
};

/// y = α·x^(1/3). α sets the fold/leaf balance — fold work ∝ z = x/y (kills, and
/// the per-segment walk over π(√z) primes), leaf work ∝ y — so the optimum drifts
/// upward with x rather than sitting at a constant. Measured optima (parabolic
/// interpolation on a 6-thread α sweep, times in s):
///
///   x      α*     α=4     at α*    gain
///   10^15  4.03   0.72    0.72       -
///   10^16  4.55   2.89    2.87      ~1%
///   10^17  6.00  14.54   11.88     18%
///   10^18  8.16  76.21   50.42     34%
///
/// Refitted 2026-07-21 after the segmented oracle (471c328..5252767) made
/// leaf-side π cheaper and pushed every optimum up: measured α* = 4.41 / 5.41 /
/// 6.55 / 7.97 / 10.58 / 12.21 at 10^15..10^20 (residuals ±0.68), with the
/// killed 10^21 anchor's bracket (α* ≥ 13.4) matching this fit's 13.45. Expect
/// another shift when the audit worklist lands; recalibrate with
/// `pi --calibrate` then, not by hand.
///
/// Previous history — refitted after the bucket ring (2c6e540) removed the
/// fold's O(nseg·π(√z)) traversal term. That term had been inflating α* at the top of the range — α was
/// the only lever that shrank z — and removing it pulled the optima down by an
/// amount that grows with x, confirming the diagnosis:
///
///   x      α* before  α* after
///   10^16    4.55       4.58     (traversal negligible here — unmoved)
///   10^17    5.72       5.34
///   10^18    7.23       6.90
///   10^19   10.34       8.65
///
/// The SLOPE survived the refit unchanged (0.6017 → 0.5980); only the intercept
/// moved. A 3-point fit over 10^16..10^18 alone had suggested 0.5023, but that
/// extrapolated to 7.92 at 10^19 against a measured 8.65 — clustered points, not a
/// real flattening. Residuals tightened from ±0.9 to ±0.34. Still first-order.
///
/// The floor matters. α* is flat at ≈4 up to 10^16 and only then climbs, which a
/// line cannot express; unclamped it reads well under 1 at 10^13, violating the
/// z < y² Legendre condition (α > 1) that B's read-off depends on. Clamped, the
/// failure mode off the bottom of the fitted range is "degrade to the old constant".
/// Untested above 10^18 — the ceiling is a guard rail, not a tuned value.
const ALPHA_A: f64 = -20.1100;
const ALPHA_B: f64 = 0.6940;
const ALPHA_LO: f64 = 4.0;
const ALPHA_HI: f64 = 24.0;

/// The built-in default α(x) at u128 width, for callers (the CLI's calibrator)
/// that want the fit's prediction without duplicating the constants.
pub fn defaultAlpha(x: u128) f64 {
    return chooseAlpha(u128, x);
}

fn chooseAlpha(comptime X: type, x: X) f64 {
    if (x < 1000) return ALPHA_LO;
    const a = ALPHA_A + ALPHA_B * @log(@as(f64, @floatFromInt(x)));
    return @min(ALPHA_HI, @max(ALPHA_LO, a));
}

fn chooseY(comptime X: type, x: X) u64 {
    const cr = icbrtG(X, x);
    const lo = cr + 1;
    const hi = isqrtG(X, x) -| 1;
    const yf = chooseAlpha(X, x) * @as(f64, @floatFromInt(cr));
    var y: u64 = if (yf >= 1.8e19) std.math.maxInt(u64) else @intFromFloat(yf);
    if (y < lo) y = lo;
    if (y > hi) y = hi;
    return y;
}

/// Dispatch: u64 path for x ≤ 2⁶⁴ (u32 primes), u128 path beyond (u64 primes).
/// Everything tunable about a run, in one place, so callers (drivers, the CLI, the
/// test suite) name what they mean instead of threading positional parameters.
/// Defaults reproduce piGourdon(x): serial, fitted α, no tracing.
pub const Config = struct {
    /// y = α·x^(1/3). null ⇒ chooseY's fitted α(x). Setting it overrides the fit.
    y: ?u64 = null,
    /// Worker threads. 1 ⇒ the serial path (not a 1-thread parallel path).
    nthreads: usize = 1,
    /// CPU to pin each worker to. null ⇒ unpinned. Length should be ≥ nthreads.
    pins: ?[]const u32 = null,
    /// Per-phase timing to stderr.
    verbose: bool = false,
    /// Sweep segment width in integers (counter bits = segw/8 bytes). null ⇒
    /// SWEEP_SEGW (32 KB bits, this laptop's L1d). Must be a multiple of 960
    /// (MASK30 phase) and ≤ 2^21 (ring-entry offset field). The cache-hierarchy
    /// knob: tune per machine, e.g. halved when running 2 threads/core.
    segw: ?usize = null,
};

/// Below this the decomposition is not well defined: chooseY's clamps collapse (y is
/// pinned to ≤ √x − 1, which is 0 at x = 1, and z = x/y then divides by zero), and
/// x* / √z stop separating. Answer directly from the oracle instead — exact, already
/// differentially tested, and trivially cheap at this size.
const DIRECT_MAX: u128 = 10_000;

fn directResult(pi: i128) GResult {
    return .{ .pi = pi, .A = 0, .B = 0, .omega = 0, .phi0 = 0, .sigma = 0, .y = 0 };
}

/// Dispatch on x: direct for tiny x, then the u64 path below 2⁶⁴ (u32 primes) and
/// the u128 path above (u64 primes).
pub fn piGourdonCfg(gpa: std.mem.Allocator, x: u128, cfg: Config) !GResult {
    if (x < 2) return directResult(0);
    if (x <= DIRECT_MAX) {
        var o = try buildPiOracle(gpa, @intCast(x));
        defer gpa.free(o.bits);
        defer gpa.free(o.pref);
        return directResult(@intCast(o.count(@intCast(x))));
    }
    if (x <= std.math.maxInt(u64)) return piGourdonV(u64, gpa, @intCast(x), cfg);
    return piGourdonV(u128, gpa, x, cfg);
}

pub fn piGourdon(gpa: std.mem.Allocator, x: u128, y_in: ?u64) !GResult {
    return piGourdonCfg(gpa, x, .{ .y = y_in });
}

/// Parallel entry: nthreads workers (pinned per `pins` if given).
pub fn piGourdonPar(gpa: std.mem.Allocator, x: u128, nthreads: usize, pins: ?[]const u32) !GResult {
    return piGourdonCfg(gpa, x, .{ .nthreads = nthreads, .pins = pins });
}

pub fn piGourdonV(comptime X: type, gpa: std.mem.Allocator, x: X, cfg: Config) !GResult {
    const y_in = cfg.y;
    const verbose = cfg.verbose;
    g_progress = verbose;
    const nthreads = cfg.nthreads;
    const pins = cfg.pins;
    const P = if (X == u64) u32 else u64; // prime element type: √x < 2³² ⇔ X = u64
    var tp = common.nowNs();
    const y = y_in orelse chooseY(X, x);
    const z: u64 = xdiv(X, x, y);
    const sqx = isqrtG(X, x);
    const x13 = icbrtG(X, x);
    const sqz = isqrt(z);
    const xstar = @max(isqrt(sqx), xdiv(X, x, y * y));

    // Explicit prime list only as far as some loop enumerates primes BY VALUE:
    // y (ω's sparse q-walk), √z (the fold), x^(1/3) (A/Σ's outer p), and A's inner
    // q ≤ √(x/p) < √(x/x*) ≤ y (since x* ≥ x/y²) — so y already covers A's q, and
    // the √(x/x*) term is belt-and-braces for the x* = x^(1/4) regime. Every π(v)
    // query above plist_max goes to the bitset oracle, which spans all of [0, √x].
    const plist_max = @max(@max(y, sqz), @max(x13, isqrtG(X, x / @as(X, @max(xstar, 1)))));
    var s = try Sieve(P).init(gpa, sqx, y, plist_max, nthreads, pins);
    defer s.deinit(gpa);
    const lap = struct {
        fn f(v: bool, tpp: *u64, name: []const u8) void {
            if (!v) return;
            const now = common.nowNs();
            std.debug.print("  [{s:>6}] {d:>8.3} s\n", .{ name, @as(f64, @floatFromInt(now - tpp.*)) / 1e9 });
            tpp.* = now;
        }
    }.f;
    lap(verbose, &tp, "sieve");

    const segw: usize = cfg.segw orelse SWEEP_SEGW;
    if (segw < 960 or segw % 960 != 0 or segw > (1 << 21)) return error.BadSegw;

    // π at segment boundaries below √x + π(√x), streamed — O(√x/segw) retained.
    const bp = try buildBoundaryPi(gpa, sqx, segw, nthreads, pins);
    const bpi = bp.bpi;
    defer gpa.free(bpi);
    const pi_sqx: i128 = @intCast(bp.total);
    lap(verbose, &tp, "bpi");

    // A/Σ phase 1 (nthreads = 1 runs the same code inline; computeTerms survives
    // only as the suite's oracle-backed differential reference).
    const t = try computeTermsPar(X, gpa, &s, x, y, sqx, x13, sqz, xstar, nthreads, pins, bpi, pi_sqx, segw);
    lap(verbose, &tp, "A/Σ");

    // ω+B fused. Phase 2: block-and-scan (nb = nthreads·8 over-partition) if parallel.
    const wb = if (nthreads > 1)
        try omegaBlockedPar(false, X, gpa, &s, x, y, z, xstar, nthreads * 8, nthreads, pins, bpi, segw)
    else
        try omegaBlocked(false, X, gpa, &s, x, y, z, xstar, 1, bpi, segw);
    const omega = wb.omega;
    const B = wb.b;
    lap(verbose, &tp, "ω+B");

    // φ₀ = Σ_{n≤y, n odd, μ(n)≠0} μ(n)·φ(x/n,1), φ(u,1)=u−⌊u/2⌋   [k=1]
    var phi0: i128 = 0;
    {
        var n: u64 = 1;
        while (n <= y) : (n += 1) {
            const lv = s.leaf[@intCast(n)];
            if (lv == 0) continue;
            if (n % 2 == 0) continue;
            const u: X = x / @as(X, n); // NOT xdiv: for x>2^64 and n∈{1,3,5}, x/n exceeds u64
            const mun: i128 = if (lv & LEAF_PLUS != 0) 1 else -1;
            phi0 += mun * @as(i128, @intCast(u - u / 2));
        }
    }
    lap(verbose, &tp, "phi0");

    // Σ closed forms
    const a = t.a;
    const bb = t.b;
    const cc = t.c;
    const dd = t.d;
    const Psx = t.P; // π(√x)
    const sig0: i128 = a - 1 + @divExact(Psx * (Psx - 1), 2) - @divExact(a * (a - 1), 2);
    const sig1: i128 = @divExact((a - bb) * (a - bb - 1), 2);
    const sig2: i128 = a * (bb - cc - @divExact(cc * (cc - 3), 2) + @divExact(dd * (dd - 3), 2));
    const sig3: i128 = @divExact(bb * (bb - 1) * (2 * bb - 1), 6) - bb - @divExact(dd * (dd - 1) * (2 * dd - 1), 6) + dd;
    const sigma = sig0 + sig1 + sig2 + sig3 + t.sig4 + t.sig5 + t.sig6;

    const pi = t.A - B + omega + phi0 + sigma;
    return .{ .pi = pi, .A = t.A, .B = B, .omega = omega, .phi0 = phi0, .sigma = sigma, .y = y };
}

pub fn main() !void {
    const gpa = std.heap.page_allocator;

    // π oracle differential check: count() and prevPrime() against an explicit
    // prime list, at EVERY v ≤ N (the oracle answers every π query in the terms,
    // so an off-by-one anywhere in the wheel indexing would silently skew π(x)).
    {
        const N: u64 = 3_000_000;
        const o = try buildPiOracle(gpa, N);
        defer gpa.free(o.bits);
        defer gpa.free(o.pref);
        const ref = try sievePrimes(u32, gpa, N);
        defer gpa.free(ref);
        var i: usize = 0;
        var bad: usize = 0;
        for (0..N + 1) |v| {
            while (i < ref.len and ref[i] <= v) i += 1;
            const want_prev: u64 = if (i == 0) 0 else ref[i - 1];
            if (o.count(v) != i or o.prevPrime(v) != want_prev) bad += 1;
        }
        std.debug.print("π oracle check (v ≤ {d}): π={d} {s} ({d} mismatches)\n\n", .{
            N, o.count(N), if (bad == 0 and o.count(N) == ref.len) "match" else "MISMATCH", bad,
        });
    }

    // Exhaustive small-x check. gourdon was only ever tested from 10^5 up, which let
    // a SIGFPE at x < 4 survive: chooseY clamps y to <= sqrt(x)-1, that is 0 at x = 1,
    // and z = x/y then divides by zero. Covers the direct path, the DIRECT_MAX
    // boundary, and the bottom of the real decomposition.
    {
        const N: u64 = 30_000;
        const ref = try buildPiOracle(gpa, N);
        defer gpa.free(ref.bits);
        defer gpa.free(ref.pref);
        var bad: usize = 0;
        var first_bad: u64 = 0;
        for (0..N + 1) |xv| {
            const got = (try piGourdonCfg(gpa, @as(u128, xv), .{})).pi;
            if (got != @as(i128, @intCast(ref.count(xv)))) {
                if (bad == 0) first_bad = @intCast(xv);
                bad += 1;
            }
        }
        std.debug.print("exhaustive x in [0,{d}]: {s}{s}\n\n", .{
            N,
            if (bad == 0) "match" else "MISMATCH",
            if (bad == 0) "" else " (see first_bad)",
        });
        if (bad != 0) std.debug.print("  {d} mismatches, first at x={d}\n", .{ bad, first_bad });
    }

    // π window differential: pwinBuild/count vs the resident oracle at EVERY v over
    // segw-aligned windows — the same machinery the ω sweep uses for C-leaves.
    {
        const SEGW: u64 = 273 * 960;
        const N: u64 = SEGW * 20; // 5.24M
        const ref = try buildPiOracle(gpa, N);
        defer gpa.free(ref.bits);
        defer gpa.free(ref.pref);
        const base = try sievePrimes(u32, gpa, isqrt(N));
        defer gpa.free(base);
        const nwin: usize = @intCast(SEGW / 240);
        const wb = try gpa.alloc(u64, nwin);
        defer gpa.free(wb);
        const wpref = try gpa.alloc(u32, nwin + 1);
        defer gpa.free(wpref);
        const pn = try gpa.alloc(u64, base.len);
        defer gpa.free(pn);
        const pwv = try gpa.alloc(u8, base.len);
        defer gpa.free(pwv);
        var win = PiWin{ .bits = wb, .pref = wpref, .lo = 0, .base = 0 };
        pwinInitCursors(u32, base, base.len, pn, pwv, 0);
        var bad: usize = 0;
        var lo: u64 = 0;
        while (lo < N) : (lo += SEGW) {
            const hi = @min(lo + SEGW, N);
            win.base = if (lo == 0) 0 else ref.count(lo - 1);
            pwinBuild(u32, &win, base, base.len, pn, pwv, lo, hi);
            var vv: u64 = lo;
            while (vv < hi) : (vv += 1) {
                if (win.count(vv) != ref.count(vv)) bad += 1;
            }
        }
        std.debug.print("π window check (v < {d}, {d} windows): {s} ({d} mismatches)\n", .{ N, N / SEGW, if (bad == 0) "match" else "MISMATCH", bad });

        // B-window differential: bwinPrev's full descending walk (refilling chunk
        // by chunk) must reproduce the oracle's prevPrime at every step.
        var bw = BWin{ .bits = wb };
        var bbad: usize = 0;
        var nsteps: usize = 0;
        var v2: u64 = N;
        while (v2 >= 2) {
            const got = bwinPrev(u32, &bw, base, base.len, v2, 0);
            const want = ref.prevPrime(v2);
            if (got != want) bbad += 1;
            nsteps += 1;
            if (want < 2) break;
            v2 = want - 1;
        }
        std.debug.print("B window check ({d} descending steps): {s} ({d} mismatches)\n\n", .{ nsteps, if (bbad == 0) "match" else "MISMATCH", bbad });
    }

    // ω differential check: naive recurrence vs O(1)-kill counter, must match exactly.
    std.debug.print("ω check (naive vs counter):\n", .{});
    for ([_]u64{ 100_000, 1_000_000, 10_000_000, 100_000_000, 1_000_000_000 }) |x| {
        const y = chooseY(u64, x);
        const z = x / y;
        const xstar = @max(isqrt(isqrt(x)), x / (y * y));
        var s = try Sieve(u32).init(gpa, isqrt(x), y, isqrt(x), 1, null);
        defer s.deinit(gpa);
        const pi = try buildPi(gpa, z);
        defer gpa.free(pi);
        const on = omegaNaive(s.primes, s.leaf, pi, x, y, xstar);
        const wb = try omegaCounter(false, u64, gpa, &s, x, y, z, xstar);
        const b_ref = try computeB(false, gpa, &s, x, y, z, isqrt(x)); // standalone B reference
        // block-consistency: nb blocks + reduction must equal the nb=1 sweep
        const wb4 = try omegaBlocked(false, u64, gpa, &s, x, y, z, xstar, 4, null, SWEEP_SEGW);
        const wb7 = try omegaBlocked(false, u64, gpa, &s, x, y, z, xstar, 7, null, SWEEP_SEGW);
        const blk_ok = wb4.omega == wb.omega and wb4.b == wb.b and wb7.omega == wb.omega and wb7.b == wb.b;
        std.debug.print("  {d:>12}  ω naive={d:>14} counter={d:>14} {s}   B fused={d} ref={d} {s}   blocks(4,7){s}\n", .{ x, on, wb.omega, if (on == wb.omega) "match" else "MISMATCH", wb.b, b_ref, if (wb.b == b_ref) "match" else "MISMATCH", if (blk_ok) "=✓" else "=MISMATCH" });
    }

    // A/Σ parallel check: 4-thread partial-sum must equal serial (phase 1 of Model-A).
    std.debug.print("\nA/Σ parallel check (serial vs 4-thread, cores 0/2/4/6):\n", .{});
    {
        const pins = [_]u32{ 0, 2, 4, 6 };
        for ([_]u64{ 1_000_000_000, 100_000_000_000, 1_000_000_000_000 }) |x| {
            const y = chooseY(u64, x);
            const z = x / y;
            const sqx = isqrt(x);
            const x13 = icbrtG(u64, x);
            const sqz = isqrt(z);
            const xstar = @max(isqrt(isqrt(x)), x / (y * y));
            var s = try Sieve(u32).init(gpa, sqx, y, sqx, 1, null);
            defer s.deinit(gpa);
            const t0 = common.nowNs();
            const ser = computeTerms(u64, &s, x, y, sqx, x13, sqz, xstar);
            const st = @as(f64, @floatFromInt(common.nowNs() - t0)) / 1e9;
            const nbpi = @as(usize, @intCast(sqx / SWEEP_SEGW)) + 2;
            const bpi = try gpa.alloc(u64, nbpi);
            defer gpa.free(bpi);
            bpi[0] = 0;
            for (1..nbpi) |k| bpi[k] = s.pio.count(@as(u64, k) * SWEEP_SEGW - 1);
            const t1 = common.nowNs();
            const par = try computeTermsPar(u64, gpa, &s, x, y, sqx, x13, sqz, xstar, 4, &pins, bpi, @intCast(s.pio.count(sqx)), SWEEP_SEGW);
            const pt = @as(f64, @floatFromInt(common.nowNs() - t1)) / 1e9;
            const ok = ser.A == par.A and ser.sig4 == par.sig4 and ser.sig5 == par.sig5 and ser.sig6 == par.sig6;
            std.debug.print("  x=10^{d}: {s}  serial {d:.4}s  4-thread {d:.4}s  ({d:.2}x)\n", .{ std.math.log10_int(x), if (ok) "match" else "MISMATCH", st, pt, st / pt });
        }
    }

    // Full parallel total: piGourdonPar(4) must equal piGourdon (serial), + speedup.
    std.debug.print("\nparallel total (serial vs 4-thread, cores 0/2/4/6):\n", .{});
    {
        const pins = [_]u32{ 0, 2, 4, 6 };
        for ([_]u64{ 1_000_000_000_000, 100_000_000_000_000 }) |x| {
            const t0 = common.nowNs();
            const ser = try piGourdon(gpa, x, null);
            const st = @as(f64, @floatFromInt(common.nowNs() - t0)) / 1e9;
            const t1 = common.nowNs();
            const par = try piGourdonPar(gpa, x, 4, &pins);
            const pt = @as(f64, @floatFromInt(common.nowNs() - t1)) / 1e9;
            std.debug.print("  x=10^{d}: pi={d} {s}  serial {d:.3}s  4-thread {d:.3}s  ({d:.2}x)\n", .{ std.math.log10_int(x), par.pi, if (ser.pi == par.pi) "match" else "MISMATCH", st, pt, st / pt });
        }
    }

    // Correctness + timing: Gourdon (counter-ω, segmented-π) vs lmo, both serial.
    std.debug.print("\ntotal (gourdon vs lmo):\n", .{});
    const xs = [_]u64{ 1_000_000_000, 10_000_000_000, 100_000_000_000, 1_000_000_000_000 };
    std.debug.print("{s:>14} {s:>16} {s:>4} {s:>11} {s:>11} {s:>8}\n", .{ "x", "pi(x)", "ok", "gourdon_s", "lmo_s", "g/lmo" });
    for (xs) |x| {
        const t0 = common.nowNs();
        const g = try piGourdon(gpa, x, null);
        const gs = @as(f64, @floatFromInt(common.nowNs() - t0)) / 1e9;

        const t1 = common.nowNs();
        const ref = try lmo.piLMO(gpa, x, null, null);
        const ls = @as(f64, @floatFromInt(common.nowNs() - t1)) / 1e9;

        const ok = g.pi == @as(i128, @intCast(ref.pi));
        std.debug.print("{d:>14} {d:>16} {s:>4} {d:>11.3} {d:>11.3} {d:>8.2}\n", .{ x, g.pi, if (ok) "y" else "NO", gs, ls, gs / ls });
    }
    std.debug.print("\nper-term profile @ 10^12:\n", .{});
    _ = try piGourdonV(u64, gpa, 1_000_000_000_000, .{ .verbose = true });

    std.debug.print("\nop-count comparison (gourdon ω/B vs lmo S2/P₂), matched x,y:\n", .{});
    for ([_]u64{ 1_000_000_000_000, 100_000_000_000_000 }) |x| {
        const y = chooseY(u64, x);
        const z = x / y;
        const xstar = @max(isqrt(isqrt(x)), x / (y * y));
        var s = try Sieve(u32).init(gpa, isqrt(x), y, isqrt(x), 1, null);
        defer s.deinit(gpa);
        std.debug.print("x=10^{d} (y={d}, z={d}):\n", .{ std.math.log10_int(x), y, z });
        std.debug.print("  gourdon(fused) ", .{});
        _ = try omegaCounter(true, u64, gpa, &s, x, y, z, xstar); // prints ωB-stats (√z fold + B queries)
        const lr = try lmo.s2AndP2FusedInstrumented(gpa, @intCast(x), y, y);
        std.debug.print("  lmo:           S2 kills={d}  S2 prefix(hard)={d}  P₂ prefix(np)={d}\n", .{ lr.kills, lr.s2q, lr.np });
    }

    // u128-path validation: force X=u128 on x<2⁶⁴; must equal the u64 path bit-for-bit.
    std.debug.print("\nu128-path check (X=u128 on x<2^64 must match u64):\n", .{});
    for ([_]u64{ 1_000_000_000_000, 100_000_000_000_000, 10_000_000_000_000_000 }) |x| {
        const g64 = try piGourdonV(u64, gpa, x, .{});
        const g128 = try piGourdonV(u128, gpa, @as(u128, x), .{});
        std.debug.print("  x=10^{d}: u64={d} u128={d} {s}\n", .{ std.math.log10_int(x), g64.pi, g128.pi, if (g64.pi == g128.pi) "match" else "MISMATCH" });
    }

    // Ceiling test: gourdon-only vs known values (u128 x reaches past 2⁶⁴).
    std.debug.print("\nceiling test (gourdon only, vs known):\n", .{});
    const cases = [_]struct { x: u128, want: i128 }{
        .{ .x = 10_000_000_000_000, .want = 346_065_536_839 },
        .{ .x = 100_000_000_000_000, .want = 3_204_941_750_802 },
        .{ .x = 1_000_000_000_000_000, .want = 29_844_570_422_669 },
        .{ .x = 10_000_000_000_000_000, .want = 279_238_341_033_925 },
    };
    for (cases) |cc| {
        const t0 = common.nowNs();
        const g = try piGourdon(gpa, cc.x, null);
        const gs = @as(f64, @floatFromInt(common.nowNs() - t0)) / 1e9;
        const ru = std.posix.getrusage(std.posix.rusage.SELF);
        const rss_mb = @divTrunc(ru.maxrss, 1024);
        std.debug.print("  {d:>20} pi={d:>19} {s:>4}  {d:>7.2} s  peakRSS={d} MB\n", .{ cc.x, g.pi, if (g.pi == cc.want) "y" else "NO", gs, rss_mb });
    }
}
