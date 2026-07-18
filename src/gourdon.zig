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
//!         pi_tab/piLE) and DR/LMO dense/sparse leaf enumeration (p>√y ⇒ m prime).
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
fn icbrt(x: u64) u64 {
    if (x == 0) return 0;
    var r: u64 = 1;
    while ((r + 1) * (r + 1) * (r + 1) <= x) r += 1;
    while (r * r * r > x) r -= 1;
    return r;
}
fn isqrtsq(x: u64) u64 { // x^(1/4)
    return isqrt(isqrt(x));
}

// --------------------------------------------------------------- small sieve
/// Primes ≤ nmax via a bit-packed Eratosthenes (nmax/8 transient bytes, not nmax·4).
fn sievePrimes(gpa: std.mem.Allocator, nmax: u64) ![]u32 {
    const N: usize = @intCast(nmax + 1);
    const nw = (N + 63) / 64;
    const bits = try gpa.alloc(u64, nw); // bit set ⇒ composite
    defer gpa.free(bits);
    @memset(bits, 0);
    const isComp = struct {
        fn f(b: []const u64, i: usize) bool {
            return b[i >> 6] & (@as(u64, 1) << @as(u6, @intCast(i & 63))) != 0;
        }
    }.f;
    var i: usize = 2;
    while (i * i < N) : (i += 1) {
        if (!isComp(bits, i)) {
            var j: usize = i * i;
            while (j < N) : (j += i) bits[j >> 6] |= @as(u64, 1) << @as(u6, @intCast(j & 63));
        }
    }
    var np: usize = 0;
    i = 2;
    while (i < N) : (i += 1) {
        if (!isComp(bits, i)) np += 1;
    }
    const primes = try gpa.alloc(u32, np);
    var k: usize = 0;
    i = 2;
    while (i < N) : (i += 1) {
        if (!isComp(bits, i)) {
            primes[k] = @intCast(i);
            k += 1;
        }
    }
    return primes;
}

/// Persistent tables: primes ≤ √x (for piLE and the term loops), and δ/μ only to y
/// (ω leaves, φ₀ — never read past y). No O(√x) factor tables.
const Sieve = struct {
    mu: []i8, // [0..y]
    delta: []u32, // smallest prime factor, [0..y]
    pi_tab: []u32, // π(v) for v ≤ y (O(1), for ω's easy-leaf shortcut)
    primes: []u32, // ascending, ≤ √x

    fn init(gpa: std.mem.Allocator, sqx: u64, y: u64) !Sieve {
        const primes = try sievePrimes(gpa, sqx);
        const Ny: usize = @intCast(y + 1);
        const delta = try gpa.alloc(u32, Ny);
        @memset(delta, 0);
        var i: usize = 2;
        while (i < Ny) : (i += 1) {
            if (delta[i] == 0) {
                var j: usize = i;
                while (j < Ny) : (j += i) {
                    if (delta[j] == 0) delta[j] = @intCast(i);
                }
            }
        }
        const mu = try gpa.alloc(i8, Ny);
        mu[0] = 0;
        if (Ny > 1) mu[1] = 1;
        i = 2;
        while (i < Ny) : (i += 1) {
            const p = delta[i];
            const m = i / p;
            mu[i] = if (m % p == 0) 0 else -mu[@intCast(m)];
        }
        const pi_tab = try gpa.alloc(u32, Ny);
        var c: u32 = 0;
        i = 0;
        while (i < Ny) : (i += 1) {
            if (i >= 2 and delta[i] == i) c += 1; // δ[i]==i ⇔ prime
            pi_tab[i] = c;
        }
        return .{ .mu = mu, .delta = delta, .pi_tab = pi_tab, .primes = primes };
    }
    fn deinit(self: *Sieve, gpa: std.mem.Allocator) void {
        gpa.free(self.mu);
        gpa.free(self.delta);
        gpa.free(self.pi_tab);
        gpa.free(self.primes);
    }
};

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

/// π(v) for v ≤ √x by binary search on the prime list (which holds every prime ≤ √x).
/// Correct ONLY for v ≤ √x — the terms that use it are all bounded there; B is not,
/// and takes the sweep instead.
fn piLE(primes: []const u32, v: u64) u64 {
    if (v < 2) return 0;
    var lo: usize = 0;
    var hi: usize = primes.len;
    while (lo < hi) {
        const mid = (lo + hi) / 2;
        if (@as(u64, primes[mid]) <= v) lo = mid + 1 else hi = mid;
    }
    return lo; // # primes ≤ v
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

/// A, Σ₄/₅/₆ and the scalars use piLE (all query points ≤ √x — proven). Only B has
/// v = x/p up to z > √x, so B alone takes a segmented sweep over its ~π(√x) points.
/// π-free terms (ω, φ₀, Σ₀–₃) are handled in piGourdonV.
fn computeTerms(s: *const Sieve, x: u64, y: u64, sqx: u64, x13: u64, sqz: u64, xstar: u64) Terms {
    const pr = s.primes;
    const a: i128 = @intCast(piLE(pr, y));
    const b: i128 = @intCast(piLE(pr, x13));
    const c: i128 = @intCast(piLE(pr, sqz));
    const d: i128 = @intCast(piLE(pr, xstar));
    const P: i128 = @intCast(piLE(pr, sqx));

    // A = Σ_{x*<p≤x^1/3} Σ_{p<q≤√(x/p)} χ(x/pq) π(x/pq),  χ = 2 if x/pq<y else 1.
    // For fixed p, q ascends ⇒ v=x/pq descends ⇒ π(v) only decreases: a monotone
    // cursor pc = π(v) replaces the per-pair binary search (O(1) amortized, no table).
    var A: i128 = 0;
    for (pr, 0..) |p32, pidx| {
        const p: u64 = p32;
        if (p <= xstar) continue;
        if (p > x13) break;
        const qhi = isqrt(x / p);
        var qi = pidx + 1;
        var pc: usize = 0; // π(v) cursor, seeded on the first q
        var seeded = false;
        while (qi < pr.len) : (qi += 1) {
            const qq: u64 = pr[qi];
            if (qq > qhi) break;
            const v = x / (p * qq); // < √x
            if (!seeded) {
                pc = piLE(pr, v);
                seeded = true;
            } else {
                while (pc > 0 and pr[pc - 1] > v) pc -= 1;
            }
            const chi: i128 = if (v < y) 2 else 1;
            A += chi * @as(i128, @intCast(pc));
        }
    }
    // Σ₄ = a·Σ_{x*<p≤√(x/y)} π(x/py)   (v ≤ x^{5/12} < √x)
    var sig4: i128 = 0;
    for (pr) |p32| {
        const p: u64 = p32;
        if (p <= xstar) continue;
        if (p > sqz) break;
        sig4 += @intCast(piLE(pr, x / (p * y)));
    }
    sig4 *= a;
    // Σ₅ = Σ_{√(x/y)<p≤x^1/3} π(x/p²)   (v < y ≤ √x)
    var sig5: i128 = 0;
    for (pr) |p32| {
        const p: u64 = p32;
        if (p <= sqz) continue;
        if (p > x13) break;
        sig5 += @intCast(piLE(pr, x / (p * p)));
    }
    // Σ₆ = −Σ_{x*<p≤x^1/3} π(√(x/p))²   (v ≤ x^{3/8} < √x)
    var sig6: i128 = 0;
    for (pr) |p32| {
        const p: u64 = p32;
        if (p <= xstar) continue;
        if (p > x13) break;
        const t: i128 = @intCast(piLE(pr, isqrt(x / p)));
        sig6 += t * t;
    }
    return .{ .A = A, .sig4 = sig4, .sig5 = sig5, .sig6 = -sig6, .a = a, .b = b, .c = c, .d = d, .P = P };
}

/// B = Σ_{y<p≤√x} π(x/p). The only term needing π(v) for v > √x (v ∈ [√x, z]).
/// Answered by a running-π segmented sweep over just [√x, z], starting at π(√x).
fn computeB(comptime INST: bool, gpa: std.mem.Allocator, s: *const Sieve, x: u64, y: u64, z: u64, sqx: u64) !i128 {
    const pr = s.primes;
    var nb: usize = 0;
    for (pr) |p32| {
        const p: u64 = p32;
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
        const p: u64 = p32;
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
fn omegaCounter(comptime INST: bool, gpa: std.mem.Allocator, s: *const Sieve, x: u64, y: u64, z: u64, xstar: u64) !struct { omega: i128, b: i128 } {
    const primes = s.primes;
    const sqx = isqrt(x);
    const sqrt_y = isqrt(y); // dense/sparse split point (DR/LMO): p≤√y dense, else sparse
    const sqz = isqrt(z); // √z: the fold bound (Legendre a′ = π(√z))
    const nay: usize = piLE(primes, y); // π(y): index nay−1 is the largest prime ≤ y
    const nax: usize = piLE(primes, xstar); // π(x*): ω leaf/stage range
    const naz: usize = piLE(primes, sqz); // π(√z): fold range
    const nsx: usize = piLE(primes, sqx); // π(√x): upper B prime
    if (naz == 0) return .{ .omega = 0, .b = 0 };
    const naz_i: i64 = @intCast(naz);
    const segw: usize = 273 * 960; // ≈256 KiB counter, MASK30-aligned (lcm(30,64)=960)
    var n_seg: u64 = 0;
    var n_mwalk: u64 = 0;
    var n_small: u64 = 0;
    var n_easy: u64 = 0;
    var n_hard: u64 = 0; // ω counter prefix() calls
    var n_kill: u64 = 0; // fold kills (to √z)
    var n_bq: u64 = 0; // B counter prefix() calls

    var ctr = try Counter3P.init(gpa, segw);
    defer ctr.deinit(gpa);
    const phi_run = try gpa.alloc(i64, nax); // ω stages only
    defer gpa.free(phi_run);
    const seg_cnt = try gpa.alloc(i64, nax);
    defer gpa.free(seg_cnt);
    const cur = try gpa.alloc(u64, nax);
    defer gpa.free(cur);
    const next = try gpa.alloc(u64, naz); // fold cursors for ALL primes ≤ √z
    defer gpa.free(next);
    const wpos = try gpa.alloc(u8, naz);
    defer gpa.free(wpos);

    @memset(phi_run, 0);
    for (0..nax) |bi| {
        cur[bi] = if (primes[bi] <= sqrt_y) y else @as(u64, nay - 1);
    }
    for (3..naz) |bi| {
        next[bi] = primes[bi]; // p·1 (p≥7 ⇒ coprime to 30), wpos=W30IDX[1]=0
        wpos[bi] = 0;
    }

    const evalPhi = struct {
        inline fn f(comptime IN: bool, vv: u64, b: usize, pp: u64, ss: *const Sieve, prs: []const u32, pr_bi: i64, ct: *const Counter3P, l: u64, sx: u64, yy: u64, cE: *u64, cH: *u64) i64 {
            var piv: i64 = -1;
            if (pp * pp > vv) {
                if (vv <= yy) piv = @intCast(ss.pi_tab[@intCast(vv)]) else if (vv <= sx) piv = @intCast(piLE(prs, vv));
            }
            if (piv >= 0) {
                if (IN) cE.* += 1;
                return 1 + @max(0, piv - @as(i64, @intCast(b)));
            }
            if (IN) cH.* += 1;
            return pr_bi + ct.prefix(@intCast(vv - l));
        }
    }.f;

    var omega: i128 = 0;
    var b_sum: i128 = 0;
    var phi_run_full: i64 = 0; // running φ(·,π(√z)) over prior segments (for B)
    var pB: usize = nsx; // B cursor: primes[pB−1] is the current largest B-prime (>y)
    var lo: u64 = 0;
    while (lo <= z) : (lo += segw) {
        const hi = @min(lo + @as(u64, segw), z + 1);
        const len: usize = @intCast(hi - lo);
        ctr.reset(len);
        if (INST) n_seg += 1;
        for (0..naz) |bi| {
            const p: u64 = primes[bi];
            if (bi < nax) {
                if (bi >= 3) seg_cnt[bi] = ctr.total; // φ(·,bi) for this segment
                if (p > 2 and p <= sqrt_y) {
                    // DENSE m-walk
                    var m: u64 = cur[bi];
                    const mlo: u64 = y / p;
                    while (m > mlo) {
                        if (INST) n_mwalk += 1;
                        const v: u64 = x / (m * p);
                        if (v >= hi) break;
                        if (v >= lo) {
                            const mm = s.mu[@intCast(m)];
                            if (mm != 0 and s.delta[@intCast(m)] > p) {
                                const phi_v: i64 = if (bi <= 2) blk: {
                                    if (INST) n_small += 1;
                                    break :blk phiSmall(v, bi);
                                } else evalPhi(INST, v, bi, p, s, primes, phi_run[bi], &ctr, lo, sqx, y, &n_easy, &n_hard);
                                omega += @as(i128, -mm) * @as(i128, phi_v);
                            }
                        }
                        m -= 1;
                    }
                    cur[bi] = m;
                } else if (p > 2) {
                    // SPARSE q-walk (p>√y): every prime q∈(p,y] is a valid leaf
                    var qc: usize = @intCast(cur[bi]);
                    while (qc > bi) {
                        if (INST) n_mwalk += 1;
                        const q: u64 = primes[qc];
                        const v: u64 = x / (p * q);
                        if (v >= hi) break;
                        if (v >= lo) omega += @as(i128, evalPhi(INST, v, bi, p, s, primes, phi_run[bi], &ctr, lo, sqx, y, &n_easy, &n_hard));
                        qc -= 1;
                    }
                    cur[bi] = qc;
                }
            }
            if (bi >= 3) { // fold p (all primes ≤ √z): kill coprime-30 multiples in [lo,hi)
                var j: u64 = next[bi];
                var wp: u8 = wpos[bi];
                while (j < hi) {
                    if (INST) n_kill += 1;
                    ctr.kill(@intCast(j - lo));
                    j += p * W30GAP[wp];
                    wp = (wp + 1) & 7;
                }
                next[bi] = j;
                wpos[bi] = wp;
            }
        }
        // Counter now holds φ(·,π(√z)) for this segment — read B's π(x/p) off it.
        while (pB > nay and x / @as(u64, primes[pB - 1]) < hi) {
            const v: u64 = x / @as(u64, primes[pB - 1]);
            if (v >= lo) {
                if (INST) n_bq += 1;
                b_sum += @as(i128, phi_run_full + ctr.prefix(@intCast(v - lo)) + naz_i - 1);
            }
            pB -= 1;
        }
        phi_run_full += ctr.total;
        for (3..nax) |bi| phi_run[bi] += seg_cnt[bi];
    }
    if (INST) std.debug.print("  ωB-stats: nax={d} naz={d} segs={d} mwalk={d} leaves(small/C/D)={d}/{d}/{d} kills={d} Bqueries={d}\n", .{ nax, naz, n_seg, n_mwalk, n_small, n_easy, n_hard, n_kill, n_bq });
    return .{ .omega = omega, .b = b_sum };
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
fn omegaNaive(primes: []const u32, mu: []const i8, delta: []const u32, pi: []const u32, x: u64, y: u64, xstar: u64) i128 {
    var omega: i128 = 0;
    for (primes, 0..) |p32, pidx| {
        const p: u64 = p32;
        if (p <= 2) continue;
        if (p > xstar) break;
        var m: u64 = y / p + 1;
        while (m <= y) : (m += 1) {
            if (mu[@intCast(m)] == 0) continue;
            if (delta[@intCast(m)] <= p) continue;
            const u = x / (p * m);
            omega += @as(i128, mu[@intCast(m)]) * @as(i128, phiRec(primes, pi, u, pidx));
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

fn chooseY(x: u64) u64 {
    const lo = icbrt(x) + 1;
    const hi = isqrt(x) -| 1;
    var y = 4 * icbrt(x);
    if (y < lo) y = lo;
    if (y > hi) y = hi;
    return y;
}

pub fn piGourdon(gpa: std.mem.Allocator, x: u64, y_in: ?u64) !GResult {
    return piGourdonV(gpa, x, y_in, false);
}

pub fn piGourdonV(gpa: std.mem.Allocator, x: u64, y_in: ?u64, verbose: bool) !GResult {
    var tp = common.nowNs();
    const y = y_in orelse chooseY(x);
    const z = x / y;
    const sqx = isqrt(x);
    const x13 = icbrt(x);
    const sqz = isqrt(z);
    const xstar = @max(isqrtsq(x), x / (y * y));

    var s = try Sieve.init(gpa, sqx, y);
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

    // A/Σ via binary-search π (all points ≤ √x).
    const t = computeTerms(&s, x, y, sqx, x13, sqz, xstar);
    lap(verbose, &tp, "A/Σ");

    // ω and B fused: one counter folded to √z serves ω's leaves and B's π(x/p).
    const wb = try omegaCounter(false, gpa, &s, x, y, z, xstar);
    const omega = wb.omega;
    const B = wb.b;
    lap(verbose, &tp, "ω+B");

    // φ₀ = Σ_{n≤y, n odd, μ(n)≠0} μ(n)·φ(x/n,1), φ(u,1)=u−⌊u/2⌋   [k=1]
    var phi0: i128 = 0;
    {
        var n: u64 = 1;
        while (n <= y) : (n += 1) {
            if (s.mu[@intCast(n)] == 0) continue;
            if (n % 2 == 0) continue;
            const u = x / n;
            phi0 += @as(i128, s.mu[@intCast(n)]) * @as(i128, @intCast(u - u / 2));
        }
    }
    lap(verbose, &tp, "phi0");

    // Σ closed forms
    const a = t.a;
    const bb = t.b;
    const cc = t.c;
    const dd = t.d;
    const P = t.P;
    const sig0: i128 = a - 1 + @divExact(P * (P - 1), 2) - @divExact(a * (a - 1), 2);
    const sig1: i128 = @divExact((a - bb) * (a - bb - 1), 2);
    const sig2: i128 = a * (bb - cc - @divExact(cc * (cc - 3), 2) + @divExact(dd * (dd - 3), 2));
    const sig3: i128 = @divExact(bb * (bb - 1) * (2 * bb - 1), 6) - bb - @divExact(dd * (dd - 1) * (2 * dd - 1), 6) + dd;
    const sigma = sig0 + sig1 + sig2 + sig3 + t.sig4 + t.sig5 + t.sig6;

    const pi = t.A - B + omega + phi0 + sigma;
    return .{ .pi = pi, .A = t.A, .B = B, .omega = omega, .phi0 = phi0, .sigma = sigma, .y = y };
}

pub fn main() !void {
    const gpa = std.heap.page_allocator;

    // ω differential check: naive recurrence vs O(1)-kill counter, must match exactly.
    std.debug.print("ω check (naive vs counter):\n", .{});
    for ([_]u64{ 100_000, 1_000_000, 10_000_000, 100_000_000, 1_000_000_000 }) |x| {
        const y = chooseY(x);
        const z = x / y;
        const xstar = @max(isqrtsq(x), x / (y * y));
        var s = try Sieve.init(gpa, isqrt(x), y);
        defer s.deinit(gpa);
        const pi = try buildPi(gpa, z);
        defer gpa.free(pi);
        const on = omegaNaive(s.primes, s.mu, s.delta, pi, x, y, xstar);
        const wb = try omegaCounter(false, gpa, &s, x, y, z, xstar);
        const b_ref = try computeB(false, gpa, &s, x, y, z, isqrt(x)); // standalone B reference
        std.debug.print("  {d:>12}  ω naive={d:>14} counter={d:>14} {s}   B fused={d} ref={d} {s}\n", .{ x, on, wb.omega, if (on == wb.omega) "match" else "MISMATCH", wb.b, b_ref, if (wb.b == b_ref) "match" else "MISMATCH" });
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
    _ = try piGourdonV(gpa, 1_000_000_000_000, null, true);

    std.debug.print("\nop-count comparison (gourdon ω/B vs lmo S2/P₂), matched x,y:\n", .{});
    for ([_]u64{ 1_000_000_000_000, 100_000_000_000_000 }) |x| {
        const y = chooseY(x);
        const z = x / y;
        const xstar = @max(isqrtsq(x), x / (y * y));
        var s = try Sieve.init(gpa, isqrt(x), y);
        defer s.deinit(gpa);
        std.debug.print("x=10^{d} (y={d}, z={d}):\n", .{ std.math.log10_int(x), y, z });
        std.debug.print("  gourdon(fused) ", .{});
        _ = try omegaCounter(true, gpa, &s, x, y, z, xstar); // prints ωB-stats (√z fold + B queries)
        const lr = try lmo.s2AndP2FusedInstrumented(gpa, @intCast(x), y, y);
        std.debug.print("  lmo:           S2 kills={d}  S2 prefix(hard)={d}  P₂ prefix(np)={d}\n", .{ lr.kills, lr.s2q, lr.np });
    }

    // Ceiling test: past the old O(z) π-table limit, gourdon-only vs known values.
    std.debug.print("\nceiling test (gourdon only, vs known):\n", .{});
    const cases = [_]struct { x: u64, want: i128 }{
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
        std.debug.print("  {d:>16} pi={d:>16} {s:>4}  {d:>7.2} s  peakRSS={d} MB\n", .{ cc.x, g.pi, if (g.pi == cc.want) "y" else "NO", gs, rss_mb });
    }
}
