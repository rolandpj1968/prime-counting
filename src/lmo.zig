//! LMO (Lagarias–Miller–Odlyzko) φ(x,a), building toward sub-linear π(x) with
//! O(x^(1/3)) memory — the path to π(10^18) in seconds.
//!
//! DERIVED FROM OUR OWN RECURSION, not from memory. φ(x,a)=φ(x,a−1)−φ(x/p_a,a−1)
//! takes primes in DECREASING index order, so with F(n,b) = μ(n)·φ(x/n,b):
//!
//!     F(n,b) = F(n,b−1) + F(n·p_b,b−1)          [μ(n·p_b) = −μ(n)]
//!
//! and the tree sums F exactly. Truncate by cutting any child with n·p_b > y:
//!
//!   • ORDINARY leaf — n ≤ y reaches b=0 ⇒ μ(n)·⌊x/n⌋. Every squarefree n ≤ y gets
//!     here (all its prefixes are ≤ n ≤ y), so S1 = Σ_{n≤y} μ(n)⌊x/n⌋. O(y), direct.
//!     Note c=0 is exactly right here — no wheel needed, there are only y of them.
//!
//!   • SPECIAL leaf — n = m·p, n > y ≥ m = n/p. Primes descend along the path, so
//!     the last one multiplied is the SMALLEST: p = P⁻(n), every prime of m above p.
//!     ⇒ S2 = Σ −μ(m)·φ(x/(m·p_b), b−1) over p_b-rough squarefree m ∈ (y/p_b, y].
//!
//! φ(x, π(y)) = S1 + S2, verified against meissel.phiOfXY.
//!
//! The b−1 index is what makes the sieve work: at step b the sieve over [1, z],
//! z = x/y, has exactly p_1..p_{b−1} removed, so sweeping b = 1..a is one monotone
//! prime-at-a-time removal, and φ(v, b−1) is a prefix count. Every special leaf has
//! v = x/(m·p_b) ≤ z because m·p_b > y — one sieve range serves them all.

const std = @import("std");
const common = @import("common.zig");
const rs = @import("rangesieve.zig");

fn icbrt(x: u64) u64 {
    if (x == 0) return 0;
    var r: u64 = @intFromFloat(std.math.pow(f64, @floatFromInt(x), 1.0 / 3.0));
    if (r == 0) r = 1;
    while (r * r * r > x) r -= 1;
    while ((r + 1) * (r + 1) * (r + 1) <= x) r += 1;
    return r;
}

/// Default knob: y = 2·x^(1/3). Measured argmin at 10^11, 10^12 and 10^13 alike
/// (6-7% better than 1.5). Unlike capped Meissel this is a SHARP interior optimum:
/// z = x/y falls as 1/α but the leaves ~ π(y)²/2 rise as α², so the curve is
/// gentle below (+36% at α=1) and brutal above (+1065% at α=16).
///
/// Caveat: the m-walk waste (see S2Result.walk) also rides α², so this argmin is a
/// LOWER bound — enumerating m properly should push it up.
pub fn defaultY(x: u64) u64 {
    return icbrt(x) * 2;
}

// ---------------------------------------------------------------- small sieves

/// μ(n) and least-prime-factor for n ≤ y, in one pass.
const SmallTables = struct {
    mu: []i8,
    lpf: []u32, // least prime factor; lpf[1] = maxInt so "lpf[m] > p" holds
    a: usize, // π(y)
    primes: []u32,

    fn init(gpa: std.mem.Allocator, y: u64) !SmallTables {
        const n: usize = @intCast(y + 1);
        const mu = try gpa.alloc(i8, n);
        @memset(mu, 1);
        const lpf = try gpa.alloc(u32, n);
        @memset(lpf, 0);
        var primes: std.ArrayList(u32) = .empty;

        var i: u64 = 2;
        while (i <= y) : (i += 1) {
            if (lpf[@intCast(i)] == 0) { // i is prime
                try primes.append(gpa, @intCast(i));
                var j = i;
                while (j <= y) : (j += i) {
                    if (lpf[@intCast(j)] == 0) lpf[@intCast(j)] = @intCast(i);
                    mu[@intCast(j)] = -mu[@intCast(j)];
                }
                var k = i * i;
                while (k <= y) : (k += i * i) mu[@intCast(k)] = 0;
            }
        }
        if (n > 1) lpf[1] = std.math.maxInt(u32);
        return .{ .mu = mu, .lpf = lpf, .a = primes.items.len, .primes = try primes.toOwnedSlice(gpa) };
    }

    fn deinit(self: *SmallTables, gpa: std.mem.Allocator) void {
        gpa.free(self.mu);
        gpa.free(self.lpf);
        gpa.free(self.primes);
    }
};

// ------------------------------------------------------------------- ordinary

pub const Foundation = struct { s1: i128, a: usize, y: u64 };

/// Ordinary leaves S1 = Σ_{n≤y} μ(n)⌊x/n⌋, and a = π(y).
pub fn ordinaryS1(gpa: std.mem.Allocator, x: u64, y: u64) !Foundation {
    var t = try SmallTables.init(gpa, y);
    defer t.deinit(gpa);
    var s1: i128 = 0;
    var n: u64 = 1;
    while (n <= y) : (n += 1) {
        const mn = t.mu[@intCast(n)];
        if (mn != 0) s1 += @as(i128, mn) * @as(i128, @intCast(x / n));
    }
    return .{ .s1 = s1, .a = t.a, .y = y };
}

// -------------------------------------------------------------------- Fenwick

/// Fenwick/BIT over [1, n] of 0/1 counts: O(log n) point update and prefix sum.
/// An implementation detail, not the idea — √-decomposition would also serve.
const Fenwick = struct {
    t: []i32,
    n: usize,

    inline fn lsb(i: usize) usize {
        return i & (~i +% 1); // two's-complement i & -i
    }

    /// All-ones init in O(n): t[i] already covers exactly lsb(i) elements.
    fn initAllOnes(gpa: std.mem.Allocator, n: usize) !Fenwick {
        const t = try gpa.alloc(i32, n + 1);
        t[0] = 0;
        for (1..n + 1) |i| t[i] = @intCast(lsb(i));
        return .{ .t = t, .n = n };
    }

    fn deinit(self: *Fenwick, gpa: std.mem.Allocator) void {
        gpa.free(self.t);
    }

    fn add(self: *Fenwick, i: usize, d: i32) void {
        var k = i;
        while (k <= self.n) : (k += lsb(k)) self.t[k] += d;
    }

    fn prefix(self: *const Fenwick, i: usize) i64 {
        var s: i64 = 0;
        var k = i;
        while (k > 0) : (k -= lsb(k)) s += self.t[k];
        return s;
    }
};

// -------------------------------------------------------------- alive-counter

/// Alive-set over one segment: **O(1) kill**, O(√S) prefix count, O(1) total.
///
/// Deliberately the opposite balance to Fenwick, because our traffic is lopsided:
/// folding the primes kills every element of [1,z] exactly once (z ≈ 1.4e9 at
/// 10^14) while the leaves only query a²/2 ≈ 2.4e7 times — ~60:1. Fenwick charges
/// O(log S) for BOTH, so it taxes the hot side to subsidise the cold one.
///
/// Two levels: a bit per element, plus an alive-count per block of `wpb` words.
/// A query sums whole-block counts then popcounts the words within its block, so
/// its cost is nblocks + wpb — minimised at wpb = √nwords.
const Counter = struct {
    bits: []u64, // 1 = alive
    cnt: []u32, // cnt[k] = alive in block k
    wpb: usize, // words per block
    nwords: usize,
    total: i64, // alive in the whole segment — makes seg_cnt O(1)

    fn init(gpa: std.mem.Allocator, seg: usize) !Counter {
        const nwords = (seg + 63) / 64;
        var wpb: usize = @intCast(common.isqrt(nwords));
        if (wpb == 0) wpb = 1;
        const nblocks = (nwords + wpb - 1) / wpb;
        return .{
            .bits = try gpa.alloc(u64, nwords),
            .cnt = try gpa.alloc(u32, nblocks),
            .wpb = wpb,
            .nwords = nwords,
            .total = 0,
        };
    }

    fn deinit(self: *Counter, gpa: std.mem.Allocator) void {
        gpa.free(self.bits);
        gpa.free(self.cnt);
    }

    /// Mark [0, len) alive, everything above dead.
    fn reset(self: *Counter, len: usize) void {
        const nw = (len + 63) / 64;
        @memset(self.bits[0..nw], ~@as(u64, 0));
        if (len % 64 != 0) self.bits[nw - 1] = (@as(u64, 1) << @as(u6, @intCast(len % 64))) - 1;
        @memset(self.bits[nw..self.nwords], 0);
        for (self.cnt, 0..) |*c, k| {
            var s: u32 = 0;
            const w0 = k * self.wpb;
            const w1 = @min(w0 + self.wpb, self.nwords);
            for (w0..w1) |w| s += @popCount(self.bits[w]);
            c.* = s;
        }
        self.total = @intCast(len);
    }

    inline fn kill(self: *Counter, i: usize) void {
        const w = i >> 6;
        const b = @as(u64, 1) << @as(u6, @intCast(i & 63));
        if (self.bits[w] & b != 0) {
            self.bits[w] &= ~b;
            self.cnt[w / self.wpb] -= 1;
            self.total -= 1;
        }
    }

    /// Alive in [0, i] inclusive.
    fn prefix(self: *const Counter, i: usize) i64 {
        const w = i >> 6;
        const blk = w / self.wpb;
        var s: i64 = 0;
        for (self.cnt[0..blk]) |c| s += c;
        for (self.bits[blk * self.wpb .. w]) |word| s += @popCount(word);
        const r: u6 = @intCast(i & 63);
        const mask: u64 = if (r == 63) ~@as(u64, 0) else (@as(u64, 1) << (r + 1)) - 1;
        s += @popCount(self.bits[w] & mask);
        return s;
    }
};

// --------------------------------------------------------------------- S2

/// `walk` counts m-candidates scanned; `leaves` how many survived. The gap is the
/// enumeration waste — we currently rescan (y/p_b, y] per b and reject most of it.
///
/// `easy` counts leaves where p_b² > v, i.e. no coprime composite is ≤ v, so
/// φ(v, b−1) = 1 + π(v) − (b−1) — a π-lookup, no sieve query needed. Equivalently
/// m·p_b³ > x, which is automatic once p_b > x^(1/3). A HARD leaf needs both
/// m > y/p and m ≤ x/p³, possible only when p < √z = x^(1/3)/√α — a window that
/// SHRINKS as α grows.
pub const S2Result = struct { s2: i128, leaves: u64, z: u64, a: usize, walk: u64 = 0, easy: u64 = 0 };

/// Special leaves S2 = Σ −μ(m)·φ(x/(m·p_b), b−1), by sweeping b = 1..a over a
/// sieve of [1, z], z = x/y, that holds exactly p_1..p_{b−1} removed at step b.
///
/// Non-segmented: the [1,z] sieve is O(x^(2/3)/α), so this is a CORRECTNESS
/// milestone, not the memory win. Segmenting [1,z] into x^(1/3) blocks (carrying
/// a running φ per b across blocks) is what buys O(x^(1/3)).
pub fn specialS2(gpa: std.mem.Allocator, x: u64, y: u64) !S2Result {
    var t = try SmallTables.init(gpa, y);
    defer t.deinit(gpa);
    const z = x / y;
    const zn: usize = @intCast(z);

    var fen = try Fenwick.initAllOnes(gpa, zn);
    defer fen.deinit(gpa);
    const alive = try gpa.alloc(bool, zn + 1);
    defer gpa.free(alive);
    @memset(alive, true);

    var s2: i128 = 0;
    var leaves: u64 = 0;

    for (t.primes) |p32| {
        const p: u64 = p32; // p = p_b; sieve currently holds exactly p_1..p_{b-1}
        // m must be p-rough squarefree with m·p > y and m ≤ y.
        var m: u64 = y / p + 1;
        while (m <= y) : (m += 1) {
            const mm = t.mu[@intCast(m)];
            if (mm == 0) continue; // squarefree only
            if (t.lpf[@intCast(m)] <= p) continue; // need P⁻(m) > p
            const v = x / (m * p); // ≤ z, since m·p > y
            s2 += @as(i128, -mm) * @as(i128, fen.prefix(@intCast(v)));
            leaves += 1;
        }
        // now fold p_b into the sieve, ready for b+1
        var j = p;
        while (j <= z) : (j += p) {
            if (alive[@intCast(j)]) {
                alive[@intCast(j)] = false;
                fen.add(@intCast(j), -1);
            }
        }
    }
    return .{ .s2 = s2, .leaves = leaves, .z = z, .a = t.a };
}

/// SEGMENTED S2 — the same sum, in O(x^(1/3)) memory instead of O(x^(2/3)).
///
/// Two ideas carry it:
///  • running φ per b: φ(v, b−1) = phi_run[b−1] + (alive in [lo, v]), where
///    phi_run[b−1] counts survivors in [1, lo) with p_1..p_{b−1} removed. Each
///    segment adds its own count to phi_run, so [1, lo) is never re-sieved.
///  • descending m-cursor per prime: segments run lo ASCENDING, and v = x/(m·p)
///    rises as m falls — so one cursor per p_b walks m down monotonically and
///    every m is touched once in TOTAL, not once per segment.
pub fn specialS2Segmented(gpa: std.mem.Allocator, x: u64, y: u64, seg: usize) !S2Result {
    var t = try SmallTables.init(gpa, y);
    defer t.deinit(gpa);
    const z = x / y;
    const a = t.a;

    var ctr = try Counter.init(gpa, seg);
    defer ctr.deinit(gpa);

    const phi_run = try gpa.alloc(i64, a); // survivors in [1, lo) per b−1
    defer gpa.free(phi_run);
    @memset(phi_run, 0);
    const seg_cnt = try gpa.alloc(i64, a); // this segment's contribution
    defer gpa.free(seg_cnt);
    const m_cur = try gpa.alloc(u64, a); // descending m cursor per prime
    defer gpa.free(m_cur);
    @memset(m_cur, y);

    var s2: i128 = 0;
    var leaves: u64 = 0;
    var walk: u64 = 0;
    var easy: u64 = 0;

    var lo: u64 = 1;
    while (lo <= z) : (lo += seg) {
        const hi = @min(lo + seg, z + 1); // segment covers [lo, hi)
        const len: usize = @intCast(hi - lo);
        ctr.reset(len);

        for (t.primes, 0..) |p32, bi| {
            const p: u64 = p32; // p = p_b, b = bi+1; segment holds primes[0..bi) removed
            seg_cnt[bi] = ctr.total; // survivors here with p_1..p_{b−1} gone — O(1)

            // walk m down while v = x/(m·p) stays inside this segment
            var m = m_cur[bi];
            const mlo = y / p; // m·p > y  ⇔  m > ⌊y/p⌋
            while (m > mlo) {
                const v = x / (m * p);
                if (v >= hi) break; // belongs to a later segment
                walk += 1;
                if (v >= lo) {
                    const mm = t.mu[@intCast(m)];
                    if (mm != 0 and t.lpf[@intCast(m)] > p) {
                        const phi_v = phi_run[bi] + ctr.prefix(@intCast(v - lo));
                        s2 += @as(i128, -mm) * @as(i128, phi_v);
                        leaves += 1;
                        if (p * p > v) easy += 1; // φ(v,b−1) = 1 + π(v) − (b−1)
                    }
                }
                m -= 1;
            }
            m_cur[bi] = m;

            // fold p_b into this segment, ready for b+1
            var j = ((lo + p - 1) / p) * p; // first multiple of p at or above lo
            while (j < hi) : (j += p) ctr.kill(@intCast(j - lo));
        }
        for (0..a) |bi| phi_run[bi] += seg_cnt[bi];
    }
    return .{ .s2 = s2, .leaves = leaves, .z = z, .a = a, .walk = walk, .easy = easy };
}

pub const PhiResult = struct { phi: i128, s1: i128, s2: i128, leaves: u64, z: u64, a: usize, y: u64, walk: u64 = 0, easy: u64 = 0 };

/// φ(x, π(y)) = S1 + S2 — the LMO decomposition end to end.
pub fn phiLMO(gpa: std.mem.Allocator, x: u64, y: u64, seg: ?usize) !PhiResult {
    const f = try ordinaryS1(gpa, x, y);
    const s = if (seg) |sz| try specialS2Segmented(gpa, x, y, sz) else try specialS2(gpa, x, y);
    return .{ .phi = f.s1 + s.s2, .s1 = f.s1, .s2 = s.s2, .leaves = s.leaves, .z = s.z, .a = s.a, .y = y, .walk = s.walk, .easy = s.easy };
}

// ---------------------------------------------------------------------- P₂

/// P₂(x,y) = Σ_{y<p≤√x} (π(x/p) − π(p) + 1) — the n ≤ x that are a product of
/// exactly two primes, both > y.
///
/// No Fenwick needed: this is the MONOTONE sweep. p ascending ⇔ x/p descending,
/// so walking segments of [1, z] upward with a running π and a cursor that walks
/// p DOWNWARD from √x visits every (p, x/p) exactly once. Within a segment the
/// v = x/p arrive in ascending order too, so one linear walk with a running count
/// answers every π(x/p) — O(z) total, no random access anywhere.
///
/// The π(p) half needs no lookup at all: primes in (y, √x] have indices a+1..A,
/// so Σ(π(p)−1) = Σ_{j=a}^{A−1} j collapses to a closed form.
pub fn p2Segmented(gpa: std.mem.Allocator, x: u64, y: u64, seg: usize) !i128 {
    const sqrt_x = common.isqrt(x);
    if (y >= sqrt_x) return 0; // no primes in (y, √x]
    const primes = try rs.basePrimes(gpa, sqrt_x);
    defer gpa.free(primes);

    var a: usize = 0; // π(y)
    for (primes) |p| {
        if (p <= y) a += 1 else break;
    }
    const A: usize = primes.len; // π(√x)
    if (A <= a) return 0;

    const z = x / y;
    const isprime = try gpa.alloc(bool, seg);
    defer gpa.free(isprime);

    var sum_pi_xp: i128 = 0;
    var running_pi: i64 = 0; // π(lo − 1)
    var idx: usize = A; // cursor: next p to process is primes[idx−1], descending

    var lo: u64 = 1;
    while (lo <= z) : (lo += seg) {
        const hi = @min(lo + seg, z + 1);
        const len: usize = @intCast(hi - lo);
        @memset(isprime[0..len], true);
        if (lo == 1) isprime[0] = false; // 1 is not prime
        for (primes) |q| {
            if (q * q >= hi) break;
            var j = @max(q * q, ((lo + q - 1) / q) * q);
            while (j < hi) : (j += q) isprime[@intCast(j - lo)] = false;
        }

        // v = x/p arrives ascending as p descends → one linear walk serves all
        var walk: u64 = lo;
        var cnt: i64 = 0; // primes in [lo, walk−1]
        while (idx > a) {
            const p = primes[idx - 1];
            const v = x / p;
            if (v >= hi) break; // this p belongs to a later segment
            while (walk <= v) : (walk += 1) {
                if (isprime[@intCast(walk - lo)]) cnt += 1;
            }
            sum_pi_xp += running_pi + cnt; // = π(x/p)
            idx -= 1;
        }
        while (walk < hi) : (walk += 1) {
            if (isprime[@intCast(walk - lo)]) cnt += 1;
        }
        running_pi += cnt;
        if (idx <= a) break; // every p consumed
    }

    // Σ_{y<p≤√x} (π(p) − 1) = Σ_{j=a}^{A−1} j, since those p have indices a+1..A
    const Ai: i128 = @intCast(A);
    const ai: i128 = @intCast(a);
    // both are products of consecutive integers, so both halves are exact
    const sub = @divExact((Ai - 1) * Ai, 2) - @divExact(ai * (ai - 1), 2);
    return sum_pi_xp - sub;
}

// ---------------------------------------------------------------------- π(x)

pub const PiResult = struct { pi: u64, phi: i128, p2: i128, y: u64, a: usize, z: u64, leaves: u64 };

/// π(x) = φ(x, a) + a − 1 − P₂(x, a), a = π(y) — LMO end to end.
pub fn piLMO(gpa: std.mem.Allocator, x: u64, y_in: ?u64, seg_in: ?usize) !PiResult {
    if (x < 2) return .{ .pi = 0, .phi = 0, .p2 = 0, .y = 0, .a = 0, .z = 0, .leaves = 0 };
    const y = y_in orelse defaultY(x);
    const seg: usize = seg_in orelse @intCast(@max(y, 1024));
    const f = try phiLMO(gpa, x, y, seg);
    const p2 = try p2Segmented(gpa, x, y, seg);
    const r = f.phi + @as(i128, @intCast(f.a)) - 1 - p2;
    return .{ .pi = @intCast(r), .phi = f.phi, .p2 = p2, .y = y, .a = f.a, .z = f.z, .leaves = f.leaves };
}
