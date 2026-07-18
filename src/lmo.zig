//! LMO (Lagarias–Miller–Odlyzko / Deléglise–Rivat) φ(x,a): sub-linear π(x) in
//! O(x^(1/3)) memory. Measured exponent 0.66 (the theoretical 2/3); π(10^14) exact
//! in ~1.2 s, ~34× the Meissel–Lehmer in meissel.zig; reaches π(10^20) past 2^64
//! (u128) in 3.4 h single-threaded. Single-thread it is ~3.8× off primecount's DR;
//! parallelism (untouched) is the biggest remaining lever. See COMBINATORIAL.md.
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
    return icbrtG(u64, x);
}

/// floor(x^(1/3)) for x of any unsigned width; the root itself always fits u64
/// (x ≤ 10²⁰ ⇒ root ≤ 4.7×10⁶). Cubes are formed in T so they don't overflow the
/// narrow type at the top of the range.
fn icbrtG(comptime T: type, x: T) u64 {
    if (x == 0) return 0;
    var r: u64 = @intFromFloat(std.math.pow(f64, @floatFromInt(x), 1.0 / 3.0));
    if (r == 0) r = 1;
    while (@as(T, r) * @as(T, r) * @as(T, r) > x) r -= 1;
    while (@as(T, r + 1) * @as(T, r + 1) * @as(T, r + 1) <= x) r += 1;
    return r;
}

/// floor(√n) for n of any unsigned width; the root always fits u64 (n ≤ 10²⁰ ⇒
/// root ≤ 10¹⁰). Squares are formed in T.
fn isqrtG(comptime T: type, n: T) u64 {
    if (n < 2) return @intCast(n);
    var x: u64 = @intFromFloat(@sqrt(@as(f64, @floatFromInt(n))));
    while (@as(T, x) * @as(T, x) > n) x -= 1;
    while (@as(T, x + 1) * @as(T, x + 1) <= n) x += 1;
    return x;
}

/// Default knob: y = 4·x^(1/3). Measured argmin at 10^11, 10^12, 10^13 and 10^14
/// alike. An interior optimum, because the two terms fight: z = x/y falls as 1/α
/// but the leaves ≈ π(y)²/2 rise as α².
///
/// It was 2 before the √y enumeration split. The old m-walk waste ALSO rode α², so
/// removing it flattened the high-α side (α=8 cost +250%, now +41%) and let the
/// optimum drift up — worth a further ~23%. α_opt has never scaled with x here:
/// flat 1.5 in capped Meissel, flat 2, now flat 4 over four powers of ten, where the
/// literature's α ~ log³x would demand a steady rise.
///
/// y is clamped to √x. The identity π(x) = φ(x,a) + a − 1 − P₂ needs every one of
/// the a primes to be ≤ x; y > x makes a count primes above x and π comes out too
/// high (π(2) returned 2 with y = 4·icbrt(2) = 4, which counts 3). So the invariant
/// is x^(1/3) ≤ y ≤ √x — the floor is the P₃ bound, the ceiling keeps p_a ≤ √x ≤ x.
/// Always satisfiable (icbrt ≤ isqrt), and it only binds below x = 4096.
pub fn defaultY(x: u128) u64 {
    return @min(icbrtG(u128, x) * 4, isqrtG(u128, x));
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

/// Ordinary leaves S1 = Σ_{n≤y} μ(n)⌊x/n⌋, and a = π(y). u64 entry point.
pub fn ordinaryS1(gpa: std.mem.Allocator, x: u64, y: u64) !Foundation {
    return ordinaryS1Gen(u64, gpa, x, y);
}

/// Generic over the value type X (u64 / u128). ⌊x/n⌋ ≤ x fits X; only the numerator
/// is wide. The result ⌊x/n⌋ can exceed u64 for n small at x = 10²⁰, so it stays X
/// until multiplied into the i128 accumulator.
pub fn ordinaryS1Gen(comptime X: type, gpa: std.mem.Allocator, x: X, y: u64) !Foundation {
    var t = try SmallTables.init(gpa, y);
    defer t.deinit(gpa);
    var s1: i128 = 0;
    var n: u64 = 1;
    while (n <= y) : (n += 1) {
        const mn = t.mu[@intCast(n)];
        if (mn != 0) s1 += @as(i128, mn) * @as(i128, @intCast(x / @as(X, n)));
    }
    return .{ .s1 = s1, .a = t.a, .y = y };
}

// -------------------------------------------------------------------- Fenwick

/// Fenwick/BIT over [1, n] of 0/1 counts: O(log n) point update and prefix sum.
///
/// KEPT ONLY as the flat specialS2's reference implementation, so that
/// flat (Fenwick) == segmented (Counter) cross-checks the Counter. It is the wrong
/// structure for the hot path — see Counter below — and measured 1.9× slower there.
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

// ------------------------------------------------------------------- mod-30 wheel

/// DR §9: "Precomputing the sieving by the first primes 2, 3, 5 saves some more time."
///
/// The profile said the fold is ~50% of cycles and that it is VOLUME, not
/// misprediction: it runs Σ_p z/p ≈ 2.76z times but only z kills land, and the
/// alive-check mispredicts at just 1.18%. The only lever is fewer iterations.
///
/// Key split: the wheel's STEPPING (which multiples to visit) is independent of the
/// array's INDEXING (bit-per-integer vs 8-per-30). All the fold win is in the
/// stepping, so we take that and leave the indexing — and hence every leaf/prefix
/// path — untouched:
///   Σ_{p≥7} (8/30)·z/p ≈ 0.46z   vs   Σ_{p≤y} z/p ≈ 2.76z   → ~6× fewer visits.
const W30 = [8]u8{ 1, 7, 11, 13, 17, 19, 23, 29 }; // residues coprime to 30
const W30GAP = [8]u8{ 6, 4, 2, 4, 2, 4, 6, 2 }; // 1→7→11→13→17→19→23→29→31; Σ = 30

// Coprimality lookups mod 30, for fast-forwarding the fold cursor to a block start.
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

/// Bit pattern of "coprime to 30" over u64 words. lcm(30, 64) = 960 bits = 15 words,
/// so a segment starting at a multiple of 960 has mask MASK30[w % 15].
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

/// φ(v, b) for b ≤ 3 in closed form — meissel.zig's base cases. Leaves with
/// p ∈ {2,3,5} need φ(v,0..2), which a 2·3·5-wheel cannot represent (it has already
/// removed them). There are only ~0.5y such leaves (~93k of 141M at 10^14).
/// b=3 is the wheel's own initial state, so the counter answers that one directly.
inline fn phiSmall(v: u64, b: usize) i64 {
    return switch (b) {
        0 => @intCast(v),
        1 => @intCast(v - v / 2),
        2 => @intCast(v - v / 2 - v / 3 + v / 6),
        else => unreachable,
    };
}

// -------------------------------------------------------------- alive-counter

/// Alive-set over one segment: **O(1) kill**, O(√S) prefix count, O(1) total.
///
/// Two levels: a bit per element, plus an alive-count per block of `wpb` words.
/// A query sums whole-block counts then popcounts the words within its block, so
/// its cost is nblocks + wpb — minimised at wpb = √nwords.
///
/// Beat Fenwick by 1.9× when adopted (α was 1.5 then, and kills outnumbered queries
/// ~60:1, so paying O(log S) on the hot kill side to subsidise the cold query side
/// was clearly wrong). NOTE the traffic has since shifted: at α=4 it is 5.4e8 kills
/// vs 1.4e8 leaves at 10^14 — only ~4:1 by count, and in COST the queries now
/// dominate the other way (1.4e8 × ~45 ≈ 6.3e9 vs 5.4e8). The Counter still wins
/// (~6.8e9 vs Fenwick's ~1.2e10), but its margin no longer comes from where this
/// comment used to claim.
///
/// SUPERSEDED by Counter3P, which is what s2AndP2Fused now uses (1.15× at 10^13,
/// 1.19× at 10^14, margin growing with x). Kept as the flat/reference counter.
///
/// The measurement corrected the attribution: I predicted ~20% from the third level
/// making the query 3·nwords^(1/3) rather than 2·√nwords (43 vs 108 at 10^14). The
/// total is indeed ~19% at 10^14 — but most of it is removing the DIVISION below, not
/// the extra level (at 10^13 the third level adds nothing over a 2-level shift).
///
/// `w / self.wpb` is the cost: wpb is a runtime divisor, and the quotient is an
/// ADDRESS, so the dependent load-modify-store sits on the critical path. Contrast
/// the fold's `((lo+p−1)/p)·p`, which measured free — there the quotient only seeded
/// a loop whose iterations were independent across primes, so ILP hid the latency.
/// Division latency is free when it feeds independent work, and expensive when it
/// feeds an address you must immediately touch.
pub const Counter = struct {
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

/// Power-of-two variants, so `w / wpb` in kill() becomes a shift. Counter's block
/// size is isqrt(nwords), a runtime divisor — a real division on the hot kill path.
/// Counter2P isolates that effect from the level count; Counter3P adds the level.
inline fn log2Floor(n: usize) u6 {
    return @intCast(63 - @clz(@as(u64, @max(n, 1))));
}

/// 2-level, power-of-two block (shift, not divide). Same shape as Counter otherwise.
pub const Counter2P = struct {
    bits: []u64,
    cnt: []u32,
    s1: u6, // block = 2^s1 words ≈ √nwords
    nwords: usize,
    total: i64,

    fn init(gpa: std.mem.Allocator, seg: usize) !Counter2P {
        const nwords = (seg + 63) / 64;
        const s1: u6 = @intCast((log2Floor(nwords) + 1) / 2);
        const nblocks = (nwords >> s1) + 1;
        return .{
            .bits = try gpa.alloc(u64, nwords),
            .cnt = try gpa.alloc(u32, nblocks),
            .s1 = s1,
            .nwords = nwords,
            .total = 0,
        };
    }
    fn deinit(self: *Counter2P, gpa: std.mem.Allocator) void {
        gpa.free(self.bits);
        gpa.free(self.cnt);
    }
    fn reset(self: *Counter2P, len: usize) void {
        const nw = (len + 63) / 64;
        @memset(self.bits[0..nw], ~@as(u64, 0));
        if (len % 64 != 0) self.bits[nw - 1] = (@as(u64, 1) << @as(u6, @intCast(len % 64))) - 1;
        @memset(self.bits[nw..self.nwords], 0);
        @memset(self.cnt, 0);
        for (self.bits, 0..) |word, w| self.cnt[w >> self.s1] += @popCount(word);
        self.total = @intCast(len);
    }
    inline fn kill(self: *Counter2P, i: usize) void {
        const w = i >> 6;
        const b = @as(u64, 1) << @as(u6, @intCast(i & 63));
        if (self.bits[w] & b != 0) {
            self.bits[w] &= ~b;
            self.cnt[w >> self.s1] -= 1;
            self.total -= 1;
        }
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

/// 3-level: bits → blocks of 2^s1 words → superblocks of 2^s2 blocks.
///
/// A query costs nsuper + 2^s2 + 2^s1, minimised at both ≈ nwords^(1/3) → 3·n^(1/3)
/// instead of the 2-level's 2·√n (43 vs 108 at 10^14). Pays one extra counter
/// decrement per kill. Worth trying because the traffic INVERTED at α=4: queries now
/// dominate cost ~12:1, the opposite of when the 2-level was chosen at α=1.5.
pub const Counter3P = struct {
    bits: []u64,
    cnt1: []u32, // per block of 2^s1 words
    cnt2: []u32, // per superblock of 2^s2 blocks
    s1: u6,
    s2: u6,
    nwords: usize,
    total: i64,

    fn init(gpa: std.mem.Allocator, seg: usize) !Counter3P {
        const nwords = (seg + 63) / 64;
        const t: u6 = @intCast(@max(1, (log2Floor(nwords) + 2) / 3)); // ≈ n^(1/3)
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
    /// Reset to the mod-30 WHEEL state — 2, 3, 5 already struck — not all-ones.
    /// Requires the segment to start at a multiple of 960 = lcm(30, 64) so the mask
    /// is word-aligned. This is φ(·, 3), so the b-loop starts at bi = 3.
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
    /// BRANCHLESS — and this only became right AFTER the mod-30 wheel.
    /// Pre-wheel it measured 0.75× (a 33% regression): the alive-check missed at just
    /// 1.18%, because p=2/3/5 dominated the fold and their already-dead patterns are
    /// short and periodic, so the predictor nailed them; paying 3 unconditional RMWs
    /// over 1.6e9 visits to remove a 98.8%-correct branch was a bad trade.
    /// The wheel deleted exactly those primes. What is left is p≥7, whose aliveness
    /// follows irregular factorizations: the miss rate jumped to 9.0% (22.5M over
    /// 2.5e8 visits, 38% of all branch misses) while visits fell 6.4×. Both terms
    /// moved, so the trade flipped.
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

// --------------------------------------------------------------------- S2

/// `walk` counts m-candidates scanned; `leaves` how many survived. The gap was the
/// enumeration waste (w/leaf ≈ 20, i.e. ~2·ln y): we rescanned (y/p_b, y] per b and
/// rejected most of it. Fixed by the √y split below — w/leaf is now ≈ 1.12, so this
/// pair is a regression guard rather than a live diagnostic.
///
/// `easy` counts leaves where p_b² > v, i.e. no coprime composite is ≤ v, so
/// φ(v, b−1) = 1 + π(v) − (b−1) — a π-lookup, no sieve query needed. Equivalently
/// m·p_b³ > x, which is automatic once p_b > x^(1/3). A HARD leaf needs both
/// m > y/p and m ≤ x/p³, possible only when p < √z = x^(1/3)/√α — a window that
/// SHRINKS as α grows. Measured: 93% easy at α=2, 99.2% at α=16, and hard/a²
/// collapses, so the whole α² leaf growth is easy leaves.
///
/// DIAGNOSTIC ONLY — we do not branch on it. Resolving easy leaves via π(v) instead
/// of the counter was tried and measured 5% SLOWER, and cannot help in principle:
/// for an easy leaf the survivors ≤ v with lpf ≥ p_b are ALREADY just 1 plus the
/// primes in [p_b, v], so ctr.prefix(v) IS φ(v,b−1). See s2AndP2Fused.
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
///  • descending cursor per prime: segments run lo ASCENDING, and v = x/(m·p) rises
///    as m falls — so one cursor per p_b walks down monotonically and every candidate
///    is touched once in TOTAL, not once per segment. Two regimes, split at √y (see
///    the cursor init below): p > √y enumerates m = q prime straight from primes[],
///    p ≤ √y walks m.
///
/// This is now piLMO's REFERENCE path, not its hot path — piLMO calls s2AndP2Fused,
/// which additionally folds P₂ into this same sweep. Kept so the two can be
/// differentially cross-checked.
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
    // Descending cursor per prime. Two regimes, split at √y:
    //  • p ≤ √y — cur = m, walking down through p-rough squarefree m ∈ (y/p, y].
    //  • p > √y — cur = a PRIME INDEX. Any p-rough m ≤ y is 1 or a prime (a
    //    composite would be ≥ lpf(m)² > y), and m=1 fails m·p > y, so m = q prime
    //    in (p, y], with m·p > y automatic since q > p > √y ≥ y/p. Enumerable
    //    straight from primes[] — no walk, no rejects. This is the a²/2 bulk.
    const sqrt_y = common.isqrt(y);
    const cur = try gpa.alloc(u64, a);
    defer gpa.free(cur);
    for (t.primes, 0..) |p32, bi| cur[bi] = if (p32 <= sqrt_y) y else a - 1;

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

            if (p <= sqrt_y) {
                // dense: walk m down while v = x/(m·p) stays inside this segment
                var m = cur[bi];
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
                cur[bi] = m;
            } else {
                // sparse: m = q prime in (p, y] only. μ(q) = −1, so −μ(m) = +1.
                var qi = cur[bi];
                while (qi > bi) {
                    const q: u64 = t.primes[@intCast(qi)];
                    const v = x / (p * q);
                    if (v >= hi) break; // belongs to a later segment
                    walk += 1;
                    if (v >= lo) {
                        const phi_v = phi_run[bi] + ctr.prefix(@intCast(v - lo));
                        s2 += @as(i128, phi_v);
                        leaves += 1;
                        if (p * p > v) easy += 1;
                    }
                    qi -= 1;
                }
                cur[bi] = qi;
            }

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
/// REFERENCE path only — piLMO uses s2AndP2Fused, which reads P₂ off the S2 counter
/// and sieves the p-ranges on the fly. This version stores every prime ≤ √x, which is
/// Θ(√x/ln x) (406 MB at 10^18); kept because it is independently verified and lets
/// the fused path be differentially cross-checked.
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

/// φ(v, b−1) for one special leaf — LMO p.555 classes (2) and (3)-easy.
///
/// If p² > v then no coprime composite is ≤ v (the smallest is p², and p = p_b is the
/// least prime not yet removed), so the survivors are exactly 1 plus the primes in
/// [p, v]:  **φ(v, b−1) = 1 + max(0, π(v) − (b−1))**. If additionally v ≤ y, that π
/// comes from a table over [1, y] — O(y) space, which is O(x^(1/3)) and affordable.
/// This one rule covers all three cheap classes:
///   • class (2), p > √(x/y): v < x/p² < y, and p > x^(1/3)/2 > x^(1/4) so p⁴ > x
///     ⇒ v ≤ x/p² < p².
///   • class (3)-easy, q ≥ x/(yp): v ≤ y, and p² > (x/y²)² > y ≥ v — which needs
///     x² > y⁵, i.e. **y ≤ x^(2/5)**. That is exactly why LMO's Truncation Rule T′
///     caps y there. α=4 satisfies it for x ≥ 4^15 ≈ 10^9.
/// Everything else (class (3)-hard, class (4)) still needs the sieve.
inline fn leafPhi(
    v: u64,
    p: u64,
    bi: usize,
    lo: u64,
    y: u64,
    pi_tab: []const u32,
    ctr: anytype,
    phi_run: []const i64,
) i64 {
    if (bi <= 2) return phiSmall(v, bi); // p ∈ {2,3,5}: the wheel cannot answer these
    if (v <= y and p * p > v) {
        const piv: i64 = pi_tab[@intCast(v)];
        return 1 + @max(0, piv - @as(i64, @intCast(bi)));
    }
    return phi_run[bi] + ctr.prefix(@intCast(v - lo));
}

/// leafPhi, plus the block-parallel bookkeeping: only the HARD path uses phi_run
/// (the running φ), so only there does the cross-block correction apply. We record
/// mu_sum[bi] += sign (the leaf's −μ(m) contribution) for exactly those leaves, so
/// the master can add prefix_φ_at_block_start[bi]·mu_sum[bi] after the fact. Closed-
/// form and π-table leaves are block-local and need no correction. For single-block
/// (serial) runs mu_sum is accumulated but never read, so the result is unchanged.
inline fn leafPhiP(
    v: u64,
    p: u64,
    bi: usize,
    lo: u64,
    y: u64,
    pi_tab: []const u32,
    ctr: anytype,
    phi_run: []const i64,
    sign: i64,
    mu_sum: []i64,
) i64 {
    if (bi <= 2) return phiSmall(v, bi);
    if (v <= y and p * p > v) {
        const piv: i64 = pi_tab[@intCast(v)];
        return 1 + @max(0, piv - @as(i64, @intCast(bi)));
    }
    mu_sum[bi] += sign; // hard leaf: this bi contributes to the cross-block φ correction
    return phi_run[bi] + ctr.prefix(@intCast(v - lo));
}

// ------------------------------------------------------------- fused sweep

/// S2 and P₂ from ONE pass over [1, z] — P₂ reads the SAME counter, for free.
///
/// After folding all a primes, the alive set in [1, z] is exactly {1} ∪ primes in
/// (y, z]: any composite with lpf > y is ≥ lpf² > y², and y² > z whenever α > 1
/// (y² = α²x^(2/3) vs z = x^(2/3)/α). So φ(v, a) = 1 + π(v) − a, i.e.
///
///     π(v) = φ(v, a) − 1 + a
///
/// and P₂'s π(x/p) is just a counter query at the end of each segment's b loop.
/// That deletes P₂'s entire separate sweep over [1, z] at zero added cost.
///
/// NOTE — what does NOT work: resolving "easy" leaves (p² > v) via π(v) instead of
/// the counter. For an easy leaf the survivors ≤ v with lpf ≥ p_b are ALREADY just
/// 1 plus the primes in [p_b, v] (a composite would be ≥ p_b² > v), so ctr.prefix(v)
/// IS φ(v,b−1) — the π-formula is a different route to the identical number, not a
/// cheaper one. Building a flat O(1)-query prefix array to serve them measured 5%
/// SLOWER: the build is O(z) but the queries are only the leaves, and z/leaves ≈ 27.
/// DR's gain must therefore come from CLUSTERING (fewer queries), not cheaper ones.
///
/// `s2` and `p2` stay separate accumulators so each can be validated on its own
/// against specialS2Segmented / p2Segmented rather than only end-to-end via π(x).
pub const FusedResult = struct {
    s2: i128,
    p2: i128,
    leaves: u64,
    easy: u64,
    walk: u64,
    z: u64,
    a: usize,
};

/// Production path: diagnostics OFF. The `walk`/`leaves`/`easy` counters cost 3.3% at
/// 10^14 and 9.3% at 10^9 (the share falls as the class-(1) binomial absorbs more
/// leaves) — `easy` alone recomputes the `v <= y and p*p > v` test that leafPhi does
/// again one line later. They are analysis instruments, not part of π(x).
pub fn s2AndP2Fused(gpa: std.mem.Allocator, x: u64, y: u64, seg: usize) !FusedResult {
    return s2AndP2FusedGen(Counter3P, false, u64, gpa, x, y, seg);
}

/// Diagnostics ON — `leaves` confirmed LMO Lemma 5.1 to 0.03% at 10^19, `walk` is the
/// regression guard on the √y split, `easy` sized the π-table class. Kept, just not paid
/// for on every run.
pub fn s2AndP2FusedInstrumented(gpa: std.mem.Allocator, x: u64, y: u64, seg: usize) !FusedResult {
    return s2AndP2FusedGen(Counter3P, true, u64, gpa, x, y, seg);
}

/// Generic over the alive-counter, the instrumentation, and the value type X (u64 /
/// u128 — the latter for x > 2⁶⁴, i.e. π(10²⁰)). Only x and the divisions with x in
/// the numerator are X-wide; every quotient (z, v = x/(m·p), x/p) fits u64 and is cast
/// down immediately, so y, a, primes, m, p, q and the whole counter stay u64.
pub fn s2AndP2FusedGen(comptime C: type, comptime INST: bool, comptime X: type, gpa: std.mem.Allocator, x: X, y: u64, seg: usize) !FusedResult {
    return s2P2Blocks(C, INST, X, gpa, x, y, seg, 1);
}

/// Everything one block sweep touches. The read-only tables/scalars have shared slice
/// headers (immutable during the sweep, so threads share them); the scratch is
/// per-thread; the blk_* outputs are written by absolute block index, so threads never
/// contend. runBlock is the single source of truth for the sweep — the serial driver
/// and the (future) threaded one both call it.
fn BlkCtx(comptime C: type, comptime X: type) type {
    return struct {
        const Self = @This();
        // read-only shared
        primes: []const u32,
        mu: []const i8,
        lpf: []const u32,
        pi_tab: []const u32,
        bp4: []const u64,
        x: X,
        y: u64,
        z: u64,
        sqrt_x: u64,
        sqrt_y: u64,
        p_cube_min: u64,
        segw: usize,
        a: usize,
        nb: usize,
        do_p2: bool,
        // per-thread scratch
        ctr: C,
        phi_run: []i64,
        seg_cnt: []i64,
        cur: []u64,
        next: []u64,
        wpos: []u8,
        pbuf: []bool,
        mu_sum: []i64,
        // per-block outputs (disjoint writes)
        blk_s2: []i128,
        blk_pi: []i128,
        blk_np: []u64,
        blk_total: []i64,
        blk_mu: []i64,
        // INST diagnostics (per-thread; summed after)
        leaves: u64 = 0,
        easy: u64 = 0,
        walk: u64 = 0,

        fn runBlock(self: *Self, comptime INST: bool, blk: usize) void {
            const x = self.x;
            const y = self.y;
            const z = self.z;
            const a = self.a;
            const segw = self.segw;
            const sqrt_y = self.sqrt_y;
            const sqrt_x = self.sqrt_x;
            const p_cube_min = self.p_cube_min;
            const do_p2 = self.do_p2;
            const pi_tab = self.pi_tab;
            const bp4 = self.bp4;
            const cur = self.cur;
            const next = self.next;
            const wpos = self.wpos;
            const phi_run = self.phi_run;
            const seg_cnt = self.seg_cnt;
            const mu_sum = self.mu_sum;
            const pbuf = self.pbuf;

            const nsegs = z / segw + 1;
            const block_lo: u64 = @min(@as(u64, blk * nsegs / self.nb) * segw, z + 1);
            const block_hi: u64 = @min(@as(u64, (blk + 1) * nsegs / self.nb) * segw, z + 1);
            if (block_lo >= block_hi) { // empty (nb > nsegs)
                self.blk_s2[blk] = 0;
                self.blk_pi[blk] = 0;
                self.blk_np[blk] = 0;
                @memset(self.blk_total[blk * (a + 1) ..][0 .. a + 1], 0);
                @memset(self.blk_mu[blk * (a + 1) ..][0 .. a + 1], 0);
                return;
            }

            // Fast-forward cursors to block_lo.
            for (self.primes, 0..) |p32, bi| {
                const p: u64 = p32;
                if (block_lo == 0) {
                    cur[bi] = if (p <= sqrt_y) y else a - 1;
                    next[bi] = p;
                    wpos[bi] = 0;
                } else {
                    const mmax: u64 = @min(y, @as(u64, @intCast(x / (@as(X, p) * @as(X, block_lo)))));
                    cur[bi] = if (p <= sqrt_y)
                        mmax
                    else if (pi_tab[@intCast(mmax)] > 0)
                        pi_tab[@intCast(mmax)] - 1
                    else
                        0;
                    var m0 = (block_lo + p - 1) / p; // ceil(block_lo/p)
                    if (m0 == 0) m0 = 1;
                    while (!COP30[@intCast(m0 % 30)]) m0 += 1;
                    next[bi] = p * m0;
                    wpos[bi] = W30IDX[@intCast(m0 % 30)];
                }
            }
            @memset(phi_run, 0);
            @memset(mu_sum, 0);

            var s2: i128 = 0;
            var sum_pi_xp: i128 = 0;
            var np: usize = 0;

            var lo: u64 = block_lo;
            while (lo < block_hi) : (lo += segw) {
                const hi = @min(lo + segw, block_hi);
                const len: usize = @intCast(hi - lo);

                self.ctr.reset(len);
                for (self.primes, 0..) |p32, bi| {
                    const p: u64 = p32;
                    if (bi >= 3) seg_cnt[bi] = self.ctr.total;

                    if (p <= sqrt_y) {
                        var m = cur[bi];
                        const mlo = y / p;
                        while (m > mlo) {
                            const v: u64 = @intCast(x / @as(X, m * p));
                            if (v >= hi) break;
                            if (INST) self.walk += 1;
                            if (v >= lo) {
                                const mm = self.mu[@intCast(m)];
                                if (mm != 0 and self.lpf[@intCast(m)] > p) {
                                    if (INST and v <= y and p * p > v) self.easy += 1;
                                    const phi_v = leafPhiP(v, p, bi, lo, y, pi_tab, &self.ctr, phi_run, -@as(i64, mm), mu_sum);
                                    s2 += @as(i128, -mm) * @as(i128, phi_v);
                                    if (INST) self.leaves += 1;
                                }
                            }
                            m -= 1;
                        }
                        cur[bi] = m;
                    } else if (p < p_cube_min) {
                        var qi = cur[bi];
                        while (qi > bi) {
                            const q: u64 = self.primes[@intCast(qi)];
                            const v: u64 = @intCast(x / @as(X, p * q));
                            if (v >= hi) break;
                            if (INST) self.walk += 1;
                            if (v >= lo) {
                                if (INST and v <= y and p * p > v) self.easy += 1;
                                s2 += @as(i128, leafPhiP(v, p, bi, lo, y, pi_tab, &self.ctr, phi_run, 1, mu_sum));
                                if (INST) self.leaves += 1;
                            }
                            qi -= 1;
                        }
                        cur[bi] = qi;
                    }

                    if (bi >= 3) {
                        var j = next[bi];
                        var wp = wpos[bi];
                        while (j < hi) {
                            self.ctr.kill(@intCast(j - lo));
                            j += p * W30GAP[wp];
                            wp = (wp + 1) & 7;
                        }
                        next[bi] = j;
                        wpos[bi] = wp;
                    }
                }

                seg_cnt[a] = self.ctr.total;
                if (do_p2) {
                    const p_hi = if (lo == 0) sqrt_x else @min(@as(u64, @intCast(x / @as(X, lo))), sqrt_x);
                    const p_lo = @max(@as(u64, @intCast(x / @as(X, hi))), y);
                    if (p_hi > p_lo) {
                        const w: usize = @intCast(p_hi - p_lo);
                        @memset(pbuf[0..w], true);
                        for (bp4) |q| {
                            if (q * q > p_hi) break;
                            var j = @max(q * q, ((p_lo + q) / q) * q);
                            while (j <= p_hi) : (j += q) pbuf[@intCast(j - p_lo - 1)] = false;
                        }
                        var k: usize = w;
                        while (k > 0) {
                            k -= 1;
                            if (!pbuf[k]) continue;
                            const v: u64 = @intCast(x / @as(X, p_lo + 1 + @as(u64, k)));
                            const phi_va = phi_run[a] + self.ctr.prefix(@intCast(v - lo));
                            sum_pi_xp += phi_va - 1 + @as(i64, @intCast(a));
                            np += 1;
                        }
                    }
                }
                for (3..@max(3, a + 1)) |bi| phi_run[bi] += seg_cnt[bi];
            }

            self.blk_s2[blk] = s2;
            self.blk_pi[blk] = sum_pi_xp;
            self.blk_np[blk] = np;
            for (0..a + 1) |bi| {
                self.blk_total[blk * (a + 1) + bi] = phi_run[bi];
                self.blk_mu[blk * (a + 1) + bi] = mu_sum[bi];
            }
        }
    };
}

/// Stitch the per-block outputs into (S2, P₂). prefix[bi] runs Σ of lower blocks'
/// block_total[bi]; each block's local φ is corrected by prefix, and P₂ by prefix[a].
fn reduceBlocks(a: usize, nb: usize, do_p2: bool, blk_s2: []const i128, blk_pi: []const i128, blk_np: []const u64, blk_total: []const i64, blk_mu: []const i64, s1_closed: i128, gpa: std.mem.Allocator) !struct { s2: i128, p2: i128 } {
    const prefix = try gpa.alloc(i64, a + 1);
    defer gpa.free(prefix);
    @memset(prefix, 0);
    var s2: i128 = s1_closed;
    var sum_pi_xp: i128 = 0;
    var np: u64 = 0;
    for (0..nb) |blk| {
        s2 += blk_s2[blk];
        sum_pi_xp += blk_pi[blk];
        var corr: i128 = 0;
        for (3..@max(3, a + 1)) |bi| corr += @as(i128, prefix[bi]) * @as(i128, blk_mu[blk * (a + 1) + bi]);
        s2 += corr;
        sum_pi_xp += @as(i128, prefix[a]) * @as(i128, @intCast(blk_np[blk]));
        np += blk_np[blk];
        for (3..@max(3, a + 1)) |bi| prefix[bi] += blk_total[blk * (a + 1) + bi];
    }
    var p2: i128 = 0;
    if (do_p2 and np > 0) {
        const Ai: i128 = @intCast(@as(u64, a) + np);
        const ai: i128 = @intCast(a);
        p2 = sum_pi_xp - (@divExact((Ai - 1) * Ai, 2) - @divExact(ai * (ai - 1), 2));
    }
    return .{ .s2 = s2, .p2 = p2 };
}

/// The fused sweep, split into `nb` contiguous blocks of [1, z] run SERIALLY here —
/// the block decomposition + cross-block φ reduction that parallelism needs, proven
/// correct against the single-block path before any threads are involved. Each block
/// sweeps its sub-range with a LOCAL running φ (phi_run init 0), records its boundary
/// totals and mu_sum, and the reduction stitches them: φ(v,b−1) in block t is
/// local + Σ_{t'<t} block_total[t'][b−1]. nb = 1 reproduces the monolithic sweep.
pub fn s2P2Blocks(comptime C: type, comptime INST: bool, comptime X: type, gpa: std.mem.Allocator, x: X, y: u64, seg: usize, nb: usize) !FusedResult {
    // Segments must start at a multiple of 960 = lcm(30, 64) for MASK30 to be
    // word-aligned, so round the segment length up to one.
    const segw = @max(@as(usize, 960), ((seg + 959) / 960) * 960);

    var t = try SmallTables.init(gpa, y);
    defer t.deinit(gpa);
    const z: u64 = @intCast(x / @as(X, y));
    const a = t.a;
    const sqrt_x = isqrtG(X, x);
    const sqrt_y = common.isqrt(y);

    // P₂ needs the primes in (y, √x] descending. We do NOT store them: for segment
    // [lo, hi) the p with x/p ∈ [lo, hi) are exactly p ∈ (⌊x/hi⌋, ⌊x/lo⌋], and as lo
    // sweeps upward those ranges TILE (y, √x] — disjoint, contiguous, each ≤ seg wide
    // — so each is sieved on the fly from base primes ≤ x^(1/4). Storing π(√x) primes
    // was Θ(√x/ln x): 406 MB at 10^18. This is ~27 KB.
    const bp4 = try rs.basePrimes(gpa, @max(common.isqrt(sqrt_x), 2));
    defer gpa.free(bp4);
    const pbuf = try gpa.alloc(bool, segw); // p-range sieve; width ≤ seg ≤ segw
    defer gpa.free(pbuf);
    const do_p2 = y < sqrt_x; // np (primes in (y, √x]) is now accumulated per block

    var ctr = try C.init(gpa, segw);
    defer ctr.deinit(gpa);

    // a+1 slots: [a] is the fully-folded state φ(·, a), which is what P₂ reads.
    const phi_run = try gpa.alloc(i64, a + 1);
    defer gpa.free(phi_run);
    @memset(phi_run, 0);
    const seg_cnt = try gpa.alloc(i64, a + 1);
    defer gpa.free(seg_cnt);
    // LMO p.555 class (1) / DR §6.1 — THE BINOMIAL.
    // If p³ ≥ x then for every q > p we have p²q > p³ ≥ x, hence x/(pq) < p. The only
    // survivor ≤ x/(pq) coprime to p_1..p_{b−1} is 1 — a surviving prime would be
    // ≥ p_b = p > x/(pq), and a coprime composite is ≥ p², bigger still. So
    // φ(x/(pq), π(p)−1) = 1 IDENTICALLY: not a cheap lookup, literally 1. With m = q
    // prime, μ(q) = −1 so −μ(m) = +1, and the entire class is one pair count:
    //
    //     S₁ = C(n₁, 2),   n₁ = #{primes p : p³ ≥ x, p ≤ y}
    //
    // ~52% of all leaves at α=4 (20.9e9 of 40.1e9 at 10^18) — constant time, no sieve,
    // no queries. This is LMO Lemma 5.1's class (1) and DR's S₁.
    //
    // p*p*p would overflow (p ≤ y = 4x^(1/3) ⇒ p³ ≤ 64x = 6.4e20 at 10^19), so the
    // threshold comes from icbrt: the least p with p³ ≥ x.
    // π(m) for m ≤ y — O(y) build, O(1) query. The table only needs to span [1, y],
    // NOT [1, z]: every leaf it serves has v ≤ y by construction. Assuming otherwise
    // is what sank the earlier flat-prefix-array attempt.
    const pi_tab = try gpa.alloc(u32, @intCast(y + 1));
    defer gpa.free(pi_tab);
    {
        var c: u32 = 0;
        for (0..@as(usize, @intCast(y + 1))) |m| {
            if (m >= 2 and t.lpf[m] == m) c += 1; // lpf[m] == m ⇔ m prime
            pi_tab[m] = c;
        }
    }

    const c3 = icbrtG(X, x);
    const p_cube_min = if (@as(X, c3) * @as(X, c3) * @as(X, c3) == x) c3 else c3 + 1;
    var n1: u64 = 0;
    for (t.primes) |p32| {
        const p: u64 = p32;
        if (p > sqrt_y and p >= p_cube_min) n1 += 1; // exactly the set skipped below
    }
    const s1_closed: i128 = @divTrunc(@as(i128, @intCast(n1)) * @as(i128, @intCast(n1 -| 1)), 2);

    // Per-block scratch (fold cursors, leaf cursor, mu_sum) + per-block outputs.
    const next = try gpa.alloc(u64, a);
    defer gpa.free(next);
    const wpos = try gpa.alloc(u8, a);
    defer gpa.free(wpos);
    const cur = try gpa.alloc(u64, a);
    defer gpa.free(cur);
    const mu_sum = try gpa.alloc(i64, a + 1);
    defer gpa.free(mu_sum);
    const blk_s2 = try gpa.alloc(i128, nb);
    defer gpa.free(blk_s2);
    const blk_pi = try gpa.alloc(i128, nb);
    defer gpa.free(blk_pi);
    const blk_np = try gpa.alloc(u64, nb);
    defer gpa.free(blk_np);
    const blk_total = try gpa.alloc(i64, nb * (a + 1));
    defer gpa.free(blk_total);
    const blk_mu = try gpa.alloc(i64, nb * (a + 1));
    defer gpa.free(blk_mu);

    var ctx = BlkCtx(C, X){
        .primes = t.primes,
        .mu = t.mu,
        .lpf = t.lpf,
        .pi_tab = pi_tab,
        .bp4 = bp4,
        .x = x,
        .y = y,
        .z = z,
        .sqrt_x = sqrt_x,
        .sqrt_y = sqrt_y,
        .p_cube_min = p_cube_min,
        .segw = segw,
        .a = a,
        .nb = nb,
        .do_p2 = do_p2,
        .ctr = ctr,
        .phi_run = phi_run,
        .seg_cnt = seg_cnt,
        .cur = cur,
        .next = next,
        .wpos = wpos,
        .pbuf = pbuf,
        .mu_sum = mu_sum,
        .blk_s2 = blk_s2,
        .blk_pi = blk_pi,
        .blk_np = blk_np,
        .blk_total = blk_total,
        .blk_mu = blk_mu,
    };
    // Serial: one block at a time. The threaded driver spawns workers that each call
    // ctx.runBlock (with their own scratch) over dispensed blocks — identical results.
    for (0..nb) |blk| ctx.runBlock(INST, blk);

    const r = try reduceBlocks(a, nb, do_p2, blk_s2, blk_pi, blk_np, blk_total, blk_mu, s1_closed, gpa);
    return .{
        .s2 = r.s2,
        .p2 = r.p2,
        .leaves = if (INST) ctx.leaves + @as(u64, @intCast(s1_closed)) else 0,
        .easy = ctx.easy,
        .walk = ctx.walk,
        .z = z,
        .a = a,
    };
}

/// Threaded fused sweep: `nthreads` workers pull blocks from a mutex-guarded counter
/// and each run BlkCtx.runBlock with their OWN scratch, sharing the read-only tables
/// (built once) and the per-block output arrays (disjoint writes). INST is always false
/// — diagnostics are a serial analysis path. Result is identical to s2P2Blocks(nb).
pub fn s2P2Parallel(comptime C: type, comptime X: type, gpa: std.mem.Allocator, x: X, y: u64, seg: usize, nb: usize, nthreads: usize, blk_ns: ?[]u64) !FusedResult {
    const segw = @max(@as(usize, 960), ((seg + 959) / 960) * 960);
    var t = try SmallTables.init(gpa, y);
    defer t.deinit(gpa);
    const z: u64 = @intCast(x / @as(X, y));
    const a = t.a;
    const sqrt_x = isqrtG(X, x);
    const sqrt_y = common.isqrt(y);
    const bp4 = try rs.basePrimes(gpa, @max(common.isqrt(sqrt_x), 2));
    defer gpa.free(bp4);
    const do_p2 = y < sqrt_x;

    const pi_tab = try gpa.alloc(u32, @intCast(y + 1));
    defer gpa.free(pi_tab);
    {
        var c: u32 = 0;
        for (0..@as(usize, @intCast(y + 1))) |m| {
            if (m >= 2 and t.lpf[m] == m) c += 1;
            pi_tab[m] = c;
        }
    }
    const c3 = icbrtG(X, x);
    const p_cube_min = if (@as(X, c3) * @as(X, c3) * @as(X, c3) == x) c3 else c3 + 1;
    var n1: u64 = 0;
    for (t.primes) |p32| {
        const p: u64 = p32;
        if (p > sqrt_y and p >= p_cube_min) n1 += 1;
    }
    const s1_closed: i128 = @divTrunc(@as(i128, @intCast(n1)) * @as(i128, @intCast(n1 -| 1)), 2);

    // shared per-block outputs
    const blk_s2 = try gpa.alloc(i128, nb);
    defer gpa.free(blk_s2);
    const blk_pi = try gpa.alloc(i128, nb);
    defer gpa.free(blk_pi);
    const blk_np = try gpa.alloc(u64, nb);
    defer gpa.free(blk_np);
    const blk_total = try gpa.alloc(i64, nb * (a + 1));
    defer gpa.free(blk_total);
    const blk_mu = try gpa.alloc(i64, nb * (a + 1));
    defer gpa.free(blk_mu);

    const Ctx = BlkCtx(C, X);
    const ctxs = try gpa.alloc(Ctx, nthreads);
    defer gpa.free(ctxs);
    var allocd: usize = 0;
    defer for (0..allocd) |i| {
        ctxs[i].ctr.deinit(gpa);
        gpa.free(ctxs[i].phi_run);
        gpa.free(ctxs[i].seg_cnt);
        gpa.free(ctxs[i].cur);
        gpa.free(ctxs[i].next);
        gpa.free(ctxs[i].wpos);
        gpa.free(ctxs[i].pbuf);
        gpa.free(ctxs[i].mu_sum);
    };
    for (0..nthreads) |i| {
        ctxs[i] = Ctx{
            .primes = t.primes,
            .mu = t.mu,
            .lpf = t.lpf,
            .pi_tab = pi_tab,
            .bp4 = bp4,
            .x = x,
            .y = y,
            .z = z,
            .sqrt_x = sqrt_x,
            .sqrt_y = sqrt_y,
            .p_cube_min = p_cube_min,
            .segw = segw,
            .a = a,
            .nb = nb,
            .do_p2 = do_p2,
            .ctr = try C.init(gpa, segw),
            .phi_run = try gpa.alloc(i64, a + 1),
            .seg_cnt = try gpa.alloc(i64, a + 1),
            .cur = try gpa.alloc(u64, a),
            .next = try gpa.alloc(u64, a),
            .wpos = try gpa.alloc(u8, a),
            .pbuf = try gpa.alloc(bool, segw),
            .mu_sum = try gpa.alloc(i64, a + 1),
            .blk_s2 = blk_s2,
            .blk_pi = blk_pi,
            .blk_np = blk_np,
            .blk_total = blk_total,
            .blk_mu = blk_mu,
        };
        allocd = i + 1;
    }

    // Atomic-counter dispenser. std.Thread.Mutex was removed in Zig 0.16's Io rework;
    // an atomic fetch-add is simpler and the contention is O(0) anyway (a block is
    // seconds of work, a fetch-add is nanoseconds). Overshoot past nb is harmless.
    const Disp = struct { next: std.atomic.Value(usize), nb: usize, blk_ns: ?[]u64 };
    var disp = Disp{ .next = std.atomic.Value(usize).init(0), .nb = nb, .blk_ns = blk_ns };
    const Worker = struct {
        fn run(ctxp: *Ctx, dp: *Disp) void {
            while (true) {
                const blk = dp.next.fetchAdd(1, .monotonic);
                if (blk >= dp.nb) break;
                if (dp.blk_ns) |ns| {
                    const t0 = common.nowNs();
                    ctxp.runBlock(false, blk);
                    ns[blk] = common.nowNs() - t0;
                } else {
                    ctxp.runBlock(false, blk);
                }
            }
        }
    };

    // Spawn nthreads-1 workers; the master runs one too. Spawn failure degrades
    // gracefully — fewer workers still drain every block.
    const threads = try gpa.alloc(std.Thread, nthreads);
    defer gpa.free(threads);
    var spawned: usize = 0;
    for (1..nthreads) |i| {
        threads[i] = std.Thread.spawn(.{}, Worker.run, .{ &ctxs[i], &disp }) catch break;
        spawned = i;
    }
    Worker.run(&ctxs[0], &disp);
    var j: usize = 1;
    while (j <= spawned) : (j += 1) threads[j].join();

    const r = try reduceBlocks(a, nb, do_p2, blk_s2, blk_pi, blk_np, blk_total, blk_mu, s1_closed, gpa);
    return .{ .s2 = r.s2, .p2 = r.p2, .leaves = 0, .easy = 0, .walk = 0, .z = z, .a = a };
}

// ---------------------------------------------------------------------- π(x)

pub const PiResult = struct { pi: u64, phi: i128, p2: i128, y: u64, a: usize, z: u64, leaves: u64 };

/// π(x) = φ(x, a) + a − 1 − P₂(x, a), a = π(y) — LMO end to end.
/// π(x) end to end. Takes u128 so x can exceed 2⁶⁴ (π(10²⁰)); dispatches to the u64
/// value type below the wall so the common case pays no u128-division tax. Everything
/// but x and its divisions stays u64 — see s2AndP2FusedGen.
pub fn piLMO(gpa: std.mem.Allocator, x: u128, y_in: ?u64, seg_in: ?usize) !PiResult {
    if (x < 2) return .{ .pi = 0, .phi = 0, .p2 = 0, .y = 0, .a = 0, .z = 0, .leaves = 0 };
    const y = y_in orelse defaultY(x);
    const seg: usize = seg_in orelse @max(y, 1024);
    if (x <= std.math.maxInt(u64)) {
        return piImpl(u64, gpa, @intCast(x), y, seg);
    } else {
        return piImpl(u128, gpa, x, y, seg);
    }
}

/// Parallel π(x): nb = nthreads·k_over blocks over the same decomposition, dispensed
/// dynamically to nthreads workers. k_over is the over-partition factor for load
/// balance (blocks are equal-segment-width for now — cost-balanced partition is next).
pub fn piLMOPar(gpa: std.mem.Allocator, x: u128, nthreads: usize, k_over: usize) !PiResult {
    if (x < 2) return .{ .pi = 0, .phi = 0, .p2 = 0, .y = 0, .a = 0, .z = 0, .leaves = 0 };
    const y = defaultY(x);
    const seg: usize = @max(y, 1024);
    const nb = @max(@as(usize, 1), nthreads * k_over);
    if (x <= std.math.maxInt(u64)) {
        return piImplPar(u64, gpa, @intCast(x), y, seg, nb, nthreads);
    } else {
        return piImplPar(u128, gpa, x, y, seg, nb, nthreads);
    }
}

fn piImplPar(comptime X: type, gpa: std.mem.Allocator, x: X, y: u64, seg: usize, nb: usize, nthreads: usize) !PiResult {
    const s1 = try ordinaryS1Gen(X, gpa, x, y);
    const f = try s2P2Parallel(Counter3P, X, gpa, x, y, seg, nb, nthreads, null);
    const phi = s1.s1 + f.s2;
    const r = phi + @as(i128, @intCast(f.a)) - 1 - f.p2;
    return .{ .pi = @intCast(r), .phi = phi, .p2 = f.p2, .y = y, .a = f.a, .z = f.z, .leaves = 0 };
}

fn piImpl(comptime X: type, gpa: std.mem.Allocator, x: X, y: u64, seg: usize) !PiResult {
    const s1 = try ordinaryS1Gen(X, gpa, x, y);
    const f = try s2AndP2FusedGen(Counter3P, false, X, gpa, x, y, seg); // S2 and P₂ share one sweep of [1,z]
    const phi = s1.s1 + f.s2;
    const r = phi + @as(i128, @intCast(f.a)) - 1 - f.p2;
    return .{ .pi = @intCast(r), .phi = phi, .p2 = f.p2, .y = y, .a = f.a, .z = f.z, .leaves = f.leaves };
}

/// icbrt, exported for analysis harnesses.
pub fn icbrtPub(x: u64) u64 {
    return icbrt(x);
}
