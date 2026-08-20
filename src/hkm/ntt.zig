//! Number-theoretic transform over the Goldilocks field, for the HKM rungs.
//!
//! Why a field and not floats: the whole point of HKM over the analytic
//! method is that every quantity is an exact integer, so there is no error
//! budget to certify. An FFT in f64 would hand that back for nothing. The
//! transform is exact or it is useless.
//!
//! p = 2^64 - 2^32 + 1 = 18446744069414584321. Two useful properties:
//!   * p - 1 = 2^32 * 3 * 5 * 17 * 257 * 65537, so there are roots of unity
//!     of every order up to 2^32 — far past any transform length we need.
//!   * 2^64 = 2^32 - 1 (mod p), so reduction of a 128-bit product is a few
//!     shifts, one 32x32 multiply and two carries. No division.
//!
//! Both the fast reduction and a u128-remainder reference are exported;
//! hkm3 checks one against the other rather than trusting the folklore.

const std = @import("std");

pub const P: u64 = 0xFFFF_FFFF_0000_0001;
const EPS: u64 = 0xFFFF_FFFF; // 2^64 mod P
const GEN: u64 = 7; // primitive root of the multiplicative group

pub fn add(a: u64, b: u64) u64 {
    const r = @addWithOverflow(a, b);
    var s = r[0];
    if (r[1] == 1) {
        s -%= P; // true sum is s + 2^64; (s + 2^64) - P == s -% P in u64
    } else if (s >= P) {
        s -= P;
    }
    return s;
}

pub fn sub(a: u64, b: u64) u64 {
    const r = @subWithOverflow(a, b);
    var s = r[0];
    if (r[1] == 1) s +%= P;
    return s;
}

/// Reference multiplication: obviously correct, uses a 128-bit division.
pub fn mulRef(a: u64, b: u64) u64 {
    return @intCast((@as(u128, a) * @as(u128, b)) % P);
}

/// Fast multiplication via the 2^64 = 2^32 - 1 folding.
pub fn mul(a: u64, b: u64) u64 {
    const x: u128 = @as(u128, a) * @as(u128, b);
    const x_lo: u64 = @truncate(x);
    const x_hi: u64 = @truncate(x >> 64);
    const hh: u64 = x_hi >> 32;
    const hl: u64 = x_hi & 0xFFFF_FFFF;
    const s0 = @subWithOverflow(x_lo, hh);
    var t0 = s0[0];
    if (s0[1] == 1) t0 -%= EPS;
    const t1 = hl *% EPS;
    const s1 = @addWithOverflow(t0, t1);
    var res = s1[0];
    if (s1[1] == 1) res +%= EPS;
    if (res >= P) res -= P;
    return res;
}

pub fn powmod(a: u64, e: u64) u64 {
    var r: u64 = 1;
    var b = a;
    var k = e;
    while (k != 0) : (k >>= 1) {
        if (k & 1 == 1) r = mul(r, b);
        b = mul(b, b);
    }
    return r;
}

pub fn inv(a: u64) u64 {
    return powmod(a, P - 2);
}

/// Field image of a signed integer, and back. Values are recovered in
/// (-P/2, P/2]; anything outside that range has wrapped and is unrecoverable,
/// which is exactly the condition hkm3 measures headroom against.
pub fn fromSigned(v: i64) u64 {
    if (v >= 0) return @as(u64, @intCast(v)) % P;
    return P - (@as(u64, @intCast(-v)) % P);
}

pub fn toSigned(v: u64) i64 {
    if (v > P / 2) return -@as(i64, @intCast(P - v));
    return @intCast(v);
}

fn bitrev(a: []u64) void {
    const n = a.len;
    var j: usize = 0;
    var i: usize = 1;
    while (i < n) : (i += 1) {
        var bit = n >> 1;
        while (j & bit != 0) : (bit >>= 1) j ^= bit;
        j |= bit;
        if (i < j) std.mem.swap(u64, &a[i], &a[j]);
    }
}

/// In-place iterative Cooley-Tukey. `a.len` must be a power of two dividing
/// 2^32 (which everything below 4 billion cells does).
pub fn transform(a: []u64, inverse: bool) void {
    const n = a.len;
    if (n <= 1) return;
    std.debug.assert(std.math.isPowerOfTwo(n));
    bitrev(a);
    var len: usize = 2;
    while (len <= n) : (len <<= 1) {
        var w = powmod(GEN, (P - 1) / @as(u64, len));
        if (inverse) w = inv(w);
        const half = len >> 1;
        var i: usize = 0;
        while (i < n) : (i += len) {
            var wn: u64 = 1;
            var j: usize = 0;
            while (j < half) : (j += 1) {
                const u = a[i + j];
                const v = mul(a[i + j + half], wn);
                a[i + j] = add(u, v);
                a[i + j + half] = sub(u, v);
                wn = mul(wn, w);
            }
        }
    }
    if (inverse) {
        const ninv = inv(@intCast(n));
        for (a) |*x| x.* = mul(x.*, ninv);
    }
}

/// Linear convolution of a and b, truncated to degree `maxdeg` inclusive.
/// The transform length is chosen to cover the FULL product degree, so no
/// term wraps around into a low coefficient — truncating afterwards is a
/// discard, not a fold. (Sizing to maxdeg+1 instead would silently alias the
/// discarded terms back into the answer.)
pub fn convolve(gpa: std.mem.Allocator, a: []const u64, b: []const u64, maxdeg: usize) ![]u64 {
    const full = a.len + b.len - 1;
    var L: usize = 1;
    while (L < full) L <<= 1;
    const fa = try gpa.alloc(u64, L);
    defer gpa.free(fa);
    const fb = try gpa.alloc(u64, L);
    defer gpa.free(fb);
    @memset(fa, 0);
    @memset(fb, 0);
    @memcpy(fa[0..a.len], a);
    @memcpy(fb[0..b.len], b);
    transform(fa, false);
    transform(fb, false);
    for (fa, fb) |*x, y| x.* = mul(x.*, y);
    transform(fa, true);
    const out_len = @min(maxdeg + 1, full);
    const out = try gpa.alloc(u64, out_len);
    @memcpy(out, fa[0..out_len]);
    return out;
}
