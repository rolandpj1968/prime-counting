//! HKM rung 7a — Lemma 21 (§4.1, "using fewer exact primes"), refereed.
//!
//! Everything so far is the t = 2 case: sieve with every prime up to sqrt N.
//! §4.1 generalises to sieving only up to N^(1/t) and recovering the rest
//! analytically, which is what buys the Otilde(N^(1/4))-space corner of their
//! trade-off curve — a different resource profile rather than a
//! constant-factor race, and the only remaining item that changes what the
//! method IS.
//!
//! Lemma 21. Let f = (1 * mu_{<= N^(1/t)}) - delta_1, i.e. the indicator of
//! integers n > 1 all of whose prime factors exceed N^(1/t). Then
//!
//!   pi(N) = sum_{k=1}^{t-1} (-1)^(k-1)/k * T_k
//!           + pi(N^(1/t))
//!           - sum_{k=2}^{t-1} (1/k) * ( pi(N^(1/k)) - pi(N^(1/t)) )
//!
//! with T_k = sum_{n<=N} f^{*k}(n), the number of ORDERED k-tuples of rough
//! integers > 1 whose product is at most N. The k-sum stops at t-1 because
//! t rough factors already exceed N.
//!
//! This rung builds none of the segmented machinery. It computes every T_k
//! and every pi(N^(1/k)) exactly by brute force and checks the identity, for
//! a range of t — the same thing R0 did for Lemma 1, and for the same reason:
//! the identity is where a sign or an off-by-one hides, and it is far cheaper
//! to find here than inside an NTT.
//!
//! Sanity, by hand, at t = 3: T_1 = pi(N) - pi(N^(1/3)) + P2 where P2 counts
//! p*q <= N with N^(1/3) < p <= q; and both entries of a T_2 pair must be
//! PRIME, since a rough composite exceeds N^(2/3) and its partner exceeds
//! N^(1/3). So T_2 = 2*P2 - (pi(sqrt N) - pi(N^(1/3))), the P2 terms cancel,
//! and what is left is exactly pi(N).
//!
//! The 1/k coefficients make the identity rational, so it is evaluated over a
//! common denominator lcm(1..t) in i128 and the exactness of the division is
//! asserted, not rounded.

const std = @import("std");
const rs = @import("rs");

fn iroot(n: u64, k: u6) u64 {
    if (k == 1) return n;
    var r: u64 = @intFromFloat(std.math.pow(f64, @floatFromInt(n), 1.0 / @as(f64, @floatFromInt(k))));
    while (r > 1 and powLe(r, k, n) == false) r -= 1;
    while (powLe(r + 1, k, n)) r += 1;
    return r;
}

/// r^k <= n, total by construction (no overflow, no saturation)
fn powLe(r: u64, k: u6, n: u64) bool {
    var acc: u64 = 1;
    var i: u6 = 0;
    while (i < k) : (i += 1) {
        if (acc > n / r) return false;
        acc *= r;
    }
    return acc <= n;
}

fn gcd(a: u64, b: u64) u64 {
    return if (b == 0) a else gcd(b, a % b);
}

/// Number of ordered k-tuples from `rough` with product <= limit.
/// pre[v] = #{rough r : r <= v}.
fn tuples(k: usize, limit: u64, rough: []const u64, pre: []const u32) u128 {
    if (k == 1) return pre[@intCast(limit)];
    var total: u128 = 0;
    for (rough) |a| {
        if (a > limit) break;
        total += tuples(k - 1, limit / a, rough, pre);
    }
    return total;
}

pub fn main(init: std.process.Init) !void {
    var args = try init.minimal.args.iterateAllocator(init.gpa);
    _ = args.skip();
    const gpa = std.heap.page_allocator;

    const usage = "Usage: hkm6 <N> [tmax]\n";
    const n_str = args.next() orelse return std.debug.print(usage, .{});
    const N = try std.fmt.parseInt(u64, n_str, 10);
    var tmax: u6 = 5;
    if (args.next()) |a| tmax = try std.fmt.parseInt(u6, a, 10);

    const base = try rs.basePrimes(gpa, std.math.sqrt(N) + 1);
    defer gpa.free(base);
    const pi_N = try rs.countInRange(u64, gpa, 2, N + 1, base);
    std.debug.print("N = {d}   pi(N) = {d}\n\n", .{ N, pi_N });
    std.debug.print("{s:>3} {s:>10} {s:>9} {s:>22} {s:>14} {s:>8}\n", .{ "t", "N^(1/t)", "|rough|", "T_1 .. T_(t-1)", "identity", "check" });

    var t: u6 = 2;
    while (t <= tmax) : (t += 1) {
        const y = iroot(N, t);
        if (y < 2) break;

        // rough = { n in (1, N] with no prime factor <= y }, plus prefix counts
        const sieve = try gpa.alloc(bool, @intCast(N + 1));
        defer gpa.free(sieve);
        @memset(sieve, true);
        sieve[0] = false;
        sieve[1] = false;
        for (base) |p| {
            if (p > y) break;
            var q: u64 = p;
            while (q <= N) : (q += p) sieve[@intCast(q)] = false;
        }
        var nrough: usize = 0;
        for (sieve) |b| nrough += @intFromBool(b);
        const rough = try gpa.alloc(u64, nrough);
        defer gpa.free(rough);
        const pre = try gpa.alloc(u32, @intCast(N + 1));
        defer gpa.free(pre);
        {
            var idx: usize = 0;
            var run: u32 = 0;
            for (sieve, 0..) |b, v| {
                if (b) {
                    rough[idx] = @intCast(v);
                    idx += 1;
                    run += 1;
                }
                pre[v] = run;
            }
        }

        // common denominator for the 1/k coefficients
        var D: u64 = 1;
        var k: u6 = 1;
        while (k <= t) : (k += 1) D = D / gcd(D, k) * k;

        var acc: i128 = 0;
        var buf: [64]u8 = undefined;
        var shown: usize = 0;
        var line: [128]u8 = undefined;
        var lw: usize = 0;
        k = 1;
        while (k <= t - 1) : (k += 1) {
            const T = tuples(k, N, rough, pre);
            const term = @as(i128, @intCast(D / k)) * @as(i128, @intCast(T));
            acc += if (k % 2 == 1) term else -term;
            if (shown < 3) {
                const printed = std.fmt.bufPrint(&buf, "{s}{d}", .{ if (shown > 0) "," else "", T }) catch buf[0..0];
                if (lw + printed.len < line.len) {
                    @memcpy(line[lw .. lw + printed.len], printed);
                    lw += printed.len;
                }
                shown += 1;
            }
        }
        const pi_y = try rs.countInRange(u64, gpa, 2, y + 1, base);
        acc += @as(i128, @intCast(D)) * @as(i128, @intCast(pi_y));
        k = 2;
        while (k <= t - 1) : (k += 1) {
            const pk = try rs.countInRange(u64, gpa, 2, iroot(N, k) + 1, base);
            acc -= @as(i128, @intCast(D / k)) * (@as(i128, @intCast(pk)) - @as(i128, @intCast(pi_y)));
        }

        const exact = @rem(acc, @as(i128, @intCast(D))) == 0;
        const val = @divTrunc(acc, @as(i128, @intCast(D)));
        std.debug.print("{d:>3} {d:>10} {d:>9} {s:>22} {d:>14} {s:>8}\n", .{ t, y, nrough, line[0..lw], val, if (exact and val == pi_N) "MATCH" else if (!exact) "NOT INTEGER" else "MISMATCH" });
        if (!exact or val != pi_N) return error.Lemma21Failed;
    }
}
