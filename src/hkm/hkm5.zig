//! HKM rung 5 — §3.2 (Newton pointwise in Fourier space, no truncation) and
//! §3.3 (partitioning primes to reduce the padding that costs).
//!
//! R4 implements §3.1: Newton in the time domain, truncating to degree K once
//! per r, which costs 2R round-trip transforms on top of the R for the E_m —
//! 3R = 36 at 1e13, and 16.1 s of the rung's 18.2 s.
//!
//! §3.2 removes the round trips. At a FIXED Fourier coordinate the recurrence
//!     r c_r = sum_{m=1..r} (-1)^(m-1) c_{r-m} e_m,  c_0 = 1
//! is a scalar recurrence in r alone: every coordinate is independent. So run
//! it pointwise, accumulate mu-hat-hat = sum_r (-1)^r c_r as you go, and take
//! ONE inverse transform at the end. R+1 transforms instead of 3R.
//! (HKM go further and observe x c'(x) = e(x) c(x), so c = exp(int e/t), which
//! computes the c_r in O(r_max log r_max) rather than O(r_max^2). At
//! r_max = 12 that is not worth a length-13 FFT, so the scalar recurrence is
//! kept — and note that c evaluated at x = -1 is exactly
//! mu-hat = exp(-sum_m E_m/m), the same identity along the other axis.)
//!
//! THE PRICE, AND THE PREDICTION. Truncation is what kept degrees bounded by
//! K. Without it, deg C_r <= r * max_p j_p, so the arrays must be padded to
//! R * max_j ~ 6K instead of 2K — a 4x longer transform to save a factor 3 in
//! transform count, and 4x as many coordinates paying the same O(R^2/2)
//! pointwise recurrence. Prediction before running: §3.2 on the whole prime
//! set LOSES to R4 at these N. The paper agrees in advance — it says only
//! that this "recovered almost the same complexity", and that its value is in
//! enabling §3.3.
//!
//! §3.3 claws the padding back by running §3.2 separately on each prime
//! interval [N^(1/2^(m+1)), N^(1/2^m)] and convolving the pieces at the end.
//! On such an interval at most ~2^(m+1) primes can multiply to <= N while
//! max_j is only K/2^m, so the required padding is ~2K for EVERY interval,
//! independent of m. Also r_max is capped by how many primes the interval
//! actually holds, which for the small-prime intervals is tiny.
//!
//! --parts 1 forces a single interval, i.e. §3.2 as written; the default
//! partitions, i.e. §3.3. Same code path, one knob.

const std = @import("std");
const rs = @import("rs");
const sg = @import("seg.zig");
const ntt = @import("ntt.zig");
const muhat = @import("muhat.zig");

fn nowNs() u64 {
    var ts: std.os.linux.timespec = undefined;
    _ = std.os.linux.clock_gettime(.MONOTONIC, &ts);
    return @as(u64, @intCast(ts.sec)) * 1_000_000_000 + @as(u64, @intCast(ts.nsec));
}
fn sec(a: u64, b: u64) f64 {
    return @as(f64, @floatFromInt(b - a)) / 1e9;
}

fn muHatSparse(gpa: std.mem.Allocator, ks: []const usize, K: usize) ![]i64 {
    const a = try gpa.alloc(i64, K + 1);
    @memset(a, 0);
    a[0] = 1;
    for (ks) |j| {
        if (j > K) continue;
        var t: usize = K + 1;
        while (t > j) {
            t -= 1;
            a[t] -= a[t - j];
        }
    }
    return a;
}

pub fn main(init: std.process.Init) !void {
    var args = try init.minimal.args.iterateAllocator(init.gpa);
    _ = args.skip();
    const gpa = std.heap.page_allocator;

    const usage = "Usage: hkm5 <N> [--parts 1] [--no-referee]\n";
    const n_str = args.next() orelse return std.debug.print(usage, .{});
    const N = try std.fmt.parseInt(u64, n_str, 10);
    var one_part = false;
    var referee = true;
    while (args.next()) |a| {
        if (std.mem.eql(u8, a, "--parts")) {
            const v = args.next() orelse return std.debug.print(usage, .{});
            one_part = std.mem.eql(u8, v, "1");
        } else if (std.mem.eql(u8, a, "--no-referee")) referee = false;
    }

    const sq = std.math.sqrt(N);
    const B = @log2(@as(f64, @floatFromInt(N)));
    const delta = B / @as(f64, @floatFromInt(sq));
    const primes = try rs.basePrimes(gpa, sq);
    defer gpa.free(primes);
    var g = try sg.build(gpa, N, delta);
    defer g.deinit(gpa);
    const K = g.K;
    const ks = try gpa.alloc(usize, primes.len);
    defer gpa.free(ks);
    for (primes, 0..) |p, i| ks[i] = g.kbar(p);

    std.debug.print("N = {d}  K = {d}  pi(sqrt N) = {d}{s}\n", .{ N, K, primes.len, if (one_part) "  (section 3.2 as written)" else "  (section 3.3 partition)" });

    const t0 = nowNs();
    var stat = muhat.Stats{};
    const acc = try muhat.fourier(gpa, ks, K, one_part, &stat);
    defer gpa.free(acc);
    const t_done = nowNs();

    std.debug.print("parts {d}   transforms {d}   sum of r_max {d}   largest transform 2^{d}   pointwise modmuls {e:.2}\n", .{ stat.parts, stat.transforms, stat.rmax_sum, std.math.log2_int(usize, @max(2, stat.L_max)), @as(f64, @floatFromInt(stat.work)) });
    std.debug.print("time {d:.2}s\n", .{sec(t0, t_done)});

    if (referee) {
        const tr = nowNs();
        const want = try muhat.sparse(gpa, ks, K);
        defer gpa.free(want);
        const ts = nowNs();
        var bad: usize = 0;
        var firstbad: usize = 0;
        for (want, acc, 0..) |w, got, i| {
            if (w != got) {
                if (bad == 0) firstbad = i;
                bad += 1;
            }
        }
        std.debug.print("referee (sparse build, {d:.2}s): {d} differing cells   {s}\n", .{ sec(tr, ts), bad, if (bad == 0) "MATCH" else "MISMATCH" });
        if (bad != 0) {
            std.debug.print("  first at k={d}: want {d} got {d}\n", .{ firstbad, want[firstbad], acc[firstbad] });
            return error.Rung5Mismatch;
        }
    }
}
