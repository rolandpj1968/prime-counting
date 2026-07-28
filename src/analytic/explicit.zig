//! Naive Riemann–von Mangoldt explicit formula for psi(x), no smoothing.
//! Referee: direct sieve sum of Lambda(n). Watch the truncated zero-sum
//! error empirically at checkpoint zero-counts.
//!
//! psi(x) = x - sum_rho x^rho/rho - ln 2pi - (1/2) ln(1 - x^-2)
//! pair(gamma) = 2 sqrt(x) (cos(gL)/2 + g sin(gL)) / (1/4 + g^2),  L = ln x

const std = @import("std");
// zig build-exe -O ReleaseFast -mcpu=native --dep rs -Mroot=explicit.zig -Mrs=../rangesieve.zig -femit-bin=explicit
const rs = @import("rs");

fn nowNs() u64 {
    var ts: std.os.linux.timespec = undefined;
    _ = std.os.linux.clock_gettime(.MONOTONIC, &ts);
    return @as(u64, @intCast(ts.sec)) * 1_000_000_000 + @as(u64, @intCast(ts.nsec));
}

const Kahan = struct {
    s: f64 = 0,
    c: f64 = 0,
    fn add(k: *Kahan, v: f64) void {
        const t = k.s + v;
        if (@abs(k.s) >= @abs(v)) {
            k.c += (k.s - t) + v;
        } else {
            k.c += (v - t) + k.s;
        }
        k.s = t;
    }
    fn val(k: *const Kahan) f64 {
        return k.s + k.c;
    }
};

fn psiDirect(gpa: std.mem.Allocator, x: u64) !f64 {
    const primes = try rs.basePrimes(gpa, x);
    defer gpa.free(primes);
    var k = Kahan{};
    for (primes) |p| k.add(@log(@as(f64, @floatFromInt(p))));
    for (primes) |p| {
        if (p * p > x) break;
        var pk = p * p;
        while (pk <= x) : (pk *= p) k.add(@log(@as(f64, @floatFromInt(p))));
    }
    return k.val();
}

pub fn main(init: std.process.Init) !void {
    var args = try init.minimal.args.iterateAllocator(init.gpa);
    _ = args.skip();
    const usage = "Usage: explicit <x> <zeros-file>\n";
    const x_str = args.next() orelse return std.debug.print(usage, .{});
    const x = try std.fmt.parseInt(u64, x_str, 10);
    const zpath = args.next() orelse return std.debug.print(usage, .{});

    const gpa = std.heap.page_allocator;

    var tio = std.Io.Threaded.init(gpa, .{});
    const data = try std.Io.Dir.cwd().readFileAlloc(tio.io(), zpath, gpa, .limited(1 << 27));
    defer gpa.free(data);

    var zeros = try std.ArrayList(f64).initCapacity(gpa, 1 << 21);
    defer zeros.deinit(gpa);
    var lines = std.mem.tokenizeScalar(u8, data, '\n');
    while (lines.next()) |line| {
        const s = std.mem.trim(u8, line, " \t\r");
        if (s.len == 0) continue;
        try zeros.append(gpa, try std.fmt.parseFloat(f64, s));
    }

    const xf: f64 = @floatFromInt(x);
    const L = @log(xf);
    const sx = @sqrt(xf);
    const ln2pi = @log(2.0 * std.math.pi);
    const tail = 0.5 * @log(1.0 - 1.0 / (xf * xf));

    const t0 = nowNs();
    const psi_true = try psiDirect(gpa, x);
    const t_sieve = @as(f64, @floatFromInt(nowNs() - t0)) / 1e9;
    std.debug.print("x = {d}  zeros = {d} (T = {d:.3})\n", .{ x, zeros.items.len, zeros.items[zeros.items.len - 1] });
    std.debug.print("psi_direct = {d:.6}   ({d:.1}s sieve)\n\n", .{ psi_true, t_sieve });
    std.debug.print("{s:>9} {s:>12} {s:>20} {s:>14} {s:>12}\n", .{ "N", "T", "psi_T", "err", "rms_pred" });

    var k = Kahan{};
    var next_cp: usize = 100;
    for (zeros.items, 0..) |g, i| {
        const th = g * L;
        k.add(2.0 * sx * (0.5 * @cos(th) + g * @sin(th)) / (0.25 + g * g));
        const n = i + 1;
        if (n == next_cp or n == zeros.items.len) {
            const psi_t = xf - k.val() - ln2pi - tail;
            // rms of the dropped tail: sqrt(x) * sqrt(sum_{gamma>T} 1/gamma^2),
            // density ln(t/2pi)/2pi  =>  integral ~ (ln(T/2pi)+1)/(2 pi T)
            const pred = sx * @sqrt((@log(g / (2.0 * std.math.pi)) + 1.0) / (2.0 * std.math.pi * g));
            std.debug.print("{d:>9} {d:>12.1} {d:>20.3} {d:>14.3} {d:>12.3}\n", .{ n, g, psi_t, psi_t - psi_true, pred });
            next_cp *= 10;
        }
    }
}
