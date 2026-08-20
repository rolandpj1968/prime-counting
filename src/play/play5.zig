//! RPJ playground
// zig run --dep rs -Mroot=play5.zig -Mrs=../rangesieve.zig -- 100

const std = @import("std");
const assert = std.debug.assert;
const rs = @import("rs");

/// floor(sqrt(n)), exact for u64.
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

fn to_f(i: u64) f64 {
    return @floatFromInt(i);
}

pub fn main(init: std.process.Init) !void {
    var args = try init.minimal.args.iterateAllocator(init.gpa);
    _ = args.skip();
    const x_str = args.next() orelse {
        std.debug.print("Usage: zig run play5.zig -- <x>\n", .{});
        return;
    };
    const x = std.fmt.parseInt(u64, x_str, 10) catch |err| {
        std.debug.print("Error: '{s}' is not a valid u64 number (Error: {any})\n", .{ x_str, err });
        return;
    };

    const y = isqrt(x);
    const z = icbrt(x);

    const gpa = std.heap.page_allocator;

    const primes2 = try rs.basePrimes(gpa, y);
    defer gpa.free(primes2);

    const primes = try rs.basePrimes(gpa, z);
    defer gpa.free(primes);

    std.debug.print("x: {d:>16} | y = x^1/2: {d:>8} | pi(y): {d:>8} | z = x^1/3: {d:>8} | pi(z): {d:>8}\n\n", .{ x, y, primes2.len, z, primes.len });

    var i: usize = 0;
    while (i < 5 and i < primes2.len) {
        const p_a = primes2[primes2.len - 1 - i];
        std.debug.print("    x / {d:>8} : {d:>8}\n", .{ p_a, x / p_a });

        var j: usize = 0;
        while (j < 1 and j < primes.len) {
            const p_b = primes[primes.len - 1 - j];
            std.debug.print("    x / {d:>8} / {d:>8} : {d:>8}\n", .{ p_a, p_b, x / p_a / p_b });

            j += 1;
        }

        i += 1;
    }

    {
        i = 0;
        const p_b = primes[primes.len - 1];
        while (i < primes2.len and x / primes2[primes2.len - 1 - i] / p_b < p_b) {
            i += 1;
        }
        std.debug.print("\n\n{d:>8} out of {d:>8} p <= x^1/2 have 2nd order phi = 1 --- {d:>6.2}%\n", .{ i, primes2.len, to_f(i) / to_f(primes2.len) * 100.0 });

        const lim = i + 5;
        while (i < lim and i < primes2.len) {
            const p_a = primes2[primes2.len - 1 - i];
            std.debug.print("    x / {d:>8} / {d:>8} : {d:>8}\n", .{ p_a, p_b, x / p_a / p_b });

            i += 1;
        }
    }

    std.debug.print("\n", .{});

    var b: usize = primes.len - 1;
    const b_lim = if (primes.len < 21) 1 else primes.len - 20;
    while (true) {
        const p_b = primes[b];
        const p_bm1 = primes[b - 1];
        const total = primes2.len - b;

        var a = primes2.len - 1;

        while (b < a) {
            const p_a = primes2[a];

            assert(x / p_a / p_b == x / p_b / p_a);
            if (x / p_a / p_b > p_bm1) {
                break;
            }

            a -= 1;
        }
        const non_1s = a - b;

        while (b < a) {
            const p_a = primes2[a];

            if (x / p_a / p_b >= p_bm1 * p_bm1) {
                break;
            }

            a -= 1;
        }
        const non_pis = a - b;

        std.debug.print("    b: {d:>8} | p_b: {d:>8} | total p_a's {d:>8} | non-1's {d:>8} - {d:>6.2}% | non-pi's {d:>8} - {d:>6.2}%\n", .{ b, p_b, total, non_1s, to_f(non_1s) / to_f(total) * 100.0, non_pis, to_f(non_pis) / to_f(total) * 100.0 });

        if (b == 0 or b <= b_lim) {
            break;
        }
        b -= 1;
    }
}
