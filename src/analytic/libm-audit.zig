//! Emit (function, argument-bits, result-bits) triples for every f64
//! transcendental the pistar certified path depends on, over the argument
//! ranges it actually uses. Bit patterns, not decimal text: the referee
//! reconstructs exact values, so no formatting loss enters the comparison.
//!
//! Referee: trigref.py (Python decimal, 80 digits, pure software — an
//! independent implementation, unlike R/numpy which re-call the same libm).
const std = @import("std");

var state: u64 = 0x2545F4914F6CDD1D;
fn nextU() u64 {
    state = state *% 6364136223846793005 +% 1442695040888963407;
    return state;
}
fn uni(lo: f64, hi: f64) f64 {
    const m: f64 = @floatFromInt(nextU() >> 11);
    return lo + (hi - lo) * (m / 9007199254740992.0);
}
fn emit(name: []const u8, arg: f64, res: f64) void {
    std.debug.print("{s} {x} {x}\n", .{ name, @as(u64, @bitCast(arg)), @as(u64, @bitCast(res)) });
}

pub fn main() void {
    const N = 3000;
    var i: usize = 0;
    while (i < N) : (i += 1) {
        // dd-reduced zero phases: theta in [-pi, pi]
        const th = uni(-std.math.pi, std.math.pi);
        emit("sin", th, @sin(th));
        emit("cos", th, @cos(th));
        // ell(): sinh/cosh at c, and sinh of sqrt(c^2 - t^2) inside the band
        const sh = uni(0.0, 35.0);
        emit("sinh", sh, std.math.sinh(sh));
        emit("cosh", sh, std.math.cosh(sh));
        // ell() out-of-band: sin of sqrt(u^2 - c^2), u up to 50c
        const big = uni(0.0, 1600.0);
        emit("sin_big", big, @sin(big));
        // exp: exp(-eps), exp(lnw - lnw0), exp(-2c), kernel damping
        const ex = uni(-70.0, 2.0);
        emit("exp", ex, @exp(ex));
        // log1p on the window offset u = log(n/x), |u| <= eps ~ 1e-7,
        // and log for L, densities
        const lp = uni(-3e-7, 3e-7);
        emit("log1p", lp, std.math.log1p(lp));
        const lg = uni(1.0, 1.0e9);
        emit("log", lg, @log(lg));
        // erfc: Gaussian kernel weight
        const ec = uni(-7.5, 7.5);
        emit("erfc", ec, erfc(ec));
        // bessI0-style: sqrt is IEEE-exact, spot-check anyway
        const sq = uni(0.0, 1.0e6);
        emit("sqrt", sq, @sqrt(sq));
    }
}
extern "c" fn erfc(x: f64) f64;
