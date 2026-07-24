//! `pi` — one binary for every π(x) implementation here, with the tuning knobs
//! exposed rather than recompiled. In the spirit of primecount's CLI.
//!
//!   pi 1e20                    π(10²⁰) with the default (fastest) algorithm
//!   pi 1e18 -a lmo -t 6        pick the algorithm and thread count
//!   pi 1e17 --alpha 6.5 -v     override the fitted α, show per-phase timing
//!   pi 1e16 --check            verify against the known π(10ⁿ) table
//!
//! Numbers accept 1e20, 10^20, 1_000_000 or plain digits.

const std = @import("std");
const common = @import("common.zig");
const gourdon = @import("gourdon.zig");
const lmo = @import("lmo.zig");
const meissel = @import("meissel.zig");

const Algo = enum { gourdon, lmo, meissel };

const Opts = struct {
    x: u128 = 0,
    algo: Algo = .gourdon,
    threads: usize = 1,
    alpha: ?f64 = null,
    y: ?u64 = null,
    pin: bool = false,
    verbose: bool = false,
    check: bool = false,
    time: bool = true,
    calibrate: bool = false,
    pin_list: ?[]const u8 = null,
    segw: ?usize = null,
    budget: f64 = 300,
    fit_a: ?f64 = null,
    fit_b: ?f64 = null,
    u128: bool = false,
    blocks: ?[]const u8 = null,
    emit: ?[]const u8 = null,
    merge: bool = false,
};

const usage =
    \\pi — combinatorial prime counting
    \\
    \\usage: pi <x> [options]
    \\
    \\  <x>                    1e20 | 10^20 | 1_000_000 | 1000000
    \\
    \\options:
    \\  -a, --algo <name>      gourdon (default) | lmo | meissel
    \\  -t, --threads <n>      worker threads (default 1; 0 = one per physical core)
    \\      --alpha <f>        override the fitted α in y = α·x^(1/3)
    \\      --y <n>            set y directly (overrides --alpha)
    \\      --pin              pin workers to cores 0,2,4,… (physical, skipping SMT)
    \\      --segw <n>         sweep segment width in integers (default 262080 =
    \\                         32 KB counter bits; multiple of 960, ≤ 2097152)
    \\      --pin-list <csv>   pin one worker per listed logical cpu (sets -t);
    \\                         e.g. 0,1 = both SMT threads of core 0
    \\  -v, --verbose          per-phase timing + periodic in-phase progress (≥60 s apart)
    \\      --check            compare against the known π(10ⁿ) table
    \\      --no-time          print only the value
    \\      --calibrate        measure α* on THIS machine at anchors 10¹⁵… and
    \\                         fit α(x) = A + B·ln x; x (if given) caps the anchors
    \\      --budget <sec>     calibration time budget (default 300)
    \\      --alpha-fit <A,B>  use a fit from a prior --calibrate run
    \\      --u128             force the u128 code path on x < 2^64 (measures the
    \\                         wide-arithmetic tax; result must be identical)
    \\      --blocks <a:b/N>   distributed task: sweep omega+B blocks [a,b) of the
    \\                         global N-grid only (with --emit; freeze y via
    \\                         --y/--alpha/--alpha-fit so all tasks agree)
    \\      --emit <file>      write the task's fragment to <file>
    \\      --merge <files..>  merge fragment files (must tile [0,N) with matching
    \\                         headers), compute A/Sigma/phi0 here, print pi
    \\  -h, --help             this text
    \\
    \\gourdon: all options. lmo: --y/--alpha serial only, --pin parallel only, u128 ok.
    \\meissel: no tuning, serial, u64 only.
    \\
;

/// Accepts 1e20, 10^20, 1_000_000 and plain digits. Exact for the forms that are
/// exact — 1e20 and 10^20 are built by repeated multiplication, not via f64, so
/// they are not subject to rounding at the top of the u128 range.
fn parseX(sraw: []const u8) !u128 {
    var buf: [64]u8 = undefined;
    var n: usize = 0;
    for (sraw) |c| {
        if (c == '_' or c == ',') continue;
        if (n == buf.len) return error.NumberTooLong;
        buf[n] = c;
        n += 1;
    }
    const s = buf[0..n];
    if (s.len == 0) return error.EmptyNumber;

    const sep = std.mem.indexOfAny(u8, s, "e^E");
    if (sep) |i| {
        // 'e' is scientific — the prefix MULTIPLIES 10^exp (1e12 = 1×10¹², 2e3 = 2000).
        // '^' is exponentiation — the prefix IS the base (10^12, 2^10).
        const pow = s[i] == '^';
        const pre_s = s[0..i];
        const exp_s = s[i + 1 ..];
        const pre: u128 = if (pre_s.len == 0) 10 else try std.fmt.parseInt(u128, pre_s, 10);
        const base: u128 = if (pow) pre else 10;
        const exp = try std.fmt.parseInt(u32, exp_s, 10);
        var v: u128 = if (pow) 1 else pre;
        var k: u32 = 0;
        while (k < exp) : (k += 1) {
            const ov = @mulWithOverflow(v, base);
            if (ov[1] != 0) return error.Overflow;
            v = ov[0];
        }
        return v;
    }
    return std.fmt.parseInt(u128, s, 10);
}

/// π(10ⁿ) for n = 0..24 (OEIS A006880 / Oliveira e Silva). Used only by --check.
const known = [_]i128{
    0,                       4,
    25,                      168,
    1229,                    9592,
    78498,                   664579,
    5761455,                 50847534,
    455052511,               4118054813,
    37607912018,             346065536839,
    3204941750802,           29844570422669,
    279238341033925,         2623557157654233,
    24739954287740860,       234057667276344607,
    2220819602560918840,     21127269486018731928,
    201467286689315906290,   1925320391606803968923,
    18435599767349200867866,
};

fn knownFor(x: u128) ?i128 {
    var v: u128 = 1;
    for (known, 0..) |k, n| {
        if (v == x) return k;
        if (n == known.len - 1) break;
        const ov = @mulWithOverflow(v, 10);
        if (ov[1] != 0) break;
        v = ov[0];
    }
    return null;
}

fn icbrt128(x: u128) u64 {
    var cr: u64 = @intFromFloat(std.math.cbrt(@as(f64, @floatFromInt(x))));
    while (cr > 1 and @as(u128, cr) * cr * cr > x) cr -= 1;
    while (@as(u128, cr + 1) * (cr + 1) * (cr + 1) <= x) cr += 1;
    return cr;
}

// ------------------------------------------------------------- α calibration
// α* is a property of the MACHINE (bandwidth, cache topology, thread count), not
// of the algorithm — so it is measured here and fitted to α(x) = A + B·ln x, and
// the run mode consumes the result via --alpha-fit. The design choices are this
// project's measured lessons, not defaults:
//   - anchors ascend 10^15, 10^16, … while the budget allows, because the fit
//     needs SPREAD: three clustered anchors once mis-measured the slope by 17%;
//   - one timing per probe: basins are flat (±25% in α ≈ 5–8% in time) and
//     run-to-run noise is 2–5%, so repeats buy nothing;
//   - a grid-edge minimum is never reported as an optimum — the grid extends
//     until the minimum is interior or a guard rail is hit (then marked "edge");
//   - every probe at one x must return the identical π: a correctness check that
//     needs no lookup table and works at any x on any machine.

const Probe = struct { a: f64, s: f64 };

fn probeLess(_: void, l: Probe, r: Probe) bool {
    return l.a < r.a;
}

/// Parabolic vertex in ln α through three bracketing probes; falls back to the
/// middle probe when the triple is not convex (flat basin / noise).
fn vertexLnA(p0: Probe, p1: Probe, p2: Probe) f64 {
    const x0 = @log(p0.a);
    const x1 = @log(p1.a);
    const x2 = @log(p2.a);
    const d = (x0 - x1) * (x0 - x2) * (x1 - x2);
    const qa = (x2 * (p1.s - p0.s) + x1 * (p0.s - p2.s) + x0 * (p2.s - p1.s)) / d;
    const qb = (x2 * x2 * (p0.s - p1.s) + x1 * x1 * (p2.s - p0.s) + x0 * x0 * (p1.s - p2.s)) / d;
    if (qa <= 0) return p1.a;
    return std.math.clamp(@exp(-qb / (2 * qa)), p0.a, p2.a);
}

const Anchor = struct { n: u32, lx: f64, astar: f64, edge: bool };

fn lsqFit(an: []const Anchor) struct { a: f64, b: f64 } {
    var sl: f64 = 0;
    var sy: f64 = 0;
    var sll: f64 = 0;
    var sly: f64 = 0;
    for (an) |q| {
        sl += q.lx;
        sy += q.astar;
        sll += q.lx * q.lx;
        sly += q.lx * q.astar;
    }
    const nn: f64 = @floatFromInt(an.len);
    const b = (nn * sly - sl * sy) / (nn * sll - sl * sl);
    return .{ .a = (sy - b * sl) / nn, .b = b };
}

const A_LO: f64 = 1.6; // z < y² needs α > 1; margin below this is never useful
const A_HI: f64 = 64.0;

fn runCalibrate(gpa: std.mem.Allocator, o: Opts, pins: ?[]const u32) !void {
    const nmax: u32 = blk: {
        if (o.x != 0) {
            const l: u32 = @intCast(std.math.log10_int(o.x));
            if (l < 17) die("--calibrate needs anchors to at least 10^17; omit x or give a larger one", .{});
            break :blk l;
        }
        break :blk 22;
    };
    std.debug.print("calibrating alpha(x): {d} thread(s), budget {d:.0} s, anchors 10^15..10^{d} as budget allows\n", .{ o.threads, o.budget, nmax });

    const t0 = common.nowNs();
    var anchors: [16]Anchor = undefined;
    var na: usize = 0;
    var est_next: f64 = 0; // projected per-probe seconds at the NEXT anchor

    var n: u32 = 15;
    while (n <= nmax and na < anchors.len) : (n += 1) {
        const el0 = @as(f64, @floatFromInt(common.nowNs() - t0)) / 1e9;
        // the first three anchors always run (a fit needs them); afterwards an
        // anchor starts only if its minimal three probes plausibly fit the budget
        if (na >= 3 and el0 + 3.0 * est_next > o.budget) break;

        var xv: u128 = 1;
        for (0..n) |_| xv *= 10;
        const cr = icbrt128(xv);

        // seed: own running fit > previous anchor + default slope > built-in default
        var a0: f64 = undefined;
        if (na >= 2) {
            const f = lsqFit(anchors[0..na]);
            a0 = f.a + f.b * @log(@as(f64, @floatFromInt(xv)));
        } else if (na == 1) {
            a0 = anchors[0].astar + 0.5980 * @log(10.0);
        } else {
            a0 = gourdon.defaultAlpha(xv);
        }
        a0 = std.math.clamp(a0, 2.2, 47.0);

        const R = 1.35;
        var probes: [12]Probe = undefined;
        var np: usize = 0;
        var piv: i128 = 0;
        var edge = false;

        // initial ascending triple, then extend past whichever edge holds the
        // minimum until it is interior (or a rail / the budget stops us)
        for ([_]f64{ a0 / R, a0, a0 * R }) |aa| {
            const ts = common.nowNs();
            const r = try gourdon.piGourdonCfg(gpa, xv, .{
                .y = @intFromFloat(aa * @as(f64, @floatFromInt(cr))),
                .nthreads = o.threads,
                .pins = pins,
            });
            const secs = @as(f64, @floatFromInt(common.nowNs() - ts)) / 1e9;
            if (np == 0) piv = r.pi else if (r.pi != piv)
                die("pi(10^{d}) changed with alpha — correctness bug, calibration aborted", .{n});
            probes[np] = .{ .a = aa, .s = secs };
            np += 1;
            std.debug.print("  10^{d}  alpha {d:>6.2}  {d:>8.2} s\n", .{ n, aa, secs });
        }
        if (knownFor(xv)) |w| if (piv != w)
            die("pi(10^{d}) = {d} disagrees with the known value {d}", .{ n, piv, w });

        while (np < probes.len) {
            std.mem.sort(Probe, probes[0..np], {}, probeLess);
            var mi: usize = 0;
            for (probes[0..np], 0..) |q, qi| {
                if (q.s < probes[mi].s) mi = qi;
            }
            var aa: f64 = 0;
            if (mi == 0) aa = probes[0].a / R else if (mi == np - 1) aa = probes[np - 1].a * R else break;
            const el = @as(f64, @floatFromInt(common.nowNs() - t0)) / 1e9;
            if (aa < A_LO or aa > A_HI or el + probes[mi].s > o.budget * 1.15) {
                edge = true;
                break;
            }
            const ts = common.nowNs();
            const r = try gourdon.piGourdonCfg(gpa, xv, .{
                .y = @intFromFloat(aa * @as(f64, @floatFromInt(cr))),
                .nthreads = o.threads,
                .pins = pins,
            });
            const secs = @as(f64, @floatFromInt(common.nowNs() - ts)) / 1e9;
            if (r.pi != piv) die("pi(10^{d}) changed with alpha — correctness bug, calibration aborted", .{n});
            probes[np] = .{ .a = aa, .s = secs };
            np += 1;
            std.debug.print("  10^{d}  alpha {d:>6.2}  {d:>8.2} s   (extending)\n", .{ n, aa, secs });
        }

        std.mem.sort(Probe, probes[0..np], {}, probeLess);
        var mi: usize = 0;
        for (probes[0..np], 0..) |q, qi| {
            if (q.s < probes[mi].s) mi = qi;
        }
        const astar = if (mi > 0 and mi < np - 1)
            vertexLnA(probes[mi - 1], probes[mi], probes[mi + 1])
        else blk: {
            edge = true;
            break :blk probes[mi].a;
        };
        anchors[na] = .{ .n = n, .lx = @log(@as(f64, @floatFromInt(xv))), .astar = astar, .edge = edge };
        na += 1;
        std.debug.print("  10^{d}  alpha* = {d:.2}{s}\n", .{ n, astar, if (edge) "  (edge — treat with suspicion)" else "" });

        var mean: f64 = 0;
        for (probes[0..np]) |q| mean += q.s;
        est_next = (mean / @as(f64, @floatFromInt(np))) * 4.6; // measured per-decade growth
    }

    if (na < 3) die("only {d} anchors fit the budget; a fit needs 3 — raise --budget", .{na});
    const f = lsqFit(anchors[0..na]);
    std.debug.print("\n{s:>6} {s:>8} {s:>8} {s:>7}\n", .{ "x", "alpha*", "fit", "resid" });
    var worst: f64 = 0;
    for (anchors[0..na]) |q| {
        const pred = f.a + f.b * q.lx;
        const res = q.astar - pred;
        if (@abs(res) > worst) worst = @abs(res);
        std.debug.print("  10^{d:<3} {d:>8.2} {d:>8.2} {s}{d:>6.2}{s}\n", .{ q.n, q.astar, pred, if (res < 0) "-" else "+", @abs(res), if (q.edge) "  (edge)" else "" });
    }
    // Err-high bias: the basin is asymmetric. Below α* the fold work blows up
    // hyperbolically (z = x/y) and terminates in the z < y² correctness wall;
    // above α* the cost is a gentle linear leaf/memory slope. Measured at 10^20:
    // −30% in α costs +16.7%, +27% costs +6.7%; at 10^18 α=1.5 was 4.8× worse
    // while α=24 was 1.7×. So the suggested fit aims one residual-RMS ABOVE the
    // symmetric least squares — cheap where the valley is flat, protective where
    // it is a cliff.
    var ssq: f64 = 0;
    for (anchors[0..na]) |q| {
        const res = q.astar - (f.a + f.b * q.lx);
        ssq += res * res;
    }
    const rms = @sqrt(ssq / @as(f64, @floatFromInt(na)));
    const bias = @max(0.3, rms);
    std.debug.print("\nalpha(x) = {d:.4} + {d:.4} * ln x   (symmetric fit; worst residual {d:.2}, {d} thread(s))\n", .{ f.a, f.b, worst, o.threads });
    std.debug.print("err-high bias +{d:.2} (rms residual; low side of the basin is the cliff)\n", .{bias});
    std.debug.print("apply with:  --alpha-fit={d:.4},{d:.4}\n", .{ f.a + bias, f.b });
    if (worst > 1.5) std.debug.print("note: residuals are large — this machine may want a per-x --alpha rather than a line\n", .{});
}

/// Yields arguments, splitting a leading "--opt=value" into "--opt" with the value
/// parked for the next eat() — so the calibrator's printed "--alpha-fit=A,B" form
/// works as well as "--alpha-fit A,B".
fn nextArg(it: *std.process.Args.Iterator, eq_val: *?[]const u8) ?[]const u8 {
    if (eq_val.*) |_| {} // an unconsumed =value means a flag got one: caught by eat/unknown
    const a = it.next() orelse return null;
    if (a.len > 2 and a[0] == '-' and a[1] == '-') {
        if (std.mem.indexOfScalar(u8, a, '=')) |e| {
            eq_val.* = a[e + 1 ..];
            return a[0..e];
        }
    }
    return a;
}

fn die(comptime fmt: []const u8, args: anytype) noreturn {
    std.debug.print("pi: " ++ fmt ++ "\n", args);
    std.process.exit(2);
}

/// Distributed task: sweep blocks [a,b) of the global N-grid, write the fragment.
fn runEmit(gpa: std.mem.Allocator, o: *const Opts, pins: ?[]const u32, y: ?u64) !void {
    const spec = o.blocks.?;
    const colon = std.mem.indexOfScalar(u8, spec, ':') orelse die("--blocks wants a:b/N", .{});
    const slash = std.mem.indexOfScalar(u8, spec, '/') orelse die("--blocks wants a:b/N", .{});
    if (slash < colon) die("--blocks wants a:b/N", .{});
    const ba = std.fmt.parseInt(usize, spec[0..colon], 10) catch die("--blocks: bad a", .{});
    const bb = std.fmt.parseInt(usize, spec[colon + 1 .. slash], 10) catch die("--blocks: bad b", .{});
    const bn = std.fmt.parseInt(usize, spec[slash + 1 ..], 10) catch die("--blocks: bad N", .{});
    if (ba >= bb or bb > bn or bn == 0) die("--blocks: need 0 <= a < b <= N", .{});
    const t0 = common.nowNs();
    const fr = gourdon.piGourdonFragment(gpa, o.x, .{
        .y = y,
        .nthreads = o.threads,
        .pins = pins,
        .verbose = o.verbose,
        .segw = o.segw,
    }, bn, ba, bb) catch |e| die("fragment failed: {s}", .{@errorName(e)});
    defer fr.frag.deinit(gpa);
    const secs = @as(f64, @floatFromInt(common.nowNs() - t0)) / 1e9;
    const nax = fr.frag.mu_sum.len;
    const buf = try gpa.alloc(u8, 4096 + 2 * nax * 24);
    defer gpa.free(buf);
    var off: usize = 0;
    off += (std.fmt.bufPrint(buf[off..], "pifrag 1\nx {d}\ny {d}\nsegw {d}\nnb {d}\nt0 {d}\nt1 {d}\nnax {d}\nomega {d}\nb {d}\nbcount {d}\ntotalfull {d}\nmusum", .{ o.x, fr.y, fr.segw, bn, ba, bb, nax, fr.frag.omega, fr.frag.b, fr.frag.bcount, fr.frag.total_full }) catch unreachable).len;
    for (fr.frag.mu_sum) |v| off += (std.fmt.bufPrint(buf[off..], " {d}", .{v}) catch unreachable).len;
    off += (std.fmt.bufPrint(buf[off..], "\ntotalsum", .{}) catch unreachable).len;
    for (fr.frag.total_sum) |v| off += (std.fmt.bufPrint(buf[off..], " {d}", .{v}) catch unreachable).len;
    off += (std.fmt.bufPrint(buf[off..], "\nend\n", .{}) catch unreachable).len;
    var tio = std.Io.Threaded.init(gpa, .{});
    defer tio.deinit();
    std.Io.Dir.cwd().writeFile(tio.io(), .{ .sub_path = o.emit.?, .data = buf[0..off] }) catch |e| die("cannot write '{s}': {s}", .{ o.emit.?, @errorName(e) });
    std.debug.print("fragment [{d},{d})/{d}  nax {d}  -> {s} ({d} bytes)  {d:.3} s\n", .{ ba, bb, bn, nax, o.emit.?, off, secs });
}

const FragFile = struct { x: u128, y: u64, segw: usize, nb: usize, t0: usize, t1: usize, fo: gourdon.FragOut };

fn expectKw(itk: anytype, want: []const u8, path: []const u8) void {
    const tok = itk.next() orelse die("'{s}': truncated (wanted {s})", .{ path, want });
    if (!std.mem.eql(u8, tok, want)) die("'{s}': expected '{s}', got '{s}'", .{ path, want, tok });
}
fn expectNum(comptime T: type, itk: anytype, path: []const u8) T {
    const tok = itk.next() orelse die("'{s}': truncated", .{path});
    return std.fmt.parseInt(T, tok, 10) catch die("'{s}': bad number '{s}'", .{ path, tok });
}

fn parseFragFile(gpa: std.mem.Allocator, io: std.Io, path: []const u8) !FragFile {
    const data = std.Io.Dir.cwd().readFileAlloc(io, path, gpa, .limited(1 << 26)) catch |e| die("cannot read '{s}': {s}", .{ path, @errorName(e) });
    defer gpa.free(data);
    var itk = std.mem.tokenizeAny(u8, data, " \n\r\t");
    expectKw(&itk, "pifrag", path);
    expectKw(&itk, "1", path);
    expectKw(&itk, "x", path);
    const x = expectNum(u128, &itk, path);
    expectKw(&itk, "y", path);
    const y = expectNum(u64, &itk, path);
    expectKw(&itk, "segw", path);
    const segw = expectNum(usize, &itk, path);
    expectKw(&itk, "nb", path);
    const nb = expectNum(usize, &itk, path);
    expectKw(&itk, "t0", path);
    const t0 = expectNum(usize, &itk, path);
    expectKw(&itk, "t1", path);
    const t1 = expectNum(usize, &itk, path);
    expectKw(&itk, "nax", path);
    const nax = expectNum(usize, &itk, path);
    if (nax == 0 or nax > (1 << 28)) die("'{s}': implausible nax {d}", .{ path, nax });
    expectKw(&itk, "omega", path);
    const omega = expectNum(i128, &itk, path);
    expectKw(&itk, "b", path);
    const b = expectNum(i128, &itk, path);
    expectKw(&itk, "bcount", path);
    const bcount = expectNum(u64, &itk, path);
    expectKw(&itk, "totalfull", path);
    const total_full = expectNum(i64, &itk, path);
    expectKw(&itk, "musum", path);
    const mu = try gpa.alloc(i64, nax);
    errdefer gpa.free(mu);
    for (mu) |*v| v.* = expectNum(i64, &itk, path);
    expectKw(&itk, "totalsum", path);
    const ts = try gpa.alloc(i64, nax);
    errdefer gpa.free(ts);
    for (ts) |*v| v.* = expectNum(i64, &itk, path);
    expectKw(&itk, "end", path);
    return .{ .x = x, .y = y, .segw = segw, .nb = nb, .t0 = t0, .t1 = t1, .fo = .{
        .omega = omega,
        .b = b,
        .bcount = bcount,
        .total_full = total_full,
        .mu_sum = mu,
        .total_sum = ts,
    } };
}

/// Merge fragments that tile [0, nb): validate headers, prove the tiling,
/// carry-walk, then compute the remaining terms here and print pi.
fn runMerge(gpa: std.mem.Allocator, o: *Opts, pins: ?[]const u32, files: []const []const u8) !void {
    if (files.len == 0) die("--merge needs fragment files", .{});
    var tio = std.Io.Threaded.init(gpa, .{});
    defer tio.deinit();
    const ffs = try gpa.alloc(FragFile, files.len);
    for (files, 0..) |p, i| ffs[i] = try parseFragFile(gpa, tio.io(), p);
    const h = ffs[0];
    for (ffs[1..]) |f| if (f.x != h.x or f.y != h.y or f.segw != h.segw or f.nb != h.nb or f.fo.mu_sum.len != h.fo.mu_sum.len)
        die("fragment headers disagree (x/y/segw/nb/nax) — fragments from different runs?", .{});
    std.mem.sort(FragFile, ffs, {}, struct {
        fn lt(_: void, a: FragFile, b: FragFile) bool {
            return a.t0 < b.t0;
        }
    }.lt);
    // tiling proof: a hole = a lost task to re-issue, an overlap = a double-spend
    var expect_t: usize = 0;
    for (ffs) |f| {
        if (f.t0 != expect_t) die("tiling broken at block {d}: next fragment covers [{d},{d})", .{ expect_t, f.t0, f.t1 });
        expect_t = f.t1;
    }
    if (expect_t != h.nb) die("tiling incomplete: covered [0,{d}) of {d} blocks", .{ expect_t, h.nb });
    const fos = try gpa.alloc(gourdon.FragOut, ffs.len);
    for (ffs, 0..) |f, i| fos[i] = f.fo;
    const omb = try gourdon.mergeFragments(gpa, h.fo.mu_sum.len, fos);
    std.debug.print("merged {d} fragments, {d} blocks — computing A/Sigma + phi0 locally\n", .{ ffs.len, h.nb });

    o.x = h.x; // for --check
    const t0 = common.nowNs();
    const r = try gourdon.piGourdonCfg(gpa, h.x, .{
        .y = h.y,
        .nthreads = o.threads,
        .pins = pins,
        .verbose = o.verbose,
        .segw = h.segw,
        .omb_in = omb,
    });
    const secs = @as(f64, @floatFromInt(common.nowNs() - t0)) / 1e9;
    if (!o.time) {
        std.debug.print("{d}\n", .{r.pi});
    } else {
        const ru = std.posix.getrusage(std.posix.rusage.SELF);
        std.debug.print("{d}\n  {d:.3} s (merge-side)   {d} thread(s)   peakRSS {d} MB\n", .{ r.pi, secs, o.threads, @divTrunc(ru.maxrss, 1024) });
    }
    if (o.check) {
        if (knownFor(o.x)) |w| {
            if (r.pi == w) {
                std.debug.print("  check: MATCH\n", .{});
            } else {
                std.debug.print("  check: MISMATCH — expected {d}\n", .{w});
                std.process.exit(1);
            }
        } else {
            std.debug.print("  check: no known value for this x\n", .{});
        }
    }
}

pub fn main(init: std.process.Init.Minimal) !void {
    const gpa = std.heap.page_allocator;
    var o = Opts{};
    var have_x = false;
    var mfiles: [128][]const u8 = undefined;
    var nmf: usize = 0;

    var it = std.process.Args.Iterator.init(init.args);
    _ = it.next(); // argv[0]
    var eq_val: ?[]const u8 = null; // value split off a --opt=value form
    while (nextArg(&it, &eq_val)) |a| {
        const eat = struct {
            fn v(iter: *std.process.Args.Iterator, ev: *?[]const u8, name: []const u8) []const u8 {
                if (ev.*) |val| {
                    ev.* = null;
                    return val;
                }
                return iter.next() orelse die("{s} needs a value", .{name});
            }
        }.v;
        if (std.mem.eql(u8, a, "-h") or std.mem.eql(u8, a, "--help")) {
            std.debug.print("{s}", .{usage});
            return;
        } else if (std.mem.eql(u8, a, "-a") or std.mem.eql(u8, a, "--algo")) {
            const v = eat(&it, &eq_val, "--algo");
            o.algo = std.meta.stringToEnum(Algo, v) orelse die("unknown algorithm '{s}' (gourdon|lmo|meissel)", .{v});
        } else if (std.mem.eql(u8, a, "-t") or std.mem.eql(u8, a, "--threads")) {
            o.threads = std.fmt.parseInt(usize, eat(&it, &eq_val, "--threads"), 10) catch die("--threads needs an integer", .{});
        } else if (std.mem.eql(u8, a, "--calibrate")) {
            o.calibrate = true;
        } else if (std.mem.eql(u8, a, "--budget")) {
            o.budget = std.fmt.parseFloat(f64, eat(&it, &eq_val, "--budget")) catch die("--budget needs seconds", .{});
        } else if (std.mem.eql(u8, a, "--alpha-fit")) {
            const v = eat(&it, &eq_val, "--alpha-fit");
            const c = std.mem.indexOfScalar(u8, v, ',') orelse die("--alpha-fit wants A,B (e.g. -17.7,0.598)", .{});
            o.fit_a = std.fmt.parseFloat(f64, v[0..c]) catch die("--alpha-fit: bad A", .{});
            o.fit_b = std.fmt.parseFloat(f64, v[c + 1 ..]) catch die("--alpha-fit: bad B", .{});
        } else if (std.mem.eql(u8, a, "--alpha")) {
            o.alpha = std.fmt.parseFloat(f64, eat(&it, &eq_val, "--alpha")) catch die("--alpha needs a number", .{});
        } else if (std.mem.eql(u8, a, "--y")) {
            o.y = std.fmt.parseInt(u64, eat(&it, &eq_val, "--y"), 10) catch die("--y needs an integer", .{});
        } else if (std.mem.eql(u8, a, "--segw")) {
            const v = eat(&it, &eq_val, "--segw");
            const n = std.fmt.parseInt(usize, v, 10) catch die("--segw needs an integer", .{});
            if (n < 960 or n % 960 != 0 or n > (1 << 21)) die("--segw must be a multiple of 960, 960..2097152", .{});
            o.segw = n;
        } else if (std.mem.eql(u8, a, "--pin-list")) {
            o.pin_list = eat(&it, &eq_val, "--pin-list");
        } else if (std.mem.eql(u8, a, "--pin")) {
            o.pin = true;
        } else if (std.mem.eql(u8, a, "-v") or std.mem.eql(u8, a, "--verbose")) {
            o.verbose = true;
        } else if (std.mem.eql(u8, a, "--u128")) {
            o.u128 = true;
        } else if (std.mem.eql(u8, a, "--blocks")) {
            o.blocks = eat(&it, &eq_val, "--blocks");
        } else if (std.mem.eql(u8, a, "--emit")) {
            o.emit = eat(&it, &eq_val, "--emit");
        } else if (std.mem.eql(u8, a, "--merge")) {
            o.merge = true;
        } else if (std.mem.eql(u8, a, "--check")) {
            o.check = true;
        } else if (std.mem.eql(u8, a, "--no-time")) {
            o.time = false;
        } else if (a.len > 0 and a[0] == '-') {
            die("unknown option '{s}' (try --help)", .{a});
        } else if (o.merge) {
            // after --merge, positionals are fragment files
            if (nmf == mfiles.len) die("--merge: too many files", .{});
            mfiles[nmf] = a;
            nmf += 1;
        } else {
            if (have_x) die("more than one x given ('{s}')", .{a});
            o.x = parseX(a) catch die("cannot parse x from '{s}'", .{a});
            have_x = true;
        }
    }
    if (!o.calibrate and !have_x and !o.merge) {
        std.debug.print("{s}", .{usage});
        std.process.exit(2);
    }

    // 0 ⇒ one thread per physical core, assuming 2 SMT siblings per core.
    const ncpu = std.Thread.getCpuCount() catch 1;
    if (o.threads == 0) o.threads = @max(1, ncpu / 2);

    var pins_buf: [256]u32 = undefined;
    var pins: ?[]const u32 = null;
    if (o.pin_list) |pl| {
        // explicit topology: comma-separated logical CPUs, one worker pinned to
        // each; thread count follows the list (SMT/cloud topologies welcome)
        var n: usize = 0;
        var itp = std.mem.tokenizeScalar(u8, pl, ',');
        while (itp.next()) |tok| {
            if (n == pins_buf.len) die("--pin-list: too many cpus", .{});
            pins_buf[n] = std.fmt.parseInt(u32, tok, 10) catch die("--pin-list: bad cpu '{s}'", .{tok});
            n += 1;
        }
        if (n == 0) die("--pin-list: empty", .{});
        o.threads = n;
        pins = pins_buf[0..n];
    } else if (o.pin) {
        const n = @min(o.threads, pins_buf.len);
        for (0..n) |i| pins_buf[i] = @intCast(i * 2); // physical cores, skipping siblings
        pins = pins_buf[0..n];
    }

    if (o.calibrate) {
        try runCalibrate(gpa, o, pins);
        return;
    }

    if (o.merge) {
        if (o.blocks != null or o.emit != null) die("--merge excludes --blocks/--emit", .{});
        if (have_x) die("--merge takes fragment files, not x (x comes from the headers)", .{});
        try runMerge(gpa, &o, pins, mfiles[0..nmf]);
        return;
    }

    // a prior calibration's fit, unless --alpha/--y override it
    if (o.alpha == null and o.y == null) if (o.fit_a) |fa| {
        o.alpha = std.math.clamp(fa + o.fit_b.? * @log(@as(f64, @floatFromInt(o.x))), 1.6, 64.0);
    };

    // --alpha is resolved to y here so the algorithms keep a single knob.
    var y: ?u64 = o.y;
    if (y == null) if (o.alpha) |al| {
        const cr = icbrt128(o.x);
        const yf = al * @as(f64, @floatFromInt(cr));
        if (yf < 1 or yf >= 1.8e19) die("--alpha {d} gives an out-of-range y", .{al});
        y = @intFromFloat(yf);
    };

    if (o.algo == .meissel and (o.alpha != null or o.y != null or o.pin))
        std.debug.print("pi: note: meissel takes no tuning parameters, ignored\n", .{});
    if (o.algo == .lmo and o.threads > 1 and (o.alpha != null or o.y != null))
        std.debug.print("pi: note: parallel lmo takes no y, --alpha/--y ignored\n", .{});
    if (o.algo != .gourdon and o.verbose)
        std.debug.print("pi: note: --verbose is gourdon-only\n", .{});

    if (o.u128 and o.algo != .gourdon) die("--u128 is gourdon-only", .{});
    if (o.u128 and o.x <= 10_000) die("--u128 needs x > 10^4 (below that the direct oracle answers)", .{});

    if (o.blocks != null or o.emit != null) {
        if (o.blocks == null or o.emit == null) die("--blocks and --emit go together", .{});
        if (o.algo != .gourdon) die("--blocks is gourdon-only", .{});
        try runEmit(gpa, &o, pins, y);
        return;
    }

    const t0 = common.nowNs();
    const gcfg: gourdon.Config = .{
        .y = y,
        .nthreads = o.threads,
        .pins = pins,
        .verbose = o.verbose,
        .segw = o.segw,
    };
    const pi: i128 = switch (o.algo) {
        // --u128 bypasses the width dispatch: same x, wide arithmetic throughout.
        .gourdon => if (o.u128 and o.x <= std.math.maxInt(u64))
            (try gourdon.piGourdonV(u128, gpa, o.x, gcfg)).pi
        else
            (try gourdon.piGourdonCfg(gpa, o.x, gcfg)).pi,
        .lmo => blk: {
            // piLMOPar takes no y; --y/--alpha therefore apply to serial lmo only.
            const r = if (o.threads > 1)
                try lmo.piLMOPar(gpa, o.x, o.threads, 8, pins)
            else
                try lmo.piLMO(gpa, o.x, y, null);
            break :blk @intCast(r.pi);
        },
        .meissel => blk: {
            if (o.x > std.math.maxInt(u64)) die("meissel is u64-only; x > 2^64 needs -a gourdon or -a lmo", .{});
            if (o.threads > 1) std.debug.print("pi: note: meissel is serial, --threads ignored\n", .{});
            break :blk @intCast(try meissel.pi(gpa, @intCast(o.x)));
        },
    };
    const secs = @as(f64, @floatFromInt(common.nowNs() - t0)) / 1e9;

    if (!o.time) {
        std.debug.print("{d}\n", .{pi});
    } else {
        const ru = std.posix.getrusage(std.posix.rusage.SELF);
        std.debug.print("{d}\n  {d:.3} s   {d} thread(s)   peakRSS {d} MB\n", .{
            pi, secs, o.threads, @divTrunc(ru.maxrss, 1024),
        });
    }

    if (o.check) {
        if (knownFor(o.x)) |w| {
            if (pi == w) {
                std.debug.print("  check: MATCH\n", .{});
            } else {
                std.debug.print("  check: MISMATCH — expected {d}\n", .{w});
                std.process.exit(1);
            }
        } else {
            std.debug.print("  check: no known value for this x (table covers 10^0..10^22)\n", .{});
        }
    }
}
