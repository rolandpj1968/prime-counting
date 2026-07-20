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
    \\  -v, --verbose          per-phase timing
    \\      --check            compare against the known π(10ⁿ) table
    \\      --no-time          print only the value
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

/// π(10ⁿ) for n = 0..22 (OEIS A006880 / Oliveira e Silva). Used only by --check.
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
    201467286689315906290,
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

fn die(comptime fmt: []const u8, args: anytype) noreturn {
    std.debug.print("pi: " ++ fmt ++ "\n", args);
    std.process.exit(2);
}

pub fn main(init: std.process.Init.Minimal) !void {
    const gpa = std.heap.page_allocator;
    var o = Opts{};
    var have_x = false;

    var it = std.process.Args.Iterator.init(init.args);
    _ = it.next(); // argv[0]
    while (it.next()) |a| {
        const eat = struct {
            fn v(iter: *std.process.Args.Iterator, name: []const u8) []const u8 {
                return iter.next() orelse die("{s} needs a value", .{name});
            }
        }.v;
        if (std.mem.eql(u8, a, "-h") or std.mem.eql(u8, a, "--help")) {
            std.debug.print("{s}", .{usage});
            return;
        } else if (std.mem.eql(u8, a, "-a") or std.mem.eql(u8, a, "--algo")) {
            const v = eat(&it, "--algo");
            o.algo = std.meta.stringToEnum(Algo, v) orelse die("unknown algorithm '{s}' (gourdon|lmo|meissel)", .{v});
        } else if (std.mem.eql(u8, a, "-t") or std.mem.eql(u8, a, "--threads")) {
            o.threads = std.fmt.parseInt(usize, eat(&it, "--threads"), 10) catch die("--threads needs an integer", .{});
        } else if (std.mem.eql(u8, a, "--alpha")) {
            o.alpha = std.fmt.parseFloat(f64, eat(&it, "--alpha")) catch die("--alpha needs a number", .{});
        } else if (std.mem.eql(u8, a, "--y")) {
            o.y = std.fmt.parseInt(u64, eat(&it, "--y"), 10) catch die("--y needs an integer", .{});
        } else if (std.mem.eql(u8, a, "--pin")) {
            o.pin = true;
        } else if (std.mem.eql(u8, a, "-v") or std.mem.eql(u8, a, "--verbose")) {
            o.verbose = true;
        } else if (std.mem.eql(u8, a, "--check")) {
            o.check = true;
        } else if (std.mem.eql(u8, a, "--no-time")) {
            o.time = false;
        } else if (a.len > 0 and a[0] == '-') {
            die("unknown option '{s}' (try --help)", .{a});
        } else {
            if (have_x) die("more than one x given ('{s}')", .{a});
            o.x = parseX(a) catch die("cannot parse x from '{s}'", .{a});
            have_x = true;
        }
    }
    if (!have_x) {
        std.debug.print("{s}", .{usage});
        std.process.exit(2);
    }

    // 0 ⇒ one thread per physical core, assuming 2 SMT siblings per core.
    const ncpu = std.Thread.getCpuCount() catch 1;
    if (o.threads == 0) o.threads = @max(1, ncpu / 2);

    var pins_buf: [256]u32 = undefined;
    var pins: ?[]const u32 = null;
    if (o.pin) {
        const n = @min(o.threads, pins_buf.len);
        for (0..n) |i| pins_buf[i] = @intCast(i * 2); // physical cores, skipping siblings
        pins = pins_buf[0..n];
    }

    // --alpha is resolved to y here so the algorithms keep a single knob.
    var y: ?u64 = o.y;
    if (y == null) if (o.alpha) |al| {
        var cr: u64 = @intFromFloat(std.math.cbrt(@as(f64, @floatFromInt(o.x))));
        while (cr > 1 and @as(u128, cr) * cr * cr > o.x) cr -= 1;
        while (@as(u128, cr + 1) * (cr + 1) * (cr + 1) <= o.x) cr += 1;
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

    const t0 = common.nowNs();
    const pi: i128 = switch (o.algo) {
        .gourdon => (try gourdon.piGourdonCfg(gpa, o.x, .{
            .y = y,
            .nthreads = o.threads,
            .pins = pins,
            .verbose = o.verbose,
        })).pi,
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
