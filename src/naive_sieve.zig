//! NaiveSieve — the textbook whole-array Sieve of Eratosthenes, written ONCE and
//! generic over a backing Store. flat []bool, hand-rolled []u64, and
//! std.DynamicBitSet differ only in how a bit is stored and read; the sieve
//! logic is identical. So we factor the logic out here and parameterise the
//! store — exactly the orthogonality you spotted.
//!
//! A Store (see src/stores/) is any namespace exposing:
//!   pub const name: []const u8;
//!   pub const State: type;
//!   pub fn footprintBytes(n) usize;
//!   pub fn init(gpa, n) !State;
//!   pub fn deinit(*State, gpa) void;
//!   pub fn clearAll(*State, n) void;              // all bits -> unset (0)
//!   pub fn get(*const State, i: u64) bool;        // is i marked composite?
//!   pub fn set(*State, i: u64) void;              // mark i composite
//!   pub fn countComposites(*const State, n) u64;  // # set bits in [0, n]
//!
//! Because Zig monomorphises and inlines, `Store.set(i)` compiles to exactly the
//! same machine code as an inline strike — the abstraction is zero-cost. We
//! already measured this: DynamicBitSet.set matched the hand-inline loop.

const std = @import("std");
const common = @import("common.zig");

pub fn NaiveSieve(comptime Store: type) type {
    return struct {
        pub const name = Store.name;
        pub const State = Store.State;

        // Storage-only operations are thin pass-throughs to the backing store
        // (inlined away at zero cost — this is just wiring).
        pub fn footprintBytes(n: u64) usize {
            return Store.footprintBytes(n);
        }

        pub fn init(gpa: std.mem.Allocator, n: u64) !State {
            return Store.init(gpa, n);
        }

        pub fn deinit(st: *State, gpa: std.mem.Allocator) void {
            Store.deinit(st, gpa);
        }

        // The actual sieve logic — written once, for every store.
        pub fn sieve(st: *State, n: u64) void {
            Store.clearAll(st, n);
            const lim = common.isqrt(n);
            var i: u64 = 2;
            while (i <= lim) : (i += 1) {
                if (!Store.get(st, i)) {
                    var j = i * i;
                    while (j <= n) : (j += i) Store.set(st, j);
                }
            }
        }

        pub fn count(st: *const State, n: u64) u64 {
            // set bits are composites in [0,n]; 0 and 1 are unset non-primes
            return n - 1 - Store.countComposites(st, n);
        }
    };
}
