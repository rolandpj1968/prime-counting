const std = @import("std");
const g = @import("gourdon.zig");
pub fn main() !void {
    const r = try g.piGourdon(std.heap.page_allocator, 100_000_000_000_000, null);
    std.debug.print("pi(10^14) = {d}\n", .{r.pi});
}
