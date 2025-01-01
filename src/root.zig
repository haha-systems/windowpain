const std = @import("std");
const testing = std.testing;
const fs = std.fs;
const json = std.json;

pub const FastaSequence = struct {
    header: []const u8,
    position: usize,
    length: usize,

    pub fn jsonStringify(self: FastaSequence, jw: anytype) !void {
        try jw.beginObject();
        try jw.objectField("header");
        try jw.write(self.header);
        try jw.objectField("position");
        try jw.write(self.position);
        try jw.objectField("length");
        try jw.write(self.length);
        try jw.endObject();
    }
};

pub const Command = enum {
    index,
    read,

    pub fn fromString(str: []const u8) !Command {
        if (std.mem.eql(u8, str, "index")) return .index;
        if (std.mem.eql(u8, str, "read")) return .read;
        return error.InvalidCommand;
    }
};

const IndexEntry = struct {
    header: []const u8,
    position: usize,
    length: usize,
};

pub fn readSequenceWindow(file: fs.File, sequence: FastaSequence, window_start: usize, window_size: usize) ![]u8 {
    if (window_start >= sequence.length) {
        return error.WindowOutOfBounds;
    }

    const actual_window_size = @min(window_size, sequence.length - window_start);
    const page_size = std.mem.page_size;

    const file_position = sequence.position + window_start;
    const aligned_offset = (file_position / page_size) * page_size;
    const offset_adjustment = file_position - aligned_offset;

    const map_size = actual_window_size + offset_adjustment;

    const ptr = try std.posix.mmap(
        null,
        map_size,
        std.posix.PROT.READ,
        .{ .TYPE = .PRIVATE },
        file.handle,
        aligned_offset,
    );
    defer std.posix.munmap(ptr[0..map_size]);

    // Count newlines in the window to allocate correct buffer size
    var newline_count: usize = 0;
    for (ptr[offset_adjustment..map_size]) |c| {
        if (c == '\n') newline_count += 1;
    }

    // Allocate buffer for sequence without newlines
    var sequence_data = try std.heap.page_allocator.alloc(u8, actual_window_size - newline_count);

    // Copy data, skipping newlines
    var write_pos: usize = 0;
    for (ptr[offset_adjustment..map_size]) |c| {
        if (c != '\n') {
            sequence_data[write_pos] = c;
            write_pos += 1;
        }
    }

    return sequence_data;
}

pub fn loadIndex(allocator: std.mem.Allocator, index_path: []const u8) !std.ArrayList(FastaSequence) {
    const index_file = try fs.cwd().openFile(index_path, .{});
    defer index_file.close();

    const index_contents = try index_file.readToEndAlloc(allocator, std.math.maxInt(usize));
    defer allocator.free(index_contents);

    var sequences = std.ArrayList(FastaSequence).init(allocator);
    errdefer sequences.deinit();

    // Parse the JSON array
    const IndexFile = []IndexEntry;
    const parsed = try json.parseFromSlice(IndexFile, allocator, index_contents, .{});
    defer parsed.deinit();

    // Convert the parsed entries to FastaSequences
    for (parsed.value) |entry| {
        const header = try allocator.dupe(u8, entry.header);
        try sequences.append(FastaSequence{
            .header = header,
            .position = entry.position,
            .length = entry.length,
        });
    }

    return sequences;
}

test "FastaSequence JSON serialization" {
    const sequence = .{
        .header = ">test sequence",
        .position = 100,
        .length = 50,
    };

    var string = std.ArrayList(u8).init(testing.allocator);
    defer string.deinit();

    try std.json.stringify(sequence, .{}, string.writer());

    try testing.expectEqualStrings("{\"header\":\">test sequence\",\"position\":100,\"length\":50}", string.items);
}

test "Command parsing" {
    try testing.expectEqual(Command.index, try Command.fromString("index"));
    try testing.expectEqual(Command.read, try Command.fromString("read"));
    try testing.expectError(error.InvalidCommand, Command.fromString("invalid"));
}

test "readSequenceWindow basic functionality" {
    // Create a temporary test file
    const test_path = "test.fasta";
    const test_content =
        \\>test sequence
        \\ACTG
        \\GTCA
        \\
    ;

    const file = try fs.cwd().createFile(test_path, .{ .read = true });
    defer {
        file.close();
        fs.cwd().deleteFile(test_path) catch {};
    }

    try file.writeAll(test_content);

    const sequence = FastaSequence{
        .header = ">test sequence",
        .position = 14, // Position after header and newline
        .length = 8, // Length of actual sequence without newlines
    };

    const window = try readSequenceWindow(file, sequence, 0, 8);
    defer std.heap.page_allocator.free(window);

    try testing.expectEqualStrings("ACTGGTCA", window);
}

test "loadIndex basic functionality" {
    // Create a temporary test index file
    const test_path = "test.json";
    const test_content =
        \\[{"header":">test1","position":0,"length":10},
        \\{"header":">test2","position":15,"length":20}]
    ;

    const file = try fs.cwd().createFile(test_path, .{ .read = true });
    defer {
        file.close();
        fs.cwd().deleteFile(test_path) catch {};
    }

    try file.writeAll(test_content);

    var sequences = try loadIndex(testing.allocator, test_path);
    defer {
        for (sequences.items) |seq| {
            testing.allocator.free(seq.header);
        }
        sequences.deinit();
    }

    try testing.expectEqual(@as(usize, 2), sequences.items.len);
    try testing.expectEqualStrings(">test1", sequences.items[0].header);
    try testing.expectEqual(@as(usize, 0), sequences.items[0].position);
    try testing.expectEqual(@as(usize, 10), sequences.items[0].length);
}
