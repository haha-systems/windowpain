const std = @import("std");
const fs = std.fs;
const print = std.debug.print;
const json = std.json;

const FastaSequence = struct {
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

const Command = enum {
    index,
    read,

    pub fn fromString(str: []const u8) !Command {
        if (std.mem.eql(u8, str, "index")) return .index;
        if (std.mem.eql(u8, str, "read")) return .read;
        return error.InvalidCommand;
    }
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

const IndexEntry = struct {
    header: []const u8,
    position: usize,
    length: usize,
};

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

pub fn scanFastaSequences(file: fs.File, allocator: std.mem.Allocator) !std.ArrayList(FastaSequence) {
    const file_size = try file.getEndPos();
    const page_size = std.mem.page_size;

    var sequences = std.ArrayList(FastaSequence).init(allocator);
    errdefer sequences.deinit();

    const chunk_size: usize = 1024 * 1024;
    var current_pos: usize = 0;
    var sequence_start_pos: ?usize = null;
    var current_sequence_length: usize = 0;
    var current_header: ?[]u8 = null;

    var line_buffer = std.ArrayList(u8).init(allocator);
    defer line_buffer.deinit();

    while (current_pos < file_size) {
        const aligned_offset = (current_pos / page_size) * page_size;
        const offset_adjustment = current_pos - aligned_offset;

        const remaining = file_size - aligned_offset;
        const adjusted_chunk_size = @min(chunk_size + offset_adjustment, remaining);

        const ptr = try std.posix.mmap(
            null,
            adjusted_chunk_size,
            std.posix.PROT.READ,
            .{ .TYPE = .PRIVATE },
            file.handle,
            aligned_offset,
        );
        defer std.posix.munmap(ptr[0..adjusted_chunk_size]);

        var i: usize = offset_adjustment;
        while (i < adjusted_chunk_size) {
            if (ptr[i] == '>' or i == adjusted_chunk_size - 1) {
                // If we were tracking a sequence, finish it
                if (sequence_start_pos) |start_pos| {
                    const sequence_length = if (ptr[i] == '>')
                        aligned_offset + i - start_pos - current_sequence_length
                    else
                        aligned_offset + i + 1 - start_pos - current_sequence_length;

                    if (current_header) |header| {
                        try sequences.append(FastaSequence{
                            .header = header,
                            .position = start_pos,
                            .length = sequence_length,
                        });
                        current_header = null;
                    }
                }

                if (ptr[i] == '>') {
                    // Start tracking new sequence
                    line_buffer.clearRetainingCapacity();

                    var j = i;
                    while (j < adjusted_chunk_size and ptr[j] != '\n') {
                        try line_buffer.append(ptr[j]);
                        j += 1;
                    }

                    if (j < adjusted_chunk_size) {
                        current_header = try allocator.dupe(u8, line_buffer.items);
                        sequence_start_pos = aligned_offset + j + 1; // Start after newline
                        current_sequence_length = 0;
                        i = j + 1;
                        continue;
                    } else {
                        current_pos = aligned_offset + i;
                        break;
                    }
                }
            } else if (ptr[i] == '\n') {
                current_sequence_length += 1; // Count newlines to subtract from sequence length
            }
            i += 1;
        }

        current_pos = aligned_offset + adjusted_chunk_size;
    }

    return sequences;
}

pub fn main() !void {
    if (.windows == @import("builtin").os.tag) {
        std.debug.print("MMap is not supported in Windows\n", .{});
        return;
    }

    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const args = try std.process.argsAlloc(allocator);
    defer std.process.argsFree(allocator, args);

    if (args.len < 2) {
        print("Usage:\n", .{});
        print("  {s} index <fasta_filename> <output.json>\n", .{args[0]});
        print("  {s} read <fasta_filename> <index.json> <sequence_index> [window_size] [window_start]\n", .{args[0]});
        return error.InvalidArguments;
    }

    const cmd = try Command.fromString(args[1]);
    switch (cmd) {
        .index => {
            if (args.len != 4) {
                print("Usage: {s} index <fasta_filename> <output.json>\n", .{args[0]});
                return error.InvalidArguments;
            }

            const input_filename = args[2];
            const output_filename = args[3];

            const file = try fs.cwd().openFile(input_filename, .{ .mode = .read_only });
            defer file.close();

            var sequences = try scanFastaSequences(file, allocator);
            defer {
                for (sequences.items) |seq| {
                    allocator.free(seq.header);
                }
                sequences.deinit();
            }

            const output_file = try fs.cwd().createFile(output_filename, .{});
            defer output_file.close();

            const writer = output_file.writer();
            try json.stringify(sequences.items, .{ .whitespace = .indent_2 }, writer);

            print("\nFASTA Index Complete:\n", .{});
            print("  Number of sequences: {d}\n", .{sequences.items.len});
            print("  Index written to: {s}\n", .{output_filename});
        },
        .read => {
            if (args.len < 5 or args.len > 7) {
                print("Usage: {s} read <fasta_filename> <index.json> <sequence_index> [window_size] [window_start]\n", .{args[0]});
                return error.InvalidArguments;
            }

            const fasta_filename = args[2];
            const index_filename = args[3];
            const sequence_index = try std.fmt.parseInt(usize, args[4], 10);

            var sequences = try loadIndex(allocator, index_filename);
            defer {
                for (sequences.items) |seq| {
                    allocator.free(seq.header);
                }
                sequences.deinit();
            }

            if (sequence_index >= sequences.items.len) {
                print("Error: Sequence index {d} is out of range (max: {d})\n", .{ sequence_index, sequences.items.len - 1 });
                return error.InvalidSequenceIndex;
            }

            const sequence = sequences.items[sequence_index];
            const window_size = if (args.len > 5)
                try std.fmt.parseInt(usize, args[5], 10)
            else
                sequence.length;
            const window_start = if (args.len > 6)
                try std.fmt.parseInt(usize, args[6], 10)
            else
                0;

            const file = try fs.cwd().openFile(fasta_filename, .{ .mode = .read_only });
            defer file.close();

            const sequence_data = try readSequenceWindow(file, sequence, window_start, window_size);
            defer std.heap.page_allocator.free(sequence_data);

            print("{s}\n", .{sequence_data});
        },
    }
}
