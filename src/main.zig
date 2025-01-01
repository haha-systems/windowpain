const root = @import("root.zig");
const std = @import("std");
const fs = std.fs;
const print = std.debug.print;
const json = std.json;
const Command = root.Command;
const FastaSequence = root.FastaSequence;
const loadIndex = root.loadIndex;
const readSequenceWindow = root.readSequenceWindow;

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
        print("  {s} index <fasta_filename> <index.json>\n", .{args[0]});
        print("  {s} read <fasta_filename> <index.json> <sequence_index> [window_size] [window_start]\n", .{args[0]});
        return error.InvalidArguments;
    }

    var timer = try std.time.Timer.start();
    const start_time = timer.lap();

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

            const end_time = timer.lap();
            const elapsed_ns = end_time - start_time;
            print("\nFASTA Index Complete:\n", .{});
            print("  Number of sequences: {d}\n", .{sequences.items.len});
            print("  Index written to: {s}\n", .{output_filename});
            print("  Time taken: {d}ms\n", .{elapsed_ns / std.time.ns_per_ms});
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

            const end_time = timer.lap();
            const elapsed_ns = end_time - start_time;
            print("\nSequence read complete:\n", .{});
            print("  Sequence length: {d}\n", .{sequence_data.len});
            print("  Time taken: {d}ms\n", .{elapsed_ns / std.time.ns_per_ms});
            print("\nSequence data:\n{s}\n", .{sequence_data});
        },
    }
}
