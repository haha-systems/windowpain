# windowpain
Windowpain is a random-access tool for reading genetic sequences from huge FASTA sequence files in milliseconds.

## Overview
Windowpain provides lightning-fast random access to FASTA sequence files by using memory mapping and efficient indexing. It allows you to read specific windows of sequence data without loading the entire file into memory.

## Features
- Memory-mapped file access for optimal performance
- JSON-based indexing of sequence positions
- Random access to sequence windows
- Handles newlines in FASTA files transparently
- Zero-copy reading of sequence data

## Installation

```bash
git clone https://github.com/haha-systems/windowpain.git
cd windowpain
zig build --release=fast
```

## Usage

```bash
./windowpain index <fasta_filename> <index.json>
./windowpain read <fasta_filename> <index.json> <sequence_index> [window_size] [window_start] [--raw]
```

## Examples

### Indexing a FASTA file and reading a sequence window

```bash
./windowpain index test.fasta test.json
./windowpain read test.fasta test.json 0 10
```

This will first index the `test.fasta` file and then read the first 10 characters of the first sequence in `test.fasta` and print it to the console.

### Using raw output to pass to another tool through `stdin`

```bash
./windowpain read test.fasta test.json 0 10 --raw | <other_tool>
```

This will read the first 10 characters of the first sequence in `test.fasta` and pass it to `other_tool`. The `other_tool` can be any tool that accepts raw sequence data on `stdin`.

> If you want to iterate over a sequence, just pass the `position += window_size` to the `read` command's `window_start` argument.

## Notes
- Windowpain uses memory mapping to read sequence data, so it may not be suitable for small files.
- The index file is created by the `index` command and should be used with the `read` command.
- The `read` command can optionally output the raw sequence data without formatting to `stdout`.