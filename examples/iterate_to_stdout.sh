#!/bin/bash

# Check if all required arguments are provided
if [ "$#" -lt 4 ]; then
    echo "Usage: $0 <fasta_file> <index_file> <sequence_index> <window_size>"
    exit 1
fi

FASTA_FILE=$1
INDEX_FILE=$2
SEQ_INDEX=$3
WINDOW_SIZE=$4
POSITION=0

# Loop until no more output is received
while true; do
    # Read the window at current position and redirect to /dev/null
    OUTPUT=$(./zig-out/bin/windowpain read "$FASTA_FILE" "$INDEX_FILE" "$SEQ_INDEX" "$WINDOW_SIZE" "$POSITION" --raw)
    
    # Check if we got any output
    if [ -z "$OUTPUT" ]; then
        break
    fi

    # Print the output to stdout
    echo "$OUTPUT"
    
    # Increment position by window size
    POSITION=$((POSITION + WINDOW_SIZE))
done