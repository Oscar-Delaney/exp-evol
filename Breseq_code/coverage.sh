#!/bin/bash

# Set the path to the directory with .gd files
ANNOTATED_GD_DIR="/home/s4528540/liv/breseq_annotated_gd"

# Set the path to the output file
OUTPUT_FILE="/home/s4528540/liv/combined_MC_data.csv"

# Print header
echo "Filename,RefGenome,Start,End,StartRange,EndRange" > "$OUTPUT_FILE"

# Extract the "MC" rows from each file and add to the output file
for GD_FILE in "$ANNOTATED_GD_DIR"/*.gd; do
    FILENAME=$(basename "$GD_FILE" .gd)
    awk -v filename="$FILENAME" '/^MC/ {print filename "," $4 "," $5 "," $6 "," $7 "," $8}' "$GD_FILE" >> "$OUTPUT_FILE"
done
