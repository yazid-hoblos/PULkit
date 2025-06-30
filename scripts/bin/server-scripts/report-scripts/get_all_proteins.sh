#!/bin/bash

# Adjust this path to the root directory containing your .h5 files
PANGENOME_DIR="Prev_pan"
OUTPUT_DIR="ppanggolin_prot_output"

mkdir -p "$OUTPUT_DIR"

# Find all .h5 files under PANGENOME_DIR and run ppanggolin for each
find "$PANGENOME_DIR" -name "*.h5" | while read -r H5FILE; do
    echo "Processing $H5FILE"
    
    # Get basename without path and extension for naming output subfolder
    PARENT_DIR=$(basename "$(dirname "$H5FILE")")
    OUT_SUBDIR="${OUTPUT_DIR}/${PARENT_DIR}"

    mkdir -p "$OUT_SUBDIR"

    ppanggolin fasta -p "$H5FILE" --output "$OUT_SUBDIR" --prot_families all -f
done

