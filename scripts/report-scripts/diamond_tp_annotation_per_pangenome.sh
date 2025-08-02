#!/bin/bash

set -euo pipefail

# Usage check
if [[ $# -lt 2 || $# -gt 3 ]]; then
    echo "Usage: $0 <input_pangenome_families_dir> <output_dir> [diamond_db (optional)]"
    exit 1
fi

INPUT_DIR="$1"        # e.g., ppanggolin_prot_output
OUTPUT_DIR="$2"       # e.g., diamond_results
DIAMOND_DB="${3:-/env/cns/proj/agc/scratch_microscope/Data/DBCAN/899/tcdb.dmnd}"  # default path

mkdir -p "$OUTPUT_DIR"

for dir in "$INPUT_DIR"/*; do
    species=$(basename "$dir")
    query="$dir/all_protein_families.faa"
    output="${OUTPUT_DIR}/${species}_tp_results.txt"

    if [[ -f "$query" ]]; then
        echo "Running diamond blastp for $species..."
        diamond blastp -q "$query" -d "$DIAMOND_DB" -o "$output" -e 1e-10
    else
        echo "Warning: Query file $query not found for $species"
    fi
done
