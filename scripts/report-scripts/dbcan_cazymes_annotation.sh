#!/bin/bash

set -euo pipefail

# Check for required arguments
if [[ $# -ne 3 ]]; then
    echo "Usage: $0 <input_pangenome_families_dir> <output_dir> <dbcan_db_dir>"
    exit 1
fi

input_dir="$1"
output_root="$2"
dbcan_db="$3"

mkdir -p "$output_root"

for dir in "$input_dir"/*; do
    species=$(basename "$dir")
    protein_file="$dir/all_protein_families.faa"
    out_dir="$output_root/$species"

    if [[ -d "$out_dir" ]]; then
        echo "Skipping $species — output already exists."
        continue
    fi

    if [[ -f "$protein_file" ]]; then
        echo "Running run_dbcan on $species..."
        run_dbcan --db_dir "$dbcan_db" "$protein_file" protein --out_dir "$out_dir"
    else
        echo "Warning: No all_protein_families.faa file found in $dir"
    fi
done
