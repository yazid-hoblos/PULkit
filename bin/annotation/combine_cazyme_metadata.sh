#!/bin/bash

set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 <metadata_dir> <dbcan_dir> <output_dir>"
    exit 1
fi

METADATA_DIR="$1"
DBCAN_DIR="$2"
OUTPUT_DIR="$3"

mkdir -p "$OUTPUT_DIR"

for meta_file in "$METADATA_DIR"/*; do
    filename=$(basename "$meta_file")

    # Remove suffix starting with _families_metadata and anything after
    species=$(echo "$filename" | sed 's/_families_metadata.*$//')

    overview_file="$DBCAN_DIR/$species/overview.txt"

    if [[ ! -f "$overview_file" ]]; then
        echo "Warning: overview file $overview_file not found, skipping"
        continue
    fi

    hmmer_ids=$(mktemp)
    overview_ids=$(mktemp)
    missing_ids=$(mktemp)

    trap 'rm -f "$hmmer_ids" "$overview_ids" "$missing_ids"' EXIT

    grep -Fw CAZyme "$meta_file" | cut -f1 | sort -u > "$hmmer_ids"
    tail -n +2 "$overview_file" | cut -f1 | sort -u > "$overview_ids"
    comm -23 "$overview_ids" "$hmmer_ids" > "$missing_ids"

    out_file="$OUTPUT_DIR/${species}_added_lines.txt"

    echo -e "families\tprotein_name\tscore" > "$out_file"

    while IFS= read -r id; do
        echo -e "${id}\tCAZyme\t0"
    done < "$hmmer_ids" >> "$out_file"

    while IFS= read -r id; do
        echo -e "${id}\tCAZyme\t0"
    done < "$missing_ids" >> "$out_file"

    trap - EXIT
    rm -f "$hmmer_ids" "$overview_ids" "$missing_ids"
done
