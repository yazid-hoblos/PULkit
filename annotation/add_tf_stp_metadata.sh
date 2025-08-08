#!/bin/bash
set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <metadata_dir> <output_dir>"
    exit 1
fi

METADATA_DIR="$1"
OUTPUT_DIR="$2"

mkdir -p "$OUTPUT_DIR"

for meta_file in "$METADATA_DIR"/*; do
    filename=$(basename "$meta_file")
    species=${filename%%_families_metadata*}  # remove suffix starting with _families_metadata

    target_file="$OUTPUT_DIR/${species}_added_lines.txt"

    if [[ ! -f "$target_file" ]]; then
        echo "Warning: target file $target_file not found, skipping"
        continue
    fi

    tmp_tf=$(mktemp)
    tmp_stp=$(mktemp)

    grep -Fw TF "$meta_file" | cut -f1 | sort -u > "$tmp_tf"
    grep -Fw STP "$meta_file" | cut -f1 | sort -u > "$tmp_stp"

    while IFS= read -r id; do
        echo -e "${id}\tTF\t0"
    done < "$tmp_tf" >> "$target_file"

    while IFS= read -r id; do
        echo -e "${id}\tSTP\t0"
    done < "$tmp_stp" >> "$target_file"

    rm -f "$tmp_tf" "$tmp_stp"
done
