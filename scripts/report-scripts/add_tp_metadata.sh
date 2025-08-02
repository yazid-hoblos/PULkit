#!/bin/bash
set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <diamond_dir> <output_dir>"
    exit 1
fi

DIAMOND_DIR="$1"
OUT_DIR="$2"

mkdir -p "$OUT_DIR"

process_file() {
    local tp_file="$1"
    local species="$2"
    local target_file="$3"

    if [[ ! -f "$target_file" ]]; then
        echo "Warning: target file $target_file not found, skipping"
        return
    fi

    local tmp_ids
    tmp_ids=$(mktemp)

    cut -f1 "$tp_file" | sort -u > "$tmp_ids"

    while IFS= read -r id; do
        echo -e "${id}\tTC\t0"
    done < "$tmp_ids" >> "$target_file"

    rm -f "$tmp_ids"
}

for tp_file in "$DIAMOND_DIR"/*_tp_results.txt; do
    [ -e "$tp_file" ] || { echo "No *_tp_results.txt files found in $DIAMOND_DIR"; break; }

    filename=$(basename "$tp_file")
    species=${filename%_tp_results.txt}
    target_file="$OUT_DIR/${species}_added_lines.txt"

    process_file "$tp_file" "$species" "$target_file"
done
