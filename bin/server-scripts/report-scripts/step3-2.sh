#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <metadata_dir> <output_dir>"
    exit 1
fi

METADATA_DIR="$1"
OUTPUT_DIR="$2"
TF_FILE="tfs"
STP_FILE="stps"

for meta_file in "$METADATA_DIR"/*; do
    filename=$(basename "$meta_file")

    # Extract species name by removing the suffix starting with _families_metadata
    species=$(echo "$filename" | sed 's/_families_metadata.*$//')

    target_file="$OUTPUT_DIR/${species}_added_lines.txt"

    if [ ! -f "$target_file" ]; then
        echo "Warning: target file $target_file not found, skipping"
        continue
    fi

    # Grep TF IDs from metadata
    grep -Ff "$TF_FILE" "$meta_file" | cut -f1 | sort | uniq > tf_ids.txt
    # Grep STP IDs from metadata
    grep -Ff "$STP_FILE" "$meta_file" | cut -f1 | sort | uniq > stp_ids.txt

    # Append TF lines
    while read -r id; do
        echo -e "${id}\tTF\t0"
    done < tf_ids.txt >> "$target_file"

    # Append STP lines
    while read -r id; do
        echo -e "${id}\tSTP\t0"
    done < stp_ids.txt >> "$target_file"

    rm tf_ids.txt stp_ids.txt
done

