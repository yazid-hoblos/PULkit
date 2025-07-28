#!/bin/bash

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <metadata_dir> <run_dir> <output_dir>"
    exit 1
fi

METADATA_DIR="$1"
RUNDIR="$2"
OUTPUT_DIR="$3"

mkdir -p "$OUTPUT_DIR"

for meta_file in "$METADATA_DIR"/*; do
    filename=$(basename "$meta_file")

    # Remove suffix starting with _families_metadata and anything after
    species=$(echo "$filename" | sed 's/_families_metadata.*$//')

    overview_file="$RUNDIR/$species/overview.txt"

    if [ ! -f "$overview_file" ]; then
        echo "Warning: overview file $overview_file not found, skipping"
        continue
    fi

    # Extract IDs from metadata file, skipping header, and keep full lines
    grep -Ff cazymes "$meta_file" | cut -f1 | sort | uniq > a_ids.txt

    # Extract IDs from overview file, skipping header, first column only
    tail -n +2 "$overview_file" | cut -f1 | sort | uniq > overview_ids.txt

    # Find IDs in overview not in metadata
    comm -23 overview_ids.txt a_ids.txt > missing_ids.txt

    out_file="$OUTPUT_DIR/${species}_added_lines.txt"

    # Write header
    echo -e "families\tprotein_name\tscore" > "$out_file"

    # Write all metadata IDs (from a_ids.txt)
    while read -r id; do
        echo -e "${id}\tCAZyme\t0"
    done < a_ids.txt >> "$out_file"

    # Write the missing overview IDs
    while read -r id; do
        echo -e "${id}\tCAZyme\t0"
    done < missing_ids.txt >> "$out_file"

    # Clean temporary files
    rm a_ids.txt overview_ids.txt missing_ids.txt
done

