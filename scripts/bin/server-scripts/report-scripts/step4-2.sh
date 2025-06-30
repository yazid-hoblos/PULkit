#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <pangenome_dir> <metadata_dir>"
    exit 1
fi

PAN_DIR="$1"
META_DIR="$2"

for meta_file in "$META_DIR"/*_added_lines.txt; do
    filename=$(basename "$meta_file")
    # Extract species name without suffix "_added_lines.txt"
    species=${filename%_added_lines.txt}

    pangenome="$PAN_DIR/$species/pangenome.h5"

    if [ ! -f "$pangenome" ]; then
        echo "Warning: Pangenome file $pangenome not found, skipping."
        continue
    fi

    echo "Running ppanggolin metadata for $species ..."

    ppanggolin metadata -p "$pangenome" --metadata "$meta_file" --source dbcan-merged-corrected --assign families -f

    echo "Done $species."
done

