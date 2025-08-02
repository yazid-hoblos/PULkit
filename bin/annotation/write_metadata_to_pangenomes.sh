#!/bin/bash
set -euo pipefail

if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]; then
    echo "Usage: $0 <pangenome_dir> <metadata_dir> [source]"
    exit 1
fi

PAN_DIR="$1"
META_DIR="$2"
SOURCE="$3"

for meta_file in "$META_DIR"/*_added_lines.txt; do
    filename=$(basename "$meta_file")
    species=${filename%_added_lines.txt}

    pangenome="$PAN_DIR/$species/pangenome.h5"

    if [[ ! -f "$pangenome" ]]; then
        echo "Warning: Pangenome file $pangenome not found, skipping."
        continue
    fi

    echo "Running ppanggolin metadata for $species with source '$SOURCE'..."

    ppanggolin metadata -p "$pangenome" --metadata "$meta_file" --source "$SOURCE" --assign families -f

    echo "Done for $species."
done
