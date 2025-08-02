#!/bin/bash

set -euo pipefail

# Usage check
if [[ $# -lt 2 || $# -gt 4 ]]; then
    echo "Usage: $0 <pangenome_base_dir> <final_output_dir> [-s source (optional)]"
    exit 1
fi

PANGENOME_BASE="$1"
FINAL_OUTPUT="$2"
TEMP_OUTPUT="temp_metadata_output"

# Handle optional -s argument
PPANGGOLIN_SOURCE_FLAG=""
if [[ $# -eq 4 && "$3" == "-s" ]]; then
    PPANGGOLIN_SOURCE_FLAG="-s $4"
fi

mkdir -p "$FINAL_OUTPUT"

for pangenome_file in "$PANGENOME_BASE"/*/pangenome.h5; do
    species_dir=$(dirname "$pangenome_file")
    species=$(basename "$species_dir")

    echo "Processing $species..."

    # Run ppanggolin write_metadata with or without -s
    if [[ -n "$PPANGGOLIN_SOURCE_FLAG" ]]; then
        ppanggolin write_metadata -p "$pangenome_file" $PPANGGOLIN_SOURCE_FLAG -o "$TEMP_OUTPUT" -f
    else
        ppanggolin write_metadata -p "$pangenome_file" -o "$TEMP_OUTPUT" -f
    fi

    # Move files with species prefix
    for file in "$TEMP_OUTPUT"/families*; do
        if [[ -f "$file" ]]; then
            mv -f "$file" "${FINAL_OUTPUT}/${species}_$(basename "$file")"
        fi
    done
done

# Clean up temp directory
rm -rf "$TEMP_OUTPUT"
