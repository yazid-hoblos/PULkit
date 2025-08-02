#!/bin/bash

set -euo pipefail

# Check input arguments
if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <pangenomes_dir> <output_dir>"
    exit 1
fi

PANGENOME_DIR="$1"
OUTPUT_DIR="$2"

mkdir -p "$OUTPUT_DIR"

# Find all .h5 files under the input directory and run `ppanggolin fasta` on each
find "$PANGENOME_DIR" -name "*.h5" | while read -r H5FILE; do
    echo "Processing $H5FILE"

    PARENT_DIR=$(basename "$(dirname "$H5FILE")")
    OUT_SUBDIR="${OUTPUT_DIR}/${PARENT_DIR}"

    mkdir -p "$OUT_SUBDIR"

    ppanggolin fasta -p "$H5FILE" --output "$OUT_SUBDIR" --prot_families all -f
done
