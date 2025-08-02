#!/bin/bash
set -euo pipefail

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <directory_with_gffs> <output_file>"
    exit 1
fi

INPUT_DIR="$1"
OUT_FILE="{$2:-ids_mapping.txt}"  # Default output file if not provided

> "$OUT_FILE"  # clear or create output file

find "$INPUT_DIR" -type f -name "genomic.gff" | while read -r gff_file; do
    grep -P "\tCDS\t" "$gff_file" | \
    awk '{
        match($0, /locus_tag=([^;]+)/, a);
        match($0, /protein_id=([^;]+)/, b);
        if (a[1] && b[1]) print a[1] "\t" b[1];
    }' >> "$OUT_FILE"
done
