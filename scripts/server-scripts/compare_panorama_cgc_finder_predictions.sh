#!/bin/bash

set -euo pipefail

if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <projection_dir> <mapping_file> <cgc_base_dir> <output_tsv>"
    exit 1
fi

PROJ_DIR="$1"
MAPPING_FILE="$2"
CGC_DIR="$3"
OUTFILE="$4"

TMP_ALPHA=$(mktemp)
TMP_BETA=$(mktemp)

echo -e "basename\tprojection_count\tmapped_count\tCGC_finder_matched\tCGC_finder_total" > "$OUTFILE"

for proj_file in "$PROJ_DIR"/*.tsv; do
    bname=$(basename "$proj_file" .tsv)

    tail -n +2 "$proj_file" | cut -f9 > "$TMP_ALPHA"
    awk 'NR==FNR {ids[$0]; next} $1 in ids' "$TMP_ALPHA" "$MAPPING_FILE" | cut -f2 > "$TMP_BETA"

    len_alpha=$(wc -l < "$TMP_ALPHA")
    len_beta=$(wc -l < "$TMP_BETA")

    CGC_FILE="${CGC_DIR}/${bname}-cgc-output/cgc_standard.out"

    if [[ -f "$CGC_FILE" ]]; then
        len_cgc=$(grep -Ff "$TMP_BETA" "$CGC_FILE" | wc -l)
        len_cgc_file=$(wc -l < "$CGC_FILE")
    else
        len_cgc=0
        len_cgc_file=0
    fi

    echo -e "${bname}\t${len_alpha}\t${len_beta}\t${len_cgc}\t${len_cgc_file}" >> "$OUTFILE"
done

rm "$TMP_ALPHA" "$TMP_BETA"
