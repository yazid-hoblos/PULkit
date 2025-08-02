#!/bin/bash

if [ -z "$1" ]; then
  echo "Usage: $0 <input_tsv_file>"
  exit 1
fi

INPUT="$1"
OUTPUT="${1%.tsv}_with_thresholds.tsv"

awk -F'\t' -v OFS='\t' '
BEGIN {
    # Define default values for missing fields
    defaults[7] = "25"
    defaults[8] = "1e-5"
    defaults[9] = "1e-5"
    defaults[10] = "0.5"
    defaults[11] = "0.5"
}
NR == 1 {
    print; next
}
NR > 1 {
    for (i = 1; i <= NF; i++) {
        if ($i == "") {
            if (i in defaults) {
                $i = defaults[i]
            } else if (i == 5) {
                val = $1
                sub(/\.hmm$/, "", val)
                $i = val
            }
        }
    }
    print
}
' "$INPUT" > "$OUTPUT"
