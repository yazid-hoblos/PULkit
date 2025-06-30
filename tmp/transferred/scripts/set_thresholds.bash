#!/bin/bash

if [ -z "$1" ]; then
  echo "Usage: $0 <input_tsv_file>"
  exit 1
fi

awk -F'\t' 'BEGIN {OFS="\t"} 
NR==1 {
    print
} 
NR>1 {
    for (i=1; i<=NF; i++) {
        if ($i == "") {
            if (i == 7) $i = "25";
            else if (i == 8 || i == 9) $i = "1e-5";
            else if (i == 10 || i == 11) $i = "0.5";
            else if (i == 5) {
                val = $1
                sub(/\.hmm$/, "", val)
                $i = val
            }
        }
    }
    print
}' "$1" > "${1%.tsv}_filled.tsv"

