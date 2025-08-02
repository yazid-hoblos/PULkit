#!/bin/bash

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <pangenomes_base_dir>"
    exit 1
fi

BASE_DIR="$1"

echo "\\begin{tabular}{|l|p{10cm}|}"
echo "\\hline"
echo "Species & Strains (GCF IDs) \\\\"
echo "\\hline"

for dir in "$BASE_DIR"/s__*; do
    species=$(basename "$dir" | sed 's/^s__//; s/_/ /g')
    strains=$(cut -f1 "$dir"/*md5* | grep '^GCF' | paste -sd ", " -)
    echo "$species & $strains \\\\"
    echo "\\hline"
done

echo "\\end{tabular}"
