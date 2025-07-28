#!/bin/bash

echo "\\begin{tabular}{|l|p{10cm}|}"
echo "\\hline"
echo "Species & Strains (GCF IDs) \\\\"
echo "\\hline"

for dir in Prev_pan/s__*; do
    species=$(basename "$dir" | sed 's/^s__//; s/_/ /g')
    strains=$(cut -f1 "$dir"/*md5* | grep '^GCF' | paste -sd ", " -)
    echo "$species & $strains \\\\"
    echo "\\hline"
done

echo "\\end{tabular}"

