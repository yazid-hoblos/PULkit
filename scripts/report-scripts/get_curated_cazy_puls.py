#!/bin/bash

if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <input_cazy_puls> <curated_puls>"
  exit 1
fi

INPUT_DIR="$1"
OUTPUT_DIR="$2"

for species_file in "$INPUT_DIR"/*.csv; do
  species_name=$(basename "$species_file" _PUL_proteins.csv)
  genes_file="${OUTPUT_DIR}/${species_name}_genes.tsv"
  output_file="${OUTPUT_DIR}/${species_name}_merged.tsv"

  if [ ! -f "$genes_file" ]; then
    echo "Genes file $genes_file not found, skipping $species_name."
    continue
  fi

  if [ -f "$output_file" ]; then
    echo "Output file $output_file already exists, skipping to avoid overwrite."
    continue
  fi

  echo "Merging: $species_file + $genes_file → $output_file"
  python3 map_cazy_old_tags.py "$genes_file" "$species_file" > "$output_file"

  if [ $? -eq 0 ]; then
    echo "Successfully merged: $species_name"
  else
    echo "Error occurred during merging for $species_name"
    rm -f "$output_file"
  fi
done
