#!/bin/bash

set -euo pipefail

# Check arguments
if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <input_pangenomes_dir> <output_root_dir>"
    exit 1
fi

input_dir="$1"
output_root="$2"

# Loop over subdirectories in input_dir
for species_dir in "$input_dir"/*; do
    species=$(basename "$species_dir")
    md5_file="$species_dir/genomes_md5sum.tsv"

    # Skip if md5 file doesn't exist
    if [[ ! -f "$md5_file" ]]; then
        echo "No MD5 file for $species, skipping..."
        continue
    fi

    # Create temporary file with genome accessions (skip header)
    accessions_file=$(mktemp)
    tail -n +2 "$md5_file" | cut -f1 > "$accessions_file"

    echo "Downloading genomes for $species..."

    # Set output directory per species
    output_dir="$output_root/$species"
    mkdir -p "$output_dir"

    # Download and unzip
    scripts/tools/datasets download genome accession --inputfile "$accessions_file" \
        --filename "$output_dir/genomes.zip" \
        --include protein,gff3

    unzip -o "$output_dir/genomes.zip" -d "$output_dir"
    rm "$output_dir/genomes.zip"
    rm "$accessions_file"

    echo "Done with $species."
done
