#!/bin/bash

set -euo pipefail

for species_dir in Prev_pan/*; do
    species=$(basename "$species_dir")
    md5_file="$species_dir/genomes_md5sum.tsv"

    # Skip if md5 file doesn't exist
    if [[ ! -f "$md5_file" ]]; then
        echo "No MD5 file for $species, skipping..."
        continue
    fi

    # Create temporary input file with genome accessions (skip header)
    accessions_file=$(mktemp)
    tail -n +2 "$md5_file" | cut -f1 > "$accessions_file"

    echo "Downloading genomes for $species..."

    # Define download and output directory
    output_dir="genomes_by_species/$species"
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

