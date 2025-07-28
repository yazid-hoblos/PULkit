#!/bin/bash

mkdir -p metadata-final

for pangenome_file in Prev_pan/*/pangenome.h5; do
    species_dir=$(dirname "$pangenome_file")
    species=$(basename "$species_dir")

    echo "Processing $species..."

    # Run ppanggolin write_metadata
    ppanggolin write_metadata -p "$pangenome_file" -s dbcan-like -o test -f

    # Move all resulting families* files, prefix with species name
    for file in test/families*; do
        if [[ -f "$file" ]]; then
            mv -f "$file" "metadata-final/${species}_$(basename "$file")"
        fi
    done
done

