#!/bin/bash

PAN_DIR="Prev_pan"
META_DIR="final_output"

for meta_file in "$META_DIR"/*_added_lines.txt; do
    filename=$(basename "$meta_file")
    # Extract species name without suffix "_added_lines.txt"
    species=${filename%_added_lines.txt}
    
    pangenome="$PAN_DIR/$species/pangenome.h5"
    
    if [ ! -f "$pangenome" ]; then
        echo "Warning: Pangenome file $pangenome not found, skipping."
        continue
    fi
    
    echo "Running ppanggolin metadata for $species ..."
    
    ppanggolin metadata -p "$pangenome" --metadata "$meta_file" --source dbcan-merged --assign families
    
    echo "Done $species."
done

