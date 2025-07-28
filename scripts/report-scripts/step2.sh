#!/bin/bash

DIAMOND_DIR="diamond_results"
FINAL_DIR="final_output"

for tp_file in "$DIAMOND_DIR"/*_tp_results.txt; do
    filename=$(basename "$tp_file")
    species=${filename%_tp_results.txt}
    
    target_file="$FINAL_DIR/${species}_added_lines.txt"
    
    if [ ! -f "$target_file" ]; then
        echo "Warning: target file $target_file not found, skipping"
        continue
    fi
    
    # Extract unique IDs from first column of diamond file
    cut -f1 "$tp_file" | sort | uniq > tp_ids.txt
    
    # Append new lines with TP instead of CAZyme
    while read -r id; do
        echo -e "${id}\tTC\tdbcan-merged\t0"
    done < tp_ids.txt >> "$target_file"
    
    rm tp_ids.txt
done

