#!/bin/bash

METADATA_DIR="metadata_output"
OUTPUT_DIR="final_output"

TF_FILE="tfs"    # file with list of IDs or patterns for TF
STP_FILE="stps"  # file with list of IDs or patterns for STP

for meta_file in "$METADATA_DIR"/*; do
    filename=$(basename "$meta_file")
    
    # Extract species prefix by removing suffix (same as before)
    species=$(echo "$filename" | sed 's/_families_metadata_from_CAZ-TF-STP.tsv$//')
    
    target_file="$OUTPUT_DIR/${species}_added_lines.txt"
    
    if [ ! -f "$target_file" ]; then
        echo "Warning: target file $target_file not found, skipping"
        continue
    fi
    
    # Grep TF IDs from metadata
    grep -Ff "$TF_FILE" "$meta_file" | cut -f1 | sort | uniq > tf_ids.txt
    # Grep STP IDs from metadata
    grep -Ff "$STP_FILE" "$meta_file" | cut -f1 | sort | uniq > stp_ids.txt
    
    # Append TF lines
    while read -r id; do
        echo -e "${id}\tTF\tdbcan-merged\t0"
    done < tf_ids.txt >> "$target_file"
    
    # Append STP lines
    while read -r id; do
        echo -e "${id}\tSTP\tdbcan-merged\t0"
    done < stp_ids.txt >> "$target_file"
    
    rm tf_ids.txt stp_ids.txt
done

