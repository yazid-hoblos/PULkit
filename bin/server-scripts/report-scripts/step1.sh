#!/bin/bash

METADATA_DIR="metadata_output"
RUNDIR="run_dbcan_output"
OUTPUT_DIR="final_output"

mkdir -p "$OUTPUT_DIR"

for meta_file in "$METADATA_DIR"/*; do
    filename=$(basename "$meta_file")
    
    # Remove suffix to get species folder name exactly
    species=$(echo "$filename" | sed 's/_families_metadata_from_CAZ-TF-STP.tsv$//')
    
    overview_file="$RUNDIR/$species/overview.txt"
    
    if [ ! -f "$overview_file" ]; then
        echo "Warning: overview file $overview_file not found, skipping"
        continue
    fi
    
    # Extract IDs from metadata file, assuming IDs in first column and skipping header
    tail -n +2 "$meta_file" | cut -f1 | sort | uniq > a_ids.txt
    
    # Extract IDs from overview file, skipping header, first column only
    tail -n +2 "$overview_file" | cut -f1 | sort | uniq > overview_ids.txt
    
    # Find IDs in overview not in metadata
    comm -23 overview_ids.txt a_ids.txt > missing_ids.txt
    
    out_file="$OUTPUT_DIR/${species}_added_lines.txt"
    
    # Write header
    echo -e "families\tprotein_name\tsource\tscore" > "$out_file"
    
    # Write the missing ID lines with fixed format
    while read -r id; do
        echo -e "${id}\tCAZyme\tdbcan-merged\t0"
    done < missing_ids.txt >> "$out_file"
    
    # Clean temporary files
    rm a_ids.txt overview_ids.txt missing_ids.txt
done

