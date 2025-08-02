#!/bin/bash

# Usage: ./create_gff_list.sh /path/to/genomes_dir > gff_list.tsv

GENOME_DIR="$1"

if [ -z "$GENOME_DIR" ]; then
    echo "Usage: $0 /path/to/genomes_dir" >&2
    exit 1
fi

for genome_path in "$GENOME_DIR"/*; do
    genome_name=$(basename "$genome_path")
    
    # Support .gff or .gff.gz
    gff_file=$(find "$genome_path" -maxdepth 1 -type f \( -name "*.gff" -o -name "*.gff.gz" \) | head -n 1)

    if [[ -n "$gff_file" ]]; then
        echo -e "${genome_name}\t${gff_file}"
    else
        echo "[WARNING] No GFF found for ${genome_name}" >&2
    fi
done

