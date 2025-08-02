#!/bin/bash

# Check for arguments
if [ -z "$1" ]; then
    echo "Usage: $0 <path_to_genomes_folder>"
    exit 1
fi

GENOMES_DIR="$1"
OUTBASE="dbcan_cgc_finder_$(basename "$GENOMES_DIR")"

mkdir -p "$OUTBASE"

for genome in "$GENOMES_DIR"/*; do
    outdir="${OUTBASE}/$(basename "$genome")-cgc-output"

    if [ ! -d "$outdir" ]; then
        echo "Running dbCAN on $genome..."
        run_dbcan easy_CGC   --db_dir dbcan-source   --input_gff "$genome/genomic.gff"   --input_raw_data "$genome/protein.faa"   --output_dir "$outdir" --mode protein   --gff_type NCBI_prok
    else
        echo "Skipping $genome — output directory already exists."
    fi
done

