#!/bin/bash

#SBATCH -J dbcan
#SBATCH -c 20
#SBATCH --mem=25GB
#SBATCH -t 24:00:00

# Check for argument
if [ -z "$1" ]; then
    echo "Usage: $0 <path_to_genome_folder>"
    exit 1
fi

GENOME_DIR="$1"
OUTBASE="easy-dbcan-results-$(basename "$GENOME_DIR")"

mkdir -p "$OUTBASE"

for genome in "$GENOME_DIR"/*; do
    outdir="${OUTBASE}/$(basename "$genome")-cgc-output"

    if [ ! -d "$outdir" ]; then
        echo "Running dbCAN on $genome..."
        run_dbcan easy_CGC   --db_dir dbcan-source   --input_gff "$genome/genomic.gff"   --input_raw_data "$genome/protein.faa"   --output_dir "$outdir" --mode protein   --gff_type NCBI_prok
    else
        echo "Skipping $genome — output directory already exists."
    fi
done

