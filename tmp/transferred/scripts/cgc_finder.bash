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
OUTBASE="dbcan-results-$(basename "$GENOME_DIR")"

mkdir -p "$OUTBASE"

for genome in "$GENOME_DIR"/*; do
    outdir="${OUTBASE}/$(basename "$genome")-cgc-output"

    if [ ! -d "$outdir" ]; then
        echo "Running dbCAN on $genome..."
        run_dbcan --db_dir /env/cns/proj/agc/scratch_microscope/Data/DBCAN/899 \
                  "$genome/protein.faa" protein \
                  -c "$genome/genomic.gff" \
                  --out_dir "$outdir" \
                  --dia_cpu 20 --hmm_cpu 20 --tf_cpu 20 --stp_cpu 20 --cgc_sig_genes tp+tf --cgc_dis 1
    else
        echo "Skipping $genome — output directory already exists."
    fi
done

