#!/bin/bash

#SBATCH -J dbcan
#SBATCH -c 20
#SBATCH --mem=25GB
#SBATCH -t 24:00:00

BASE_DIR="genomes_by_species"
DBCAN_DB="dbcan-source"
OUTBASE="easy-dbcan-results"

mkdir -p "$OUTBASE"

for species_dir in "$BASE_DIR"/*; do
    species_name=$(basename "$species_dir")
    species_outdir="${OUTBASE}/${species_name}"
    mkdir -p "$species_outdir"

    for genome_dir in "$species_dir"/ncbi_dataset/data/*; do
        if [ -d "$genome_dir" ]; then
            genome_id=$(basename "$genome_dir")
            outdir="${species_outdir}/${genome_id}-cgc-output"

            if [ ! -d "$outdir" ]; then
                echo "Running dbCAN on ${species_name}/${genome_id}..."
                run_dbcan easy_CGC \
                    --db_dir "$DBCAN_DB" \
                    --input_gff "${genome_dir}/genomic.gff" \
                    --input_raw_data "${genome_dir}/protein.faa" \
                    --output_dir "$outdir" \
                    --mode protein \
                    --gff_type NCBI_prok
            else
                echo "Skipping ${species_name}/${genome_id} — output already exists."
            fi
        fi
    done
done

