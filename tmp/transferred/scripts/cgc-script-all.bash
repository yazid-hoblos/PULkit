#!/bin/bash

#SBATCH -J dbcan
#SBATCH -c 20
#SBATCH --mem 25GB
#SBATCH -t 24:00:00

for genome in systems_detection/Bt_systems/s__Bacteroides_thetaiotaomicron/Bt_genomes/*; do
    outdir="dbcan_predictions/Bt_tp-tf/$(basename "$genome")-cgc-output"
    
    if [ ! -d "$outdir" ]; then
        echo "Running dbCAN on $genome..."
        run_dbcan --db_dir /env/cns/proj/agc/scratch_microscope/Data/DBCAN/899 "$genome/protein.faa" protein -c "$genome/genomic.gff" --out_dir "$outdir" --dia_cpu 20 --hmm_cpu 20 --tf_cpu 20 --stp_cpu 20 --cgc_sig_genes tp+tf
    else
        echo "Skipping $genome — output directory already exists."
    fi
done

