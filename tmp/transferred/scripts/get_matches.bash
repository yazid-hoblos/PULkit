#!/bin/bash

for g in panorama-systems/detected_Bt_systems/Bt_t1_s/s__Bacteroides_thetaiotaomicron/CAZ-TF-TC-STP-simple/projection/*; do
    tail -n +2 "$g" | cut -f9 > alpha
    awk 'NR==FNR{ids[$0]; next} $1 in ids' alpha uniq_ids_mapping | cut -f2 > beta
    b=$(basename "$g" .tsv)

    len_alpha=$(wc -l < alpha)
    len_beta=$(wc -l < beta)

    cgc_file="dbcan-results/${b}-cgc-output/cgc_standard.out"

    if [[ -f "$cgc_file" ]]; then
        len_cgc=$(grep -Ff beta "$cgc_file" | wc -l)
        len_cgc_file=$(wc -l < "$cgc_file")
    else
        len_cgc=0
        len_cgc_file=0
    fi

    echo -e "${b}\t${len_alpha}\t${len_beta}\t${len_cgc}\t${len_cgc_file}"
done


