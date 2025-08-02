#!/bin/bash

set -euo pipefail

# Check arguments
if [[ $# -ne 3 ]]; then
    echo "Usage: $0 <input_extracted_cazy_puls> <genomes_dir> <output_dir>"
    exit 1
fi

input_dir="$1"
genomes_dir="$2"
output_dir="$3"

mkdir -p "$output_dir"

for species_file in "$input_dir"/*.csv; do
  species_name=$(basename "$species_file" _cazy_puls.csv)
  gff_dir="$genomes_dir/s__${species_name}/ncbi_dataset/data/"

  if [ ! -d "$gff_dir" ]; then
    echo "Directory $gff_dir does not exist. Skipping $species_name."
    continue
  fi

  cut -d, -f4 "$species_file" | tr ';' '\n' > accessions.tmp
  dos2unix accessions.tmp  

  awk -v OFS='\t' '
    FNR==NR { ids[$1]; next }
    $3 == "gene" {
      old = ""; locus = "";
      if (match($9, /old_locus_tag=([^;]+)/, a)) old = a[1];
      if (match($9, /locus_tag=([^;]+)/, b)) locus = b[1];
      if (old in ids) print $1, $4, $5, $7, locus, old;
    }
  ' accessions.tmp "$gff_dir"/*/*.gff > "$output_dir/${species_name}_genes.tsv"

  echo "Processed $species_name -> $output_dir/${species_name}_genes.tsv"
done

rm accessions.tmp
