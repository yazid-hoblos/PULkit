output_dir="cazy_curated"
mkdir -p "$output_dir"

for species_file in Prev_CAZy/*.csv; do
  species_name=$(basename "$species_file" _PUL_proteins.csv)
  gff_dir="genomes_by_species/s__${species_name}//ncbi_dataset/data/"

  if [ ! -d "$gff_dir" ]; then
    echo "Directory $gff_dir does not exist. Skipping $species_name."
    continue
  fi

  cut -d, -f4 "$species_file" | tr ';' '\n' > a
  dos2unix a  

  awk -v OFS='\t' '
    FNR==NR { ids[$1]; next }
    $3 == "gene" {
      old = ""; locus = "";
      if (match($9, /old_locus_tag=([^;]+)/, a)) old = a[1];
      if (match($9, /locus_tag=([^;]+)/, b)) locus = b[1];
      if (old in ids) print $1, $4, $5, $7, locus, old;
    }
  ' a "$gff_dir"/*/*.gff > "$output_dir/${species_name}_genes.tsv"

  echo "Processed $species_name -> $output_dir/${species_name}_genes.tsv"
done

rm a

