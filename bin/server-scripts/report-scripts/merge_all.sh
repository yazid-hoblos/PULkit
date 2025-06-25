for species_file in Prev_CAZy/*.csv; do
  species_name=$(basename "$species_file" _PUL_proteins.csv)
  genes_file="cazy_curated/${species_name}_genes.tsv"
  output_file="cazy_curated/${species_name}_merged.tsv"

  if [ ! -f "$genes_file" ]; then
    echo "Genes file $genes_file not found, skipping $species_name."
    continue
  fi

  if [ -f "$output_file" ]; then
    echo "Output file $output_file already exists, skipping to avoid overwrite."
    continue
  fi

  python synced/merge.py "$genes_file" "$species_file" > "$output_file"
  echo "Merged for $species_name saved to $output_file"
done

