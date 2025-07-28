mkdir -p run_dbcan_output

for dir in ppanggolin_prot_output/*; do
    species=$(basename "$dir")
    protein_file="$dir/all_protein_families.faa"
    out_dir="run_dbcan_output/$species"

    if [[ -d "$out_dir" ]]; then
        echo "Skipping $species — output already exists."
        continue
    fi

    if [[ -f "$protein_file" ]]; then
        echo "Running run_dbcan on $species..."
        run_dbcan --db_dir /env/cns/proj/agc/scratch_microscope/Data/DBCAN/899 "$protein_file" protein --out_dir "$out_dir"
    else
        echo "Warning: No all_protein_families.faa file found in $dir"
    fi
done

