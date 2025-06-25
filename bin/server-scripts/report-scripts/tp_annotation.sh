mkdir -p diamond_results

for dir in ppanggolin_prot_output/*; do
    species=$(basename "$dir")
    query="$dir/all_protein_families.faa"
    output="diamond_results/${species}_tp_results.txt"

    if [[ -f "$query" ]]; then
        echo "Running diamond blastp for $species..."
        diamond blastp -q "$query" -d /env/cns/proj/agc/scratch_microscope/Data/DBCAN/899/tcdb.dmnd -o "$output" -e 1e-10
    else
        echo "Warning: Query file $query not found for $species"
    fi
done

