import os
import pandas as pd

def load_cgc_finder(file_path):
    df = pd.read_csv(file_path, sep='\t')
    df = df[['system number', 'gene.ID', 'contig', 'start', 'stopstrand']]
    df = df.rename(columns={'system number': 'system_number', 'gene.ID': 'gene_id', 'stopstrand': 'stop'})
    df['start'] = df['start'].astype(int)
    df['stop'] = df['stop'].astype(int)
    return df

def load_cazy_curated(file_path):
    df = pd.read_csv(file_path, sep='\t', header=None,
                     names=['pulsystem', 'gene_id', 'contig', 'start', 'stop', 'strand', 'locus_tag'])
    df['start'] = df['start'].astype(int)
    df['stop'] = df['stop'].astype(int)
    return df

def overlap(row1, row2, threshold=0.1):
    if row1['contig'] != row2['contig']:
        return False
    start = max(row1['start'], row2['start'])
    stop = min(row1['stop'], row2['stop'])
    overlap_len = max(0, stop - start)
    union_len = max(row1['stop'], row2['stop']) - min(row1['start'], row2['start'])
    return overlap_len / union_len >= threshold

def compare_predictions(cgc_df, cazy_df):
    matched = []
    unmatched_cgc = []
    unmatched_cazy = cazy_df.copy()

    for idx_cgc, row_cgc in cgc_df.iterrows():
        found = False
        for idx_cazy, row_cazy in cazy_df.iterrows():
            if overlap(row_cgc, row_cazy):
                matched.append((row_cgc.to_dict(), row_cazy.to_dict()))
                unmatched_cazy = unmatched_cazy.drop(idx_cazy)
                found = True
                break
        if not found:
            unmatched_cgc.append(row_cgc.to_dict())

    return matched, unmatched_cgc, unmatched_cazy.to_dict(orient='records')

def main():
    cazy_path = "cazy_curated/Prevotella_denticola_merged.tsv"
    cazy_df = load_cazy_curated(cazy_path)

    projection_dir = "all_prev_systems/Prevotella_denticola/dbcan-merged-corrected/projection/"
    files = [f for f in os.listdir(projection_dir) if f.endswith('.tsv')]

    for file in files:
        genome_id = file.replace('.tsv', '')  # e.g. 'GCF_000193395.1'
        cgc_path = os.path.join(projection_dir, file)
        cgc_df = load_cgc_finder(cgc_path)

        # Filter CAZy rows by genome_id in 'contig' column (assuming genome ID is part of contig name)
        filtered_cazy_df = cazy_df[cazy_df['contig'].str.contains(genome_id)]

        print(f"\nComparing genome: {genome_id}")
        print(f"CGC-Finder genes: {len(cgc_df)}, CAZy genes: {len(filtered_cazy_df)}")

        matched, unmatched_cgc, unmatched_cazy = compare_predictions(cgc_df, filtered_cazy_df)

        print(f"Matched genes: {len(matched)}")
        print(f"CGC-only genes: {len(unmatched_cgc)}")
        print(f"CAZy-only genes: {len(unmatched_cazy)}")

        # Optionally: Save or further process these lists per genome

if __name__ == "__main__":
    main()

