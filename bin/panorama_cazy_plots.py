import os
import pandas as pd
import matplotlib.pyplot as plt

def load_cgc_finder(file_path):
    df = pd.read_csv(file_path, sep='\t')
    df = df[['system number', 'gene.ID', 'contig', 'start', 'stop', 'strand']]
    df = df.rename(columns={'system number': 'system_number', 'gene.ID': 'gene_id'})
    df['start'] = df['start'].astype(int)
    df['stop'] = df['stop'].astype(int)
    return df

def load_cazy(file_path):
    df = pd.read_csv(file_path, sep='\t', header=None,
                     names=['system_number', 'gene_id', 'contig', 'start', 'stop', 'strand', 'locus_tag'])
    df = df.rename(columns={'system_number': 'cazy_system'})
    df['start'] = df['start'].astype(int)
    df['stop'] = df['stop'].astype(int)
    return df[['cazy_system', 'contig', 'start', 'stop']]

def genes_overlap(row1, row2):
    if row1['contig'] != row2['contig']:
        return False
    return not (row1['stop'] < row2['start'] or row2['stop'] < row1['start'])

def match_genes(cgc_df, cazy_df):
    matched_pairs = []
    cazy_by_contig = {contig: group for contig, group in cazy_df.groupby('contig')}
    for idx_cgc, row_cgc in cgc_df.iterrows():
        contig = row_cgc['contig']
        if contig not in cazy_by_contig:
            continue
        cazy_group = cazy_by_contig[contig]
        for idx_cazy, row_cazy in cazy_group.iterrows():
            if genes_overlap(row_cgc, row_cazy):
                matched_pairs.append({
                    'cgc_system': row_cgc['system_number'],
                    'cazy_system': row_cazy['cazy_system'],
                })
                break
    return pd.DataFrame(matched_pairs)

def plot_scatter(match_df, cgc_df, cazy_df, title, output_path):
    cgc_sizes = cgc_df.groupby('system_number').size().rename('cgc_total_genes')
    cazy_sizes = cazy_df.groupby('cazy_system').size().rename('cazy_total_genes')

    matched_counts = match_df.groupby(['cgc_system', 'cazy_system']).size().rename('matched_genes').reset_index()
    merged = matched_counts.merge(cgc_sizes, left_on='cgc_system', right_index=True)
    merged = merged.merge(cazy_sizes, left_on='cazy_system', right_index=True)

    merged['overlap_pct_cgc'] = merged['matched_genes'] / merged['cgc_total_genes']
    merged['overlap_pct_cazy'] = merged['matched_genes'] / merged['cazy_total_genes']

    plt.figure(figsize=(10, 8))
    scatter = plt.scatter(merged['cgc_total_genes'], merged['cazy_total_genes'],
                          s=merged['matched_genes'] * 20,
                          c=merged['overlap_pct_cgc'], cmap='viridis', alpha=0.7, edgecolors='k')
    plt.colorbar(scatter, label='% Overlap relative to PANORAMA system size')
    plt.xlabel('PANORAMA System Size (total gene families)')
    plt.ylabel('CAZy System Size (total genes)')
    plt.title(title)
    plt.tight_layout()

    plt.savefig(output_path)
    plt.close()

def main():
    cgc_dir = "all_prev_systems_tp_only/Prevotella_denticola/dbcan-merged-corrected/projection"
    cazy_file = "cazy_curated/Prevotella_denticola_merged.tsv"
    output_dir = "panorama_cazy_plots"
    os.makedirs(output_dir, exist_ok=True)

    print("Loading CAZy curated data...")
    cazy_df = load_cazy(cazy_file)
    print(f"Loaded CAZy data with {len(cazy_df)} genes.")

    print("Loading all CGC files from projection directory...")
    cgc_files = [os.path.join(cgc_dir, f) for f in os.listdir(cgc_dir) if f.endswith('.tsv')]
    print(f"Found {len(cgc_files)} CGC files.")

    all_matches = []
    all_cgc_dfs = []
    all_cazy_dfs = []

    for cgc_file in cgc_files:
        species_name = os.path.splitext(os.path.basename(cgc_file))[0]
        print(f"\nProcessing species: {species_name}")

        cgc_df = load_cgc_finder(cgc_file)
        species_contigs = cgc_df['contig'].unique()
        cazy_species_df = cazy_df[cazy_df['contig'].isin(species_contigs)]

        if cazy_species_df.empty:
            print(f"No CAZy genes found for species {species_name}. Skipping.")
            continue

        match_df = match_genes(cgc_df, cazy_species_df)
        if match_df.empty:
            print(f"No matched genes found for {species_name}.")
            continue

        species_plot_path = os.path.join(output_dir, f"{species_name}_plot.png")
        plot_scatter(match_df, cgc_df, cazy_species_df,
                     f'System Size vs Matched Genes (PANORAMA vs CAZy) - {species_name}', species_plot_path)
        print(f"Saved plot for {species_name} at {species_plot_path}")

        all_matches.append(match_df)
        all_cgc_dfs.append(cgc_df)
        all_cazy_dfs.append(cazy_species_df)

    if all_matches:
        combined_matches = pd.concat(all_matches, ignore_index=True)
        combined_cgc = pd.concat(all_cgc_dfs, ignore_index=True)
        combined_cazy = pd.concat(all_cazy_dfs, ignore_index=True)

        combined_plot_path = os.path.join(output_dir, "combined_plot.png")
        plot_scatter(combined_matches, combined_cgc, combined_cazy,
                     'Combined System Size vs Matched Genes (All Species)', combined_plot_path)
        print(f"Saved combined plot at {combined_plot_path}")
    else:
        print("No matches found across all species.")

if __name__ == "__main__":
    main()
