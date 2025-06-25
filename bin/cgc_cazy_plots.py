import os
import pandas as pd
import matplotlib.pyplot as plt

def load_cgc_standard(file_path):
    df = pd.read_csv(file_path, sep='\t')
    df = df.rename(columns={
        'CGC#': 'cgc_system',
        'Contig ID': 'contig',
        'Gene Start': 'start',
        'Gene Stop': 'stop',
    })
    df['start'] = df['start'].astype(int)
    df['stop'] = df['stop'].astype(int)
    return df[['cgc_system', 'contig', 'start', 'stop']]

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
                    'cgc_system': row_cgc['cgc_system'],
                    'cazy_system': row_cazy['cazy_system'],
                })
                break
    return pd.DataFrame(matched_pairs)

def plot_scatter(match_df, cgc_df, cazy_df, title, output_path):
    cgc_sizes = cgc_df.groupby('cgc_system').size().rename('cgc_total_genes')
    cazy_sizes = cazy_df.groupby('cazy_system').size().rename('cazy_total_genes')

    matched_counts = match_df.groupby(['cgc_system', 'cazy_system']).size().rename('matched_genes').reset_index()
    merged = matched_counts.merge(cgc_sizes, on='cgc_system').merge(cazy_sizes, on='cazy_system')

    merged['overlap_pct_cgc'] = merged['matched_genes'] / merged['cgc_total_genes']
    merged['overlap_pct_cazy'] = merged['matched_genes'] / merged['cazy_total_genes']

    plt.figure(figsize=(10, 8))
    scatter = plt.scatter(merged['cgc_total_genes'], merged['cazy_total_genes'],
                          s=merged['matched_genes'] * 20,
                          c=merged['overlap_pct_cgc'], cmap='viridis', alpha=0.7, edgecolors='k')
    plt.colorbar(scatter, label='% Overlap relative to CGC system size')
    plt.xlabel('CGC-Finder System Size (total genes)')
    plt.ylabel('CAZy System Size (total genes)')
    plt.title(title)
    # plt.grid(True)
    plt.tight_layout()

    plt.savefig(output_path)
    plt.close()  # close figure to free memory

def main():
    dbcan_results_dir = "dbcan-results-Pd_genomes"
    cazy_file = "cazy_curated/Prevotella_denticola_merged.tsv"
    output_dir = "plots"
    os.makedirs(output_dir, exist_ok=True)

    cazy_df = load_cazy(cazy_file)
    print(f"Loaded CAZy data with {len(cazy_df)} genes.")

    species_files = []
    for root, dirs, files in os.walk(dbcan_results_dir):
        for file in files:
            if file.endswith("cgc_standard.out"):
                species_files.append(os.path.join(root, file))
    print(f"Found {len(species_files)} CGC files (species) to process.")

    all_matches = []
    all_cgc_dfs = []
    all_cazy_dfs = []

    for cgc_file in species_files:
        species_name = os.path.basename(os.path.dirname(cgc_file))
        print(f"\nProcessing species: {species_name}")

        cgc_df = load_cgc_standard(cgc_file)
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
                     f'System Size vs Matched Genes (dbCAN vs CAZy) - {species_name}', species_plot_path)
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
