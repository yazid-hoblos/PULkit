#!/usr/bin/env python3
"""
Generate Venn diagrams and heatmaps comparing gene predictions across PANORAMA, CAZy, and dbCAN.

Usage:
    python venn_heatmaps_all3.py [-h] --panorama DIR --cazy FILE --dbcan DIR [-o OUTPUT_DIR]

Arguments:
    --panorama DIR    Directory containing PANORAMA CGC TSV files
    --cazy FILE      Path to CAZy curated TSV file
    --dbcan DIR      Directory containing dbCAN results
    -o/--output      Output directory for plots (default: systems_comparison/)

Example:
    python venn_heatmaps_all3.py --panorama projection_dir --cazy cazy_curated.tsv --dbcan dbcan_results -o plots/
"""


import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
import seaborn as sns
import argparse

def load_projections(file_path):
    df = pd.read_csv(file_path, sep='\t')
    df = df.rename(columns={'system number': 'cgc_system', 'gene.ID': 'gene_id'})
    df['start'] = df['start'].astype(int)
    df['stop'] = df['stop'].astype(int)
    return df[['cgc_system', 'contig', 'start', 'stop', 'gene_id']]

def load_cazy_curated(file_path):
    df = pd.read_csv(file_path, sep='\t', header=None,
                     names=['cazy_system', 'gene_id', 'contig', 'start', 'stop', 'strand', 'locus_tag'])
    df['start'] = df['start'].astype(int)
    df['stop'] = df['stop'].astype(int)
    return df[['cazy_system', 'contig', 'start', 'stop', 'gene_id']]

def load_dbcan(file_path):
    df = pd.read_csv(file_path, sep='\t')
    df = df.rename(columns={'CGC#': 'dbcan_system', 'Contig ID': 'contig', 'Gene Start': 'start', 'Gene Stop': 'stop'})
    df['start'] = df['start'].astype(int)
    df['stop'] = df['stop'].astype(int)
    df = df.rename(columns={'Protein ID': 'gene_id'})
    return df[['dbcan_system', 'contig', 'start', 'stop', 'gene_id']]

def merge_datasets(cgc_df, cazy_df, dbcan_df):
    merged = pd.merge(cgc_df, cazy_df, on=['contig', 'start', 'stop'], how='outer')
    merged = pd.merge(merged, dbcan_df, on=['contig', 'start', 'stop'], how='outer')
    merged = merged.rename(columns={
        'cgc_system': 'PANORAMA',
        'cazy_system': 'CAZy',
        'dbcan_system': 'dbCAN'
    })
    return merged

def gene_overlap_counts(merged_df):
    # cgc_present = merged_df['cgc_system'].notnull()
    # cazy_present = merged_df['cazy_system'].notnull()
    # dbcan_present = merged_df['dbcan_system'].notnull()

    # Sets of gene ids per method for proper Venn counting (unique genes)
    # Using (contig, start, stop, gene_id) tuples to identify unique genes
    def gene_set(df, col_prefix):
        return set(
            df[df[f'{col_prefix}'].notnull()][['contig', 'start', 'stop', 'gene_id']]
            .apply(tuple, axis=1)
        )

    cgc_genes = gene_set(merged_df, 'PANORAMA')
    cazy_genes = gene_set(merged_df, 'CAZy')
    dbcan_genes = gene_set(merged_df, 'dbCAN')

    only_cgc = cgc_genes - cazy_genes - dbcan_genes
    only_cazy = cazy_genes - cgc_genes - dbcan_genes
    only_dbcan = dbcan_genes - cgc_genes - cazy_genes

    cgc_cazy = (cgc_genes & cazy_genes) - dbcan_genes
    cgc_dbcan = (cgc_genes & dbcan_genes) - cazy_genes
    cazy_dbcan = (cazy_genes & dbcan_genes) - cgc_genes

    all_three = cgc_genes & cazy_genes & dbcan_genes

    return {
        'PANORAMA only': len(only_cgc),
        'CAZy only': len(only_cazy),
        'dbCAN only': len(only_dbcan),
        'PANORAMA & CAZY': len(cgc_cazy),
        'PANORAMA & dbCAN': len(cgc_dbcan),
        'CAZy & dbCAN': len(cazy_dbcan),
        'All three': len(all_three)
    }

def plot_venn(overlap_counts, output_path):
    plt.figure(figsize=(7,7))
    venn3(subsets=(
        overlap_counts['PANORAMA only'],
        overlap_counts['CAZy only'],
        overlap_counts['PANORAMA & CAZY'],
        overlap_counts['dbCAN only'],
        overlap_counts['PANORAMA & dbCAN'],
        overlap_counts['CAZy & dbCAN'],
        overlap_counts['All three']),
        set_labels=('PANORAMA', 'CAZy', 'dbCAN'))
    plt.title("Component Gene Predictions Overlap Across Methods")
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

def system_gene_matches(merged_df, system_col1, system_col2):
    systems = merged_df[system_col1].dropna().unique()
    records = []

    for sys in systems:
        genes_in_sys = merged_df[merged_df[system_col1] == sys]
        total_genes = len(genes_in_sys)
        matched_genes = genes_in_sys[genes_in_sys[system_col2].notnull()]
        frac_matched = len(matched_genes) / total_genes if total_genes > 0 else 0
        records.append({'system': sys, 'total_genes': total_genes, 'matched_genes': len(matched_genes), 'frac_matched': frac_matched})
    return pd.DataFrame(records)

def plot_system_overlap_heatmap(merged_df, sys_col1, sys_col2, output_path):
    matched = merged_df.dropna(subset=[sys_col1, sys_col2])
    matrix = matched.groupby([sys_col1, sys_col2]).size().unstack(fill_value=0)
    matrix_frac = matrix.div(matrix.sum(axis=1), axis=0).fillna(0)

    plt.figure(figsize=(10,8))
    sns.heatmap(matrix_frac, cmap='Blues', cbar_kws={'label': 'Fraction of genes matched'})
    plt.xlabel(sys_col2)
    plt.ylabel(sys_col1)
    plt.title(f'System overlap: {sys_col1} vs {sys_col2}')
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def main():
    """Parse arguments and generate visualizations."""
    parser = argparse.ArgumentParser(description="Compare gene predictions across PANORAMA, CAZy, and dbCAN")
    parser.add_argument("--panorama", required=True, help="Directory containing PANORAMA CGC TSV files")
    parser.add_argument("--cazy", required=True, help="Path to CAZy curated TSV file")
    parser.add_argument("--dbcan", required=True, help="Directory containing dbCAN results")
    parser.add_argument("-o", "--output", default="systems_comparison", help="Output directory for plots")
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output, exist_ok=True)
    
    print("Loading all CGC files from PANORAMA directory...")
    projection_files = glob.glob(os.path.join(args.panorama, "*.tsv"))
    projection_dfs = []
    for f in projection_files:
        print(f"Loading {f}")
        projection_dfs.append(load_projections(f))
    projection_df = pd.concat(projection_dfs, ignore_index=True)

    print("Loading CAZy curated data...")
    cazy_df = load_cazy_curated(args.cazy)

    print("Loading all dbCAN files...")
    dbcan_files = glob.glob(os.path.join(args.dbcan, "*-cgc-output", "cgc_standard.out"))
    dbcan_dfs = []
    for f in dbcan_files:
        print(f"Loading {f}")
        dbcan_dfs.append(load_dbcan(f))
    dbcan_df = pd.concat(dbcan_dfs, ignore_index=True)

    print("Merging datasets across all genomes...")
    merged = merge_datasets(projection_df, cazy_df, dbcan_df)
    print(f"Total unique genes (union): {len(merged)}")

    print("Calculating gene overlap counts...")
    overlaps = gene_overlap_counts(merged)
    print("Gene overlaps:")
    for k, v in overlaps.items():
        print(f"{k}: {v}")

    venn_path = os.path.join(args.output, "venn_diagram.png")
    print(f"Saving Venn diagram to {venn_path}")
    plot_venn(overlaps, venn_path)

    pairs = [
        ('PANORAMA', 'CAZy'),
        ('PANORAMA', 'dbCAN'),
        ('CAZy', 'dbCAN')
    ]

    for sys1, sys2 in pairs:
        print(f"\nAnalyzing system overlap {sys1} vs {sys2}...")
        df_match = system_gene_matches(merged, sys1, sys2)
        print(df_match.head())

        heatmap_path = os.path.join(args.output, f'heatmap_{sys1}_vs_{sys2}.png')
        print(f"Saving heatmap to {heatmap_path}")
        plt.figure(figsize=(10,8))
        plot_system_overlap_heatmap(merged, sys1, sys2, heatmap_path)

if __name__ == "__main__":
    main()

# panorama_dir = "all_prev_systems_tp_only/Prevotella_denticola/dbcan-merged-corrected/projection"
# cazy_path = "cazy_curated/Prevotella_denticola_merged.tsv"
# dbcan_dir = "dbcan-results-Pd_genomes"
