#!/usr/bin/env python3
''' This script compares predictions from CGC-Finder with curated CAZy data.
It loads CGC-Finder predictions and curated CAZy data, checks for overlaps,
and generates a summary of matched and unmatched genes.
'''
import os
import pandas as pd
import argparse

def load_cgc_finder(file_path):
    df = pd.read_csv(file_path, sep='\t')
    df = df[['system number', 'gene.ID', 'contig', 'start', 'stop', 'strand']]
    df = df.rename(columns={'system number': 'system_number', 'gene.ID': 'gene_id'})
    df['start'] = df['start'].astype(int)
    df['stop'] = df['stop'].astype(int)
    return df

def load_cazy_curated(file_path):
    df = pd.read_csv(file_path, sep='\t', header=None,
                     names=['system_number', 'gene_id', 'contig', 'start', 'stop', 'strand', 'locus_tag'])
    df['start'] = df['start'].astype(int)
    df['stop'] = df['stop'].astype(int)
    return df

def overlap(row1, row2, threshold=0.9):
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
            # if overlap(row_cgc, row_cazy):
            if (row_cgc['contig'] == row_cazy['contig'] and
                row_cgc['start'] == row_cazy['start'] and
                row_cgc['stop'] == row_cazy['stop']):
                matched.append((row_cgc.to_dict(), row_cazy.to_dict()))
                unmatched_cazy = unmatched_cazy.drop(idx_cazy)
                found = True
                break
        if not found:
            unmatched_cgc.append(row_cgc.to_dict())

    return matched, unmatched_cgc, unmatched_cazy.to_dict(orient='records')


def plot_comparison_bar(summary_df):
    import matplotlib.ticker as mticker
    import matplotlib.pyplot as plt
    
    colors = ['#66c2a5', '#fc8d62', '#8da0cb']  # soft teal, coral, lavender
    fig, ax = plt.subplots(figsize=(14, 7))
    
    bars = summary_df[['Both', 'PANORAMA only', 'CAZy only']].plot(
        kind='bar',
        ax=ax,
        color=colors,
        alpha=0.9,
        edgecolor='black',
        width=0.75,
        legend=False
    )
    
    ax.yaxis.grid(True, linestyle='--', alpha=0.5)
    ax.set_axisbelow(True)
    
    ax.set_title("PANORAMA Predictions vs Empirical PULs", fontsize=18, weight='bold')
    ax.set_xlabel("Species", fontsize=14)
    ax.set_ylabel("Number of PUL Component Genes", fontsize=14)
    
    ax.set_xticks(range(len(summary_df)))
    ax.set_xticklabels(['VPI-5482', '7330'], rotation=45, ha='right', fontsize=12) #summary_df.index
    ax.yaxis.set_major_locator(mticker.MaxNLocator(integer=True))
    
    ax.legend(['Both', 'PANORAMA only', 'Empirical only'], title='Prediction Source', fontsize=12, title_fontsize=13, loc='upper right')
    
    for p in ax.patches:
        height = p.get_height()
        if height > 0:
            ax.annotate(f'{int(height)}',
                        (p.get_x() + p.get_width() / 2, height),
                        ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig("empirical.png", dpi=800)
    plt.show()
    

def main():
    parser = argparse.ArgumentParser(description="Compare CGC-Finder predictions with curated CAZy data.")
    parser.add_argument('--panorama_path', type=str, required=True,
                        help='Path to CGC-Finder predictions file (TSV format)')
    parser.add_argument('--cazy_path', type=str, required=True,
                        help='Path to curated CAZy data file (TSV format)')
    args = parser.parse_args()

    cazy_df = load_cazy_curated(args.cazy_path)
    
    projection_dir = os.path.dirname(args.panorama_path)
    files = [f for f in os.listdir(projection_dir) if f.endswith('.tsv')]
    summary = []
    for file in files:
        genome_id = file.replace('.tsv', '')  
        cgc_path = os.path.join(projection_dir, file)
        cgc_df = load_cgc_finder(cgc_path)
        relevant_contigs = cgc_df['contig'].unique()

        # Filter CAZy rows by genome_id in 'contig' column (assuming genome ID is part of contig name)
        filtered_cazy_df = cazy_df[cazy_df['contig'].isin(relevant_contigs)]
        if filtered_cazy_df.empty:
            print(f"No CAZy data found for genome: {genome_id}")
            continue
        print(f"\nComparing genome: {genome_id}")
        print(f"CGC-Finder genes: {len(cgc_df)}, CAZy genes: {len(filtered_cazy_df)}")

        matched, unmatched_cgc, unmatched_cazy = compare_predictions(cgc_df, filtered_cazy_df)

        print(f"Matched genes: {len(matched)}")
        print(f"CGC-only genes: {len(unmatched_cgc)}")
        print(f"CAZy-only genes: {len(unmatched_cazy)}")
        
        summary.append({
        'genome': genome_id,
        'Both': len(matched),
        'PANORAMA only': len(unmatched_cgc),
        'CAZy only': len(unmatched_cazy),
        })
        
        # cgc_systems = cgc_df['system_number'].unique()
        # cazy_systems = filtered_cazy_df['system_number'].unique()

        # # For each pair of systems, compare genes and count matched
        # for cgc_sys in cgc_systems:
        #     cgc_group = cgc_df[cgc_df['system_number'] == cgc_sys]
        #     for cazy_sys in cazy_systems:
        #         cazy_group = filtered_cazy_df[filtered_cazy_df['system_number'] == cazy_sys]

        #         matched, _, _ = compare_predictions(cgc_group, cazy_group)

        #         if len(matched) > 0:
        #             summary.append({
        #                 'genome': genome_id,
        #                 'cgc_system': cgc_sys,
        #                 'cazy_system': cazy_sys,
        #                 'cgc_total_genes': len(cgc_group),
        #                 'cazy_total_genes': len(cazy_group),
        #                 'matched_genes': len(matched),
        #             })

    # summary_df = pd.DataFrame(summary)
    # summary_df.to_csv("system_gene_pair_overlap_summary.tsv", sep='\t', index=False)
    # print(summary_df)
    
    summary_df = pd.DataFrame(summary)
    summary_df.set_index('genome', inplace=True)
    summary_df.sort_index(inplace=True)
    
    plot_comparison_bar(summary_df)


if __name__ == "__main__":
    main()

