#!/usr/bin/env python3
"""
Generate comparison plots between PANORAMA, dbCAN, and CAZy annotations.

Usage:
    python comparison_plots.py [-h] [--panorama DIR] [--dbcan DIR] [--cazy FILE] [-o OUTPUT]
    
At least two tools must be specified for comparison.
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
from tqdm import tqdm
from itertools import combinations

def load_panorama(file_path):
    """Load PANORAMA projection file."""
    df = pd.read_csv(file_path, sep='\t')
    df = df.rename(columns={'system number': 'panorama_system'})
    return df[['panorama_system', 'contig', 'start', 'stop']]

def load_dbcan(file_path):
    """Load dbCAN CGC standard output file."""
    df = pd.read_csv(file_path, sep='\t')
    df = df.rename(columns={
        'CGC#': 'dbcan_system',
        'Contig ID': 'contig',
        'Gene Start': 'start',
        'Gene Stop': 'stop',
    })
    df['start'] = df['start'].astype(int)
    df['stop'] = df['stop'].astype(int)
    return df[['dbcan_system', 'contig', 'start', 'stop']]

def load_cazy(file_path):
    """Load CAZy curated file."""
    df = pd.read_csv(file_path, sep='\t', header=None,
                     names=['system_number', 'gene_id', 'contig', 'start', 'stop', 'strand', 'locus_tag'])
    df = df.rename(columns={'system_number': 'cazy_system'})
    df['start'] = df['start'].astype(int)
    df['stop'] = df['stop'].astype(int)
    return df[['cazy_system', 'contig', 'start', 'stop']]

def genes_overlap(row1, row2):
    """Check if two genes overlap."""
    if row1['contig'] != row2['contig']:
        return False
    return not (row1['stop'] < row2['start'] or row2['stop'] < row1['start'])

def match_genes(df1, df2, tool1, tool2):
    """Match genes between two datasets."""
    matched_pairs = []
    df2_by_contig = {contig: group for contig, group in df2.groupby('contig')}
    
    for _, row1 in df1.iterrows():
        contig = row1['contig']
        if contig not in df2_by_contig:
            continue
        group2 = df2_by_contig[contig]
        for _, row2 in group2.iterrows():
            if genes_overlap(row1, row2):
                matched_pairs.append({
                    f'{tool1}_system': row1[f'{tool1}_system'],
                    f'{tool2}_system': row2[f'{tool2}_system'],
                })
                break
    return pd.DataFrame(matched_pairs)

def plot_comparison(df1, df2, match_df, tool1, tool2, title, output_path):
    """Create scatter plot comparing system sizes between two tools."""
    sizes1 = df1.groupby(f'{tool1}_system').size().rename(f'{tool1}_total_genes')
    sizes2 = df2.groupby(f'{tool2}_system').size().rename(f'{tool2}_total_genes')

    matched_counts = match_df.groupby([f'{tool1}_system', f'{tool2}_system']).size().rename('matched_genes').reset_index()
    merged = matched_counts.merge(sizes1, left_on=f'{tool1}_system', right_index=True)\
                         .merge(sizes2, left_on=f'{tool2}_system', right_index=True)

    plt.figure(figsize=(10, 8))
    scatter = plt.scatter(merged[f'{tool1}_total_genes'], merged[f'{tool2}_total_genes'],
                         s=merged['matched_genes'] * 20,
                         c=merged['matched_genes'] / merged[f'{tool1}_total_genes'],
                         cmap='viridis', alpha=0.7, edgecolors='k')
    plt.colorbar(scatter, label=f'% Overlap relative to {tool1} system size')
    plt.xlabel(f'{tool1} System Size (total genes)')
    plt.ylabel(f'{tool2} System Size (total genes)')
    plt.title(title)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def process_genome_files(tool1_files, tool2_files, tool1, tool2, output_dir):
    """Process files for a specific genome and create comparison plots."""
    results = []
    
    # Handle CAZy's single file case
    if tool2 == 'cazy' and len(tool2_files) == 1:
        cazy_df = globals()[f'load_{tool2}'](tool2_files[0])
        
        for f1 in tqdm(tool1_files, desc=f"Processing {tool1} files"):
            genome_name = os.path.splitext(os.path.basename(os.path.dirname(f1)))[0]
            df1 = globals()[f'load_{tool1}'](f1)
            
            # Filter CAZy data for this genome's contigs
            species_contigs = df1['contig'].unique()
            df2 = cazy_df[cazy_df['contig'].isin(species_contigs)]
            
            if df2.empty:
                print(f"No CAZy genes found for genome {genome_name}. Skipping.")
                continue
            
            match_df = match_genes(df1, df2, tool1, tool2)
            
            if not match_df.empty:
                genome_output = os.path.join(output_dir, f"{genome_name}.png")
                plot_comparison(df1, df2, match_df,
                              tool1, tool2,
                              f'System Size Comparison: {tool1} vs {tool2}\n{genome_name}',
                              genome_output)
                results.extend([df1, df2, match_df])
                print(f"Created plot for {genome_name}")
    
    # Handle tool vs tool comparison (both have per-genome files)
    else:
        for f1 in tqdm(tool1_files, desc=f"Processing {tool1} files"):
            genome_name = os.path.splitext(os.path.basename(os.path.dirname(f1)))[0]
            df1 = globals()[f'load_{tool1}'](f1)
            
            # Find matching file for tool2
            matching_f2 = [f2 for f2 in tool2_files if genome_name in f2]
            if not matching_f2:
                continue
                
            df2 = globals()[f'load_{tool2}'](matching_f2[0])
            match_df = match_genes(df1, df2, tool1, tool2)
            
            if not match_df.empty:
                genome_output = os.path.join(output_dir, f"{genome_name}.png")
                plot_comparison(df1, df2, match_df,
                              tool1, tool2,
                              f'System Size Comparison: {tool1} vs {tool2}\n{genome_name}',
                              genome_output)
                results.extend([df1, df2, match_df])
                print(f"Created plot for {genome_name}")
    
    return results

def main():
    parser = argparse.ArgumentParser(description="Generate comparison plots between PANORAMA, dbCAN, and CAZy")
    parser.add_argument("--panorama", help="Directory containing PANORAMA projection files")
    parser.add_argument("--dbcan", help="Directory containing dbCAN results")
    parser.add_argument("--cazy", help="Path to CAZy curated TSV file")
    parser.add_argument("-o", "--output", default="comparison_plots", 
                       help="Output directory for plots")
    
    args = parser.parse_args()
    
    # Check that at least two tools are specified
    tools = {'panorama': args.panorama, 'dbcan': args.dbcan, 'cazy': args.cazy}
    active_tools = {k: v for k, v in tools.items() if v is not None}
    if len(active_tools) < 2:
        parser.error("At least two tools must be specified for comparison")
    
    # Create main output directory
    os.makedirs(args.output, exist_ok=True)
    
    # Load and organize files
    files = {}
    if args.cazy:
        print("Loading CAZy data...")
        files['cazy'] = [args.cazy]
    
    if args.panorama:
        print("Loading PANORAMA files...")
        files['panorama'] = [os.path.join(args.panorama, f) 
                           for f in os.listdir(args.panorama) 
                           if f.endswith('.tsv')]
    
    if args.dbcan:
        print("Loading dbCAN files...")
        files['dbcan'] = []
        for root, _, filenames in os.walk(args.dbcan):
            for f in filenames:
                if f in ["cgc_standard.out", "cgc_standard_out.tsv"]:
                    files['dbcan'].append(os.path.join(root, f))
    
    # Process each tool combination
    for tool1, tool2 in combinations(active_tools.keys(), 2):
        print(f"\nProcessing {tool1} vs {tool2} comparisons...")
        
        # Create directory for this comparison
        comparison_dir = os.path.join(args.output, f"{tool1}_vs_{tool2}")
        os.makedirs(comparison_dir, exist_ok=True)
        
        # Create per-genome plots
        genome_dir = os.path.join(comparison_dir, "per_genome")
        os.makedirs(genome_dir, exist_ok=True)
        
        results = process_genome_files(files[tool1], files[tool2], 
                                     tool1, tool2, genome_dir)
        
        if results:
            # Create combined plot
            df1 = pd.concat([r for r in results[::3]], ignore_index=True)
            df2 = pd.concat([r for r in results[1::3]], ignore_index=True)
            match_df = pd.concat([r for r in results[2::3]], ignore_index=True)
            
            combined_output = os.path.join(comparison_dir, "combined_comparison.png")
            plot_comparison(df1, df2, match_df,
                          tool1, tool2,
                          f'Combined System Size Comparison: {tool1} vs {tool2}',
                          combined_output)
            print(f"Created combined plot for {tool1} vs {tool2}")

if __name__ == "__main__":
    main()