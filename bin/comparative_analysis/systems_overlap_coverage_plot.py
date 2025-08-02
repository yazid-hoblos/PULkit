#!/usr/bin/env python3
"""
Generate coverage plots comparing system predictions between PANORAMA, CAZy, and dbCAN.

Usage:
    python systems_coverage_plot.py [-h] [--panorama DIR] [--dbcan DIR] [--cazy FILE] [-o OUTPUT]
    
At least two tools must be specified for comparison.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
from tqdm import tqdm

def load_panorama(file_path):
    """Load PANORAMA projection file."""
    df = pd.read_csv(file_path, sep='\t')
    df =df.rename(columns={'system number': 'system_number'})
    return df[['system_number', 'contig', 'start', 'stop']]

def load_dbcan(file_path):
    """Load dbCAN CGC standard output file."""
    df = pd.read_csv(file_path, sep='\t')
    df = df.rename(columns={
        'CGC#': 'system_number',
        'Contig ID': 'contig',
        'Gene Start': 'start',
        'Gene Stop': 'stop',
    })
    return df[['system_number', 'contig', 'start', 'stop']]

def load_cazy(file_path):
    """Load CAZy curated file."""
    df = pd.read_csv(file_path, sep='\t', header=None,
                     names=['system_number', 'gene_id', 'contig', 'start', 'stop', 'strand', 'locus_tag'])
    return df[['system_number', 'contig', 'start', 'stop']]

def genes_overlap(row1, row2):
    """Check if two genes overlap."""
    if row1['contig'] != row2['contig']:
        return False
    return not (row1['stop'] < row2['start'] or row2['stop'] < row1['start'])


def extract_genome_id(filepath):
    """Extract genome ID from filepath."""
    # Get the basename or directory name depending on the path
    if 'cgc_standard_out' in filepath or 'cgc_standard.out' in filepath:
        # For dbCAN files: path/to/GCF_000191765.1-cgc-output/cgc_standard_out.tsv
        dirname = os.path.basename(os.path.dirname(filepath))
        return dirname.replace('-cgc-output', '')
    else:
        # For PANORAMA files: GCF_000191765.1.tsv
        basename = os.path.basename(filepath)
        return basename.replace('.tsv', '')

def load_data_with_genome_ids(args):
    """Load data from all tools with proper genome IDs."""
    data = {}
    
    if args.cazy:
        print("Loading CAZy data...")
        data['cazy'] = load_cazy(args.cazy)
    
    if args.panorama:
        print("Loading PANORAMA data...")
        panorama_data = []
        for f in os.listdir(args.panorama):
            if f.endswith('.tsv'):
                genome_id = extract_genome_id(f)
                df = load_panorama(os.path.join(args.panorama, f))
                df['genome'] = genome_id
                panorama_data.append(df)
        data['panorama'] = pd.concat(panorama_data)
    
    if args.dbcan:
        print("Loading dbCAN data...")
        dbcan_data = []
        for root, _, files in os.walk(args.dbcan):
            for f in files:
                if f in ["cgc_standard.out", "cgc_standard_out.tsv"]:
                    filepath = os.path.join(root, f)
                    genome_id = extract_genome_id(filepath)
                    df = load_dbcan(filepath)
                    df['genome'] = genome_id
                    dbcan_data.append(df)
        data['dbcan'] = pd.concat(dbcan_data)
    
    return data


def compute_system_overlap(df1, df2):
    """Compute system overlap between two datasets."""
    overlap_data = []
    df2_by_contig = {contig: group for contig, group in df2.groupby('contig')}
    
    for sys1, group1 in tqdm(df1.groupby('system_number')):
        total_genes1 = len(group1)
        matched_genes = set()
        
        for _, gene1 in group1.iterrows():
            contig = gene1['contig']
            if contig not in df2_by_contig:
                continue
                
            group2 = df2_by_contig[contig]
            for _, gene2 in group2.iterrows():
                if genes_overlap(gene1, gene2):
                    matched_genes.add(gene2['system_number'])
                    
        for sys2 in matched_genes:
            overlap_data.append({
                'system1': sys1,
                'system2': sys2,
                'matched_genes': sum(1 for _, gene1 in group1.iterrows()
                                   if any(genes_overlap(gene1, gene2) 
                                        for _, gene2 in df2[df2['system_number'] == sys2].iterrows())),
                'total_genes1': total_genes1,
                'total_genes2': len(df2[df2['system_number'] == sys2])
            })
    
    return pd.DataFrame(overlap_data)


def compute_system_overlap_by_genome(df1, df2, tool1, tool2):
    """Compute system overlap between two datasets, separated by genome."""
    overlap_data = []
    
    if tool1 == 'cazy':
        tool1, tool2 = tool2, tool1  # Ensure cazy is always df1 for consistency
        df1, df2 = df2, df1  # Swap dataframes
    
        df2 = df2.copy()  # Create a copy to avoid SettingWithCopyWarning
        df2['genome'] = ''  # Initialize empty genome column
        
        # Create a mapping of contigs to genomes from df1
        contig_to_genome = pd.Series(df1['genome'].values, index=df1['contig']).to_dict()
        
        # Map contigs in df2 to their corresponding genomes
        df2['genome'] = df2['contig'].map(contig_to_genome)
        
        # Remove rows where we couldn't identify the genome
        df2 = df2[df2['genome'] != '']
    
    for genome in tqdm(df1['genome'].unique(), desc="Processing genomes"):
        df1_genome = df1[df1['genome'] == genome]
        df2_genome = df2[df2['genome'] == genome]
        
        if df2_genome.empty:
            continue
            
        df2_by_contig = {contig: group for contig, group in df2_genome.groupby('contig')}
        
        for sys1, group1 in df1_genome.groupby('system_number'):
            total_genes1 = len(group1)
            matched_genes = set()
            
            for _, gene1 in group1.iterrows():
                contig = gene1['contig']
                if contig not in df2_by_contig:
                    continue
                    
                group2 = df2_by_contig[contig]
                for _, gene2 in group2.iterrows():
                    if genes_overlap(gene1, gene2):
                        matched_genes.add(gene2['system_number'])
                        
            for sys2 in matched_genes:
                overlap_data.append({
                    'genome': genome,
                    'system1': sys1,
                    'system2': sys2,
                    'matched_genes': sum(1 for _, gene1 in group1.iterrows()
                                       if any(genes_overlap(gene1, gene2) 
                                            for _, gene2 in df2_genome[df2_genome['system_number'] == sys2].iterrows())),
                    'total_genes1': total_genes1,
                    'total_genes2': len(df2_genome[df2_genome['system_number'] == sys2])
                })
    
    return pd.DataFrame(overlap_data)

# def compute_system_overlap_by_genome(df1, df2, tool1, tool2):
#     """Compute system overlap between two datasets, separated by genome."""
#     overlap_data = []
#     # Extract genome name from contig ID (assuming format: genome_contig)
#     df1['genome'] = df1['contig'].str.split('_').str[0]
#     df2['genome'] = df2['contig'].str.split('_').str[0]
    
#     for genome in tqdm(df1['genome'].unique(), desc="Processing genomes"):
#         df1_genome = df1[df1['genome'] == genome]
#         df2_genome = df2[df2['genome'] == genome]
        
#         if df2_genome.empty:
#             continue
            
#         df2_by_contig = {contig: group for contig, group in df2_genome.groupby('contig')}
        
#         for sys1, group1 in df1_genome.groupby('system_number'):
#             total_genes1 = len(group1)
#             matched_genes = set()
            
#             for _, gene1 in group1.iterrows():
#                 contig = gene1['contig']
#                 if contig not in df2_by_contig:
#                     continue
                    
#                 group2 = df2_by_contig[contig]
#                 for _, gene2 in group2.iterrows():
#                     if genes_overlap(gene1, gene2):
#                         matched_genes.add(gene2['system_number'])
                        
#             for sys2 in matched_genes:
#                 overlap_data.append({
#                     'genome': genome,
#                     'system1': sys1,
#                     'system2': sys2,
#                     'matched_genes': sum(1 for _, gene1 in group1.iterrows()
#                                        if any(genes_overlap(gene1, gene2) 
#                                             for _, gene2 in df2_genome[df2_genome['system_number'] == sys2].iterrows())),
#                     'total_genes1': total_genes1,
#                     'total_genes2': len(df2_genome[df2_genome['system_number'] == sys2])
#                 })
    
#     return pd.DataFrame(overlap_data)

def generate_coverage_plots_by_genome(overlap_df, tool1, tool2, output_dir):
    """Generate coverage plots for system overlap, both per-genome and combined."""
    thresholds = np.arange(0, 1.01, 0.05)
    overlap_df['coverage'] = overlap_df['matched_genes'] / overlap_df['total_genes2']
    
    genomes = overlap_df['genome'].unique()
    coverage_counts = pd.DataFrame(index=thresholds, columns=genomes)
    
    #seperate plot per genome 
    
    # Create per-genome directory
    # genome_dir = os.path.join(output_dir, "per_genome")
    # os.makedirs(genome_dir, exist_ok=True)

    # for genome in genomes:
    #     genome_df = overlap_df[overlap_df['genome'] == genome]
    #     for thresh in thresholds:
    #         passing = genome_df[genome_df['coverage'] >= thresh]
    #         coverage_counts.at[thresh, genome] = len(passing['system2'].unique())
        
    #     # Generate individual genome plot
    #     plt.figure(figsize=(10, 6))
    #     plt.plot(thresholds, coverage_counts[genome], marker='o')
    #     plt.xlabel(f"Coverage Threshold of {tool2} System")
    #     plt.ylabel(f"Number of {tool2} Systems Covered by {tool1}")
    #     plt.title(f"Coverage of {tool2} Systems by {tool1} Predictions\n{genome}")
    #     plt.grid(True)
    #     plt.tight_layout()
    #     plt.savefig(os.path.join(genome_dir, f"coverage_{genome}.png"), dpi=300, bbox_inches='tight')
    #     plt.close()
    
    # Combined plot
    plt.figure(figsize=(12, 7))
    for genome in genomes:
        genome_df = overlap_df[overlap_df['genome'] == genome]
        for thresh in thresholds:
            passing = genome_df[genome_df['coverage'] >= thresh]
            coverage_counts.at[thresh, genome] = len(passing['system2'].unique())
        plt.plot(coverage_counts.index, coverage_counts[genome], 
                marker='o', label=genome, alpha=0.3)
    
    plt.xlabel(f"Coverage Threshold of {tool2} System")
    plt.ylabel(f"Number of {tool2} Systems Covered by {tool1}")
    plt.title(f"Coverage of {tool2} Systems by {tool1} Predictions\nAll Genomes")
    plt.legend(title="Genome", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"per_genome_coverage_combined.png"), dpi=300, bbox_inches='tight')
    plt.close()
    
    #aggregate all genomes for a single plot
    
    # all_genomes_covered = []
    # for thresh in thresholds:
    #     passing = overlap_df[overlap_df['coverage'] >= thresh]
    #     all_genomes_covered.append(len(passing['system2'].unique()))
    
    # plt.plot(thresholds, all_genomes_covered, 
    #         marker='o', label='All Genomes', 
    #         color='red', linewidth=2, alpha=1.0)
    
    # plt.xlabel(f"Coverage Threshold of {tool2} System")
    # plt.ylabel(f"Number of {tool2} Systems Covered by {tool1}")
    # plt.title(f"Coverage of {tool2} Systems by {tool1} Predictions")
    # plt.legend(title="Genome", bbox_to_anchor=(1.05, 1), loc='upper left')
    # plt.grid(True)
    # plt.tight_layout()
    # plt.savefig(os.path.join(output_dir, f"coverage_all.png"), dpi=300, bbox_inches='tight')
    # plt.close()


def generate_coverage_plot(overlap_df, tool1, tool2, output_path):
    """Generate coverage plot for system overlap."""
    thresholds = np.arange(0, 1.01, 0.05)
    overlap_df['coverage'] = overlap_df['matched_genes'] / overlap_df['total_genes2']
    
    systems_covered = []
    for thresh in thresholds:
        passing = overlap_df[overlap_df['coverage'] >= thresh]
        systems_covered.append(len(passing['system2'].unique()))
    
    plt.figure(figsize=(10, 6))
    plt.plot(thresholds, systems_covered, marker='o')
    plt.xlabel(f"Coverage Threshold of {tool2} System")
    plt.ylabel(f"Number of {tool2} Systems Covered by {tool1}")
    plt.title(f"Coverage of {tool2} Systems by {tool1} Predictions")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Generate coverage plots between prediction tools")
    parser.add_argument("--panorama", help="Directory containing PANORAMA files")
    parser.add_argument("--dbcan", help="Directory containing dbCAN results")
    parser.add_argument("--cazy", help="Path to CAZy curated file")
    parser.add_argument("-o", "--output", default="coverage_plots",
                       help="Output directory for plots")
    
    args = parser.parse_args()
    
    # Validate input
    tools = {'panorama': args.panorama, 'dbcan': args.dbcan, 'cazy': args.cazy}
    active_tools = {k: v for k, v in tools.items() if v is not None}
    if len(active_tools) < 2:
        parser.error("At least two tools must be specified")
    
    # Create output directory
    os.makedirs(args.output, exist_ok=True)
    
    data = load_data_with_genome_ids(args)
    
    for tool1 in active_tools:
        for tool2 in active_tools:
            if tool1 >= tool2:  # Skip self-comparisons and duplicates
                continue
            
            print(f"\nComputing overlap between {tool1} and {tool2}...")
            tool_dir = os.path.join(args.output, f"{tool1}_vs_{tool2}")
            os.makedirs(tool_dir, exist_ok=True)
            
            overlap_df_by_genome = compute_system_overlap_by_genome(data[tool1], data[tool2], tool1, tool2)
            if not overlap_df_by_genome.empty:
                generate_coverage_plots_by_genome(overlap_df_by_genome, tool1, tool2, tool_dir)
                print(f"Generated plots in: {tool_dir}")
                
            overlap_df = compute_system_overlap(data[tool1], data[tool2])
            
            if not overlap_df.empty:
                output_path = os.path.join(args.output, f"coverage_{tool1}_vs_{tool2}.png")
                generate_coverage_plot(overlap_df, tool1, tool2, output_path)
                print(f"Generated plot: {output_path}")   
                         

if __name__ == "__main__":
    main()