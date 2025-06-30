#!/usr/bin/env python3
"""
Generate similarity heatmaps for PUL systems based on model and context gene families.

Usage:
    python post-processing.py [-h] --systems FILE [-o OUTPUT]
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import argparse
import os
from tqdm import tqdm

def overlap_min_similarity(set1, set2):
    """Calculate similarity as intersection over minimum size."""
    intersection = set1 & set2
    min_len = min(len(set1), len(set2))
    return len(intersection) / min_len if min_len > 0 else 0

def generate_similarity_heatmaps(systems_file, output_dir="similarity_plots"):
    """Generate and save similarity heatmaps for model and context gene families."""
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Load data
    print(f"Loading systems from: {systems_file}")
    df = pd.read_csv(systems_file, sep='\t')

    # Convert to sets
    model_sets = [set(str(x).split(', ')) for x in df['model_GF']]
    context_sets = [set(str(x).split(', ')) for x in df['context_GF']]

    n = len(df)
    print(f"Total systems: {n}")
    model_sim = np.zeros((n, n))
    context_sim = np.zeros((n, n))

    # Compute similarity matrices
    for i in tqdm(range(n), desc="Computing similarities", unit="system"):
        # print(f"Processing system {i+1}/{n}")
        for j in range(n):
            model_sim[i, j] = overlap_min_similarity(model_sets[i], model_sets[j])
            context_sim[i, j] = overlap_min_similarity(context_sets[i], context_sets[j])

    # Generate clustermaps
    model_cluster = sns.clustermap(model_sim, cmap='Blues', xticklabels=False, 
                                 yticklabels=False, method='average', metric='euclidean')
    plt.close(model_cluster.fig)

    context_cluster = sns.clustermap(context_sim, cmap='Oranges', xticklabels=False, 
                                   yticklabels=False, method='average', metric='euclidean')
    plt.close(context_cluster.fig)

    # Get reordered indices
    model_idx = model_cluster.dendrogram_row.reordered_ind
    context_idx = context_cluster.dendrogram_row.reordered_ind

    # Reorder matrices
    model_reordered = model_sim[np.ix_(model_idx, model_idx)]
    context_reordered = context_sim[np.ix_(context_idx, context_idx)]

    # Plot combined heatmaps
    fig, axes = plt.subplots(1, 2, figsize=(16, 8))

    sns.heatmap(model_reordered, cmap='Blues', ax=axes[0], cbar=True, 
                xticklabels=False, yticklabels=False)
    axes[0].set_title("Model_GF Similarity (Clustered)")

    sns.heatmap(context_reordered, cmap='Oranges', ax=axes[1], cbar=True, 
                xticklabels=False, yticklabels=False)
    axes[1].set_title("Context_GF Similarity (Clustered)")

    plt.tight_layout()
    output_path = os.path.join(output_dir, "similarity_heatmaps.png")
    plt.savefig(output_path, dpi=500, bbox_inches='tight')
    print(f"Saved heatmaps to: {output_path}")
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Generate similarity heatmaps for PUL systems")
    parser.add_argument("--systems", required=True, help="Path to systems.tsv file")
    parser.add_argument("-o", "--output", default="similarity_plots", 
                       help="Output directory for plots (default: similarity_plots)")
    
    args = parser.parse_args()
    generate_similarity_heatmaps(args.systems, args.output)

if __name__ == "__main__":
    main()