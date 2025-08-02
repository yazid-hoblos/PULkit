#!/usr/bin/env python3
"""
Generate similarity heatmaps for PUL systems based on subset relationships
between model and context gene families from two system files.

Usage:
    python post-processing.py --systems1 FILE1 --systems2 FILE2 [-o OUTPUT]
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import argparse
import os
from tqdm import tqdm
from scipy.cluster.hierarchy import linkage, leaves_list

def subset_similarity(set1, set2):
    """Return 1.0 if set1 is subset of set2 or vice versa, else 0.0."""
    return 1.0 if set1.issubset(set2) or set2.issubset(set1) else 0.0

def reorder_matrix(sim_matrix):
    """Reorder rows and columns of the matrix using hierarchical clustering."""
    if sim_matrix.shape[0] < 2 or sim_matrix.shape[1] < 2:
        return sim_matrix, np.arange(sim_matrix.shape[0]), np.arange(sim_matrix.shape[1])
    
    row_linkage = linkage(sim_matrix, method='average', metric='euclidean')
    col_linkage = linkage(sim_matrix.T, method='average', metric='euclidean')
    row_order = leaves_list(row_linkage)
    col_order = leaves_list(col_linkage)
    reordered = sim_matrix[np.ix_(row_order, col_order)]
    return reordered, row_order, col_order

def generate_similarity_heatmaps(file1, file2, output_dir="similarity_plots"):
    """Generate and save similarity heatmaps comparing systems in file1 vs file2 using subset matching."""
    os.makedirs(output_dir, exist_ok=True)

    print(f"Loading systems from: {file1}")
    df1 = pd.read_csv(file1, sep='\t')
    print(f"Loading systems from: {file2}")
    df2 = pd.read_csv(file2, sep='\t')

    # Robust splitting of gene family strings
    model_sets1 = [set(map(str.strip, str(x).split(','))) for x in df1['model_GF']]
    context_sets1 = [set(map(str.strip, str(x).split(','))) for x in df1['context_GF']]
    model_sets2 = [set(map(str.strip, str(x).split(','))) for x in df2['model_GF']]
    context_sets2 = [set(map(str.strip, str(x).split(','))) for x in df2['context_GF']]

    n1 = len(df1)
    n2 = len(df2)

    model_sim = np.zeros((n1, n2))
    context_sim = np.zeros((n1, n2))

    print("Computing subset-based similarity matrices...")
    for i in tqdm(range(n1), desc="Computing similarities", unit="system"):
        for j in range(n2):
            model_sim[i, j] = subset_similarity(model_sets1[i], model_sets2[j])
            context_sim[i, j] = subset_similarity(context_sets1[i], context_sets2[j])

    print("Clustering model_GF similarity matrix...")
    model_reordered, _, _ = reorder_matrix(model_sim)

    print("Clustering context_GF similarity matrix...")
    context_reordered, _, _ = reorder_matrix(context_sim)

    fig, axes = plt.subplots(1, 2, figsize=(18, 10))

    sns.heatmap(model_reordered, cmap='Blues', ax=axes[0], cbar=True,
                xticklabels=False, yticklabels=False, vmin=0, vmax=1,
                linewidths=0.1, linecolor='gray')
    axes[0].set_title("Model_GF Subset Match (Clustered)")

    sns.heatmap(context_reordered, cmap='Oranges', ax=axes[1], cbar=True,
                xticklabels=False, yticklabels=False, vmin=0, vmax=1,
                linewidths=0.1, linecolor='gray')
    axes[1].set_title("Context_GF Subset Match (Clustered)")

    plt.tight_layout()
    output_path = os.path.join(output_dir, "subset_similarity_heatmaps.png")
    plt.show()
    print(f"Saved subset-based heatmaps to: {output_path}")
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Compare PUL systems from two files using subset relationships")
    parser.add_argument("--systems1", required=True, help="Path to first systems.tsv file")
    parser.add_argument("--systems2", required=True, help="Path to second systems.tsv file")
    parser.add_argument("-o", "--output", default="similarity_plots", help="Output directory")

    args = parser.parse_args()
    generate_similarity_heatmaps(args.systems1, args.systems2, args.output)

if __name__ == "__main__":
    main()
