#!/usr/bin/env python3
"""
Generate UMAP visualization and DBSCAN clustering for PUL systems.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
# from sklearn.manifold import TSNE
from sklearn.cluster import DBSCAN
import seaborn as sns
import umap.umap_ as umap
import os
from tqdm import tqdm

def compute_similarity_matrix(gf_sets):
    """Compute pairwise similarity matrix for gene family sets."""
    n = len(gf_sets)
    similarity = np.zeros((n, n))
    for i in tqdm(range(n), desc="Computing similarities"):
        for j in range(n):
            inter = len(gf_sets[i].intersection(gf_sets[j]))
            min_len = min(len(gf_sets[i]), len(gf_sets[j]))
            similarity[i, j] = inter / min_len if min_len > 0 else 0
    return similarity

def process_and_plot(df, output_dir="pul_clustering"):
    """Process and generate plots for both model and context GFs."""
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Process both types of gene families
    gf_types = {
        'model_GF': "Model Gene Families",
        'context_GF': "Context Gene Families"
    }
    
    # Create subplots
    fig, axes = plt.subplots(1, 2, figsize=(16, 8))
    
    for idx, (col, title) in enumerate(gf_types.items()):
        print(f"\nProcessing {title}...")
        
        # Parse gene families
        gf_sets = [set(str(x).split(", ")) for x in df[col]]
        
        # Compute similarity and distance
        similarity = compute_similarity_matrix(gf_sets)
        distance = 1 - similarity
        
        # UMAP dimensionality reduction
        embedding = umap.UMAP(metric="precomputed", 
                            random_state=42).fit_transform(distance)
        
        # DBSCAN clustering
        clustering = DBSCAN(eps=0.4, min_samples=5, 
                          metric="precomputed").fit(distance)
        
        # Plot on corresponding subplot
        sns.scatterplot(x=embedding[:, 0], y=embedding[:, 1], 
                       hue=clustering.labels_, palette='tab10', 
                       s=100, ax=axes[idx])
        axes[idx].set_title(f"PUL System Clustering ({title})")
        axes[idx].set_xlabel("UMAP Dim 1")
        axes[idx].set_ylabel("UMAP Dim 2")
        axes[idx].legend(title="Cluster")
        axes[idx].grid(True)
        
        # Print clustering statistics
        n_clusters = len(set(clustering.labels_)) - (1 if -1 in clustering.labels_ else 0)
        n_noise = list(clustering.labels_).count(-1)
        print(f"Number of clusters: {n_clusters}")
        print(f"Number of noise points: {n_noise}")
    
    plt.tight_layout()
    output_path = os.path.join(output_dir, "pul_systems_clustering.png")
    # plt.savefig(output_path, dpi=500, bbox_inches='tight')
    plt.show()
    print(f"\nPlot saved to: {output_path}")
    plt.close()

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Generate UMAP visualization for PUL systems")
    parser.add_argument("--systems", required=True, help="Path to systems.tsv file")
    parser.add_argument("-o", "--output", default="pul_clustering", 
                       help="Output directory for plots")
    
    args = parser.parse_args()
    
    # Load data
    print(f"Loading systems from: {args.systems}")
    df = pd.read_csv(args.systems, sep="\t")
    
    # Process and generate plots
    process_and_plot(df, args.output)

if __name__ == "__main__":
    main()