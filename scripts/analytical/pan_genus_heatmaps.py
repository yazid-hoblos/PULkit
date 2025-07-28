#!/usr/bin/env python3
"""
Pan-Genus System Similarity Heatmap Generator

This script generates a clustered heatmap showing similarity between bacterial systems
based on shared gene family clusters. It reads gene family cluster mappings and 
system data to compute pairwise similarities and visualize them as a heatmap.
"""

import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
from pathlib import Path


def load_cluster_mappings(clusters_file):
    """Load gene family to cluster mappings from TSV file."""
    cluster_map = {}
    with open(clusters_file) as f:
        for line in f:
            cluster_id, gf_id = line.strip().split()
            cluster_map[gf_id] = cluster_id
    
    print(f"Loaded {len(cluster_map)} gene family to cluster mappings")
    return cluster_map


def compute_similarity(set1, set2):
    """Compute similarity between two sets of clusters."""
    if not set1 or not set2:
        return 0
    intersection = set1 & set2
    return len(intersection) / min(len(set1), len(set2)) / 2


def load_systems_data(root_dir, cluster_map):
    """Load and process systems data from all species directories."""
    all_systems = []
    system_labels = []
    systems_filename = "dbcan-merged-corrected/systems.tsv"

    for species_dir in sorted(os.listdir(root_dir)):
        system_path = os.path.join(root_dir, species_dir, systems_filename)
        if not os.path.exists(system_path):
            print(f"Skipping {system_path} (not found)")
            continue

        df = pd.read_csv(system_path, sep='\t')
        for idx, row in df.iterrows():
            gfs = str(row.get("model_GF", ""))
            gf_list = [gf.strip() for gf in gfs.split(",") if gf.strip()]
            clusters = set(cluster_map.get(gf) for gf in gf_list if gf in cluster_map)
            
            if clusters:
                all_systems.append(clusters)
                system_labels.append(f"{species_dir}#{idx}")

    print(f"Processed {len(all_systems)} systems with clusters")
    return all_systems, system_labels


def compute_similarity_matrix(all_systems):
    """Compute pairwise similarity matrix for all systems."""
    n = len(all_systems)
    sim_matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(n):
            sim_matrix[i, j] = compute_similarity(all_systems[i], all_systems[j])

    return sim_matrix


def create_clustered_heatmap(similarity_matrix, output_file, figsize=(14, 14)):
    """Create and save a clustered heatmap."""
    # Convert similarity to distance
    distance_matrix = similarity_matrix.copy()
    np.fill_diagonal(distance_matrix, 0)
    
    # Create clustered heatmap
    sns.set(style="white")
    clust = sns.clustermap(distance_matrix,
                          figsize=figsize,
                          cmap="viridis",
                          xticklabels=False,
                          yticklabels=False,
                          metric="euclidean",
                          method="average")

    clust.ax_heatmap.set_title("Clustered System Distance Heatmap (1 - Cluster Similarity)")
    plt.tight_layout()
    
    # Save the plot
    clust.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Heatmap saved to: {output_file}")
    return clust


def main():
    """Main function to generate the heatmap."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument("clusters_file",
                       help="Path to gene family clusters TSV file")
    parser.add_argument("root_directory", 
                       help="Root directory containing species subdirectories")
    parser.add_argument("-o", "--output", default="report/system_cluster_heatmap.png",
                       help="Output file path for the heatmap (default: report/system_cluster_heatmap.png)")
    
    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.exists(args.clusters_file):
        print(f"Error: Clusters file '{args.clusters_file}' not found")
        return 1
    
    if not os.path.exists(args.root_directory):
        print(f"Error: Root directory '{args.root_directory}' not found")
        return 1
    
    # Create output directory
    output_dir = os.path.dirname(args.output)
    if output_dir:
        Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Load data
    print("Loading cluster mappings...")
    cluster_map = load_cluster_mappings(args.clusters_file)
    
    print("Loading systems data...")
    all_systems, system_labels = load_systems_data(args.root_directory, cluster_map)
    
    if not all_systems:
        print("Error: No systems found. Check your input paths.")
        return 1
    
    # Compute similarity matrix
    print("Computing similarity matrix...")
    similarity_matrix = compute_similarity_matrix(all_systems)
    
    # Create heatmap
    print("Creating clustered heatmap...")
    create_clustered_heatmap(similarity_matrix, args.output, tuple(args.figsize))
    
    print("Analysis complete!")
    return 0


if __name__ == "__main__":
    exit(main())