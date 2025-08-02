#!/usr/bin/env python3
"""
UMAP Analysis of Gene Family Clusters

This script generates a UMAP projection of bacterial systems based on gene family clusters.
It reads gene family cluster mappings, computes similarity between systems based on shared clusters,
and visualizes the results using UMAP with options for outlier filtering and analysis.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import umap
import seaborn as sns
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


def load_systems_data(root_dir, systems_filename, cluster_map):
    """Load and process systems data from all species directories."""
    all_systems = []
    system_labels = []
    species_labels = []

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
                species_labels.append(species_dir)

    print(f"Processed {len(all_systems)} systems with clusters")
    return all_systems, system_labels, species_labels


def compute_similarity_matrix(all_systems):
    """Compute pairwise similarity matrix for all systems."""
    n = len(all_systems)
    sim_matrix = np.zeros((n, n))
    
    for i in range(n):
        for j in range(i, n):
            sim = compute_similarity(all_systems[i], all_systems[j])
            sim_matrix[i, j] = sim
            sim_matrix[j, i] = sim
    
    return 1 - sim_matrix  # Convert to distance matrix


def perform_umap(distance_matrix, random_state=42):
    """Perform UMAP dimensionality reduction."""
    reducer = umap.UMAP(metric="precomputed", random_state=random_state)
    return reducer.fit_transform(distance_matrix)


def identify_outliers(embedding, method="statistical", threshold_factor=2, percentile=95):
    """Identify outliers based on distance from center."""
    center_x, center_y = np.mean(embedding[:, 0]), np.mean(embedding[:, 1])
    distances = np.sqrt((embedding[:, 0] - center_x)**2 + (embedding[:, 1] - center_y)**2)
    
    if method == "statistical":
        threshold = np.mean(distances) + threshold_factor * np.std(distances)
    elif method == "percentile":
        threshold = np.percentile(distances, percentile)
    else:
        raise ValueError("Method must be 'statistical' or 'percentile'")
    
    return distances > threshold, distances, threshold


def analyze_outliers(embedding, all_systems, system_labels, species_labels, n_top=5):
    """Analyze and report top outliers."""
    center_x, center_y = np.mean(embedding[:, 0]), np.mean(embedding[:, 1])
    distances = np.sqrt((embedding[:, 0] - center_x)**2 + (embedding[:, 1] - center_y)**2)
    
    top_indices = np.argsort(distances)[-n_top:][::-1]
    
    print(f"\n=== TOP {n_top} OUTLIER SYSTEMS ===")
    print(f"Center coordinates: ({center_x:.2f}, {center_y:.2f})")
    
    for i, idx in enumerate(top_indices, 1):
        species = species_labels[idx]
        system = system_labels[idx]
        distance = distances[idx]
        x, y = embedding[idx]
        clusters = sorted(list(all_systems[idx]))
        
        print(f"\n{i}. Distance: {distance:.2f} | Coords: ({x:.2f}, {y:.2f})")
        print(f"   Species: {species}")
        print(f"   System: {system}")
        print(f"   Gene clusters ({len(clusters)}): {clusters}")
    
    outlier_species = [species_labels[idx] for idx in top_indices]
    unique_clusters = len(set().union(*[all_systems[idx] for idx in top_indices]))
    
    print(f"\n=== SUMMARY ===")
    print(f"Species represented: {set(outlier_species)}")
    print(f"Unique gene clusters across top outliers: {unique_clusters}")


def create_umap_plot(embedding, species_labels, title, output_file, show_plot=False):
    """Create and save UMAP plot."""
    plt.figure(figsize=(14, 10))
    
    palette = sns.color_palette("tab20", n_colors=len(set(species_labels)))
    species_to_color = {sp: palette[i] for i, sp in enumerate(sorted(set(species_labels)))}
    colors = [species_to_color[sp] for sp in species_labels]
    
    plt.scatter(embedding[:, 0], embedding[:, 1], c=colors, s=15, alpha=0.7)
    plt.title(title)
    plt.xlabel("UMAP 1")
    plt.ylabel("UMAP 2")
    plt.grid(True)
    
    # Legend
    handles = [plt.Line2D([0], [0], marker='o', color='w', label=sp,
                         markerfacecolor=species_to_color[sp], markersize=8)
               for sp in species_to_color]
    plt.legend(handles=handles, title="Species", bbox_to_anchor=(1.05, 1), loc="upper left")
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    
    if show_plot:
        plt.show()
    else:
        plt.close()
    
    print(f"Plot saved: {output_file}")


def create_comparison_plot(embedding_all, embedding_filtered, species_all, species_filtered, 
                          output_file, show_plot=False):
    """Create side-by-side comparison plot."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))
    
    # Original plot
    palette_all = sns.color_palette("tab20", n_colors=len(set(species_all)))
    species_to_color_all = {sp: palette_all[i] for i, sp in enumerate(sorted(set(species_all)))}
    colors_all = [species_to_color_all[sp] for sp in species_all]
    
    ax1.scatter(embedding_all[:, 0], embedding_all[:, 1], c=colors_all, s=15, alpha=0.7)
    ax1.set_title("Original (with outliers)")
    ax1.set_xlabel("UMAP 1")
    ax1.set_ylabel("UMAP 2")
    ax1.grid(True)
    
    # Filtered plot
    palette_filtered = sns.color_palette("tab20", n_colors=len(set(species_filtered)))
    species_to_color_filtered = {sp: palette_filtered[i] for i, sp in enumerate(sorted(set(species_filtered)))}
    colors_filtered = [species_to_color_filtered[sp] for sp in species_filtered]
    
    ax2.scatter(embedding_filtered[:, 0], embedding_filtered[:, 1], c=colors_filtered, s=15, alpha=0.7)
    ax2.set_title("Filtered (outliers removed)")
    ax2.set_xlabel("UMAP 1")
    ax2.set_ylabel("UMAP 2")
    ax2.grid(True)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    
    if show_plot:
        plt.show()
    else:
        plt.close()
    
    print(f"Comparison plot saved: {output_file}")


def main():
    """Main function to run the UMAP analysis."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument("--clusters", default="all_species_clusters.tsv",
                       help="Path to gene family clusters TSV file (default: all_species_clusters.tsv)")
    parser.add_argument("--panorama", default="all_prev_systems_tp_only",
                       help="Root directory containing species subdirectories (default: all_prev_systems_tp_only)")
    parser.add_argument("--systems-file", default="dbcan-merged-corrected/systems.tsv",
                       help="Relative path to systems file within each species directory")
    parser.add_argument("--output-dir", default="report",
                       help="Output directory for plots (default: report)")
    parser.add_argument("--analyze-outliers", type=int, default=5,
                       help="Number of top outliers to analyze (default: 5)")
    parser.add_argument("--show-plots", action="store_true",
                       help="Display plots interactively")
    
    args = parser.parse_args()
    
    # Create output directory
    Path(args.output_dir).mkdir(exist_ok=True)
    
    # Load data
    print("Loading cluster mappings...")
    cluster_map = load_cluster_mappings(args.clusters)
    
    print("Loading systems data...")
    all_systems, system_labels, species_labels = load_systems_data(
        args.root_dir, args.systems_file, cluster_map)
    
    if not all_systems:
        print("No systems found. Check your input paths.")
        return
    
    # Compute similarity and perform UMAP
    print("Computing similarity matrix...")
    distance_matrix = compute_similarity_matrix(all_systems)
    
    print("Performing UMAP...")
    embedding = perform_umap(distance_matrix, args.random_seed)
    
    # Analyze outliers
    if args.analyze_outliers > 0:
        analyze_outliers(embedding, all_systems, system_labels, species_labels, args.analyze_outliers)
    
    # Create original plot
    output_file = os.path.join(args.output_dir, "umap_original.png")
    create_umap_plot(embedding, species_labels, 
                    "UMAP Projection of Systems by Cluster Similarity", 
                    output_file, args.show_plots)
    
    # Identify and filter outliers
    outlier_mask, distances, threshold = identify_outliers(
        embedding, args.outlier_method, args.outlier_threshold, args.outlier_percentile)
    
    print(f"\nRemoving {np.sum(outlier_mask)} outliers (threshold: {threshold:.2f})")
    print(f"Plotting {np.sum(~outlier_mask)} systems")
    
    # Create filtered data
    embedding_filtered = embedding[~outlier_mask]
    species_filtered = [species_labels[i] for i in range(len(species_labels)) if not outlier_mask[i]]
    
    # Create filtered plot
    output_file_filtered = os.path.join(args.output_dir, "umap_filtered.png")
    create_umap_plot(embedding_filtered, species_filtered,
                    "UMAP Projection of Systems (Outliers Removed)",
                    output_file_filtered, args.show_plots)
    
    # Create comparison plot
    if not args.skip_comparison:
        comparison_file = os.path.join(args.output_dir, "umap_comparison.png")
        create_comparison_plot(embedding, embedding_filtered, species_labels, species_filtered,
                             comparison_file, args.show_plots)


if __name__ == "__main__":
    main()