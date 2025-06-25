#!/usr/bin/env python3

"""
Perform PCA analysis on systems.tsv gene family data and generate visualization.

Usage:
    python systems_pca_cluseting.py [-h] --input DIR [-o OUTPUT] [-n N_CLUSTERS]

Arguments:
    --systems DIR      Directory containing systems.tsv file
    -o/--output     Output file path for plot (default: pca_systems.png)
    -n/--clusters   Number of clusters for K-means (default: 3)
"""

import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import MultiLabelBinarizer
from sklearn.cluster import KMeans
import seaborn as sns

def find_optimal_clusters(X, output_dir, max_clusters=10):
    """Determine optimal number of clusters using elbow method."""
    inertias = []
    
    for k in range(1, max_clusters + 1):
        kmeans = KMeans(n_clusters=k, random_state=42)
        kmeans.fit(X)
        inertias.append(kmeans.inertia_)
    
    # Calculate the rate of change
    diffs = np.diff(inertias)
    rates_of_change = np.diff(diffs)
    
    # Find the elbow point (where rate of change is maximum)
    optimal_clusters = np.argmin(rates_of_change) + 2
    
    # Plot elbow curve
    plt.figure(figsize=(8, 4))
    plt.plot(range(1, max_clusters + 1), inertias, 'bo-')
    plt.axvline(x=optimal_clusters, color='r', linestyle='--', label=f'Optimal k={optimal_clusters}')
    plt.xlabel('Number of Clusters (k)')
    plt.ylabel('Inertia')
    plt.title('Elbow Method for Optimal k')
    plt.legend()
    plt.grid(True)
    plt.savefig(os.path.join(output_dir, "elbow_method.png"), dpi=300, bbox_inches='tight')
    plt.close()
    
    return optimal_clusters

def run_pca_analysis(input_dir, n_clusters=None, output_path="pca_systems.png"):
    # Load systems.tsv
    systems_file = os.path.join(input_dir, "systems.tsv")
    if not os.path.exists(systems_file):
        raise FileNotFoundError(f"systems.tsv not found in {input_dir}")
    
    os.makedirs(output_path, exist_ok=True)
    
    df = pd.read_csv(systems_file, sep="\t")

    # Extract model_GF sets
    gf_sets = [set(str(x).split(", ")) for x in df["model_GF"]]

    # Convert to binary matrix
    mlb = MultiLabelBinarizer()
    X = mlb.fit_transform(gf_sets)

    # PCA
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X)

    # Determine number of clusters if not specified
    if n_clusters is None:
        n_clusters = find_optimal_clusters(X_pca, output_path)
        print(f"Automatically determined optimal number of clusters: {n_clusters}")

    # Clustering
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    labels = kmeans.fit_predict(X_pca)

    # Plot
    plt.figure(figsize=(8,6))
    scatter = sns.scatterplot(x=X_pca[:,0], y=X_pca[:,1], 
                            hue=labels, palette='tab10', s=100)
    plt.title(f"PCA of Systems (Model_GF) - {n_clusters} clusters")
    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% var)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% var)")
    plt.legend(title="Cluster")
    plt.grid(True)
    plt.tight_layout()
    
    # Create output directory if it doesn't exist
    plt.savefig(os.path.join(output_path,"pca_plot.png"), dpi=500, bbox_inches='tight')
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Perform PCA analysis on systems.tsv gene family data")
    parser.add_argument("--systems", required=True, help="Directory containing systems.tsv file")
    parser.add_argument("-o", "--output", default="pca_systems", 
                       help="Output file path for plot")
    parser.add_argument("-n", "--clusters", type=int, default=None, 
                       help="Number of clusters (optional, will auto-determine if not specified)")
    
    args = parser.parse_args()
    run_pca_analysis(args.systems, args.clusters, args.output)

if __name__ == "__main__":
    main()