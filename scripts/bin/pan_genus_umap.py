import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import umap
import seaborn as sns

# === Paths ===
clusters_tsv = "all_species_clusters.tsv"
root_dir = "all_prev_systems_tp_only"
systems_filename = "dbcan-merged-corrected/systems.tsv"

# === Read cluster file ===
cluster_map = {}
with open(clusters_tsv) as f:
    for line in f:
        cluster_id, gf_id = line.strip().split()
        cluster_map[gf_id] = cluster_id

print(f"Loaded {len(cluster_map)} gene family to cluster mappings")

def compute_similarity(set1, set2):
    if not set1 or not set2:
        return 0
    intersection = set1 & set2
    return len(intersection) / min(len(set1), len(set2))/2

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

# Compute similarity matrix
n = len(all_systems)
sim_matrix = np.zeros((n, n))
for i in range(n):
    for j in range(i, n):
        sim = compute_similarity(all_systems[i], all_systems[j])
        sim_matrix[i, j] = sim
        sim_matrix[j, i] = sim

# Convert similarity to distance (1 - similarity)
dist_matrix = 1 - sim_matrix

# UMAP dimensionality reduction
reducer = umap.UMAP(metric="precomputed", random_state=42)
embedding = reducer.fit_transform(dist_matrix)

# Plot
plt.figure(figsize=(14, 10))
palette = sns.color_palette("tab20", n_colors=len(set(species_labels)))
species_to_color = {sp: palette[i] for i, sp in enumerate(sorted(set(species_labels)))}

colors = [species_to_color[sp] for sp in species_labels]

plt.scatter(embedding[:, 0], embedding[:, 1], c=colors, s=15, alpha=0.7)
plt.title("UMAP projection of Prevotella Systems by Cluster Similarity")
plt.xlabel("UMAP 1")
plt.ylabel("UMAP 2")
plt.grid(True)

# Legend for species
handles = [plt.Line2D([0], [0], marker='o', color='w', label=sp, markerfacecolor=species_to_color[sp], markersize=8)
           for sp in species_to_color]
plt.legend(handles=handles, title="Species", bbox_to_anchor=(1.05, 1), loc="upper left")

plt.tight_layout()
plt.savefig("report/prevotella_systems_umap.png", dpi=300)
plt.show()
print("UMAP plot saved to: report/prevotella_systems_umap.png")

# Add this code after your UMAP embedding generation and before plotting

# Statistical outlier detection using distance from center
# center_x, center_y = np.mean(embedding[:, 0]), np.mean(embedding[:, 1])
# distances_from_center = np.sqrt((embedding[:, 0] - center_x)**2 + (embedding[:, 1] - center_y)**2)

# # Define outliers as points beyond 2 standard deviations from mean distance
# threshold = np.mean(distances_from_center) + 2 * np.std(distances_from_center)
# outlier_indices = np.where(distances_from_center > threshold)[0]

# print(f"\n=== OUTLIER SYSTEMS (>{threshold:.2f} units from center) ===")
# print(f"Found {len(outlier_indices)} outlier systems:")

# for idx in outlier_indices:
#     species = species_labels[idx]
#     system = system_labels[idx]
#     distance = distances_from_center[idx]
#     x, y = embedding[idx]
    
#     # Show the actual gene clusters for this system
#     clusters = all_systems[idx]
#     cluster_list = sorted(list(clusters))
    
#     print(f"\nDistance: {distance:.2f} | Coords: ({x:.2f}, {y:.2f})")
#     print(f"Species: {species}")
#     print(f"System: {system}")
#     print(f"Gene clusters: {cluster_list}")

# # Summary by species
# outlier_species = [species_labels[idx] for idx in outlier_indices]
# species_counts = {}
# for sp in outlier_species:
#     species_counts[sp] = species_counts.get(sp, 0) + 1

# print(f"\n=== OUTLIER SUMMARY BY SPECIES ===")
# for species, count in sorted(species_counts.items()):
#     print(f"{species}: {count} outlier systems")

# Add this code after your UMAP embedding generation and before plotting

# Find top 5 outliers based on distance from center
center_x, center_y = np.mean(embedding[:, 0]), np.mean(embedding[:, 1])
distances_from_center = np.sqrt((embedding[:, 0] - center_x)**2 + (embedding[:, 1] - center_y)**2)

# Get indices of top 5 most distant points
top_5_outlier_indices = np.argsort(distances_from_center)[-5:][::-1]  # Sort descending

print(f"\n=== TOP 5 OUTLIER SYSTEMS (most distant from center) ===")
print(f"Center coordinates: ({center_x:.2f}, {center_y:.2f})")

for i, idx in enumerate(top_5_outlier_indices, 1):
    species = species_labels[idx]
    system = system_labels[idx]
    distance = distances_from_center[idx]
    x, y = embedding[idx]
    
    # Show the actual gene clusters for this system
    clusters = all_systems[idx]
    cluster_list = sorted(list(clusters))
    
    print(f"\n{i}. Distance from center: {distance:.2f}")
    print(f"   Coordinates: ({x:.2f}, {y:.2f})")
    print(f"   Species: {species}")
    print(f"   System: {system}")
    print(f"   Gene clusters ({len(clusters)} total): {cluster_list}")

# Quick summary
outlier_species = [species_labels[idx] for idx in top_5_outlier_indices]
print(f"\n=== SUMMARY ===")
print(f"Species represented: {set(outlier_species)}")
print(f"Total unique gene clusters across top 5 outliers: {len(set().union(*[all_systems[idx] for idx in top_5_outlier_indices]))}")



# Add this code after your UMAP embedding generation to create a plot without outliers

# Calculate outliers
center_x, center_y = np.mean(embedding[:, 0]), np.mean(embedding[:, 1])
distances_from_center = np.sqrt((embedding[:, 0] - center_x)**2 + (embedding[:, 1] - center_y)**2)

# Define outlier threshold (adjust as needed)
# Method 1: Statistical threshold (2 standard deviations)
threshold = np.mean(distances_from_center) + 2 * np.std(distances_from_center)

# Method 2: Percentile threshold (top 5% most distant)
# threshold = np.percentile(distances_from_center, 95)

outlier_mask = distances_from_center > threshold
non_outlier_mask = ~outlier_mask

print(f"Removing {np.sum(outlier_mask)} outliers out of {len(embedding)} total systems")
print(f"Plotting {np.sum(non_outlier_mask)} systems")

# Create filtered data
embedding_filtered = embedding[non_outlier_mask]
species_labels_filtered = [species_labels[i] for i in range(len(species_labels)) if non_outlier_mask[i]]

# Plot without outliers
plt.figure(figsize=(14, 10))
palette = sns.color_palette("tab20", n_colors=len(set(species_labels_filtered)))
species_to_color = {sp: palette[i] for i, sp in enumerate(sorted(set(species_labels_filtered)))}
colors = [species_to_color[sp] for sp in species_labels_filtered]

plt.scatter(embedding_filtered[:, 0], embedding_filtered[:, 1], c=colors, s=15, alpha=0.7)
plt.title("UMAP projection of Prevotella Systems by Cluster Similarity (Outliers Removed)")
plt.xlabel("UMAP 1")
plt.ylabel("UMAP 2")
plt.grid(True)

# Legend for species
handles = [plt.Line2D([0], [0], marker='o', color='w', label=sp, 
                     markerfacecolor=species_to_color[sp], markersize=8)
           for sp in species_to_color]
plt.legend(handles=handles, title="Species", bbox_to_anchor=(1.05, 1), loc="upper left")
plt.tight_layout()
plt.savefig("report/prevotella_systems_umap_no_outliers.png", dpi=300)
plt.show()

# Optional: Create a side-by-side comparison
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))

# Original plot
palette_all = sns.color_palette("tab20", n_colors=len(set(species_labels)))
species_to_color_all = {sp: palette_all[i] for i, sp in enumerate(sorted(set(species_labels)))}
colors_all = [species_to_color_all[sp] for sp in species_labels]

ax1.scatter(embedding[:, 0], embedding[:, 1], c=colors_all, s=15, alpha=0.7)
ax1.set_title("Original (with outliers)")
ax1.set_xlabel("UMAP 1")
ax1.set_ylabel("UMAP 2")
ax1.grid(True)

# Filtered plot
ax2.scatter(embedding_filtered[:, 0], embedding_filtered[:, 1], c=colors, s=15, alpha=0.7)
ax2.set_title("Filtered (outliers removed)")
ax2.set_xlabel("UMAP 1")
ax2.set_ylabel("UMAP 2")
ax2.grid(True)

plt.tight_layout()
# plt.savefig("report/prevotella_systems_comparison.png", dpi=300)
plt.show()

print("Plots saved:")
print("- Single filtered plot: report/prevotella_systems_umap_no_outliers.png")
print("- Side-by-side comparison: report/prevotella_systems_comparison.png")