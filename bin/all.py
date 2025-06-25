import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# === Paths ===
clusters_tsv = "all_species_clusters.tsv"
root_dir = "all_prev_systems_tp_only"
systems_filename = "dbcan-merged-corrected/systems.tsv"

# === Read cluster file ===
# Map each gene family ID -> cluster ID
cluster_map = {}
with open(clusters_tsv) as f:
    for line in f:
        cluster_id, gf_id = line.strip().split()
        cluster_map[gf_id] = cluster_id

print(f"Loaded {len(cluster_map)} gene family to cluster mappings")

# === Similarity function ===
def compute_similarity(set1, set2):
    if not set1 or not set2:
        return 0
    intersection = set1 & set2
    return len(intersection) / min(len(set1), len(set2))/2

# === Load all systems and convert model_GF to cluster sets ===
all_systems = []
system_labels = []

for species_dir in sorted(os.listdir(root_dir)):
    system_path = os.path.join(root_dir, species_dir, systems_filename)
    if not os.path.exists(system_path):
        print(f"Skipping {system_path} (not found)")
        continue

    df = pd.read_csv(system_path, sep='\t')
    for idx, row in df.iterrows():
        gfs = str(row.get("model_GF", ""))
        # Convert gene families to clusters
        gf_list = [gf.strip() for gf in gfs.split(",") if gf.strip()]
        clusters = set(cluster_map.get(gf) for gf in gf_list if gf in cluster_map)
        if clusters:
            all_systems.append(clusters)
            system_labels.append(f"{species_dir}#{idx}")

print(f"Processed {len(all_systems)} systems with clusters")

# === Compute similarity matrix ===
n = len(all_systems)
sim_matrix = np.zeros((n, n))

for i in range(n):
    for j in range(n):
        sim_matrix[i, j] = compute_similarity(all_systems[i], all_systems[j])

# === Distance matrix for clustering ===
dist_matrix = sim_matrix
np.fill_diagonal(dist_matrix, 0)

# === Plot clustered heatmap ===
sns.set(style="white")
clust = sns.clustermap(dist_matrix,
                       figsize=(14, 14),
                       cmap="viridis",
                       xticklabels=False,
                       yticklabels=False,
                       metric="euclidean",
                       method="average")

clust.ax_heatmap.set_title("Clustered System Distance Heatmap (1 - Cluster Similarity)")
plt.tight_layout()
os.makedirs("report", exist_ok=True)
clust.savefig("report/prevotella_system_cluster_heatmap.png", dpi=300)
# plt.close()
# clust.show()
print("Heatmap saved to: report/prevotella_system_cluster_heatmap.png")
