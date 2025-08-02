import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from collections import defaultdict


if len(sys.argv) < 2:
    print("Usage: python script.py <input_systems_file.tsv>")
    sys.exit(1)

input_file = sys.argv[1]

# Read only columns 1 and 5 (0-based index: columns 1 and 5)
df = pd.read_csv(input_file, sep="\t", usecols=["system number", "model_GF"])
# Parse gene families into sets
systems = {
    row["system number"]: set(map(str.strip, str(row["model_GF"]).split(',')))
    for _, row in df.iterrows()
}
# Create similarity matrix (Modified Jaccard Index)
system_ids = list(systems.keys())
similarity_matrix = pd.DataFrame(index=system_ids, columns=system_ids, dtype=float)

for i in system_ids:
    for j in system_ids:
        inter = systems[i].intersection(systems[j])
        union = systems[i].union(systems[j])
        similarity = len(inter) / min(len(systems[i]), len(systems[j])) if min(len(systems[i]), len(systems[j])) > 0 else 0
        similarity_matrix.loc[i, j] = similarity

similarity_matrix.index = similarity_matrix.index.astype(str)
similarity_matrix.columns = similarity_matrix.columns.astype(str)

# Plot heatmap
sns.clustermap(
    similarity_matrix.astype(float),
    annot=False,
    fmt=".2f",
    cmap="YlGnBu",
    figsize=(10, 10),
    metric="euclidean",
    method="average",
    cbar_kws={"label": "Similarity"},  
    xticklabels=True,
    yticklabels=True)

plt.title("Overlap Jaccard")
plt.xlabel("Systems")
plt.ylabel("Systems")
plt.tight_layout()
plt.show()

# Identify high similarity pairs
threshold = 0.7  # 70%
high_similarity_pairs = []

for i, sys1 in enumerate(similarity_matrix.index):
    for j, sys2 in enumerate(similarity_matrix.columns):
        if j > i:
            similarity = similarity_matrix.loc[sys1, sys2]
            if similarity > threshold:
                high_similarity_pairs.append((sys1, sys2, similarity))

# Display results
for pair in high_similarity_pairs:
    print(f"System {pair[0]} - System {pair[1]}: Similarity = {pair[2]:.2f}")


# Hierarchical clustering
distance_matrix = 1 - similarity_matrix
condensed_distance = squareform(distance_matrix.values)
Z = linkage(condensed_distance, method='average')

# Assign cluster labels
cut_threshold = 0.3
cluster_labels = fcluster(Z, cut_threshold, criterion='distance')

# Group systems by cluster
clustered_systems = defaultdict(list)
for system_id, cluster_id in zip(similarity_matrix.index, cluster_labels):
    clustered_systems[cluster_id].append(system_id)

# Build supersystems
supersystems = {}
for cluster_id, system_ids in clustered_systems.items():
    gene_sets = [systems[int(sys_id)] for sys_id in system_ids]
    core_genes = set.intersection(*gene_sets)
    accessory_genes = set.union(*gene_sets) - core_genes

    supersystems[cluster_id] = {
        'core': core_genes,
        'accessory': accessory_genes,
        'members': system_ids
    }

# Display supersystems
print("\n === Supersystems === \n")
for cluster_id, data in supersystems.items():
    if len(data['members']) < 2:
        continue
    print(f"Cluster {cluster_id}:")
    print(f"  Members: {', '.join(data['members'])}")
    print(f"  Core genes: {', '.join(data['core']) if data['core'] else 'None'}")
    print(f"  Accessory genes: {', '.join(data['accessory']) if data['accessory'] else 'None'}")
    print()
