import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load TSV file
# Replace 'your_file.tsv' with your actual path
df = pd.read_csv("super_systems", sep="\t", header=None, names=["system_id", "families"])

# Parse gene families into sets
systems = {
    row["system_id"]: set(map(str.strip, row["families"].split(',')))
    for _, row in df.iterrows()
}

# Create similarity matrix (Jaccard Index)
system_ids = list(systems.keys())
similarity_matrix = pd.DataFrame(index=system_ids, columns=system_ids, dtype=float)

for i in system_ids:
    for j in system_ids:
        inter = systems[i].intersection(systems[j])
        union = systems[i].union(systems[j])
        similarity = len(inter) / min(len(systems[i]), len(systems[j]))
        similarity_matrix.loc[i, j] = similarity


similarity_matrix.index = similarity_matrix.index.astype(str)
similarity_matrix.columns = similarity_matrix.columns.astype(str)


# Plot heatmap
plt.figure(figsize=(10, 8))
# sns.clustermap(similarity_matrix.astype(float), annot=False, square=True, cmap="viridis", linewidths=0.5, linecolor='black', cbar_kws={"label": "Jaccard Index"})
g=sns.clustermap(
    similarity_matrix.astype(float),
    annot=False,               # Show similarity values
    fmt=".2f",                # Format values
    cmap="YlGnBu",            # Color map
    figsize=(10, 10),         # Size of the figure
    metric="euclidean",       # Distance metric
    method="average",         # Linkage method
    cbar_kws={"label": "Similarity"},
    xticklabels=True,  # Show x-axis labels
    yticklabels=True,  # Show y-axis labels
)



plt.title("System Gene Family Similarity (Jaccard Index)")
plt.xlabel("System")
plt.ylabel("System")
plt.tight_layout()
plt.show()



# import pandas as pd
# from scipy.cluster.hierarchy import linkage, fcluster
# from scipy.spatial.distance import squareform

# # Assume `similarity_matrix` is already computed as in your earlier code
# # Use a distance matrix (1 - similarity) for clustering
# distance_matrix = 1 - similarity_matrix.astype(float)

# # Convert to condensed distance matrix for linkage
# condensed_dist = squareform(distance_matrix.values, checks=False)

# # Perform hierarchical clustering
# Z = linkage(condensed_dist, method='average')  # or 'ward', 'single', etc.

# threshold = 0.5  # distance threshold (1 - similarity)
# cluster_labels = fcluster(Z, threshold, criterion='distance')

# # Map cluster labels back to system IDs
# system_ids = similarity_matrix.index.tolist()
# clusters = {}
# for sys_id, label in zip(system_ids, cluster_labels):
#     clusters.setdefault(label, []).append(sys_id)

# # Display resulting clusters
# for cluster_id, members in clusters.items():
#     print(f"Cluster {cluster_id}: {members}")



threshold = 0.7  # 80%
high_similarity_pairs = []

# Iterate over upper triangle to avoid duplicates
for i, sys1 in enumerate(similarity_matrix.index):
    for j, sys2 in enumerate(similarity_matrix.columns):
        if j > i:
            similarity = similarity_matrix.loc[sys1, sys2]
            if similarity > threshold:
                high_similarity_pairs.append((sys1, sys2, similarity))

# Display results
for pair in high_similarity_pairs:
    print(f"System {pair[0]} - System {pair[1]}: Similarity = {pair[2]:.2f}")



from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

# Convert similarity matrix to distance (1 - similarity)
distance_matrix = 1 - similarity_matrix
condensed_distance = squareform(distance_matrix.values)

# Hierarchical clustering
Z = linkage(condensed_distance, method='average')  # or 'ward', 'complete'

# Cut tree at desired threshold
threshold = 0.3  # Max distance for systems to be in the same cluster (i.e. 80% similarity)
cluster_labels = fcluster(Z, threshold, criterion='distance')

# Group systems by cluster
from collections import defaultdict
clustered_systems = defaultdict(list)
for system_id, cluster_id in zip(similarity_matrix.index, cluster_labels):
    clustered_systems[cluster_id].append(system_id)

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

print("Supersystems:")
for cluster_id, data in supersystems.items():
    if len(data['members']) < 2:
        continue
    print(f"Cluster {cluster_id}:")
    print(f"  Members: {', '.join(data['members'])}")
    print(f"  Core genes: {', '.join(data['core']) if data['core'] else 'None'}")
    print(f"  Accessory genes: {', '.join(data['accessory']) if data['accessory'] else 'None'}")
    print()