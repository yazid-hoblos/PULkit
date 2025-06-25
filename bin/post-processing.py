import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Load data
df = pd.read_csv("Bt_sys/s__Bacteroides_thetaiotaomicron/dbcan-merged/systems.tsv", sep='\t')
# df = pd.read_csv("systems.tsv", sep='\t')


# Convert to sets
model_sets = [set(str(x).split(', ')) for x in df['model_GF']]
context_sets = [set(str(x).split(', ')) for x in df['context_GF']]

n = len(df)
print(f"Total systems: {n}")
model_sim = np.zeros((n, n))
context_sim = np.zeros((n, n))

# Custom similarity: intersection over min size
def overlap_min_similarity(set1, set2):
    intersection = set1 & set2
    min_len = min(len(set1), len(set2))
    return len(intersection) / min_len if min_len > 0 else 0

# Compute similarity matrices
for i in range(n):
    print(f"Processing system {i+1}")
    for j in range(n):
        model_sim[i, j] = overlap_min_similarity(model_sets[i], model_sets[j])
        context_sim[i, j] = overlap_min_similarity(context_sets[i], context_sets[j])
        # print(f"Similarity ({i+1}, {j+1}): Model={model_sim[i, j]:.3f}, Context={context_sim[i, j]:.3f}")

# Generate clustermaps but do NOT show, just get the reordered data and indices
model_cluster = sns.clustermap(model_sim, cmap='Blues', xticklabels=False, yticklabels=False, method='average', metric='euclidean')
plt.close(model_cluster.fig)  # Close to prevent auto display

context_cluster = sns.clustermap(context_sim, cmap='Oranges', xticklabels=False, yticklabels=False, method='average', metric='euclidean')
plt.close(context_cluster.fig)  # Close to prevent auto display

# Get reordered indices from dendrogram leaves
model_idx = model_cluster.dendrogram_row.reordered_ind
context_idx = context_cluster.dendrogram_row.reordered_ind

# Reorder matrices accordingly
model_reordered = model_sim[np.ix_(model_idx, model_idx)]
context_reordered = context_sim[np.ix_(context_idx, context_idx)]

# Plot both reordered heatmaps side by side in one figure
fig, axes = plt.subplots(1, 2, figsize=(16, 8))

sns.heatmap(model_reordered, cmap='Blues', ax=axes[0], cbar=True, xticklabels=False, yticklabels=False)
axes[0].set_title("Model_GF Similarity (Clustered)")

sns.heatmap(context_reordered, cmap='Oranges', ax=axes[1], cbar=True, xticklabels=False, yticklabels=False)
axes[1].set_title("Context_GF Similarity (Clustered)")

plt.tight_layout()
# plt.savefig("report/Bt_pul_similarity_clustered_combined.png", dpi=500)
plt.show()
