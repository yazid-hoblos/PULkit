import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from sklearn.cluster import DBSCAN
from sklearn.metrics import pairwise_distances
import seaborn as sns

# Load and parse the systems.tsv file
df = pd.read_csv("all_prev_systems_tp_only/Prevotella_denticola/dbcan-merged-corrected/systems.tsv", sep="\t")
gf_sets = [set(str(x).split(", ")) for x in df["model_GF"]]  # or use context_GF

# Compute pairwise similarity
n = len(gf_sets)
similarity = np.zeros((n, n))
for i in range(n):
    for j in range(n):
        inter = len(gf_sets[i].intersection(gf_sets[j]))
        min_len = min(len(gf_sets[i]), len(gf_sets[j]))
        similarity[i, j] = inter / min_len if min_len > 0 else 0

# Convert similarity to distance
distance = 1 - similarity

# Dimensionality reduction
tsne = TSNE(metric="precomputed", random_state=42)
# embedding = tsne.fit_transform(distance)
import umap.umap_ as umap
embedding = umap.UMAP(metric="precomputed", random_state=42).fit_transform(distance)
# Clustering (e.g. DBSCAN)
clustering = DBSCAN(eps=0.3, min_samples=2, metric="precomputed").fit(distance)

# Plot
plt.figure(figsize=(8, 6))
sns.scatterplot(x=embedding[:, 0], y=embedding[:, 1], hue=clustering.labels_, palette='tab10', s=100)
plt.title("PUL System Clustering (Model_GF)")
plt.xlabel("TSNE Dim 1")
plt.ylabel("TSNE Dim 2")
plt.legend(title="Cluster")
plt.grid(True)
plt.tight_layout()
# plt.show()
plt.savefig("report/pul_systems_tsne_clustering.png", dpi=500)
