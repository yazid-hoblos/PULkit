import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import MultiLabelBinarizer
from sklearn.cluster import KMeans
import seaborn as sns

# Load systems.tsv
df = pd.read_csv("all_prev_systems_tp_only/Prevotella_denticola/dbcan-merged-corrected/systems.tsv", sep="\t")

# Extract model_GF sets
gf_sets = [set(str(x).split(", ")) for x in df["model_GF"]]

# Convert to binary matrix
mlb = MultiLabelBinarizer()
X = mlb.fit_transform(gf_sets)

# PCA
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X)

# (Optional) Clustering
kmeans = KMeans(n_clusters=3, random_state=42)
labels = kmeans.fit_predict(X_pca)

# Plot
plt.figure(figsize=(8,6))
sns.scatterplot(x=X_pca[:,0], y=X_pca[:,1], hue=labels, palette='tab10', s=100)
plt.title("PCA of Systems (Model_GF)")
plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% var)")
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% var)")
plt.legend(title="Cluster")
plt.grid(True)
plt.tight_layout()
plt.show()
