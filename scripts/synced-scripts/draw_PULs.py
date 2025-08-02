import pandas as pd
import matplotlib.pyplot as plt

# Load the extracted data
df = pd.read_csv("PUL_gene_metadata.csv")

# Count number of components (genes) per PUL
component_counts = df["PUL_ID"].value_counts().sort_index()

# Plot histogram of how many PULs have N components
plt.figure(figsize=(10, 6))
component_counts.value_counts().sort_index().plot(kind='bar', color='skyblue', edgecolor='black')

plt.title("Frequency of Number of Components per PUL")
plt.xlabel("Number of Components")
plt.ylabel("Number of PULs")
plt.xticks(rotation=0)
plt.grid(axis='y', linestyle='--', alpha=0.7)

plt.tight_layout()
plt.show()
