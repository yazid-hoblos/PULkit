import pandas as pd

# Load full dataset
df = pd.read_csv("PUL_gene_metadata.csv")

# Count number of components per PUL
component_counts = df["PUL_ID"].value_counts()

# Find PUL IDs with 4 or fewer components
small_puls = component_counts[component_counts <= 4].index

# Filter original dataframe for those PULs
filtered_df = df[df["PUL_ID"].isin(small_puls)]

# Save or inspect
filtered_df.to_csv("PULs_with_4_or_less_components.csv", index=False)
print(filtered_df.head())
