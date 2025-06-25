import pandas as pd
import matplotlib.pyplot as plt

# Load TSV (adjust path to your file)
df = pd.read_csv("all_prev_systems_tp_only/Prevotella_denticola/dbcan-merged-corrected/systems.tsv", sep='\t')

# Inspect columns
print(df.columns)

# Convert system number to int for plotting
df['system number'] = df['system number'].astype(int)

# Extract system number and completeness
systems = df['system number']
completeness = df['completeness']

# Plot
plt.figure(figsize=(10, 5))
plt.plot(systems, completeness, marker='o', linestyle='-', color='teal')
plt.xlabel("System Number")
plt.ylabel("Degree of Similarity (Completeness)")
plt.title("Degree of Similarity between PUL Systems")
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.show()
