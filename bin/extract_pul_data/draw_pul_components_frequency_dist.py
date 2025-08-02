import pandas as pd
import matplotlib.pyplot as plt
import sys 

if len(sys.argv) != 2:
    print("Usage: python draw_PULs.py <dbcan_pul_components.csv>")
    sys.exit(1)

puls_file = sys.argv[1]

# Load the extracted data
df = pd.read_csv(puls_file, sep=',')

# Count number of components (genes) per PUL
component_counts = df["PUL_ID"].value_counts().sort_index()

# Count how many PULs have N components
distribution = component_counts.value_counts().sort_index()

# Plot frequency distribution of number of components per PUL
plt.figure(figsize=(10, 6))
bars = plt.bar(distribution.index, distribution.values, color="#4A90E2", edgecolor='black', width=0.7)

# Add value labels above bars
for bar in bars:
    height = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2, height + 0.5, f'{int(height)}', ha='center', va='bottom', fontsize=8)

plt.title("Frequency of Number of Components per PUL", fontsize=14, weight='bold')
plt.xlabel("Number of Components", fontsize=12)
plt.ylabel("Number of PULs", fontsize=12)
plt.xticks(distribution.index)
plt.grid(axis='y', linestyle='--', alpha=0.5)

plt.tight_layout()
plt.show()