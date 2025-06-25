import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load your data file
df = pd.read_csv("files/system_gene_pair_overlap_summary.tsv", sep='\t')

# Calculate coverage of CAZy system by matched genes
df['coverage'] = df['matched_genes'] / df['cazy_total_genes']

# Define thresholds to test
thresholds = np.arange(0, 1.01, 0.05)  # from 0 to 1 in steps of 0.05

# Prepare dataframe to hold counts
genomes = df['genome'].unique()
coverage_counts = pd.DataFrame(index=thresholds, columns=genomes)

for genome in genomes:
    genome_df = df[df['genome'] == genome]
    for thresh in thresholds:
        # Select rows where coverage >= threshold
        passing = genome_df[genome_df['coverage'] >= thresh]

        # Count unique CAZy systems per genome by using the combined key (genome, cazy_system)
        # Here it's within one genome, so just unique cazy_system strings within that genome
        coverage_counts.at[thresh, genome] = passing['cazy_system'].nunique()

# Convert counts to int
coverage_counts = coverage_counts.fillna(0).astype(int)

# Plotting
plt.figure(figsize=(12, 7))
for genome in genomes:
    plt.plot(coverage_counts.index, coverage_counts[genome], marker='o', label=genome, alpha=0.3)

plt.xlabel("Coverage Threshold of CAZy System by PANORAMA System")
plt.ylabel("Number of CAZy Systems Predicted by PANORAMA")
plt.title("Coverage of CAZy Systems by PANORAMA Predictions")
plt.legend(title="Genome ID", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True)
plt.tight_layout()
plt.savefig("report/coverage_plot.png", dpi=500)
plt.show()
