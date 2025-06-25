import pandas as pd
import matplotlib.pyplot as plt

def compute_pul_size_freq(path, label_prefix):
    df = pd.read_csv(path, sep="\t", header=None)
    df.columns = ["PUL_label", "GeneID", "Contig", "Start", "Stop", "Strand", "Locus"]
    
    # Extract PUL number (works for both formats)
    df['PUL_num'] = df['PUL_label'].str.extract(r'(?:Predicted|Literature-derived) PUL\s*(\d+)')
    df['PUL_id'] = df['PUL_num'] + "_" + df['Contig']
    
    # Compute sizes
    pul_sizes = df.groupby('PUL_id').size()
    size_freq = pul_sizes.value_counts().sort_index()
    return size_freq

# Load both datasets
file1 = "cazy_curated/Prevotella_denticola_merged.tsv"
file2 = "Bt_empirical_PULs.tsv"

freq1 = compute_pul_size_freq(file1, "Predicted")
freq2 = compute_pul_size_freq(file2, "Literature-derived")

# Combine both distributions into a DataFrame for plotting
all_sizes = sorted(set(freq1.index).union(freq2.index))
df_plot = pd.DataFrame({
    'Predicted': freq1.reindex(all_sizes, fill_value=0),
    'Literature-derived': freq2.reindex(all_sizes, fill_value=0)
}, index=all_sizes)

# Plot
plt.figure(figsize=(10,6))
df_plot.plot(kind='bar', width=0.8, color=['skyblue', 'salmon'], edgecolor='black')
plt.title("Frequency of Unique PUL Sizes (Predicted vs Literature-derived)")
plt.xlabel("Number of Genes per PUL")
plt.ylabel("Number of PULs")
plt.xticks(rotation=0)
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.legend(title="Dataset")
plt.tight_layout()
# plt.show()
plt.savefig("latex/figures/report/pul_size_distribution_comparison.png", dpi=500)
