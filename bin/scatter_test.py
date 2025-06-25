import matplotlib.pyplot as plt
import pandas as pd

def plot_scatter(summary_file):
    df = pd.read_csv(summary_file, sep='\t')

    # Calculate % overlaps relative to system sizes
    df['overlap_pct_cgc'] = df['matched_genes'] / df['cgc_total_genes']
    df['overlap_pct_cazy'] = df['matched_genes'] / df['cazy_total_genes']

    plt.figure(figsize=(10, 8))
    scatter = plt.scatter(df['cgc_total_genes'], df['cazy_total_genes'],
                          s=df['matched_genes']*10,  # scale point size by matched genes
                          c=df['overlap_pct_cgc'], cmap='viridis', alpha=0.7)
    plt.colorbar(scatter, label='% Overlap relative to CGC system size')
    plt.xlabel('CGC System Size (Total Genes)')
    plt.ylabel('CAZy System Size (Total Genes)')
    plt.title('System Size vs Matched Genes')
    plt.savefig("system_size_vs_matched_genes.png")

if __name__ == "__main__":
    plot_scatter("system_gene_pair_overlap_summary.tsv")
