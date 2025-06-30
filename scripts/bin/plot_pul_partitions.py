import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import argparse

def plot_summary_partition_and_counts(systems_file, output_dir="partition_summary_plots"):
    os.makedirs(output_dir, exist_ok=True)
    df = pd.read_csv(systems_file, sep="\t")

    required_cols = {'partition', 'strict', 'extended', 'split'}
    if not required_cols.issubset(df.columns):
        raise ValueError(f"Missing required columns in file: {required_cols}")

    species_name = Path(systems_file).parent.parent.name.replace("_", " ")
    ### --- Plot 1: Frequency of Partition Types ---
    plt.figure(figsize=(8, 5))
    sns.countplot(data=df, x='partition', order=df['partition'].value_counts().index, palette='Set2')
    plt.title(f"Partition Type Frequency – {species_name}", fontsize=13, weight='bold')
    plt.xlabel("Partition Type")
    plt.ylabel("Number of Systems")
    plt.xticks(rotation=45)
    plt.tight_layout()
    out_path1 = os.path.join(output_dir, f"{species_name}_partition_frequency.png")
    plt.savefig(out_path1, dpi=300)
    plt.show()
    plt.close()
    print(f"Saved: {out_path1}")

    ### --- Plot 2: Mean Strict / Extended / Split per system ---
    mean_values = df[['strict', 'extended', 'split']].mean().reset_index()
    mean_values.columns = ['Type', 'Mean Count']

    plt.figure(figsize=(6, 4))
    sns.barplot(data=mean_values, x='Type', y='Mean Count', palette='pastel')
    plt.title(f"Average Gene Type Count – {species_name}", fontsize=13, weight='bold')
    plt.ylabel("Mean Count per System")
    plt.tight_layout()
    out_path2 = os.path.join(output_dir, f"{species_name}_mean_unit_counts.png")
    plt.savefig(out_path2, dpi=300)
    plt.show()
    plt.close()
    print(f"Saved: {out_path2}")

def main():
    parser = argparse.ArgumentParser(description="Plot partition and unit types from a single systems.tsv file")
    parser.add_argument("--systems", required=True, help="Path to systems.tsv file")
    parser.add_argument("-o", "--output", default="partition_plots", help="Output directory for plots")
    args = parser.parse_args()

    plot_summary_partition_and_counts(args.systems, args.output)

if __name__ == "__main__":
    main()