#!/usr/bin/env python3
"""
PUL Size Distribution Analysis

This script analyzes the size distribution of Polysaccharide Utilization Loci (PULs)
from TSV files and creates comparative bar plots showing the frequency of different
PUL sizes across datasets.
"""

import pandas as pd
import matplotlib.pyplot as plt
import argparse
from pathlib import Path


def compute_pul_size_freq(file_path):
    """
    Compute PUL size frequency from a TSV file.
    
    Args:
        file_path: Path to TSV file with PUL data
        
    Returns:
        Series with PUL sizes as index and frequencies as values
    """
    df = pd.read_csv(file_path, sep="\t", header=None)
    df.columns = ["PUL_label", "GeneID", "Contig", "Start", "Stop", "Strand", "Locus"]
    
    # Extract PUL number (works for both predicted and literature-derived formats)
    df['PUL_num'] = df['PUL_label'].str.extract(r'(?:Predicted|Literature-derived) PUL\s*(\d+)')
    df['PUL_id'] = df['PUL_num'] + "_" + df['Contig']
    
    # Compute sizes (number of genes per PUL)
    pul_sizes = df.groupby('PUL_id').size()
    size_freq = pul_sizes.value_counts().sort_index()
    
    return size_freq


def create_comparison_plot(datasets, labels, output_file, title=None, figsize=(10, 6), 
                          colors=None, show_plot=False):
    """
    Create a comparative bar plot of PUL size distributions.
    
    Args:
        datasets: List of frequency Series
        labels: List of dataset labels
        output_file: Output file path
        title: Plot title
        figsize: Figure size tuple
        colors: List of colors for each dataset
        show_plot: Whether to display the plot
    """
    # Get all unique sizes across datasets
    all_sizes = sorted(set().union(*[freq.index for freq in datasets]))
    
    # Create DataFrame for plotting
    df_plot = pd.DataFrame(index=all_sizes)
    for dataset, label in zip(datasets, labels):
        df_plot[label] = dataset.reindex(all_sizes, fill_value=0)
    
    # Set default colors if not provided
    if colors is None:
        colors = ['skyblue', 'salmon', 'lightgreen', 'orange', 'purple'][:len(labels)]
    
    # Create plot
    plt.figure(figsize=figsize)
    df_plot.plot(kind='bar', width=0.8, color=colors[:len(labels)], 
                edgecolor='black', ax=plt.gca())
    
    # Customize plot
    if title is None:
        title = "PUL Size Distribution Comparison"
    plt.title(title)
    plt.xlabel("Number of Genes per PUL")
    plt.ylabel("Number of PULs")
    plt.xticks(rotation=0)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.legend(title="Dataset")
    plt.tight_layout()
    
    # Save plot
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    
    if show_plot:
        plt.show()
    else:
        plt.close()
    
    print(f"Plot saved: {output_file}")


def print_summary(datasets, labels):
    """Print summary statistics for each dataset."""
    print("\n=== PUL SIZE SUMMARY ===")
    for dataset, label in zip(datasets, labels):
        total_puls = dataset.sum()
        min_size = dataset.index.min()
        max_size = dataset.index.max()
        mean_size = sum(size * freq for size, freq in dataset.items()) / total_puls
        
        print(f"\n{label}:")
        print(f"  Total PULs: {total_puls}")
        print(f"  Size range: {min_size}-{max_size} genes")
        print(f"  Mean size: {mean_size:.1f} genes")
        print(f"  Most common size: {dataset.idxmax()} genes ({dataset.max()} PULs)")


def main():
    """Main function to analyze PUL size distributions."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument("files", nargs='+',
                       help="Input TSV files containing PUL data")
    parser.add_argument("-l", "--labels", nargs='+',
                       help="Labels for each dataset (default: use filenames)")
    parser.add_argument("-o", "--output", default="report/pul_size_distribution.png",
                       help="Output file path (default: report/pul_size_distribution.png)")
    parser.add_argument("--show", action="store_true",
                       help="Display the plot interactively")
    parser.add_argument("--summary", action="store_true",
                       help="Print summary statistics")
    
    args = parser.parse_args()
    
    # Validate inputs
    if len(args.files) < 1:
        print("Error: At least one input file is required")
        return 1
    
    for file_path in args.files:
        if not Path(file_path).exists():
            print(f"Error: File '{file_path}' not found")
            return 1
    
    # Create output directory
    output_dir = Path(args.output).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate labels if not provided
    if args.labels:
        if len(args.labels) != len(args.files):
            print("Error: Number of labels must match number of files")
            return 1
        labels = args.labels
    else:
        labels = [Path(f).stem for f in args.files]
    
    # Load and process datasets
    print(f"Processing {len(args.files)} datasets...")
    datasets = []
    for file_path, label in zip(args.files, labels):
        print(f"Loading: {file_path} (as '{label}')")
        freq = compute_pul_size_freq(file_path)
        datasets.append(freq)
    
    # Create comparison plot
    create_comparison_plot(
        datasets=datasets,
        labels=labels,
        output_file=args.output,
        title=args.title,
        figsize=tuple(args.figsize),
        colors=args.colors,
        show_plot=args.show
    )
    
    # Print summary if requested
    if args.summary:
        print_summary(datasets, labels)
    
    print("Analysis complete!")
    return 0


if __name__ == "__main__":
    exit(main())