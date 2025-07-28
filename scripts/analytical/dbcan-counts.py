#!/usr/bin/env python3
"""
dbCAN System Count Analysis

This script analyzes dbCAN CGC (Carbohydrate-active enzyme Gene Cluster) results
to count and visualize the number of systems per species and genome.
It processes directories containing CGC files and generates bar plots.
"""

import os
import argparse
import matplotlib.pyplot as plt
import numpy as np
from statistics import mean, median
from pathlib import Path


def load_cgc_data(base_path, species_filter=None):
    """
    Load CGC data from directory structure.
    
    Args:
        base_path: Base directory containing species subdirectories
        species_filter: Optional filter for species names (e.g., "s__Prevotella_")
    
    Returns:
        Tuple of lists: (species_names, total_systems, mean_systems, median_systems, genome_counts)
    """
    species_names = []
    total_systems = []
    mean_systems = []
    median_systems = []
    genome_counts_per_species = []
    
    for species_dir in sorted(os.listdir(base_path)):
        if species_filter and not species_dir.startswith(species_filter):
            continue
        
        species_path = os.path.join(base_path, species_dir)
        if not os.path.isdir(species_path):
            continue
        
        genome_counts = []
        for genome_dir in os.listdir(species_path):
            genome_path = os.path.join(species_path, genome_dir)
            if not os.path.isdir(genome_path):
                continue
            
            try:
                cgc_file = next(f for f in os.listdir(genome_path) if f.startswith("cgc_standard"))
                cgc_path = os.path.join(genome_path, cgc_file)
                
                with open(cgc_path) as f:
                    system_ids = set(line.split('\t')[0].strip() for line in f if line.strip())
                    genome_counts.append(len(system_ids))
            except StopIteration:
                continue
        
        if not genome_counts:
            continue  # skip species with no valid genomes
        
        # Clean up species names
        species_label = species_dir.replace("s__", "").replace("Prevotella_", "P. ")
        species_names.append(species_label)
        total_systems.append(sum(genome_counts))
        mean_systems.append(mean(genome_counts))
        median_systems.append(median(genome_counts))
        genome_counts_per_species.append(len(genome_counts))
    
    print(f"Processed {len(species_names)} species")
    return species_names, total_systems, mean_systems, median_systems, genome_counts_per_species


def create_bar_plot(species_names, values, genome_counts, title, ylabel, output_file, 
                    figsize=(16, 7), show_genome_counts=True):
    """Create and save a bar plot."""
    # Create labels with genome counts if requested
    if show_genome_counts:
        labels = [f"{name} ({count})" for name, count in zip(species_names, genome_counts)]
        xlabel = "Species (number of genomes)"
    else:
        labels = species_names
        xlabel = "Species"
    
    x = np.arange(len(species_names))
    
    plt.figure(figsize=figsize)
    plt.bar(x, values, color='steelblue', alpha=0.8)
    
    plt.xticks(x, labels, rotation=45, ha='right')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(axis='y', linestyle='--', alpha=0.5)
    plt.tight_layout()
    
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Plot saved: {output_file}")


def create_comparison_plot(species_names, total_systems, mean_systems, median_systems, 
                          genome_counts, output_file, figsize=(18, 8)):
    """Create a comparison plot with total, mean, and median values."""
    labels = [f"{name} ({count})" for name, count in zip(species_names, genome_counts)]
    x = np.arange(len(species_names))
    width = 0.25
    
    plt.figure(figsize=figsize)
    
    plt.bar(x - width, total_systems, width, label='Total Systems', alpha=0.8)
    plt.bar(x, mean_systems, width, label='Mean Systems per Genome', alpha=0.8)
    plt.bar(x + width, median_systems, width, label='Median Systems per Genome', alpha=0.8)
    
    plt.xticks(x, labels, rotation=45, ha='right')
    plt.xlabel("Species (number of genomes)")
    plt.ylabel("Number of Systems")
    plt.title("System Counts Comparison per Species")
    plt.legend()
    plt.grid(axis='y', linestyle='--', alpha=0.5)
    plt.tight_layout()
    
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Comparison plot saved: {output_file}")


def print_summary(species_names, total_systems, mean_systems, median_systems, genome_counts):
    """Print summary statistics."""
    print("\n=== SUMMARY STATISTICS ===")
    print(f"Total species analyzed: {len(species_names)}")
    print(f"Total genomes: {sum(genome_counts)}")
    print(f"Total systems across all species: {sum(total_systems)}")
    print(f"Average systems per species: {mean(total_systems):.1f}")
    print(f"Average systems per genome: {mean(mean_systems):.1f}")
    
    print("\n=== TOP 5 SPECIES BY TOTAL SYSTEMS ===")
    sorted_indices = sorted(range(len(total_systems)), key=lambda i: total_systems[i], reverse=True)
    for i, idx in enumerate(sorted_indices[:5], 1):
        print(f"{i}. {species_names[idx]}: {total_systems[idx]} total ({genome_counts[idx]} genomes)")


def main():
    """Main function to run the dbCAN count analysis."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument("base_path",
                       help="Base directory containing species subdirectories with CGC results")
    parser.add_argument("-f", "--filter", default="s__Prevotella_",
                       help="Species filter prefix (default: s__Prevotella_)")
    parser.add_argument("-o", "--output-dir", default="report",
                       help="Output directory for plots (default: report)")
    parser.add_argument("--summary", action="store_true",
                       help="Print summary statistics")
    
    args = parser.parse_args()
    
    # Validate input
    if not os.path.exists(args.base_path):
        print(f"Error: Base path '{args.base_path}' not found")
        return 1
    
    # Create output directory
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    
    # Load data
    print(f"Loading CGC data from: {args.base_path}")
    print(f"Filtering species with prefix: {args.filter}")
    
    species_names, total_systems, mean_systems, median_systems, genome_counts = load_cgc_data(
        args.base_path, args.filter)
    
    if not species_names:
        print("No species found matching the filter criteria")
        return 1
    
    # Create plots
    figsize = tuple(args.figsize)
    show_counts = not args.no_genome_counts
    
    # Total systems plot
    total_output = os.path.join(args.output_dir, "systems_counts_total.png")
    create_bar_plot(species_names, total_systems, genome_counts,
                   f"Total Number of Systems per Species (filter: {args.filter})",
                   "Total Systems", total_output, figsize, show_counts)
    
    # Mean systems plot
    mean_output = os.path.join(args.output_dir, "systems_counts_mean.png")
    create_bar_plot(species_names, mean_systems, genome_counts,
                   f"Mean Systems per Genome by Species (filter: {args.filter})",
                   "Mean Systems per Genome", mean_output, figsize, show_counts)
    
    # Median systems plot
    median_output = os.path.join(args.output_dir, "systems_counts_median.png")
    create_bar_plot(species_names, median_systems, genome_counts,
                   f"Median Systems per Genome by Species (filter: {args.filter})",
                   "Median Systems per Genome", median_output, figsize, show_counts)
    
    # Comparison plot
    if not args.skip_comparison:
        comparison_output = os.path.join(args.output_dir, "systems_counts_comparison.png")
        create_comparison_plot(species_names, total_systems, mean_systems, median_systems,
                             genome_counts, comparison_output, (18, 8))
    
    # Print summary
    if args.summary:
        print_summary(species_names, total_systems, mean_systems, median_systems, genome_counts)
    
    print("Analysis complete!")
    return 0


if __name__ == "__main__":
    exit(main())