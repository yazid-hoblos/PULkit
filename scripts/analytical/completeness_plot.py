#!/usr/bin/env python3
"""
Generate completeness plot for PUL systems.

Usage:
    python completeness.py [-h] --systems FILE [-o OUTPUT]
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import argparse
from pathlib import Path

def plot_completeness(systems_file, output_dir="completeness_plots"):
    """Generate completeness plot from systems.tsv file."""
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Load data
    print(f"Loading systems from: {systems_file}")
    df = pd.read_csv(systems_file, sep='\t')
    
    # Ensure required columns exist
    required_cols = ['system number', 'completeness']
    if not all(col in df.columns for col in required_cols):
        raise ValueError(f"Missing required columns. File must contain: {required_cols}")
    
    # Convert system number to int and sort
    df['system number'] = df['system number'].astype(int)
    df = df.sort_values('system number')
    
    # Extract data for plotting
    systems = df['system number']
    completeness = df['completeness']
    
    # Create plot
    plt.figure(figsize=(12, 6))
    plt.plot(systems, completeness, marker='o', linestyle='-', color='teal', 
             markersize=6, linewidth=2, alpha=0.7)
    
    # Customize plot
    plt.xlabel("System Number", fontsize=12)
    plt.ylabel("Completeness", fontsize=12)
    genome_name = Path(systems_file).parent.parent.name.replace('Prevotella_', 'P. ')
    plt.title(f"PUL System Completeness - {genome_name}", fontsize=14)
    
    # Add grid and style
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tick_params(axis='both', which='major', labelsize=10)
    
    # Add mean line
    mean_completeness = completeness.mean()
    plt.axhline(y=mean_completeness, color='red', linestyle='--', alpha=0.5,
                label=f'Mean: {mean_completeness:.2f}')
    
    # Add statistics annotation
    stats_text = (f"Total Systems: {len(systems)}\n"
                 f"Mean: {mean_completeness:.2f}\n"
                 f"Median: {completeness.median():.2f}\n"
                 f"Max: {completeness.max():.2f}\n"
                 f"Min: {completeness.min():.2f}")
    plt.text(0.02, 0.98, stats_text, transform=plt.gca().transAxes, 
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.legend()
    plt.tight_layout()
    
    # Save plot
    output_path = os.path.join(output_dir, f"completeness_{genome_name}.png")
    # plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved plot to: {output_path}")
    plt.show()
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Generate completeness plot for PUL systems")
    parser.add_argument("--systems", required=True, help="Path to systems.tsv file")
    parser.add_argument("-o", "--output", default="completeness_plots", 
                       help="Output directory for plots")
    
    args = parser.parse_args()
    plot_completeness(args.systems, args.output)

if __name__ == "__main__":
    main()