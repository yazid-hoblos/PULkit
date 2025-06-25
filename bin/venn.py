#!/usr/bin/env python3
"""
Generate a Venn diagram comparing intersections between three sets of gene families.

Usage:
    python venn.py [-h] set1 set2 set3 [-o OUTPUT] [-l LABELS]

Arguments:
    set1        Path to first input file containing gene families (one per line)
    set2        Path to second input file containing gene families
    set3        Path to third input file containing gene families
    -o/--output Output path for the Venn diagram image (default: venn_diagram.png)
    -l/--labels Labels for the three sets (comma-separated)

Example:
    python venn.py cazy_families.txt pfam_domains.txt egg_nog.txt -o my_venn.png -l "CAZy,Pfam,EggNOG"
"""

import argparse
from matplotlib import pyplot as plt
from matplotlib_venn import venn3

def load_set(filename):
    """Load a set of gene families from a file, one entry per line."""
    with open(filename) as f:
        return set(line.strip() for line in f if line.strip())

def create_venn_diagram(files, labels=None, output="venn_diagram.png"):
    """
    Create a Venn diagram from three input files of gene families.
    
    Args:
        files: List of three input file paths
        labels: List of three labels for the sets
        output: Output file path for the diagram
    """
    # Load sets from files
    sets = [load_set(f) for f in files]
    
    # Compute regions for venn diagram
    only_1 = sets[0] - sets[1] - sets[2]
    only_2 = sets[1] - sets[0] - sets[2]
    only_3 = sets[2] - sets[0] - sets[1]

    one_two = (sets[0] & sets[1]) - sets[2]
    one_three = (sets[0] & sets[2]) - sets[1]
    two_three = (sets[1] & sets[2]) - sets[0]

    all_sets = sets[0] & sets[1] & sets[2]

    # Create and save Venn diagram
    venn3(subsets=(
        len(only_1),
        len(only_2),
        len(one_two),
        len(only_3),
        len(one_three),
        len(two_three),
        len(all_sets)
    ), set_labels=labels)

    plt.title("Gene Family Set Intersections")
    plt.savefig(output, dpi=500, bbox_inches='tight')
    plt.close()

"""Parse arguments and create the Venn diagram."""
parser = argparse.ArgumentParser(description="Generate a Venn diagram comparing three sets of gene families.")
parser.add_argument("set1", help="First input file containing gene families")
parser.add_argument("set2", help="Second input file containing gene families")
parser.add_argument("set3", help="Third input file containing gene families")
parser.add_argument("-o", "--output", default="venn_diagram.png",
                    help="Output path for the Venn diagram (default: venn_diagram.png)")
parser.add_argument("-l", "--labels", 
                    help="Labels for the three sets (comma-separated)")

args = parser.parse_args()

labels = args.labels.split(",") if args.labels else None
create_venn_diagram([args.set1, args.set2, args.set3], labels, args.output)