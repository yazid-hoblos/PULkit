import os
import matplotlib.pyplot as plt
import numpy as np
from statistics import mean, median

base_path = "/home/yaz/remote/report/easy-dbcan-results"
species_names = []
total_systems = []
mean_systems = []
median_systems = []
genome_counts_per_species = []

for species_dir in sorted(os.listdir(base_path)):
    if not species_dir.startswith("s__Prevotella_"):
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

    species_label = species_dir.replace("s__", "").replace("Prevotella_", "P. ")
    species_names.append(species_label)
    total_systems.append(sum(genome_counts))
    mean_systems.append(mean(genome_counts))
    median_systems.append(median(genome_counts))
    genome_counts_per_species.append(len(genome_counts))  # Store the number of genomes

# Create labels with genome counts
species_labels_with_counts = [
    f"{name} ({genomes})" for name, genomes in zip(species_names, genome_counts_per_species)
]

# Plotting grouped bars
x = np.arange(len(species_names))
width = 0.25

plt.figure(figsize=(16, 7))
plt.bar(x, total_systems, width, label='Total Systems')

plt.xticks(x, species_labels_with_counts, rotation=45, ha='right')
plt.xlabel("Species (number of genomes)")
plt.ylabel("Number of Systems")
plt.title("Total Number of Systems per Prevotella Species")
plt.legend()
plt.tight_layout()
plt.grid(axis='y', linestyle='--', alpha=0.5)
plt.savefig('prevotella_systems_counts_total.png', dpi=500)
# plt.show()
