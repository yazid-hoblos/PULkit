import matplotlib.pyplot as plt

# Input data: (directory, genome count)
data = [
    ("easy-dbcan-results/s__Prevotella_bryantii", 16),
    ("easy-dbcan-results/s__Prevotella_caecicola", 33),
    ("easy-dbcan-results/s__Prevotella_copri_B", 63),
    ("easy-dbcan-results/s__Prevotella_copri_E", 17),
    ("easy-dbcan-results/s__Prevotella_denticola", 20),
    ("easy-dbcan-results/s__Prevotella_histicola", 20),
    ("easy-dbcan-results/s__Prevotella_hominis", 15),
    ("easy-dbcan-results/s__Prevotella_intermedia", 35),
    ("easy-dbcan-results/s__Prevotella_melaninogenica", 34),
    ("easy-dbcan-results/s__Prevotella_merdae", 29),
    ("easy-dbcan-results/s__Prevotella_nigrescens", 16),
    ("easy-dbcan-results/s__Prevotella_pallens", 35),
    ("easy-dbcan-results/s__Prevotella_pectinovora", 22),
    ("easy-dbcan-results/s__Prevotella_rodentium", 291),
    ("easy-dbcan-results/s__Prevotella_sp000436035", 15),
    ("easy-dbcan-results/s__Prevotella_sp002409785", 19),
    ("easy-dbcan-results/s__Prevotella_sp900313215", 55),
    ("easy-dbcan-results/s__Prevotella_sp900545525", 151),
    ("easy-dbcan-results/s__Prevotella_stercorea", 42)
]

# Extract cleaned species names and genome counts
species_dirs, genome_counts = zip(*data)
species_names = [d.split('/')[-1].replace("s__Prevotella_", "") for d in species_dirs]
species_labels = [f"{name} ({count})" for name, count in zip(species_names, genome_counts)]

# First method total runtime (in minutes)
total_runtime_min_method1 = 27 * 60 + 14   # 27 min 14 sec = 27.2333 min
total_genomes = sum(genome_counts)
minutes_per_genome_method1 = total_runtime_min_method1 / total_genomes
estimated_runtimes_method1 = [count * minutes_per_genome_method1 for count in genome_counts]

# Second method: PANORAM total runtime (in seconds)
total_runtime_sec_method2 = 3519 / 60 + 20  # 3519 seconds = 58.65 minutes
seconds_per_genome_method2 = total_runtime_sec_method2 / total_genomes
estimated_runtimes_method2 = [count * seconds_per_genome_method2 for count in genome_counts]

# Plot
plt.figure(figsize=(12, 6))
plt.plot(species_labels, estimated_runtimes_method1, marker='o', linestyle='-', color='royalblue', label='dbCAN (27h14m)')
plt.plot(species_labels, estimated_runtimes_method2, marker='s', linestyle='--', color='orange', label='PANORAMA (1h18m)')

plt.xticks(rotation=45, ha='right', fontsize=8)
plt.xlabel("Species (Genome Count)")
plt.ylabel("Runtime")
plt.title("Runtime per Species")
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend()
plt.tight_layout()
plt.savefig("latex/figures/report/runtime_per_species.png", dpi=500)
plt.show()
