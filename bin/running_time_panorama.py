import matplotlib.pyplot as plt
import numpy as np

# Data and species labels from before
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

species_dirs, genome_counts = zip(*data)
species_names = [d.split('/')[-1].replace("s__Prevotella_", "") for d in species_dirs]
species_labels = [f"{name} ({count})" for name, count in zip(species_names, genome_counts)]
total_genomes = sum(genome_counts)

# PANORAM times (minutes)
systems_detection = 3
annotation = 15
systems_projection = 58
total_panoram_time = systems_detection + annotation + systems_projection  # 76 minutes total

# Calculate minutes per genome for each phase (linear scaling assumption)
min_per_genome_detection = systems_detection / total_genomes
min_per_genome_annotation = annotation / total_genomes
min_per_genome_projection = systems_projection / total_genomes

# Estimated runtimes per species for each phase
detection_times = [count * min_per_genome_detection for count in genome_counts]
annotation_times = [count * min_per_genome_annotation for count in genome_counts]
projection_times = [count * min_per_genome_projection for count in genome_counts]

# Plot stacked bar chart
fig, ax = plt.subplots(figsize=(14, 7))

bar_width = 0.7
x = np.arange(len(species_labels))

p1 = ax.bar(x, detection_times, bar_width, label='Systems Detection (3 min total)', color='skyblue')
p2 = ax.bar(x, annotation_times, bar_width, bottom=detection_times, label='Annotation (15 min total)', color='lightgreen')

# Bottom for projection bars = detection + annotation times
bottom_projection = [d + a for d, a in zip(detection_times, annotation_times)]
p3 = ax.bar(x, projection_times, bar_width, bottom=bottom_projection, label='Systems Projection (58 min total)', color='salmon')

ax.set_xticks(x)
ax.set_xticklabels(species_labels, rotation=45, ha='right', fontsize=8)
ax.set_ylabel('Runtime (minutes)')
ax.set_xlabel('Species (Genome Count)')
ax.set_title('PANORAM Runtime per Species (Split by Phase)')
ax.legend()
ax.grid(True, linestyle='--', alpha=0.5)

plt.tight_layout()
plt.savefig("report/panoram_runtime_per_species.png", dpi=500)
plt.show()
