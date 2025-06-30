import matplotlib.pyplot as plt

# Species names (extracted and cleaned)
species = [
    "P. bryantii", "P. caecicola", "P. copri B", "P. copri E",
    "P. denticola", "P. histicola", "P. hominis", "P. intermedia",
    "P. melaninogenica", "P. merdae", "P. nigrescens", "P. pallens",
    "P. pectinovora", "P. rodentium", "P. sp000436035", "P. sp002409785",
    "P. sp900313215", "P. sp900545525", "P. stercorea", 'B. thetaiotaomicron'
]

# Genome counts from `ls $f | wc -l`
genome_counts = [
    16, 33, 63, 17, 20, 20, 15, 35, 34, 29, 16, 35, 22, 291, 15, 19, 55, 151, 42, 327
]

# Systems counts from `wc -l systems.tsv`
system_counts = [
    41, 16, 48, 41, 16, 15, 60, 44, 43, 53, 20, 10, 37, 33, 81, 38, 119, 33, 36, 540
]

# Sort data by species for consistent plotting
# species, genome_counts, system_counts = zip(*sorted(zip(species, genome_counts, system_counts)))
# sort by genome_counts
species, genome_counts, system_counts = zip(*sorted(zip(species, genome_counts, system_counts), key=lambda x: x[1]))

# Plotting
x = range(len(species))
width = 0.4

plt.figure(figsize=(14, 7))
plt.bar(x, genome_counts, width=width, label='Genomes', align='center')
plt.bar([i + width for i in x], system_counts, width=width, label='Systems', align='center')

plt.xlabel('Species')
plt.ylabel('Count')
plt.title('Number of Genomes and Systems per Prevotella Species')
plt.xticks([i + width/2 for i in x], species, rotation=70)
plt.legend()
plt.tight_layout()
plt.grid(axis='y', linestyle='--', alpha=0.5)
# plt.savefig('prevotella_genomes_systems_counts.png', dpi=500)
plt.show()
