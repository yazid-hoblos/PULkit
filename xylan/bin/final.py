import os, sys
from collections import defaultdict
import re

def infer_family(annotation):
    # annotation = annotation.lower()
    if "hypothetical" in annotation:
        return "Hypothetical"
    # if "Sus" in annotation or "Ton" in annotation:
    #     return "SusD"
    return "Other"

def load_protein_categories(category_files):
    protein_to_category = defaultdict(list)
    for category, filename in category_files.items():
        filepath = filename
        if not os.path.isfile(filepath):
            print(f"Warning: category file not found: {filepath}")
            continue
        with open(filepath, 'r') as f:
            for line in f:
                line=line.strip()
                if not line or '|' not in line:
                    continue
                _, protein_id = line.split('|', 1)
                protein_to_category[protein_id.strip()].append(category)
    reduced_protein_to_category = {}
    for protein_id, categories in protein_to_category.items():
        if len(categories) == 1:
            reduced_protein_to_category[protein_id] = next(iter(categories))
        elif 'CAZyme' in categories:
            reduced_protein_to_category[protein_id] = 'CAZyme'
        elif 'Transporter' in categories:
            reduced_protein_to_category[protein_id] = 'Transporter'
        elif 'STP' in categories:
            reduced_protein_to_category[protein_id] = 'STP'
        elif 'TF' in categories:
            reduced_protein_to_category[protein_id] = 'TF'
        elif 'Peptidase' in categories:
            reduced_protein_to_category[protein_id] = 'Peptidase'
        elif 'Sulfatase' in categories:
            reduced_protein_to_category[protein_id] = 'Sulfatase'
        else:
            reduced_protein_to_category[protein_id] = next(iter(categories))
        # if 'TF' in categories:
            # reduced_protein_to_category[protein_id] = 'TF'
    return reduced_protein_to_category

def load_protein_annotations(annotation_file):
    protein_annotation_map = {}
    with open(annotation_file, 'r') as f:
        for line in f:
            parts = line.strip().split(',')
            if len(parts) < 11:
                continue
            if parts[2].strip() != "-":
                protein_id = parts[2].strip()
            elif parts[1].strip() != "-":
                protein_id = parts[1].strip()
            else:
                protein_id = parts[3].strip()
            annotation = parts[9].strip()
            protein_annotation_map[protein_id] = annotation
    return protein_annotation_map

def group_proteins_by_pul(annotation_file):
    pul_proteins = defaultdict(list)
    with open(annotation_file, 'r') as f:
        for line in f:
            parts = line.strip().split(',')
            if len(parts) < 11:
                continue
            pul_id = parts[0].strip()
            if parts[2].strip() != "-":
                protein_id = parts[2].strip()
            elif parts[1].strip() != "-":
                protein_id = parts[1].strip()
            else:
                protein_id = parts[3].strip()
            pul_proteins[pul_id].append(protein_id)
    return pul_proteins

if __name__ == "__main__":
    category_files = {
        "CAZyme": "cazymes",
        "TF": "tfs",
        "STP": "stps",
        "Peptidase": "peptidases",
        "Transporter": "tcs",
        "Sulfatase": "sulfatases"
    }
    annotation_file = sys.argv[1] if len(sys.argv) > 1 else "xylan_pul_data.csv"

    protein_to_category = load_protein_categories(category_files)
    protein_annotations = load_protein_annotations(annotation_file)
    pul_proteins = group_proteins_by_pul(annotation_file)

    enriched_pul_data = defaultdict(list)
    for pul_id, proteins in pul_proteins.items():
        for pid in proteins:
            if pid not in protein_to_category:
                annotation = protein_annotations[pid]
                category = infer_family(annotation)
            else:
                category = protein_to_category[pid]

            enriched_pul_data[pul_id].append({
                "protein_id": pid,
                "category": category
            })

    # Example printout
    for pul_id, components in enriched_pul_data.items():
        print(f"{pul_id}:")
        for comp in components:
            print(f"  {comp['protein_id']} | {comp['category']}")

exit()
## PATTERNS

from collections import Counter

# Store PUL category sequences
pul_category_sequences = {}

# Sort PULs by number of components (descending)
sorted_puls = sorted(enriched_pul_data.items(), key=lambda x: -len(x[1]))
pul_category_sequences = {}
for pul_id, components in sorted_puls:
    categories = [comp["category"] for comp in components]
    pul_category_sequences[pul_id] = categories


# Print individual PUL category sequences
print("\nCategory sequences per PUL:")
for pul_id, seq in pul_category_sequences.items():
    print(f"{pul_id}: {' → '.join(seq)}")

# Count common n-grams across all PULs
def get_ngrams(seq, n):
    return [tuple(seq[i:i+n]) for i in range(len(seq)-n+1)]

ngram_counter = Counter()

for seq in pul_category_sequences.values():
    for n in range(2, 5):  # 2-grams to 4-grams
        ngrams = get_ngrams(seq, n)
        ngram_counter.update(ngrams)

# Show top consensus patterns
print("\nTop consensus patterns (n-grams):")
for pattern, count in ngram_counter.most_common(100):
    print(f"{' → '.join(pattern)} ({count} occurrences)")


## Overall consensus pattern

# Step 1: Get all n-grams
def get_ngrams(seq, n):
    return [tuple(seq[i:i+n]) for i in range(len(seq)-n+1)]

ngram_counter = Counter()
for seq in pul_category_sequences.values():
    for n in range(2, 6):  # 2-grams to 5-grams
        ngram_counter.update(get_ngrams(seq, n))

# Step 2: Sort from longest and most frequent to shortest
sorted_ngrams = sorted(ngram_counter.items(), key=lambda x: (-len(x[0]), -x[1]))

# Step 3: Greedily build one consensus sequence
consensus = []
for pattern, _ in sorted_ngrams:
    if not consensus:
        consensus.extend(pattern)
    else:
        # Try to extend consensus if pattern overlaps
        for i in range(1, len(pattern)):
            if consensus[-i:] == list(pattern[:i]):
                consensus.extend(pattern[i:])
                break
        else:
            continue  # skip if no overlap

# Remove duplicates while preserving order
seen = set()
final_consensus = []
for cat in consensus:
    if cat not in seen:
        final_consensus.append(cat)
        seen.add(cat)

# Output final consensus sequence
print("\nConsensus category sequence:")
print(" → ".join(final_consensus))


## PLOT

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from collections import defaultdict

# Assume pul_category_sequences is already sorted and available

# Assign each category a unique color
categories = sorted({cat for seq in pul_category_sequences.values() for cat in seq})
color_map = {cat: plt.cm.tab20(i) for i, cat in enumerate(categories)}

# Plot
fig, ax = plt.subplots(figsize=(12, len(pul_category_sequences) * 0.5))

for i, (pul_id, seq) in enumerate(pul_category_sequences.items()):
    for j, category in enumerate(seq):
        ax.barh(i, 1, left=j, color=color_map[category], edgecolor='black')

# Customize ticks
ax.set_yticks(range(len(pul_category_sequences)))
ax.set_yticklabels(list(pul_category_sequences.keys()))
ax.set_xlabel("Gene Position")
ax.set_title("PUL Category Sequences (Sorted by Length)")

# Legend
handles = [mpatches.Patch(color=color_map[cat], label=cat) for cat in categories]
plt.legend(handles=handles, bbox_to_anchor=(1.05, 1), loc='upper left', title="Category")

plt.tight_layout()
plt.show()

import matplotlib.cm as cm
custom_order = ["CAZyme", "Transporter", "SusD", "STP", "TF", "Sulfatase", "Peptidase", "Hypothetical", "Other"]

# Filter actual categories used in data (preserve order from custom_order)
used_categories = [cat for cat in custom_order if any(cat in seq for seq in pul_category_sequences.values())]

# Assign consistent colors
cmap = cm.get_cmap("Set3", len(used_categories))
color_map = {
    "CAZyme": "#1f77b4",       # blue
    "Transporter": "#ff7f0e",  # orange
    "SusD": "#2ca02c",         # green
    "STP": "#d62728",          # red
    "TF": "#9467bd",           # purple
    "Sulfatase": "#8c564b",    # brown
    "Peptidase": "#e377c2",    # pink
    "Hypothetical": "#7f7f7f", # gray
    "Other": "#17becf"         # cyan
}

# Sort PULs from largest to smallest by sequence length
sorted_puls = sorted(pul_category_sequences.items(), key=lambda x: len(x[1]), reverse=True)

# Plot
fig, ax = plt.subplots(figsize=(12, len(sorted_puls) * 0.5))

for i, (pul_id, seq) in enumerate(sorted_puls):
    # Sort each PUL's gene categories according to custom_order
    sorted_seq = sorted(seq, key=lambda x: custom_order.index(x) if x in custom_order else 999)
    for j, category in enumerate(sorted_seq):
        ax.barh(i, 1, left=j, color=color_map.get(category, "gray"), edgecolor='black')

# Customize ticks
ax.set_yticks(range(len(sorted_puls)))
ax.set_yticklabels([pul_id for pul_id, _ in sorted_puls])
ax.set_xlabel("Gene Position (Grouped by Custom Category Order)")
ax.set_title("PUL Category Sequences (Grouped by Custom Order, Largest to Smallest)")
ax.invert_yaxis()

# Legend
handles = [mpatches.Patch(color=color_map[cat], label=cat) for cat in used_categories]
plt.legend(handles=handles, title="Category", bbox_to_anchor=(1.05, 1), loc='upper left')

plt.tight_layout()
plt.show()


## Disregard order - convert to sets

from collections import defaultdict, Counter

# Assume enriched_pul_data is already populated
# Structure: { PUL_ID: [ { protein_id, category, annotation }, ... ] }

# Step 1: Convert each PUL to a set of categories
pul_category_sets = {}
for pul_id, components in enriched_pul_data.items():
    category_set = frozenset(comp["category"] for comp in components)
    pul_category_sets[pul_id] = category_set

# Step 2: Count how often each category set appears
category_set_counter = Counter(pul_category_sets.values())

# Step 3: Print the most common sets
print("\nMost common category sets across PULs:")
for i, (cat_set, count) in enumerate(category_set_counter.most_common(), 1):
    categories_sorted = sorted(cat_set)
    print(f"{i}. {categories_sorted} -> {count} PUL(s)")

# Optional: Print which PULs belong to each category set
print("\nCategory set to PUL mapping:")
for cat_set, count in category_set_counter.most_common():
    matching_puls = [pul_id for pul_id, s in pul_category_sets.items() if s == cat_set]
    print(f"{sorted(cat_set)}: {matching_puls}")
    
    
## New PLOT

# import matplotlib.pyplot as plt
# import matplotlib.patches as mpatches
# from collections import defaultdict

# pul_ids = list(pul_category_sequences.keys())
# max_len = max(len(seq) for seq in pul_category_sequences.values())

# # Get all unique categories
# all_categories = sorted(set(cat for seq in pul_category_sequences.values() for cat in seq))

# cmap = plt.cm.get_cmap("tab20", len(all_categories))
# color_map = {cat: cmap(i) for i, cat in enumerate(all_categories)}

# fig, ax = plt.subplots(figsize=(max_len * 1.2, len(pul_ids) * 0.8))

# # For each gene position, group PULs by category at that position
# for gene_pos in range(max_len):
#     # Group puls by category at this position
#     category_to_puls = defaultdict(list)
#     for pul_id in pul_ids:
#         seq = pul_category_sequences[pul_id]
#         if gene_pos < len(seq):
#             category_to_puls[seq[gene_pos]].append(pul_id)
#         else:
#             # no gene at this position for this PUL, treat separately if needed
#             pass

#     # Assign vertical positions within this column by stacking categories
#     y_offset = 0
#     for cat in all_categories:
#         pul_list = category_to_puls.get(cat, [])
#         for i, pul_id in enumerate(pul_list):
#             y = y_offset + i
#             ax.scatter(gene_pos, y, s=600, color=color_map[cat], edgecolors='k', marker='s')
#         y_offset += len(pul_list)

# # Customize axes
# ax.set_xlim(-0.5, max_len - 0.5)
# ax.set_xlabel("Gene Position")
# ax.set_title("PUL Categories Grouped by Category at Each Gene Position")

# # Remove y-axis labels since vertical position is per column grouping, not PULs directly
# ax.set_yticks([])

# # X ticks
# ax.set_xticks(range(max_len))

# # Legend
# patches = [mpatches.Patch(color=color_map[cat], label=cat) for cat in all_categories]
# ax.legend(handles=patches, title="Category", bbox_to_anchor=(1.05, 1), loc='upper left')

# plt.tight_layout()
# plt.show()




### SAVE

# import matplotlib.pyplot as plt
# import matplotlib.patches as mpatches
# from collections import defaultdict

# pul_ids = list(pul_category_sequences.keys())
# max_len = max(len(seq) for seq in pul_category_sequences.values())
# all_categories = sorted(set(cat for seq in pul_category_sequences.values() for cat in seq))

# cmap = plt.cm.get_cmap("tab20", len(all_categories))
# color_map = {cat: cmap(i) for i, cat in enumerate(all_categories)}

# fig_width = min(max_len * 1.2, 36)  # cap width to 30 inches max for visibility
# fig_height = max(len(pul_ids) * 0.7, 12)  # minimum height for clarity
# fig, ax = plt.subplots(figsize=(fig_width, fig_height))

# marker_size = 400  # larger squares for visibility

# for gene_pos in range(max_len):
#     category_to_puls = defaultdict(list)
#     for pul_id in pul_ids:
#         seq = pul_category_sequences[pul_id]
#         if gene_pos < len(seq):
#             category_to_puls[seq[gene_pos]].append(pul_id)

#     y_offset = 0
#     for cat in all_categories:
#         pul_list = category_to_puls.get(cat, [])
#         for i, pul_id in enumerate(pul_list):
#             y = y_offset + i
#             ax.scatter(gene_pos, y, s=marker_size, color=color_map[cat], edgecolors='k', marker='s')
#         y_offset += len(pul_list)

# ax.set_xlim(-0.5, max_len - 0.5)
# ax.set_xlabel("Gene Position")
# ax.set_title("PUL Categories Grouped by Category at Each Gene Position")

# # Optionally add horizontal grid lines every ~10 units for readability
# ax.yaxis.grid(True, linestyle='--', alpha=0.3)

# ax.set_yticks([])  # no y labels since vertical positions shift per gene pos

# ax.set_xticks(range(max_len))

# patches = [mpatches.Patch(color=color_map[cat], label=cat) for cat in all_categories]
# ax.legend(handles=patches, title="Category", bbox_to_anchor=(1.05, 1), loc='upper left')

# plt.tight_layout()
# plt.savefig("puls_by_category.png", dpi=100)
# plt.show()
