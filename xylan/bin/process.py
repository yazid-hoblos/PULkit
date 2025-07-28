import os
import sys
from collections import defaultdict


def infer_family(annotation):
    if "hypothetical" in annotation:
        return "Hypothetical"
    if "SusD" in annotation or "RagB" in annotation:
        return "SusD"
    return "Other"


def load_protein_categories_from_directory(data_dir, suscd_file=None):
    file_mapping = {
        "CAZyme": os.path.join(data_dir, "overview.tsv"),
        "TF": os.path.join(data_dir, "TF_hmm_results.tsv"),
        "STP": os.path.join(data_dir, "STP_hmm_results.tsv"),
        "Peptidase": os.path.join(data_dir, "diamond.out.peptidase"),
        "Transporter": os.path.join(data_dir, "diamond.out.tc"),
        "Sulfatase": os.path.join(data_dir, "diamond.out.sulfatlas")
    }

    protein_to_category = defaultdict(list)

    # CAZyme from overview.tsv (column 1)
    if os.path.isfile(file_mapping["CAZyme"]):
        with open(file_mapping["CAZyme"]) as f:
            for line in f:
                if not line.strip():
                    continue
                fields = line.strip().split("\t")
                if fields:
                    # protein_to_category[fields[0]].append("CAZyme")
                    if fields[2].strip() != '-':
                        family = fields[2].strip().split("+")[0].split("_")[0].split("(")[0]
                    elif fields[3].strip() != '-':
                        family = fields[3].strip().split("+")[0].split("_")[0].split("(")[0]
                    else:
                        family = fields[4].strip().split("+")[0].split("_")[0].split("(")[0]
                    protein_to_category[fields[0]].append(family)

    # TFs and STPs from column 3
    for category in ["TF", "STP"]:
        if os.path.isfile(file_mapping[category]):
            with open(file_mapping[category]) as f:
                for line in f:
                    if not line.strip():
                        continue
                    fields = line.strip().split("\t")
                    if len(fields) >= 3:
                        if not category in protein_to_category[fields[2]]:
                            protein_to_category[fields[2]].append(category)

    # Peptidases, Transporters, Sulfatases from column 3
    for category in ["Peptidase", "Transporter", "Sulfatase"]:
        if os.path.isfile(file_mapping[category]):
            with open(file_mapping[category]) as f:
                for line in f:
                    if not line.strip():
                        continue
                    fields = line.strip().split("\t")
                    if len(fields) >= 3:
                        if not category in protein_to_category[fields[2]]:
                            protein_to_category[fields[2]].append(category)

    if suscd_file:
        suscd_categories = load_suscd_categories(suscd_file)
        for pid, cat in suscd_categories.items():
            protein_to_category[pid] = [cat]

    # Prioritize category assignment
    reduced_protein_to_category = {}
    for protein_id, categories in protein_to_category.items():
        reduced_protein_to_category[protein_id] = protein_to_category[protein_id][0] 
    #     # if 'TF' in categories:
    #     #     print(f"Warning: {protein_id} has multiple categories: {categories}")
    #     if 'CAZyme' in categories:
    #         reduced_protein_to_category[protein_id] = 'CAZyme'
    #     elif 'Transporter' in categories:
    #         reduced_protein_to_category[protein_id] = 'Transporter'
    #     elif 'STP' in categories and 'TF' in categories:
    #         reduced_protein_to_category[protein_id] = 'TF-STP'
    #     else:
    #         reduced_protein_to_category[protein_id] = categories[0]
    #     if len(categories) == 1:
    #         reduced_protein_to_category[protein_id] = categories[0]
    #     elif 'CAZyme' in categories:
    #         reduced_protein_to_category[protein_id] = 'CAZyme'
    #     elif 'Transporter' in categories:
    #         reduced_protein_to_category[protein_id] = 'Transporter'
    #     elif 'STP' in categories:
    #         reduced_protein_to_category[protein_id] = 'STP'
    #     elif 'TF' in categories:
    #         reduced_protein_to_category[protein_id] = 'TF'
    #     elif 'Peptidase' in categories:
    #         reduced_protein_to_category[protein_id] = 'Peptidase'
    #     elif 'Sulfatase' in categories:
    #         reduced_protein_to_category[protein_id] = 'Sulfatase'
    #     else:
    #         reduced_protein_to_category[protein_id] = categories[0]

    return reduced_protein_to_category


def load_suscd_categories(suscd_file):
    protein_to_suscd = defaultdict(list)
    if os.path.isfile(suscd_file):
        with open(suscd_file) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.strip().split()  # splits by ANY whitespace
                # if len(fields) >= 3:
                domain_name = fields[0]      # e.g., SusC or SusD HMM name
                protein_id = fields[2]       # target name (query protein)
                if "SusD" in domain_name:
                    protein_to_suscd[protein_id] = "SusD"
                else:
                    protein_to_suscd[protein_id] = "SusC"
    return protein_to_suscd


def load_protein_annotations(annotation_file):
    protein_annotation_map = {}
    with open(annotation_file, 'r') as f:
        for line in f:
            parts = line.strip().split(',')
            if len(parts) < 11:
                continue
            protein_id = parts[0].strip() + '|' + parts[3].strip()
            annotation = parts[9].strip()
            protein_annotation_map[protein_id] = annotation
    return protein_annotation_map


def group_proteins_by_pul(annotation_file):
    pul_proteins = defaultdict(list)
    with open(annotation_file, 'r') as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split(',')
            if len(parts) < 11:
                continue
            pul_id = parts[0].strip()
            protein_id = pul_id + '|' + parts[3].strip()
            pul_proteins[pul_id].append(protein_id)
    return pul_proteins


data_dir = sys.argv[1] if len(sys.argv) > 1 else "cgc_easy_xylan_updated"
annotation_file = sys.argv[2] if len(sys.argv) > 2 else "xylan_pul_data.csv"

protein_to_category = load_protein_categories_from_directory(data_dir, suscd_file=sys.argv[3] if len(sys.argv) > 3 else None)
protein_annotations = load_protein_annotations(annotation_file)
pul_proteins = group_proteins_by_pul(annotation_file)

unannotated = set()

enriched_pul_data = defaultdict(list)
for pul_id, proteins in pul_proteins.items():
    for pid in proteins:
        if pid not in protein_to_category:
            unannotated.add(pid)
            annotation = protein_annotations[pid]
            category = infer_family(annotation)
        else:
            category = protein_to_category[pid]
        enriched_pul_data[pul_id].append({"protein_id": pid, "category": category})

# for pul_id, components in enriched_pul_data.items():
#     print(f"{pul_id}:")
#     for comp in components:
#         print(f"  {comp['protein_id'].split('|')[1]} | {comp['category']}")
    

from collections import Counter

pul_category_sequences = {}

# Sort PULs by number of components (descending)
sorted_puls = sorted(enriched_pul_data.items(), key=lambda x: -len(x[1]))
pul_category_sequences = {}
for pul_id, components in sorted_puls:
    categories = [comp["category"] for comp in components]
    pul_category_sequences[pul_id] = categories

# save pul_category_sequences to a file
output_file = "xylan_pul_category_sequences_with_susCD.txt"
with open(output_file, 'w') as f:
    for pul_id, seq in pul_category_sequences.items():
        seq = [cat for cat in seq if cat not in ["Other", "Hypothetical", 'TF', 'STP', 'Transporter', 'Peptidase', 'Sulfatase']]
        f.write(f"{pul_id}: {' → '.join(seq)}\n")

# Print individual PUL category sequences
print("\nCategory sequences per PUL:")

all_cazyme_families = []
for pul_id, seq in pul_category_sequences.items():
    if pul_id == 'PUL0415' or pul_id == 'PUL0345':
        print(f"{pul_id}: {' → '.join(seq)}")
    # remove Other and Hypothetical categories
    seq = [cat for cat in seq if cat not in ["Other", "Hypothetical", 'TF', 'STP', 'Transporter', 'Peptidase', 'Sulfatase', 'SusD', 'SusC']]
    all_cazyme_families.extend(set(seq))
    # print(f"{pul_id}: {' → '.join(seq)}")
print("Most common CAZyme families:")
cazyme_counter = Counter(all_cazyme_families)
for family, count in cazyme_counter.most_common(20):
    print(f"{family}: {count} occurrences")
print("Total unique CAZyme families:", len(cazyme_counter))
    
print(len(list(unannotated)), "unannotated proteins")
print(len(list(set(protein_to_category.keys()))), "annotated proteins")

# for protein in unannotated:
#     if infer_family(protein_annotations[protein]) == "Other":
#         print(f"{protein} -  {protein_annotations[protein]}")


# Consensus patterns (n-grams)
def get_ngrams(seq, n):
    return [tuple(seq[i:i+n]) for i in range(len(seq)-n+1)]

ngram_counter = Counter()

for seq in pul_category_sequences.values():
    for n in range(2, 5):  # 2-grams to 4-grams
        ngrams = get_ngrams(seq, n)
        ngram_counter.update(ngrams)

# Show top consensus patterns
# print("\nTop consensus patterns (n-grams):")
# for pattern, count in ngram_counter.most_common(100):
#     print(f"{' → '.join(pattern)} ({count} occurrences)")

# # Step 2: Sort from longest and most frequent to shortest
# sorted_ngrams = sorted(ngram_counter.items(), key=lambda x: (-len(x[0]), -x[1]))

# # Step 3: Greedily build one consensus sequence
# consensus = []
# for pattern, _ in sorted_ngrams:
#     if not consensus:
#         consensus.extend(pattern)
#     else:
#         # Try to extend consensus if pattern overlaps
#         for i in range(1, len(pattern)):
#             if consensus[-i:] == list(pattern[:i]):
#                 consensus.extend(pattern[i:])
#                 break
#         else:
#             continue  # skip if no overlap

# # Remove duplicates while preserving order
# seen = set()
# final_consensus = []
# for cat in consensus:
#     if cat not in seen:
#         final_consensus.append(cat)
#         seen.add(cat)

# # Output final consensus sequence
# print("\nConsensus category sequence:")
# print(" → ".join(final_consensus))
# print(" → ".join(consensus))