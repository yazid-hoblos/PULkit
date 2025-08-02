"""Loads dbCAN and susCD annotation files to infer PUL component families."""

import os
import argparse
from collections import defaultdict
from collections import Counter
import re
from plot_pul_components import plot_pul_categories


def infer_family(annotation):
    """
    Classify proteins into functional categories using regex patterns.
    Returns tuple of (category, original_annotation).
    """
    if not annotation or annotation.strip() == '':
        return ("Unknown", annotation)
    
    # Convert to lowercase for case-insensitive matching
    ann_lower = annotation.lower()
    
    # 1. Sugar/Carbohydrate Metabolism (most common in dataset)
    if re.search(r'(xylose\s+isomerase|xylulose\s+kinase|xylulokinase|mannonate\s+dehydratase|glucuronate\s+isomerase| \
                 uronate\s+isomerase|aldose\s+1-epimerase|transaldolase|transketolase|fructose.*aldolase|ribulokinase| \
                 l-arabinose\s+isomerase|d-mannonate\s+hydrolase|d-mannonate\s+oxidoreductase|mannonate\s+dehydrogenase| \
                 hexuronic\s+acid\s+isomerase)', ann_lower):
        return ("Sugar_metabolism", annotation)
    
    # 2. Mobile Genetic Elements
    if re.search(r'(transposase|integrase|is[0-9]+|insertion\s+sequence)', ann_lower):
        return ("Mobile_elements", annotation)
    
    # 3. Oxidoreductases
    if re.search(r'(sdr\s+family\s+oxidoreductase|dehydrogenase|oxidoreductase|reductase)', ann_lower):
        return ("Oxidoreductases", annotation)
    
    # 4. Ribosomal Proteins
    if re.search(r'(ribosomal\s+protein|50s\s+ribosomal|30s\s+ribosomal)', ann_lower):
        return ("Ribosomal", annotation)
    
    # 5. tRNA Processing
    if re.search(r'(trna.*biosynthesis|2-thiocytidine\s+biosynthesis|trna.*thiouridylase|aminoacyl-trna\s+hydrolase)', ann_lower):
        return ("RNA_processing", annotation)
    
    # 6. Cell Adhesion/Surface
    if re.search(r'(fasciclin|beta-ig-h3)', ann_lower):
        return ("Cell_adhesion", annotation)
    
    # # 7. Vitamin Biosynthesis
    # if re.search(r'(thiamine\s+biosynthesis|apbe)', ann_lower):
    #     return ("Vitamin_biosynthesis", annotation)
    
    # # 8. Nucleotide Metabolism
    # if re.search(r'(nucleoside\s+phosphorylase|nucleoside\s+hydrolase|purine\s+nucleoside)', ann_lower):
    #     return ("Nucleotide_metabolism", annotation)
    
    # 9. Kinases 
    if re.search(r'(kinase|phosphofructokinase)', ann_lower):
        return ("Kinases", annotation)
    
    # 10. Hydrolases
    if re.search(r'(hydrolase|sgnh.*hydrolase|gdsl.*hydrolase|exonuclease|endonuclease)', ann_lower):
        return ("Hydrolases", annotation)
    
    # 11. Transport/Surface Proteins
    if re.search(r'(sugar-binding\s+protein|transport\s+system|outer\s+membrane|efflux\s+protein|lipoprotein|sorting\s+domain|glycan-binding\s+surface)', ann_lower):
        return ("Transport/Surface", annotation)
    
    # # 12. DNA Repair/Processing
    # if re.search(r'(excinuclease|dna-binding\s+protein|hu\s+family)', ann_lower):
    #     return ("DNA Processing", annotation)
    
    # 13. General Metabolic Enzymes
    if re.search(r'(aldolase|epimerase|decarboxylase|synthetase|transferase|isomerase)', ann_lower):
        return ("General_metabolism", annotation)
    
    # # 14. Stress Response/Chaperones
    # if re.search(r'(hsp20|alpha\s+crystallin|heat\s+shock)', ann_lower):
    #     return ("Stress_response", annotation)
    
    # 15. Conserved Unknown Function (DUF domains)
    if re.search(r'(duf[0-9]+|domain.*unknown\s+function)', ann_lower):
        return ("DUF", annotation)
    
    # 16. Peptidases/Proteases
    if re.search(r'(peptidase|protease|deformylase)', ann_lower):
        return ("Peptidase", annotation)
    
    # 17. Phosphatases
    if re.search(r'(phosphatase)', ann_lower):
        return ("Phosphatase", annotation)
    
    # 18. Operon/Regulon 
    if re.search(r'\b(operon|regulon|[A-Z][a-z]{2}[A-Z0-9])\b', annotation):
        return ("Operon/Regulon", annotation)
    
    # Uncomment if needed:
    # if "SusD" in annotation or "RagB" in annotation:
    #     return ("SusD", annotation)
    
    # Default categories
    if "hypothetical" in ann_lower:
        return ("Hypothetical", annotation)
    
    # Everything else
    return ("Other", annotation)


# Loads dbCAN annotations
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
                    if fields[2].strip() != '-':
                        family = fields[2].strip().split("+")[0].split("_")[0].split("(")[0]
                    elif fields[3].strip() != '-':
                        family = fields[3].strip().split("+")[0].split("_")[0].split("(")[0]
                    else:
                        family = fields[4].strip().split("+")[0].split("_")[0].split("(")[0]
                    protein_to_category[fields[0]].append(('CAZyme', family))

    # TFs and STPs from column 3
    for category in ["TF", "STP"]:
        if os.path.isfile(file_mapping[category]):
            with open(file_mapping[category]) as f:
                for line in f:
                    if not line.strip():
                        continue
                    fields = line.strip().split("\t")
                    if len(fields) >= 3:
                        # if not category in protein_to_category[fields[2]]:
                        protein_to_category[fields[2]].append((category, fields[0]))

    # Peptidases, Transporters, Sulfatases from column 3
    for category in ["Peptidase", "Transporter", "Sulfatase"]:
        if os.path.isfile(file_mapping[category]):
            with open(file_mapping[category]) as f:
                for line in f:
                    if not line.strip():
                        continue
                    fields = line.strip().split("\t")
                    if len(fields) >= 3:
                        # if not category in protein_to_category[fields[2]]:
                        protein_to_category[fields[2]].append((category, fields[0]))

    if suscd_file:
        suscd_categories = load_suscd_categories(suscd_file)
        for pid, cat in suscd_categories.items():
            protein_to_category[pid] = [cat]

    # Prioritize category assignment (in case of multiple annotation per protein)
    # order commented below: CAZyme > Transporter > STP/TF 
    reduced_protein_to_category = {}
    for protein_id, categories in protein_to_category.items():
        reduced_protein_to_category[protein_id] = categories[0]  # Default to first category 
    #     if 'CAZyme' in categories:
    #         reduced_protein_to_category[protein_id] = 'CAZyme'
    #     elif 'Transporter' in categories:
    #         reduced_protein_to_category[protein_id] = 'Transporter'
    #     elif 'STP' in categories and 'TF' in categories:
    #         reduced_protein_to_category[protein_id] = 'TF-STP'
    #     else:
    #         reduced_protein_to_category[protein_id] = categories[0]

    return reduced_protein_to_category


# Loads SusCD annotations
def load_suscd_categories(suscd_file):
    protein_to_suscd = defaultdict(list)
    if os.path.isfile(suscd_file):
        with open(suscd_file) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.strip().split()  # splits by ANY whitespace
                # if len(fields) >= 3:
                domain_name = fields[0]      
                protein_id = fields[2]       # target name (query protein)
                if "SusD" in domain_name:
                    protein_to_suscd[protein_id] = ("Transporter", "SusD")
                else:
                    protein_to_suscd[protein_id] = ("Transporter", "SusC")
    return protein_to_suscd


# Loads the output of add_genbank_records.py
def load_protein_annotations(annotation_file):
    protein_annotation_map = {}
    with open(annotation_file, 'r') as f:
        for line in f:
            if "PUL_ID" in line or not line.strip():
                continue
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
        for line in f:
            if "PUL_ID" in line or not line.strip():
                continue
            parts = line.strip().split(',')
            if len(parts) < 11:
                continue
            pul_id = parts[0].strip()
            protein_id = pul_id + '|' + parts[3].strip()
            pul_proteins[pul_id].append(protein_id)
    return pul_proteins


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    # Positional arguments
    parser.add_argument("dbcan_dir",
                       help="dbCAN easy CGC output directory")
    parser.add_argument("annotation_file", 
                       help="Output file of add_genbank_records.py")
    parser.add_argument("--suscd",
                       help="SusCD annotation file (optional)")

    args = parser.parse_args()

    # create output directory if it does not exist and overwrite if it does
    output_dir = "pul_component_patterns"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    else:
        for file in os.listdir(output_dir):
            file_path = os.path.join(output_dir, file)
            if os.path.isfile(file_path):
                os.remove(file_path)
    
    dbcan_dir = args.dbcan_dir
    annotation_file = args.annotation_file
    suscd_file= args.suscd 

    protein_to_category = load_protein_categories_from_directory(dbcan_dir, suscd_file=suscd_file)
    protein_annotations = load_protein_annotations(annotation_file)
    pul_proteins = group_proteins_by_pul(annotation_file)

    unannotated = set()
    enriched_pul_data = defaultdict(list)

    for pul_id, proteins in pul_proteins.items():
        for pid in proteins:
            if pid in protein_to_category:
                category = protein_to_category[pid]
            else:
                unannotated.add(pid)
                annotation = protein_annotations[pid]
                category = infer_family(annotation)
            enriched_pul_data[pul_id].append({"protein_id": pid, "category": category})

    output_file = os.path.join(output_dir, "detailed_components.txt")
    with open(output_file, 'w') as f:
        for pul_id, components in enriched_pul_data.items():
            f.write(f"=== {pul_id} ===\n")
            for comp in components:
                type_cat, category = comp['category']
                protein_id = comp['protein_id'].split('|')[1]
                f.write(f"  {protein_id:<15} | {type_cat:<18} | {category}\n")

    pul_category_sequences = {}
    # Sort PULs by number of components (descending)
    sorted_puls = sorted(enriched_pul_data.items(), key=lambda x: -len(x[1]))
    pul_category_sequences = {}
    for pul_id, components in sorted_puls:
        categories = [comp["category"] for comp in components]
        pul_category_sequences[pul_id] = categories

    # save pul_category_sequences to a file
    output_file = os.path.join(output_dir, "category_sequences.txt")
    with open(output_file, 'w') as f:
        for pul_id, seq in pul_category_sequences.items():
            f.write(f"{pul_id}: {' → '.join([s[1] if 'Sus' in s[1] else s[0] for s in seq])}\n")

    # Print summary statistics
    print("=== Summary Statistics ===")
    print(len(pul_proteins), "total PULs found")
    print(len(list(unannotated)), "total unannotated proteins")
    print(len(list(set(protein_to_category.keys()))), "total annotated proteins\n")

    print("=== Most common CAZyme families ===")
    all_cazyme_families = []
    for pul_id, seq in pul_category_sequences.items():
        cazyme_fams = [s[1] for s in seq if s[0] == 'CAZyme']
        all_cazyme_families.extend(cazyme_fams)

    cazyme_counter = Counter(all_cazyme_families)
    for family, count in cazyme_counter.most_common(20):
        print(f"{family}: {count} occurrences")
    print("Total unique CAZyme families:", len(cazyme_counter), "\n")
    
    # To examine the genbank annotations of unannotated proteins left as "Other"
    # for protein in unannotated:
    #     if infer_family(protein_annotations[protein])[0] in ("Other"):
    #         print(f"{protein} -  {protein_annotations[protein]}")


    pul_category_sequences_reduced = {pul_id: [s[1] if 'Sus' in s[1] else s[0] for s in seq] for pul_id, seq in pul_category_sequences.items()}
    # Consensus patterns (n-grams)
    def get_ngrams(seq, n):
        return [tuple(seq[i:i+n]) for i in range(len(seq)-n+1)]

    ngram_counter = Counter()

    for seq in pul_category_sequences_reduced.values():
        for n in range(2, 5):  # 2-grams to 4-grams
            ngrams = get_ngrams(seq, n)
            ngram_counter.update(ngrams)

    # Show top consensus patterns
    print("=== Top consensus patterns (n-grams) ===")
    for pattern, count in ngram_counter.most_common(20):
        print(f"{' → '.join(pattern)} ({count} occurrences)")

    # # Sort from longest and most frequent to shortest
    # sorted_ngrams = sorted(ngram_counter.items(), key=lambda x: (-len(x[0]), -x[1]))

    # # Greedily build one consensus sequence
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

    org_file = 'data/xylan_pul_org.csv'
    plot_pul_categories(pul_category_sequences_reduced, 'all_categories_considered.png', pul_to_organism_file=org_file)
    custom_order = ["CAZyme", "Transporter", "SusC", "SusD", "STP", "TF-STP", "Sulfatase", "Peptidase", "Hypothetical", "Other"] # all TFs were found to be STPs as well
    plot_pul_categories(pul_category_sequences_reduced, 'custom_ordered_plot.png', pul_to_organism_file=org_file, custom_order=custom_order)



if __name__ == "__main__":
    main()