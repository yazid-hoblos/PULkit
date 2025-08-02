'''Infers annotation through simple regex and prints the corresponding gene components per PUL'''

import re
from collections import defaultdict
import json
import sys

def infer_family(annotation):
    annotation = annotation.lower()

    # Explicit family matches
    if match := re.search(r'glyco(?:s|z)id[ae]? hydrolase(?: family)? (\d+)', annotation):
        return f"GH{match.group(1)}"
    if match := re.search(r'family (\d+).*glyco(?:s|z)id[ae]? hydrolase', annotation):
        return f"GH{match.group(1)}"
    if match := re.search(r'polysaccharide lyase(?: family)? (\d+)', annotation):
        return f"PL{match.group(1)}"
    if match := re.search(r'carbohydrate esterase(?: family)? (\d+)', annotation):
        return f"CE{match.group(1)}"
    if match := re.search(r'(?:cbm|carbohydrate[- ]binding module)(?: family)? (\d+)', annotation):
        return f"CBM{match.group(1)}"

    # Known enzyme activity
    if "xylanase" in annotation or "xyn" in annotation:
        return "Xylanase"
    if "xylosidase" in annotation or "xyloside" in annotation or "beta-xylosidase" in annotation:
        return "Xylosidase"
    if "glucanase" in annotation:
        return "Glucanase"
    if "glucosidase" in annotation:
        return "Glucosidase"
    if "epimerase" in annotation:
        return "Epimerase"
    if "isomerase" in annotation:
        return "Isomerase"
    if "esterase" in annotation:
        return "Esterase"
    if "transaldolase" in annotation:
        return "Transaldolase"
    if "transglutaminase" in annotation:
        return "Transglutaminase"
    if "dehydratase" in annotation:
        return "Dehydratase"
    if "decarboxylase" in annotation:
        return "Decarboxylase"
    if "galactosidase" in annotation:
        return "Galactosidase"
    if "endo-1" in annotation:
        return "Endo-1"

    # Sus system
    if "susc" in annotation or "raga" in annotation or "tonb-linked outer membrane protein" in annotation:
        return "SusC"
    if "susd" in annotation or "ragb" in annotation:
        return "SusD"

    # Transporters
    if "abc transporter" in annotation or "abc-type transporter" in annotation:
        return "ABC Transporter"
    if "permease" in annotation or "transport protein" in annotation or "transport system" in annotation:
        return "Transporter"
    if "solute-binding" in annotation or "lipo-binding" in annotation:
        return "Binding Protein"
    if "tonb" in annotation and "receptor" in annotation:
        return "TonB Receptor"
    if "fluoride ion transporter" in annotation:
        return "Ion Transporter"

    # Regulators and sensors
    if "two-component" in annotation or "sensor kinase" in annotation:
        return "Two-Component System"
    if "regulator" in annotation:
        return "Regulator"
    if "sigma" in annotation:
        return "Sigma Factor"

    # Housekeeping and other functional proteins
    if "tRNA" in annotation:
        return "tRNA-related"
    if "transposase" in annotation:
        return "Transposase"
    if "enolase" in annotation:
        return "Enolase"
    if "atpase" in annotation:
        return "ATPase"
    if "anaphase" in annotation:
        return "Cell Cycle"
    if "vitamin b12" in annotation or "btu" in annotation:
        return "Vitamin B12 Transport"

    # Default/hypothetical
    if "hypothetical" in annotation:
        return "Hypothetical"
    return "Other"


def parse_pul_file(filename):
    pul_dict = defaultdict(list)

    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split(',')
            if len(parts) < 11:
                continue  # Skip malformed lines

            pul_id = parts[0]
            locus_tag = parts[2]
            protein_id = parts[3]
            annotation = parts[9]
            family = infer_family(annotation)

            component = {
                "locus_tag": locus_tag,
                "protein_id": protein_id,
                "function": annotation,
                "family": family
            }

            pul_dict[pul_id].append(component)

    return pul_dict


if len(sys.argv) < 2:
    print("Usage: python infer_family_from_annotation.py <pul_file>")
    print("Input is expected to have same format as add_genbank_records.py output.")
    sys.exit(1)
    
pul_file = sys.argv[1]
pul_data = parse_pul_file(pul_file)

# Pretty print JSON
# print(json.dumps(pul_data, indent=2))

for pul_id, components in pul_data.items():
    families = [comp["family"] for comp in components]
    print(pul_id, "->", families)

