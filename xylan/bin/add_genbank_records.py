import pandas as pd
from Bio import Entrez, SeqIO
from time import sleep
import sys

import ssl
ssl._create_default_https_context = ssl._create_unverified_context

Entrez.email = "yazidhoblos5@gmail.com"

df = pd.read_csv(sys.argv[1])

annotations = []
aa_sequences = []

for pid in df["Protein_ID"]:
    try:
        if pd.isna(pid) or pid.strip() == "":
            annotations.append("")
            aa_sequences.append("")
            continue

        handle = Entrez.efetch(db="protein", id=pid.strip(), rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()

        # Extract annotation (product name)
        annotation = record.description
        for feature in record.features:
            if feature.type == "CDS":
                if "product" in feature.qualifiers:
                    annotation = feature.qualifiers["product"][0].replace(",", ';')
                    break

        annotations.append(annotation)

        # Extract amino acid sequence from the record (protein sequence)
        aa_seq = str(record.seq)
        aa_sequences.append(aa_seq)

        sleep(0.1)

    except Exception as e:
        print(f"Error with {pid}: {e}")
        annotations.append("")
        aa_sequences.append("")

    print(f"Processed {pid} - Annotation: {annotations[-1]} - AA length: {len(aa_sequences[-1])}")

df["Annotation"] = annotations
df["Amino_Acid_Sequence"] = aa_sequences

df.to_csv("xylan_gene_info_corrected_with_seq.csv", index=False)
print("✅ Done. Saved enriched CSV to PUL_gene_metadata_annotated.csv")
