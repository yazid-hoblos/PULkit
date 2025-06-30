import pandas as pd
from Bio import Entrez, SeqIO
from time import sleep

import ssl
from urllib.request import urlopen
ssl._create_default_https_context = ssl._create_unverified_context


# Your email (required by NCBI)
Entrez.email = "yazidhoblos5@gmail.com"

# Load your PUL gene CSV
df = pd.read_csv("PUL_gene_metadata.csv")

# Placeholder for annotations
annotations = []

# Iterate through protein IDs
for pid in df["Protein_ID"]:
    try:
        # Skip if missing or NaN
        if pd.isna(pid) or pid.strip() == "":
            annotations.append("")
            continue

        # Fetch record from NCBI
        handle = Entrez.efetch(db="protein", id=pid.strip(), rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")

        # Try to extract from CDS or product annotations
        annotation = record.description  # fallback
        for feature in record.features:
            if feature.type == "CDS":
                if "product" in feature.qualifiers:
                    annotation = feature.qualifiers["product"][0]
                    break

        annotations.append(annotation)

        sleep(0.5)  # Be nice to NCBI servers

    except Exception as e:
        print(f"Error with {pid}: {e}")
        annotations.append("")
    print(f"Processed {pid} - Annotation: {annotations[-1]}")

# Add annotations to DataFrame
df["Annotation"] = annotations

# Save new CSV
df.to_csv("PUL_gene_metadata_annotated.csv", index=False)
print("✅ Done. Saved enriched CSV to PUL_gene_metadata_annotated.csv")
