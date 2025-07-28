import requests
from bs4 import BeautifulSoup
import csv

base_url = "https://aca.unl.edu/dbCAN_PUL/dbCAN_PUL/CGC?clusterid=PUL{:04d}"
results = []

for i in range(1, 674):  # PUL0001 to PUL0673
    pul_id = f"PUL{i:04d}"
    url = base_url.format(i)
    response = requests.get(url)
    soup = BeautifulSoup(response.text, "html.parser")
    print(soup.prettify())
    exit(0)

    # Try parsing <div id="message"> layout first
    gene_blocks = soup.find_all("div", id=lambda x: x and x.startswith("message"))
    print(gene_blocks)
    exit(0)

    if gene_blocks:
        for block in gene_blocks:
            data = {
                "PUL_ID": pul_id,
                "Gene_Name": "",
                "Locus_Tag": "",
                "Protein_ID": "",
                "Protein_Link": "",
                "Gene_Position": "",
                "GenBank_Contig": "",
                "Contig_Link": "",
                "EC_Number": ""
            }

            strong_tags = block.find_all("strong")
            for strong in strong_tags:
                label = strong.text.strip()
                # print(f"Processing {pul_id} - {label}")
                if "Gene name" in label:
                    if strong.next_sibling:
                        data["Gene_Name"] = strong.next_sibling.strip()

                elif "contig" in label:
                    a_tag = strong.find_next("a")
                    if a_tag:
                        data["GenBank_Contig"] = a_tag.text.strip()
                        data["Contig_Link"] = a_tag.get("href", "")


                elif "Locus tag" in label:
                    a_tag = strong.find_next("a")
                    if a_tag:
                        data["Locus_Tag"] = a_tag.text.strip()
                        data["Protein_Link"] = a_tag.get("href", "")

                elif "Protein ID" in label:
                    a_tag = strong.find_next("a")
                    if a_tag:
                        data["Protein_ID"] = a_tag.text.strip()
                        data["Protein_Link"] = a_tag.get("href", "")

                elif "Location and strand" in label:
                    if strong.next_sibling:
                        data["Gene_Position"] = strong.next_sibling.strip()
                    
                elif "EC number" in label:
                    ec_number = "-"
                    a_tag = strong.find_next("a")
                    if a_tag and all(char.isdigit() or char=='.' for char in a_tag.text) and len(a_tag.text) >= 7:
                        ec_number = a_tag.text.strip()
                    data["EC_Number"] = ec_number
                    # a_tag = strong.find_next("a")
                    # if a_tag:
                    #     data["EC Number"] = a_tag.text.strip()

            results.append(data)

    print(f"Processed {pul_id} - Found {len(gene_blocks)} gene blocks")
    # exit(0)

# Save results to CSV
fieldnames = [
    "PUL_ID", "Gene_Name", "Locus_Tag", "Protein_ID", 
    "Protein_Link", "Gene_Position", "GenBank_Contig", 
    "Contig_Link", "EC_Number"
]

with open("PUL_gene_metadata.csv", "w", newline='', encoding='utf-8') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(results)

print("✅ Done. Extracted data saved to PUL_gene_metadata.csv")
