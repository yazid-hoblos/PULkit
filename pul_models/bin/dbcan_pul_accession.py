import requests
from bs4 import BeautifulSoup
import csv

def extract_general_info(soup):
    table_rows = soup.find_all("tr")
    info = {
        "PUL_ID": "-",
        "PubMed": "-",
        "Characterization_Method": "-",
        "Genomic_Accession": "-",
        "Nucleotide_Position_Range": "-",
        "Substrate": "-",
        "Loci": "-",
        "Species": "-",
        "Degradation_or_Biosynthesis": "-"
    }

    for row in table_rows:
        cells = row.find_all("td")
        if len(cells) != 2:
            continue
        label = cells[0].text.strip()
        value_cell = cells[1]

        if label == "PubMed":
            a_tag = value_cell.find("a")
            info["PubMed"] = a_tag.text.strip() if a_tag else value_cell.text.strip()

        elif label == "Characterization method":
            info["Characterization_Method"] = value_cell.text.strip().replace(',', ';')

        elif label == "Genomic accession number":
            info["Genomic_Accession"] = value_cell.text.strip()

        elif label == "Nucleotide position range":
            info["Nucleotide_Position_Range"] = value_cell.text.strip()

        elif label == "Substrate":
            info["Substrate"] = value_cell.text.strip()

        elif label == "Loci":
            info["Loci"] = value_cell.text.strip()

        elif label == "Species":
            info["Species"] = value_cell.text.strip()

        elif label == "Degradation or Biosynthesis":
            info["Degradation_or_Biosynthesis"] = value_cell.text.strip()
    return info


base_url = "https://aca.unl.edu/dbCAN_PUL/dbCAN_PUL/CGC?clusterid=PUL{:04d}"

pul_metadata_list = []
gene_info_list = []

for i in range(1, 797):  # PUL0001 to PUL0796 (latest in dbCAN-PUL last I checked mid-July)
    pul_id = f"PUL{i:04d}"
    url = base_url.format(i)
    response = requests.get(url)
    soup = BeautifulSoup(response.text, "html.parser")

    # Extract general info once per PUL
    general_info = extract_general_info(soup)
    general_info["PUL_ID"] = pul_id
    pul_metadata_list.append(general_info)

    # Extract component genes information
    gene_tab_div = soup.find("div", {"data-tab": "2", "class": "ui attached tab segment"})
    if gene_tab_div:
        gene_table = gene_tab_div.find("table", class_="ui celled table")
        if gene_table:
            rows = gene_table.find_all("tr")
            length = len(rows[1:])
            for row in rows[1:]:  # Skip header row
                cells = row.find_all("td")
                if len(cells) < 6:
                    continue
                gene_data = {
                    "PUL_ID": pul_id,
                    "Gene_Name": cells[0].text.strip(),
                    "Locus_Tag": cells[1].text.strip(),
                    "Protein_ID": cells[2].text.strip(),
                    "Gene_Position": cells[3].text.strip(),
                    "GenBank_Contig": cells[4].text.strip(),
                    "EC_Number": cells[5].text.strip().replace(',', ';')
                }
                
                gene_info_list.append(gene_data)
                
    print(f"Processed {pul_id} - Found {length} gene entries")
    
    # Find gene blocks (alternative way to extract gene info)
    # gene_blocks = soup.find_all("div", id=lambda x: x and x.startswith("message"))

    # if gene_blocks:
    #     for block in gene_blocks:
    #         data = {
    #             "PUL_ID": pul_id,
    #             "Gene_Name": "",
    #             "Locus_Tag": "",
    #             "Protein_ID": "",
    #             "Protein_Link": "",
    #             "Gene_Position": "",
    #             "GenBank_Contig": "",
    #             "Contig_Link": "",
    #             "EC_Number": ""
    #         }

    #         strong_tags = block.find_all("strong")
    #         for strong in strong_tags:
    #             label = strong.text.strip()
    #             if "Gene name" in label and strong.next_sibling:
    #                 data["Gene_Name"] = strong.next_sibling.strip()

    #             elif "contig" in label:
    #                 a_tag = strong.find_next("a")
    #                 if a_tag:
    #                     data["GenBank_Contig"] = a_tag.text.strip()
    #                     data["Contig_Link"] = a_tag.get("href", "")

    #             elif "Locus tag" in label:
    #                 a_tag = strong.find_next("a")
    #                 if a_tag:
    #                     data["Locus_Tag"] = a_tag.text.strip()
    #                     data["Protein_Link"] = a_tag.get("href", "")

    #             elif "Protein ID" in label:
    #                 a_tag = strong.find_next("a")
    #                 if a_tag:
    #                     data["Protein_ID"] = a_tag.text.strip()
    #                     data["Protein_Link"] = a_tag.get("href", "")

    #             elif "Location and strand" in label and strong.next_sibling:
    #                 data["Gene_Position"] = strong.next_sibling.strip()

    #             elif "EC number" in label:
    #                 ec_number = "-"
    #                 a_tag = strong.find_next("a")
    #                 if a_tag and all(c.isdigit() or c == '.' for c in a_tag.text) and len(a_tag.text) >= 7:
    #                     ec_number = a_tag.text.strip()
    #                 data["EC_Number"] = ec_number

    #         gene_info_list.append(data)

    # print(f"Processed {pul_id} - Found {len(gene_blocks)} gene blocks")

# Save PUL metadata into CSV file (one row per PUL)
metadata_fields = [
    "PUL_ID", "PubMed", "Characterization_Method", "Genomic_Accession",
    "Nucleotide_Position_Range", "Substrate", "Loci", "Species", "Degradation_or_Biosynthesis"
]

with open("dbcan_pul_metadata.csv", "w", newline='', encoding='utf-8') as f:
    writer = csv.DictWriter(f, fieldnames=metadata_fields)
    writer.writeheader()
    writer.writerows(pul_metadata_list)

# Save gene info into CSV file (one row per gene)
gene_fields = [
    "PUL_ID", "Gene_Name", "Locus_Tag", "Protein_ID",
    "Protein_Link", "Gene_Position", "GenBank_Contig",
    "Contig_Link", "EC_Number"
]

with open("dbcan_pul_gene_info.csv", "w", newline='', encoding='utf-8') as f:
    writer = csv.DictWriter(f, fieldnames=gene_fields)
    writer.writeheader()
    writer.writerows(gene_info_list)

print("✅ Done. Metadata saved to dbcan_pul_metadata.csv and gene info saved to dbcan_pul_gene_info.csv")