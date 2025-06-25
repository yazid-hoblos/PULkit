import requests
from bs4 import BeautifulSoup
import urllib.parse
import csv
import time
import sys

if len(sys.argv) != 2:
    print("Usage: python script.py '<Species name>'")
    sys.exit(1)

species_input = sys.argv[1]
species_input = species_input.replace("_", " ")
headers = {
    "User-Agent": "Mozilla/5.0 (compatible; Python script)"
}

base_url = "https://www.cazy.org/PULDB/"
species_encoded = urllib.parse.quote_plus(species_input)
species_url = f"{base_url}index.php?sp_name={species_encoded}"

print(f"Fetching PULDB data for: {species_input}")
print(f"URL: {species_url}")

# headers = {
#     "User-Agent": "Mozilla/5.0 (compatible; Python script)"
# }

# base_url = "https://www.cazy.org/PULDB/"

# species_input = input("Enter species name (e.g., Prevotella denticola): ").strip()
# species_encoded = urllib.parse.quote_plus(species_input)
# species_url = f"{base_url}index.php?sp_name={species_encoded}"

def write_to_csv(csv_filename, pul_list):
    csv_filename = f"{species_input.replace(' ', '_')}_PUL_proteins.csv"
    with open(csv_filename, mode='w', newline='', encoding='utf-8') as csvfile:
        writer = csv.writer(csvfile)
        # Header row
        writer.writerow(["PUL", "Species", "NCBI_Taxon_ID", "Protein_IDs"])

        for pul_name, species, pul_url in pul_list:
            print(f"Processing {pul_name} ({species}) → {pul_url}")
            proteins, taxon_id = get_pul_proteins(pul_url)
            # time.sleep(1)
            # Join all protein IDs with semicolon (or comma if you prefer)
            protein_str = ";".join(proteins)
            writer.writerow([pul_name, species, taxon_id,protein_str])

    print(f"\n✅ CSV file created: {csv_filename}")

def get_pul_proteins(pul_url):
    proteins = []
    response = requests.get(pul_url, headers=headers)
    if response.status_code == 200:
        soup = BeautifulSoup(response.content, "html.parser")

        taxon_id = "N/A"
        for a_tag in soup.find_all("a", href=True):
            href = a_tag["href"]
            if "wwwtax.cgi?id=" in href:
                taxon_id = a_tag.get_text(strip=True)
                break

        for a_tag in soup.find_all("a", href=True):
            href = a_tag["href"]
            if "prot=" in href:
                protein_name = a_tag.get_text(strip=True)
                proteins.append(protein_name) 
    else:
        print(f"Failed to fetch PUL page: {pul_url}")
    return proteins, taxon_id



pul_list = []

response = requests.get(species_url, headers=headers)
taxon_id = "N/A"
if response.status_code == 200:
    soup = BeautifulSoup(response.content, "html.parser")
    tables = soup.find_all("table")
    if len(tables) >= 3:
        table = tables[2]
        rows = table.find_all("tr")[1:]
        for row in rows:
            cells = row.find_all("td")
            if len(cells) >= 2:
                species = cells[0].get_text(strip=True)
                pul_tag = cells[1].find("a", href=True)
                if pul_tag:
                    pul_name = pul_tag.get_text(strip=True)
                    pul_url = base_url + pul_tag["href"]
                    pul_list.append((pul_name, species,pul_url))
    else:
        print("⚠️ Could not find PUL table for this species.")
else:
    print(f"❌ Failed to fetch page for species: {species_input} (status code {response.status_code})")

if pul_list:
    print(f"Found {len(pul_list)} PULs for species: {species_input}")
    write_to_csv(f"{species_input.replace(' ', '_')}_PUL_proteins.csv", pul_list)


