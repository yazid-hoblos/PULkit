import requests
import sys

def get_gcf_ids(taxon_id):
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    key="ffbc985055058237daf2ef7c859bdb3d6d08"
    params = {
        "db": "assembly",
        "term": f"txid{taxon_id}[Organism:exp]",
        "retmode": "json",
        "retmax": 100  # you can increase this if needed
        ,"api_key": key
    }
    response = requests.get(url, params=params)
    response.raise_for_status()
    ids = response.json()["esearchresult"].get("idlist", [])
    if not ids:
        return []

    summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    response = requests.get(summary_url, params={
        "db": "assembly",
        "id": ",".join(ids),
        "retmode": "json"
    })
    response.raise_for_status()
    summaries = response.json()["result"]
    return [summaries[i]["assemblyaccession"] for i in summaries if i != "uids"]

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py taxon_ids.txt")
        sys.exit(1)

    filename = sys.argv[1]
    with open(filename, 'r') as f:
        taxon_ids = [line.strip() for line in f if line.strip()]

    for taxon_id in taxon_ids:
        try:
            gcf_ids = get_gcf_ids(taxon_id)
            for gcf in gcf_ids:
                print(f"{taxon_id}\t{gcf}")
        except Exception as e:
            print(f"Error fetching taxon {taxon_id}: {e}", file=sys.stderr)
