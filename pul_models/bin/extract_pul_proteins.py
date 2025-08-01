''' Creates a FASTA file of all sequences in the input CSV file'''

import csv
import sys

if len(sys.argv) != 3:
    print("Usage: python extract.py <input_file> <output_fasta>")
    print("Input file is output of add_genbank_records.py")
    sys.exit(1)

input_file = sys.argv[1]
output_fasta = sys.argv[2]

with open(input_file, newline='') as infile, open(output_fasta, 'w') as outfile:
    reader = csv.reader(infile, delimiter=',') 
    for row in reader:
        if 'ID' in row[0]:  # skip header row3
            continue
        pul_id = row[0]
        protein_id = row[3]
        if protein_id == "-":
            continue
        seq = row[10].strip()
        header = f">{pul_id}|{protein_id}"
        outfile.write(f"{header}\n{seq}\n")
