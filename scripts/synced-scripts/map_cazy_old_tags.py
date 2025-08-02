import csv
import sys 


if len(sys.argv) != 3:
    print("Usage: python merge.py <mappings_file> <cazy_data_file>")
    sys.exit(1)

file1 = sys.argv[1]
file2 = sys.argv[2]

# Load mappings into a dictionary
mapping = {}
with open(file1) as f:
    for line in f:
        contig, start, end, strand, locus_tag, old_tag = line.strip().split("\t")
        mapping[old_tag] = (contig, start, end, strand, locus_tag)

# Process PUL file
with open(file2) as f:
    reader = csv.reader(f)
    for row in reader:
        pul_id = row[0].strip()
        old_tags = row[3].split(";")
        for old_tag in old_tags:
            if old_tag in mapping:
                contig, start, end, strand, locus_tag = mapping[old_tag]
                print(f"{pul_id}\t{old_tag}\t{contig}\t{start}\t{end}\t{strand}\t{locus_tag}")
            else:
                print(f"{pul_id}\t{old_tag}\tNOT_FOUND", file=sys.stderr)
