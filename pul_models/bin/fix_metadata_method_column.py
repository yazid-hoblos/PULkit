'''Fix method column issue in metadata file; should be fixed if you re-run dbcan_pul_accession.py'''

import csv
import sys

if len(sys.argv) != 3:
    print("Usage: python extract.py <input_file> <output_fasta>")
    print("Input is the metadata output file of dbcan_pul_accession.py")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

with open(input_file, newline='', encoding='utf-8') as infile, \
     open(output_file, 'w', newline='', encoding='utf-8') as outfile:

    # Read CSV with default delimiter ',' and quotechar '"'
    reader = csv.reader(infile, delimiter=',', quotechar='"')
    writer = csv.writer(outfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

    for row in reader:
        # Replace commas with semicolons inside each field
        new_row = [field.replace(',', ';') for field in row]
        writer.writerow(new_row)
