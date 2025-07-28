import csv

input_file = "data/xylan_gene_info_corrected_with_seq.csv"
output_file = "fixed_metadata2.csv"

with open(input_file, newline='', encoding='utf-8') as infile, \
     open(output_file, 'w', newline='', encoding='utf-8') as outfile:

    # Read CSV with default delimiter ',' and quotechar '"'
    reader = csv.reader(infile, delimiter=',', quotechar='"')
    writer = csv.writer(outfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

    for row in reader:
        # Replace commas with semicolons inside each field
        new_row = [field.replace(',', ';') for field in row]
        writer.writerow(new_row)
