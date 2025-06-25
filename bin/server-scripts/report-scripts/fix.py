import csv

cazymes = set()
with open('cazymes') as f:
    for line in f:
        cazymes.add(line.strip())

with open('HMM-Databases/complete_with_stp_corrected.updated_25.tsv') as infile, \
     open('HMM-Databases/complete_with_stp_corrected.updated_25.modified.tsv', 'w', newline='') as outfile:
    
    reader = csv.reader(infile, delimiter='\t')
    writer = csv.writer(outfile, delimiter='\t')

    for row in reader:
        if row[0] == 'name':  # header
            writer.writerow(row)
            continue
        protein_name = row[4]
        if protein_name not in cazymes:
            row[7] = '1e-4'  # eval_threshold
            row[8] = '1e-4'  # ieval_threshold
        writer.writerow(row)

