'''Takes projection file as input and removes systems that are matched by same GF but have no GF only corresponding to them.'''

import sys

if len(sys.argv) != 2:
    print("Usage: python script.py <input_projection_file>")
    sys.exit(1)
    
input_file = sys.argv[1]
output_file = f"{input_file[:-4]}_reduced.tsv"

with open(input_file) as f:
    rows = [line.strip().split('\t') for line in f if line.strip()]

# Identify all system ids
individuals = set(row[0].strip() for row in rows if ',' not in row[0] if row[6].strip() != '') # rows corresponding to model GFs and a single system

# Process each row and update first column if needed
cleaned_rows = []
for row in rows:
    original = row[0].strip()
    if ',' in original:
        parts = [p.strip() for p in original.split(', ')]
        present_parts = [p for p in parts if p in individuals]
        if present_parts:
            row[0] = ','.join(present_parts)  # update first column
        else:
            row[0] = original  # leave it unchanged if nothing matches
    cleaned_rows.append(row)

# Write output
with open(output_file, 'w') as f:
    for row in cleaned_rows:
        f.write('\t'.join(row) + '\n')

print(f"Modified file written to: {output_file}")
