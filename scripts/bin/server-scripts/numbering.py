# take input file as argument
import sys
import os
import csv


# This script processes a TSV file to clean up the first column by removing
# any values that do not match a predefined set of individuals.
if len(sys.argv) != 2:
    print("Usage: python clean_tsv.py <input_file>")
    sys.exit(1)
input_file = sys.argv[1]
output_file = "cleaned.tsv"

# Step 1: Read all rows, split by tabs
with open(input_file) as f:
    rows = [line.strip().split('\t') for line in f if line.strip()]

# Step 2: Identify all individual values in first column

# print([row[6].strip() for row in rows if row])  # Debug: print all first column values
individuals = set(row[0].strip() for row in rows if ',' not in row[0] if row[6].strip() != '') 
print(f"Identified {individuals}")
# Step 3: Process each row and update first column if needed
cleaned_rows = []
for row in rows:
    original = row[0].strip()
    if ',' in original:
        # print(f"Processing row: {original}")
        parts = [p.strip() for p in original.split(', ')]
        # print(f"Split parts: {parts}")
        present_parts = [p for p in parts if p in individuals]
        # print(f"Present parts: {present_parts}")
        if present_parts:
            row[0] = ','.join(present_parts)  # update first column
            # print(f"Updated row: {row[0]}")
        else:
            row[0] = original  # leave it unchanged if nothing matches
    # print(f"Final row: {row}")
    cleaned_rows.append(row)

# Step 4: Write output
with open(output_file, 'w') as f:
    for row in cleaned_rows:
        f.write('\t'.join(row) + '\n')

print(f"Modified file written to: {output_file}")
