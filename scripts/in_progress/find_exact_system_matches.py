import csv
import argparse

def read_systems(file_path):
    systems = []
    with open(file_path, 'r', encoding='utf-8') as f:
        reader = csv.reader(f, delimiter='\t')
        headers = next(reader)  # Skip header
        for row in reader:
            if len(row) < 6:
                continue  # Skip malformed rows
            system_num = row[0]
            model_GF = set(f.strip() for f in row[4].split(',') if f.strip())
            context_GF = set(f.strip() for f in row[5].split(',') if f.strip())
            systems.append({
                "system_num": system_num,
                "model_GF": model_GF,
                "context_GF": context_GF
            })
    return systems

def find_matching_systems(systems1, systems2):
    matches = []
    for sys1 in systems1:
        for sys2 in systems2:
            if sys1['model_GF'] == sys2['model_GF'] and sys1['context_GF'] == sys2['context_GF']:
                matches.append((sys1['system_num'], sys2['system_num']))
    return matches

def main():
    parser = argparse.ArgumentParser(description="Compare systems from two files based on model_GF and context_GF.")
    parser.add_argument("file1", help="Path to the first TSV file")
    parser.add_argument("file2", help="Path to the second TSV file")
    args = parser.parse_args()

    systems1 = read_systems(args.file1)
    systems2 = read_systems(args.file2)

    matches = find_matching_systems(systems1, systems2)

    print(f"\nFound {len(matches)} matching systems:")
    for s1, s2 in matches:
        print(f"System {s1} in {args.file1} matches System {s2} in {args.file2}")

if __name__ == "__main__":
    main()

