import csv
import argparse

def read_systems(file_path):
    systems = []
    with open(file_path, 'r', encoding='utf-8') as f:
        reader = csv.reader(f, delimiter='\t')
        headers = next(reader)
        for row in reader:
            if len(row) < 6:
                continue
            system_num = row[0]
            model_GF = set(f.strip() for f in row[4].split(',') if f.strip())
            context_GF = set(f.strip() for f in row[5].split(',') if f.strip())
            systems.append({
                "system_num": system_num,
                "model_GF": model_GF,
                "context_GF": context_GF
            })
    return systems

def compare_systems(systems1, systems2):
    matched_1 = set()
    matched_2 = set()
    matches = []

    for i, sys1 in enumerate(systems1):
        for j, sys2 in enumerate(systems2):
            if sys1['model_GF'] == sys2['model_GF'] and sys1['context_GF'] == sys2['context_GF']:
                matches.append((sys1['system_num'], sys2['system_num']))
                matched_1.add(i)
                matched_2.add(j)
                break  # assume one-to-one matching

    unmatched_1 = [systems1[i]['system_num'] for i in range(len(systems1)) if i not in matched_1]
    unmatched_2 = [systems2[i]['system_num'] for i in range(len(systems2)) if i not in matched_2]

    return matches, unmatched_1, unmatched_2

def main():
    parser = argparse.ArgumentParser(description="Compare systems from two files based on model_GF and context_GF.")
    parser.add_argument("file1", help="Path to the first TSV file")
    parser.add_argument("file2", help="Path to the second TSV file")
    args = parser.parse_args()

    systems1 = read_systems(args.file1)
    systems2 = read_systems(args.file2)

    matches, unmatched1, unmatched2 = compare_systems(systems1, systems2)

    print(f"\n✅ Found {len(matches)} matching systems:")
    for s1, s2 in matches:
        print(f"  - File1 System {s1} <==> File2 System {s2}")

    print(f"\n❌ {len(unmatched1)} unmatched systems in {args.file1}:")
    for s in unmatched1:
        print(f"  - System {s}")

    print(f"\n❌ {len(unmatched2)} unmatched systems in {args.file2}:")
    for s in unmatched2:
        print(f"  - System {s}")

if __name__ == "__main__":
    main()

