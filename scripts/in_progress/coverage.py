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

def compute_best_coverage(source_systems, target_systems):
    results = []
    for src in source_systems:
        best_model_cov = 0
        best_context_cov = 0
        best_match_num = None

        for tgt in target_systems:
            if len(src['model_GF']) == 0 or len(src['context_GF']) == 0:
                continue  # Avoid divide by zero

            model_overlap = len(src['model_GF'] & tgt['model_GF'])
            context_overlap = len(src['context_GF'] & tgt['context_GF'])

            model_cov = model_overlap / len(src['model_GF']) * 100
            context_cov = context_overlap / len(src['context_GF']) * 100

            avg_cov = (model_cov + context_cov) / 2

            if avg_cov > (best_model_cov + best_context_cov) / 2:
                best_model_cov = model_cov
                best_context_cov = context_cov
                best_match_num = tgt['system_num']

        results.append({
            "system_num": src['system_num'],
            "best_match": best_match_num,
            "model_GF_coverage": round(best_model_cov, 2),
            "context_GF_coverage": round(best_context_cov, 2)
        })
    return results

def print_coverage_report(file_label, coverage_data):
    print(f"\n🔎 Coverage for systems in {file_label}:")
    for item in coverage_data:
        match = item["best_match"] or "None"
        print(f"  - System {item['system_num']} best matches System {match}")
        print(f"    Model_GF coverage:   {item['model_GF_coverage']}%")
        print(f"    Context_GF coverage: {item['context_GF_coverage']}%")

def main():
    parser = argparse.ArgumentParser(description="Compute system match coverage between two TSV files.")
    parser.add_argument("file1", help="Path to the first TSV file")
    parser.add_argument("file2", help="Path to the second TSV file")
    args = parser.parse_args()

    systems1 = read_systems(args.file1)
    systems2 = read_systems(args.file2)

    cov1 = compute_best_coverage(systems1, systems2)
    cov2 = compute_best_coverage(systems2, systems1)

    print_coverage_report(args.file1, cov1)
    print_coverage_report(args.file2, cov2)

if __name__ == "__main__":
    main()

