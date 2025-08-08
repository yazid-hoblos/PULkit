'''Identifies best system matches based on model GF overlaps between two pangenome systems files.'''

import pandas as pd
import sys

def gf_to_set(gf_str):
    if pd.isna(gf_str):
        return set()
    return set(item.strip() for item in gf_str.split(','))

def jaccard_similarity(set1, set2):
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    return intersection / union if union > 0 else 0

def find_best_matches(file1, file2, output_file):
    df1 = pd.read_csv(file1, sep='\t')
    df2 = pd.read_csv(file2, sep='\t')

    df1 = df1[['system number', 'model_GF', 'context_GF']].rename(columns={'system number': 'system_id_1'})
    df2 = df2[['system number', 'model_GF', 'context_GF']].rename(columns={'system number': 'system_id_2'})

    df1['model_GF_set'] = df1['model_GF'].apply(gf_to_set)
    df1['context_GF_set'] = df1['context_GF'].apply(gf_to_set)
    df2['model_GF_set'] = df2['model_GF'].apply(gf_to_set)
    df2['context_GF_set'] = df2['context_GF'].apply(gf_to_set)

    results = []
    unmatched_systems = []
    
    for _, row1 in df1.iterrows():
        best_match = ''
        best_score = context_best_score = 0

        for _, row2 in df2.iterrows():
            # Calculate similarity (both model and context)
            model_score = jaccard_similarity(row1['model_GF_set'], row2['model_GF_set'])

            if model_score > best_score:
                best_score = model_score
                best_match = row2['system_id_2']
                context_best_score = jaccard_similarity(row1['context_GF_set'], row2['context_GF_set'])

        if best_match:
            results.append({
            'system_id': row1['system_id_1'],
            'best_match': str(best_match),
            'overlap_jaccard': round(best_score, 2),
            'context_overlap_jaccard': round(context_best_score, 2)}) 
        else:
            unmatched_systems.append(row1['system_id_1'])

    results.sort(key=lambda x: (x['overlap_jaccard'], x['context_overlap_jaccard']), reverse=True)
    results_df = pd.DataFrame(results)
    results_df.to_csv(output_file, index=False, sep='\t')
    print(results_df.to_string(index=False))
    print(f'\nUnmatched systems: {unmatched_systems}\n')
    print(f"Results saved to {output_file} ({len(results_df)} rows)")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python match_best_systems.py <file1.tsv> <file2.tsv> [output_file.tsv]")
        sys.exit(1)

    output_file = sys.argv[3] if len(sys.argv) > 3 else 'comparison_results.csv'
    find_best_matches(sys.argv[1], sys.argv[2], output_file)
