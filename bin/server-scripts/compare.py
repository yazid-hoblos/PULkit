import pandas as pd
import sys

def gf_to_set(gf_str):
    if pd.isna(gf_str):
        return set()
    return set(item.strip() for item in gf_str.split(','))

def jaccard_similarity(set1, set2):
    # if not set1 and not set2:
    #     return 1.0
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    return intersection / union if union > 0 else 0

def find_best_matches(file1, file2):
    df1 = pd.read_csv(file1, sep='\t')
    df2 = pd.read_csv(file2, sep='\t')

    df1 = df1[['system number', 'model_GF', 'context_GF']].rename(columns={'system number': 'system_id_1'})
    df2 = df2[['system number', 'model_GF', 'context_GF']].rename(columns={'system number': 'system_id_2'})

    df1['model_GF_set'] = df1['model_GF'].apply(gf_to_set)
    df1['context_GF_set'] = df1['context_GF'].apply(gf_to_set)
    df2['model_GF_set'] = df2['model_GF'].apply(gf_to_set)
    df2['context_GF_set'] = df2['context_GF'].apply(gf_to_set)

    results = []

    for _, row1 in df1.iterrows():
        best_match = None
        best_score = -1

        for _, row2 in df2.iterrows():
            # Calculate similarity (you can change this to context_GF_jaccard or average)
            # score = jaccard_similarity(row1['model_GF_set'], row2['model_GF_set'])
            score = jaccard_similarity(row1['context_GF_set'], row2['context_GF_set'])
            # score = (model_score + context_score) / 2

            if score > best_score:
                best_score = score
                best_match = row2['system_id_2']

        results.append({
            'system_id_1': row1['system_id_1'],
            'best_match_system_id_2': best_match,
            'best_model_GF_jaccard': best_score
        })

    results_df = pd.DataFrame(results)
    print(results_df)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python match_best_systems.py <file1.tsv> <file2.tsv>")
        sys.exit(1)

    find_best_matches(sys.argv[1], sys.argv[2])
