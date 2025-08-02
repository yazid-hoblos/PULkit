import pandas as pd
from collections import Counter


def parse_dbcan_pul_data(pul_file, substrate_mapping_file):
    "Parse the dbCAN-PUL Database and substrate mapping file"
    puls = pd.read_csv(pul_file)
    substrate_mapping = pd.read_csv(substrate_mapping_file, skiprows=1)    
    puls['PULID'] = puls['PULID'].astype(str)
    substrate_mapping['PUL ID'] = substrate_mapping['PUL ID'].astype(str)
    
    merged = pd.merge(
        puls, 
        substrate_mapping[['PUL ID', 'updated_substrate (04/14/2023)']], 
        left_on='PULID', 
        right_on='PUL ID', 
        how='left')
    # Use updated substrate where available, fall back to substrate_final
    merged['substrate'] = merged['updated_substrate (04/14/2023)'].fillna(merged['substrate_final'])
    return merged

def extract_pul_features(puls):
    """Extract features from each PUL in the database."""
    pul_features = []
    
    for idx, pul in puls.iterrows():
        pul_id = pul['PULID']
        organism = pul['organism_name']
        substrate = pul['substrate']
        
        cazymes_str = pul['cazymes_predicted_dbcan'] if pd.notnull(pul['cazymes_predicted_dbcan']) else ""

        # Parse CAZymes
        cazymes = []
        if cazymes_str and isinstance(cazymes_str, str):
            for caz in cazymes_str.split(','):
                if '|' in caz:
                    # For cases like "GH30|GH30_8", take the first part
                    caz = caz.split('|')[0].strip() 
                if '_' in caz:
                    caz = caz.split('_')[0].strip()
                
                if caz.startswith(('GH', 'GT', 'CE', 'PL', 'CBM')): #AA
                    cazymes.append(caz)
        
        gene_count = int(pul['num_genes']) if pd.notnull(pul['num_genes']) else 0
        cazyme_count = int(pul['num_cazymes']) if pd.notnull(pul['num_cazymes']) else 0
        
        # Extract genomic range to calculate PUL size
        size = 0
        if pd.notnull(pul['nucleotide_position_range']):
            try:
                start, end = map(int, pul['nucleotide_position_range'].split('-'))
                size = end - start
            except:
                size = 0
        
        feature_dict = {
            'PUL_ID': pul_id,
            'Organism': organism,
            'Substrate': substrate,
            'CAZymes': cazymes,
            'CAZyme_Count': cazyme_count,
            'Gene_Count': gene_count,
            'Size': size,
            'Verification': pul['verification_final'] if pd.notnull(pul['verification_final']) else ""
        }
        
        cazyme_counts = Counter(cazymes)
        for family, count in cazyme_counts.items():
            feature_dict[f'CAZyme_{family}'] = count
        
        pul_features.append(feature_dict)
        pul_features_db = pd.DataFrame(pul_features)
        
        # pul_features_db.to_csv('pul_features.csv', index=True)
        
    return pul_features_db


def group_puls_by_substrate(pul_features, min_group_size=5):
    """Group PULs by substrate"""
    substrate_groups = {}
    pul_features['Substrate_Clean'] = pul_features['Substrate'].str.lower().str.strip()
    
    substrate_counts = pul_features['Substrate_Clean'].value_counts()
    
    # Only consider substrates with enough examples
    valid_substrates = substrate_counts[substrate_counts >= min_group_size].index
    
    for substrate in valid_substrates:
        group = pul_features[pul_features['Substrate_Clean'] == substrate]
        
        # Get common CAZyme families associated with the substrate
        all_cazymes = []
        for cazyme_list in group['CAZymes']:
            all_cazymes.extend(cazyme_list)
        
        common_cazymes = Counter(all_cazymes).most_common() 
        
        # Store group information
        substrate_groups[substrate] = {
            'count': len(group),
            'pul_ids': group['PUL_ID'].tolist(),
            'common_cazymes': common_cazymes
        }
        
        print(f"Substrate: {substrate}, Count: {len(group)}, PUL IDs: {', '.join(group['PUL_ID'].tolist())}") 
        print(f"Common CAZymes: {common_cazymes[:5]}")
        print("-" * 50)
    
    return substrate_groups


def identify_pul_signatures(pul_features, substrate_groups, t=0.5):
    """Identify signature patterns for each substrate group"""
    '''t: threshold for min percentage of occurence of CAZymes in a substrate group'''
    signatures = {}
    
    for substrate, group_info in substrate_groups.items():
        pul_ids = group_info['pul_ids']
        substrate_puls = pul_features[pul_features['PUL_ID'].isin(pul_ids)]
        
        cazyme_combinations = []
        for _, pul in substrate_puls.iterrows():
            cazymes = pul['CAZymes']

            if len(cazymes) < 2:
                continue
                
            for i in range(len(cazymes)):
                for j in range(i+1, len(cazymes)):
                    combination = f"{cazymes[i]}+{cazymes[j]}"
                    cazyme_combinations.append(combination)
        
        common_combinations = Counter(cazyme_combinations)
        
        threshold = len(substrate_puls) * t
        
        # Filter combinations that meet the threshold
        significant_combinations = [combo for combo, count in common_combinations.items() 
                                  if count >= threshold]
        
        # Get individual significant CAZymes
        significant_cazymes = [caz for caz, count in group_info['common_cazymes'] if count >= threshold]
        
        signatures[substrate] = {
            'cazyme_signatures': significant_cazymes,
            'combination_signatures': significant_combinations}
        
        print(f"Substrate: {substrate}")
        print(f"Signature CAZymes: {', '.join(significant_cazymes[:5])}")
        print(f"Signature Combinations: {', '.join(significant_combinations[:5])}")
        print("-" * 50)
    
    return signatures


# def compare_results(rf_results, rule_results):
#     """Visualize the evaluation results of the classifiers"""
#     import matplotlib.pyplot as plt
#     import seaborn as sns
    
#     rf_df = pd.DataFrame(rf_results).T.reset_index()
#     rf_df.columns = ['Substrate', 'RF_Accuracy', 'RF_Correct', 'RF_Total']
    
#     rule_df = pd.DataFrame(rule_results).T.reset_index()
#     rule_df.columns = ['Substrate', 'Rule_Accuracy', 'Rule_Correct', 'Rule_Total']
    
#     merged_df = pd.merge(rf_df, rule_df, on='Substrate')
    
#     plt.figure(figsize=(12, 6))
    
#     sns.barplot(data=merged_df, x='Substrate', y='RF_Accuracy', color='blue', label='RF Accuracy', alpha=0.6)
#     sns.barplot(data=merged_df, x='Substrate', y='Rule_Accuracy', color='red', label='Rule Accuracy', alpha=0.7)
    
#     plt.xticks(rotation=45)
#     plt.ylabel('Accuracy')
#     plt.title('Comparison of RF and Rule-based Classifier Accuracies')
#     plt.legend()
    
#     plt.tight_layout()
#     plt.show()


def compare_predictions(rf_results, rule_results):
    """Visualize and compare the predictions of the rule-based and RF classifiers"""
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    import numpy as np
    
    rf_df = pd.DataFrame(rf_results).T.reset_index()
    rf_df.columns = ['Substrate', 'RF_Accuracy', 'RF_Correct', 'RF_Total']
    
    rule_df = pd.DataFrame(rule_results).T.reset_index()
    rule_df.columns = ['Substrate', 'Rule_Accuracy', 'Rule_Correct', 'Rule_Total']
    
    merged_df = pd.merge(rf_df, rule_df, on='Substrate')
    
    merged_df['Accuracy_Diff'] = merged_df['RF_Accuracy'] - merged_df['Rule_Accuracy']
    merged_df['Better_Method'] = merged_df['Accuracy_Diff'].apply(
        lambda x: 'Random Forest' if x > 0 else ('Rule-based' if x < 0 else 'Equal'))
    
    merged_df = merged_df.sort_values('Substrate')
    
    plt.figure(figsize=(14, 6))
    gs = GridSpec(1, 2)
    
    # 1. Side-by-side bar chart 
    ax1 = plt.subplot(gs[0, 0])
    substrates = merged_df['Substrate']
    x = np.arange(len(substrates))
    width = 0.35
    
    ax1.bar(x - width/2, merged_df['RF_Accuracy'], width, label='Random Forest', color='#4C72B0')
    ax1.bar(x + width/2, merged_df['Rule_Accuracy'], width, label='Rule-based', color='#C44E52')
    
    for i, v in enumerate(merged_df['RF_Accuracy']):
        ax1.text(i - width/2, v + 0.01, f'{v:.2f}', ha='center', va='bottom', fontsize=8)
    
    for i, v in enumerate(merged_df['Rule_Accuracy']):
        ax1.text(i + width/2, v + 0.01, f'{v:.2f}', ha='center', va='bottom', fontsize=8)
    
    ax1.set_xticks(x)
    ax1.set_xticklabels(substrates, rotation=45, ha='right')
    ax1.set_ylim(0, 1.1)
    ax1.set_ylabel('Accuracy')
    ax1.set_title('Side-by-Side Comparison of Classifier Accuracies')
    ax1.legend()
    ax1.grid(axis='y', linestyle='--', alpha=0.7)
    
    # 2. Accuracy difference plot
    ax2 = plt.subplot(gs[0, 1])
    bars = ax2.bar(x, merged_df['Accuracy_Diff'], color=[
        '#4C72B0' if diff > 0 else '#C44E52' for diff in merged_df['Accuracy_Diff']
    ])
    
    ax2.axhline(y=0, color='black', linestyle='-', alpha=0.3)
    
    for i, bar in enumerate(bars):
        height = bar.get_height()
        if height < 0:
            va = 'top'
            offset = -0.01
        else:
            va = 'bottom'
            offset = 0.01
        ax2.text(bar.get_x() + bar.get_width()/2, height + offset,
                f'{merged_df["Accuracy_Diff"].iloc[i]:.2f}',
                ha='center', va=va, fontsize=8)
    
    ax2.set_xticks(x)
    ax2.set_xticklabels(substrates, rotation=45, ha='right')
    ax2.set_ylabel('RF Accuracy - Rule Accuracy')
    ax2.set_title('Performance Difference (RF - Rule)')
    ax2.grid(axis='y', linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.4)
    plt.show()
    # plt.savefig('classifiers_comparison.png', dpi=500)
    

def compare_confusion_matrices(rf_results, rule_results):
    import matplotlib.pyplot as plt
    
    rf_df = pd.DataFrame(rf_results).T.reset_index()
    rf_df.columns = ['Substrate', 'RF_Accuracy', 'RF_Correct', 'RF_Total']
    
    rule_df = pd.DataFrame(rule_results).T.reset_index()
    rule_df.columns = ['Substrate', 'Rule_Accuracy', 'Rule_Correct', 'Rule_Total']
    
    merged_df = pd.merge(rf_df, rule_df, on='Substrate')
    
    merged_df['Accuracy_Diff'] = merged_df['RF_Accuracy'] - merged_df['Rule_Accuracy']
    merged_df['Better_Method'] = merged_df['Accuracy_Diff'].apply(
        lambda x: 'Random Forest' if x > 0 else ('Rule-based' if x < 0 else 'Equal'))
    
    merged_df = merged_df.sort_values('Substrate')
    
    plt.figure(figsize=(12, 6))
    
    table_data = merged_df[['Substrate', 'RF_Accuracy', 'Rule_Accuracy', 
                           'RF_Correct', 'RF_Total', 'Rule_Correct', 'Rule_Total', 'Better_Method']]
    
    cell_text = []
    for i, row in table_data.iterrows():
        cell_text.append([
            row['Substrate'],
            f"{row['RF_Accuracy']:.3f}",
            f"{row['Rule_Accuracy']:.3f}",
            f"{int(row['RF_Correct'])}/{int(row['RF_Total'])}",
            f"{int(row['Rule_Correct'])}/{int(row['Rule_Total'])}",
            row['Better_Method']])
    
    column_labels = ['Substrate', 'RF Accuracy', 'Rule Accuracy', 
                     'RF (Correct/Total)', 'Rule (Correct/Total)', 'Better Method']
    
    table = plt.table(cellText=cell_text, colLabels=column_labels, 
                      loc='center', cellLoc='center')
    
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.5)
    
    for i, row in enumerate(cell_text):
        if row[5] == 'Random Forest':
            table[(i+1, 5)].set_facecolor('#DCE5F5')  # Light blue
        elif row[5] == 'Rule-based':
            table[(i+1, 5)].set_facecolor('#F5DCDC')  # Light red
    
    for i in range(len(column_labels)):
        table[(0, i)].set_facecolor('#F0F0F0')
        table[(0, i)].set_text_props(weight='bold')
    
    plt.axis('off')
    plt.show()
    # plt.savefig('confusion_matrices_comparison.png', dpi=900)