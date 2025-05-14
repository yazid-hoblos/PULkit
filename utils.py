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
                
                if caz.startswith(('GH', 'GT', 'CE', 'PL', 'CBM')):
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

