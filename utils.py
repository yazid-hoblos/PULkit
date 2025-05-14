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
        
        if pd.notnull(pul['cazymes_predicted_dbcan']):
            cazymes_str = pul['cazymes_predicted_dbcan']
        else:
            cazymes_str = pul['cazymes_predicted_dbCAN2'] if pd.notnull(pul['cazymes_predicted_dbCAN2']) else ""
        
        # Parse CAZymes
        cazymes = []
        if cazymes_str and isinstance(cazymes_str, str):
            for caz in cazymes_str.split(','):
                if '|' in caz:
                    # For cases like "GH30|GH30_8", take the first part
                    caz = caz.split('|')[0].strip()
                
                if '_' in caz:
                    caz = caz.split('_')[0].strip()
                
                if caz.startswith(('GH', 'GT', 'CE', 'PL', 'AA', 'CBM')):
                    cazymes.append(caz)
        
        gene_count = pul['num_genes'] if pd.notnull(pul['num_genes']) else 0
        cazyme_count = pul['num_cazymes'] if pd.notnull(pul['num_cazymes']) else 0
        
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
    
    return pd.DataFrame(pul_features)

puls = parse_dbcan_pul_data('dbCAN-PUL_v5.csv', 'dbCAN-PUL.substrate.mapping.csv')
pul_features = extract_pul_features(puls)

print(f"\nTotal PULs processed: {len(pul_features)}")
print(f"Unique substrates: {pul_features['Substrate'].nunique()}")
