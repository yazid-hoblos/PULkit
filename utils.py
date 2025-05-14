import pandas as pd

# Parse the dbCAN-PUL Database and substrate mapping file
def parse_dbcan_pul_data(pul_file, substrate_mapping_file):
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