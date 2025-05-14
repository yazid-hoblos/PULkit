from utils import *

puls = parse_dbcan_pul_data('../dbCAN-PUL_v5.csv', '../dbCAN-PUL.substrate.mapping.csv')
pul_features = extract_pul_features(puls)

print(f"\nTotal PULs processed: {len(pul_features)}")
print(f"Unique substrates: {pul_features['Substrate'].nunique()}")
