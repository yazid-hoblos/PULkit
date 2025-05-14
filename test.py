from utils import *
from PUL_rules import *
from rf_classifier import *

puls = parse_dbcan_pul_data('../dbCAN-PUL_v5.csv', '../dbCAN-PUL.substrate.mapping.csv')
pul_features = extract_pul_features(puls)

print(f"\nTotal PULs processed: {len(pul_features)}")
print(f"Unique substrates: {pul_features['Substrate'].nunique()}")

substrate_groups = group_puls_by_substrate(pul_features, min_group_size=10)

signatures = identify_pul_signatures(pul_features, substrate_groups)

rule_sets = create_pul_rule_sets(signatures, pul_features, substrate_groups)

print("\nTraining machine learning classifier...")
clf, le, family_to_idx, X_test, y_test = prepare_ml_classifier(pul_features, substrate_groups)

print("\nPredicting PUL types using rule-based approach...")
rule_predictions = predict_pul_type_with_rules(pul_features, rule_sets)
print(rule_predictions)



    
