import pandas as pd

def create_pul_rule_sets(signatures, pul_features, substrate_groups):
    """Create rule sets for each PUL type based on signatures."""
    rule_sets = {}
    
    for substrate, sig in signatures.items():
        pul_ids = substrate_groups[substrate]['pul_ids']
        sub_puls = pul_features[pul_features['PUL_ID'].isin(pul_ids)]
        
        # get metrics for this substrate group
        avg_cazyme_count = sub_puls['CAZyme_Count'].mean()
        median_size = sub_puls['Size'].median()
        
        rule_sets[substrate] = {
            'required_cazymes': sig['cazyme_signatures'][:3] if sig['cazyme_signatures'] else [],
            'cazyme_combinations': sig['combination_signatures'][:3] if sig['combination_signatures'] else [],
            'min_cazymes': max(1, int(avg_cazyme_count * 0.5)),  # At least 50% of average
            'typical_size': median_size}
        
        print(f"Rule set for {substrate}:")
        print(f"Required CAZymes (any of): {', '.join(rule_sets[substrate]['required_cazymes'])}")
        print(f"Signature combinations (any of): {', '.join(rule_sets[substrate]['cazyme_combinations'])}")
        print(f"Minimum CAZyme count: {rule_sets[substrate]['min_cazymes']}")
        print(f"Typical size: {rule_sets[substrate]['typical_size']}")
        print("-" * 50)
    
    return rule_sets


def predict_pul_type_with_rules(pul_features, rule_sets):
    """Predict PUL type based on substrate-based rule models"""
    
    predictions = []
    
    for _, pul in pul_features.iterrows():
        pul_id = pul['PUL_ID']
        cazymes = pul['CAZymes']
        substrate_clean = pul['Substrate_Clean'] if 'Substrate_Clean' in pul else pul['Substrate'].lower().strip()
        
        scores = {}
        for substrate, rules in rule_sets.items():
            score = 0
            
            for req_caz in rules['required_cazymes']:
                if req_caz in cazymes:
                    score += 2
            
            for combo in rules['cazyme_combinations']:
                caz1, caz2 = combo.split('+')
                if caz1 in cazymes and caz2 in cazymes:
                    score += 3
            
            if len(cazymes) >= rules['min_cazymes']:
                score += 1
            
            scores[substrate] = score
        
        if scores:
            best_match = max(scores.items(), key=lambda x: x[1])
            if best_match[1] > 2:  # minimum threshold
                predictions.append({
                    'PUL_ID': pul_id,
                    'Predicted_Type': best_match[0],
                    'Confidence_Score': best_match[1],
                    'True_Type': substrate_clean
                })
            else:
                predictions.append({
                    'PUL_ID': pul_id,
                    'Predicted_Type': 'Unknown',
                    'Confidence_Score': 0,
                    'True_Type': substrate_clean
                })
        else:
            predictions.append({
                'PUL_ID': pul_id,
                'Predicted_Type': 'Unknown',
                'Confidence_Score': 0,
                'True_Type': substrate_clean
            })
    
    return pd.DataFrame(predictions)
