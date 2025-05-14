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