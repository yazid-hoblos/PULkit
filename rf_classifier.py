import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import classification_report
import pandas as pd

def prepare_ml_classifier(pul_features, substrate_groups):
    """Prepare and train a machine learning classifier for PUL substrate prediction"""

    valid_substrates = list(substrate_groups.keys())
    valid_puls = pul_features[pul_features['Substrate_Clean'].isin(valid_substrates)].copy()
    
    # Encode substrates
    le = LabelEncoder()
    valid_puls['Substrate_Encoded'] = le.fit_transform(valid_puls['Substrate_Clean'])
    
    #  Feature matrix
    all_cazyme_families = set()
    for cazyme_list in valid_puls['CAZymes']:
        all_cazyme_families.update(cazyme_list)
    
    # One-hot encoding matrix for CAZymes
    X = np.zeros((len(valid_puls), len(all_cazyme_families)))
    family_to_idx = {family: i for i, family in enumerate(sorted(all_cazyme_families))}
        
    for i, cazyme_list in enumerate(valid_puls['CAZymes']):
        for family in cazyme_list:
            if family in family_to_idx:  
                X[i, family_to_idx[family]] = 1
    
    # Add other features
    if 'CAZyme_Count' in valid_puls and 'Gene_Count' in valid_puls:
        additional_features = np.array(valid_puls[['CAZyme_Count', 'Gene_Count']])
        X = np.hstack((X, additional_features))
    
    y = valid_puls['Substrate_Encoded'].values
    
    # Split into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    
    # Train
    clf = RandomForestClassifier(n_estimators=100, random_state=42)
    clf.fit(X_train, y_train)
    
    # Evaluate
    y_pred = clf.predict(X_test)
    print("Classification Report:")
    print(classification_report(y_test, y_pred, target_names=le.inverse_transform(np.unique(y))))
    
    # check feature importances
    if len(clf.feature_importances_) >= len(family_to_idx): # remove additional features (i.e. counts)
        importances = clf.feature_importances_[:len(family_to_idx)]
        indices = np.argsort(importances)[::-1]
        
        print("\nMost important CAZyme families:")
        for i in range(min(10, len(indices))):
            family = list(family_to_idx.keys())[indices[i]]
            print(f"{family}: {importances[indices[i]]:.4f}")
    
    return clf, le, family_to_idx, X_test, y_test


def predict_pul_type_with_rf(clf, le, family_to_idx, pul_features):
    """Predict PUL type using the trained Random Forest classifier"""
    
    all_cazyme_families = set(family_to_idx.keys())
    X = np.zeros((len(pul_features), len(all_cazyme_families)))
    
    for i, cazyme_list in enumerate(pul_features['CAZymes']):
        for family in cazyme_list:
            if family in family_to_idx:
                X[i, family_to_idx[family]] = 1
    
    # Add other features
    if 'CAZyme_Count' in pul_features and 'Gene_Count' in pul_features:
        additional_features = np.array(pul_features[['CAZyme_Count', 'Gene_Count']])
        X = np.hstack((X, additional_features))
    
    # Predict
    y_pred = clf.predict(X)
    
    # Decode predictions
    pul_features['Predicted_Type'] = le.inverse_transform(y_pred)
    pul_features['True_Type'] = pul_features['Substrate_Clean'].apply(lambda x: x if x in le.classes_ else 'Unknown')
    
    return pul_features


def evaluate_classification(predictions):
    """Evaluate the classification performance"""
    # Filter out unknowns 
    valid_predictions = predictions
    #predictions[(predictions['Predicted_Type'] != 'Unknown') & 
    #                              (predictions['True_Type'] != 'Unknown')] 
    
    # print the count of unknowns
    unknown_count = len(predictions[predictions['Predicted_Type'] == 'Unknown'])
    print(f"Number of unknown predictions: {unknown_count}")
    
    # overall accuracy
    correct = (valid_predictions['Predicted_Type'] == valid_predictions['True_Type']).sum()
    total = len(valid_predictions)
    accuracy = correct / total if total > 0 else 0
    
    print(f"Overall accuracy: {accuracy:.4f} ({correct}/{total})")
    
    # per-substrate metrics
    results = {}
    for substrate in valid_predictions['True_Type'].unique():
        substrate_puls = valid_predictions[valid_predictions['True_Type'] == substrate]
        substrate_correct = (substrate_puls['Predicted_Type'] == substrate_puls['True_Type']).sum()
        substrate_total = len(substrate_puls)
        substrate_accuracy = substrate_correct / substrate_total if substrate_total > 0 else 0
        
        results[substrate] = {
            'accuracy': substrate_accuracy,
            'correct': substrate_correct,
            'total': substrate_total
        }
        
        print(f"Substrate {substrate}: {substrate_accuracy:.4f} ({substrate_correct}/{substrate_total})")
    
    # confusion matrix
    confusion = pd.crosstab(
        valid_predictions['True_Type'], 
        valid_predictions['Predicted_Type'], 
        rownames=['True'], 
        colnames=['Predicted']
    )
    
    print("\nConfusion Matrix:")
    print(confusion)
    
    return results, confusion