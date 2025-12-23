import numpy as np
from typing import List, Dict, Any
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, mean_squared_error

# --- Multi-Modal Fusion & Prediction Heads ---
def prepare_model_inputs(features: List[Dict[str, Any]]):
    # Flatten fingerprints and descriptors, concatenate context features
    X = []
    y_synergy = []
    y_effect_size = []
    for feat in features:
        # Example: concatenate fp1, fp2, desc1, desc2, synergy_score, benefit_risk_ratio
        x_vec = np.concatenate([
            feat['fp1'], feat['fp2'],
            np.array(list(feat['desc1'].values())),
            np.array(list(feat['desc2'].values())),
            np.array([
                feat.get('synergy_score', 0.0),
                feat.get('ai_predicted', 0.0),
                feat.get('benefit_risk_ratio', 1.0)
            ])
        ])
        X.append(x_vec)
        # Synergy/antagonism classification (1=synergy, 0=antagonism/additive/unknown)
        y_synergy.append(1 if feat.get('interaction_effect', '').lower() == 'synergistic' else 0)
        # Effect size regression (map qualitative to numeric)
        effect_map = {'large': 1.0, 'medium': 0.7, 'small': 0.4, 'additive': 0.5, 'antagonistic': 0.0, 'unknown': 0.5}
        y_effect_size.append(effect_map.get(feat.get('effect_size', '').lower(), 0.5))
    return np.array(X), np.array(y_synergy), np.array(y_effect_size)

def train_and_benchmark_nextgen_model(X, y_synergy, y_effect_size):
    # Split for classification (synergy) and regression (effect size)
    X_train, X_test, y_train, y_test = train_test_split(X, y_synergy, test_size=0.3, random_state=42)
    clf = RandomForestClassifier(n_estimators=100, random_state=42)
    clf.fit(X_train, y_train)
    y_pred = clf.predict(X_test)
    print("Synergy/Antagonism Classification Report:")
    print(classification_report(y_test, y_pred))

    # Regression for effect size
    X_train_r, X_test_r, y_train_r, y_test_r = train_test_split(X, y_effect_size, test_size=0.3, random_state=42)
    reg = RandomForestRegressor(n_estimators=100, random_state=42)
    reg.fit(X_train_r, y_train_r)
    y_pred_r = reg.predict(X_test_r)
    print("Effect Size Regression MSE:", mean_squared_error(y_test_r, y_pred_r))

# Example usage:
# features = extract_features_for_model(...)
# X, y_synergy, y_effect_size = prepare_model_inputs(features)
# train_and_benchmark_nextgen_model(X, y_synergy, y_effect_size)
