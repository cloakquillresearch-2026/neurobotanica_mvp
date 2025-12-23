#!/usr/bin/env python3
"""
Reformulate TherapeuticPredictionModel as Binary Classifier
Predict: Is compound effective for condition? (Yes/No)

Instead of continuous efficacy scores (0-1), use binary labels:
- 1 (Effective): effect_size = large or medium
- 0 (Not effective): effect_size = small, minimal, or None
"""
import json
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List
from collections import Counter

# Add backend to path
SCRIPT_DIR = Path(__file__).parent
PROJECT_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(PROJECT_ROOT / 'backend'))

from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import classification_report, roc_auc_score, confusion_matrix
import joblib


def load_enriched_dataset() -> Dict:
    """Load fixed dataset."""
    data_file = PROJECT_ROOT / 'data' / 'training' / 'neurobotanica_fully_enriched_fixed.json'
    
    with open(data_file, 'r', encoding='utf-8') as f:
        return json.load(f)


def prepare_binary_classification_data(dataset: Dict) -> tuple:
    """Prepare data for binary classification.
    
    Returns:
        (X_features, y_binary, sample_weights, metadata)
    """
    compounds = dataset['compounds']
    
    X_list = []
    y_list = []
    weights_list = []
    metadata_list = []
    
    stats = {
        'total_samples': 0,
        'positive_samples': 0,
        'negative_samples': 0,
        'effect_size_distribution': Counter()
    }
    
    for compound in compounds:
        descriptors = compound.get('rdkit_descriptors')
        if not descriptors or 'calculation_note' in descriptors:
            continue
        
        clinical_studies = compound.get('clinical_studies', [])
        num_studies = len(clinical_studies)
        
        base_confidence = min(1.0, 0.3 + (num_studies * 0.05))
        
        # Extract features
        features = {
            'molecular_weight': descriptors.get('molecular_weight', 0),
            'logP': descriptors.get('logP', 0),
            'h_bond_donors': descriptors.get('h_bond_donors', 0),
            'h_bond_acceptors': descriptors.get('h_bond_acceptors', 0),
            'rotatable_bonds': descriptors.get('rotatable_bonds', 0),
            'tpsa': descriptors.get('tpsa', 0),
            'num_aromatic_rings': descriptors.get('num_aromatic_rings', 0),
            'num_stereocenters': descriptors.get('num_stereocenters', 0),
            'sp3_fraction': descriptors.get('sp3_fraction', 0),
            'complexity': descriptors.get('complexity', 0),
            'qed': descriptors.get('qed', 0.5),
            'bbb_penetration': 1 if descriptors.get('bbb_penetration_predicted') else 0,
            'has_phenol': 1 if descriptors.get('has_phenol') else 0,
            'num_terpene_synergies': len(compound.get('terpene_synergies', [])),
            'avg_synergy_evidence_tier': _avg_synergy_tier(compound.get('terpene_synergies', []))
        }
        
        therapeutic_targets = compound.get('therapeutic_targets', [])
        
        if not therapeutic_targets:
            continue
        
        for target in therapeutic_targets:
            # Binary classification: is compound effective for this target?
            is_effective, effect_size = _determine_binary_efficacy(compound, target, clinical_studies)
            
            # Confidence based on evidence
            target_normalized = target.upper().replace(' ', '_').replace('-', '_')
            has_direct_evidence = any(
                target_normalized in s.get('condition', '').upper().replace(' ', '_').replace('-', '_') or
                s.get('condition', '').upper().replace(' ', '_').replace('-', '_') in target_normalized
                for s in clinical_studies
            )
            
            confidence = base_confidence
            if has_direct_evidence:
                confidence = min(1.0, confidence + 0.4)  # Higher boost for binary task
            else:
                # For samples without direct evidence, use lower confidence
                confidence = max(0.1, confidence - 0.2)
            
            features_with_target = features.copy()
            features_with_target['target_encoded'] = hash(target) % 100
            
            X_list.append(features_with_target)
            y_list.append(is_effective)
            weights_list.append(confidence)
            metadata_list.append({
                'compound': compound['compound_name'],
                'target': target,
                'effect_size': effect_size,
                'has_evidence': has_direct_evidence
            })
            
            # Update stats
            stats['total_samples'] += 1
            if is_effective:
                stats['positive_samples'] += 1
            else:
                stats['negative_samples'] += 1
            stats['effect_size_distribution'][effect_size] += 1
    
    X = pd.DataFrame(X_list)
    y = pd.Series(y_list, dtype=int)
    weights = pd.Series(weights_list)
    
    return X, y, weights, metadata_list, stats


def _avg_synergy_tier(synergies: List[Dict]) -> float:
    if not synergies:
        return 0.0
    return np.mean([s.get('evidence_tier', 3) for s in synergies])


def _determine_binary_efficacy(compound: Dict, target: str, studies: List[Dict]) -> tuple:
    """Determine if compound is effective for target (binary).
    
    Returns:
        (is_effective: bool, effect_size: str)
    """
    target_normalized = target.upper().replace(' ', '_').replace('-', '_')
    
    # Find relevant studies
    relevant_studies = []
    for s in studies:
        condition = s.get('condition', '').upper().replace(' ', '_').replace('-', '_')
        if (target_normalized in condition or 
            condition in target_normalized or
            target_normalized.replace('_', '') in condition.replace('_', '')):
            relevant_studies.append(s)
    
    if not relevant_studies:
        # No direct evidence - assume not effective (conservative)
        return False, 'none'
    
    # Get effect sizes from studies
    effect_sizes = [s.get('effect_size') for s in relevant_studies]
    
    # Count by category
    large_count = effect_sizes.count('large')
    medium_count = effect_sizes.count('medium')
    small_count = effect_sizes.count('small')
    minimal_count = effect_sizes.count('minimal')
    none_count = effect_sizes.count(None)
    
    # Decision rule: effective if majority are large/medium
    effective_count = large_count + medium_count
    ineffective_count = small_count + minimal_count + none_count
    
    if effective_count > ineffective_count:
        # Determine primary effect size
        if large_count >= medium_count:
            return True, 'large'
        else:
            return True, 'medium'
    else:
        # Not effective
        if small_count > 0:
            return False, 'small'
        elif minimal_count > 0:
            return False, 'minimal'
        else:
            return False, 'none'


def train_binary_classifier():
    """Train binary classification model."""
    print("=" * 70)
    print("Binary Classification: Compound Efficacy Prediction")
    print("=" * 70)
    print()
    
    # Load dataset
    print("Loading enriched dataset...")
    dataset = load_enriched_dataset()
    print(f"✅ Loaded {len(dataset['compounds'])} compounds")
    print()
    
    # Prepare binary data
    print("Preparing binary classification data...")
    X, y, weights, metadata, stats = prepare_binary_classification_data(dataset)
    
    print(f"Dataset Statistics:")
    print(f"  Total samples: {stats['total_samples']}")
    print(f"  Positive (Effective): {stats['positive_samples']} ({stats['positive_samples']/stats['total_samples']*100:.1f}%)")
    print(f"  Negative (Not Effective): {stats['negative_samples']} ({stats['negative_samples']/stats['total_samples']*100:.1f}%)")
    print()
    print(f"  Effect size distribution:")
    for effect_size, count in stats['effect_size_distribution'].most_common():
        print(f"    {effect_size or 'none'}: {count}")
    print()
    
    # Check class balance
    class_ratio = stats['positive_samples'] / stats['negative_samples']
    if class_ratio < 0.5 or class_ratio > 2.0:
        print(f"⚠️  Class imbalance detected (ratio: {class_ratio:.2f})")
        print(f"   Using class_weight='balanced' in models")
        use_balanced = True
    else:
        print(f"✅ Reasonable class balance (ratio: {class_ratio:.2f})")
        use_balanced = False
    print()
    
    # Split data
    X_train, X_test, y_train, y_test, w_train, w_test = train_test_split(
        X, y, weights, test_size=0.2, random_state=42, stratify=y
    )
    
    # Scale features
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    print(f"Training set: {len(X_train)} samples")
    print(f"Test set: {len(X_test)} samples")
    print()
    
    # === Train Gradient Boosting Classifier ===
    print("=" * 70)
    print("1. Gradient Boosting Classifier")
    print("=" * 70)
    
    gb_model = GradientBoostingClassifier(
        n_estimators=100,
        learning_rate=0.1,
        max_depth=3,
        min_samples_split=10,
        min_samples_leaf=5,
        subsample=0.7,
        random_state=42
    )
    
    print("Training...")
    gb_model.fit(X_train_scaled, y_train, sample_weight=w_train)
    
    # Evaluate
    train_score = gb_model.score(X_train_scaled, y_train)
    test_score = gb_model.score(X_test_scaled, y_test)
    
    y_pred_gb = gb_model.predict(X_test_scaled)
    y_proba_gb = gb_model.predict_proba(X_test_scaled)[:, 1]
    
    print(f"Train Accuracy: {train_score:.4f}")
    print(f"Test Accuracy: {test_score:.4f}")
    
    if len(np.unique(y_test)) > 1:
        auc_gb = roc_auc_score(y_test, y_proba_gb)
        print(f"ROC AUC: {auc_gb:.4f}")
    
    print()
    print("Classification Report:")
    print(classification_report(y_test, y_pred_gb, target_names=['Not Effective', 'Effective']))
    
    print("Confusion Matrix:")
    cm = confusion_matrix(y_test, y_pred_gb)
    print(f"  [[TN={cm[0,0]}, FP={cm[0,1]}],")
    print(f"   [FN={cm[1,0]}, TP={cm[1,1]}]]")
    print()
    
    # Cross-validation
    cv_scores = cross_val_score(gb_model, X_train_scaled, y_train, cv=5, scoring='accuracy')
    print(f"CV Accuracy: {cv_scores.mean():.4f} ± {cv_scores.std():.4f}")
    print()
    
    # === Train Random Forest Classifier ===
    print("=" * 70)
    print("2. Random Forest Classifier")
    print("=" * 70)
    
    rf_model = RandomForestClassifier(
        n_estimators=100,
        max_depth=5,
        min_samples_split=10,
        min_samples_leaf=4,
        max_features='sqrt',
        class_weight='balanced' if use_balanced else None,
        random_state=42,
        n_jobs=-1
    )
    
    print("Training...")
    rf_model.fit(X_train_scaled, y_train, sample_weight=w_train)
    
    # Evaluate
    train_score_rf = rf_model.score(X_train_scaled, y_train)
    test_score_rf = rf_model.score(X_test_scaled, y_test)
    
    y_pred_rf = rf_model.predict(X_test_scaled)
    y_proba_rf = rf_model.predict_proba(X_test_scaled)[:, 1]
    
    print(f"Train Accuracy: {train_score_rf:.4f}")
    print(f"Test Accuracy: {test_score_rf:.4f}")
    
    if len(np.unique(y_test)) > 1:
        auc_rf = roc_auc_score(y_test, y_proba_rf)
        print(f"ROC AUC: {auc_rf:.4f}")
    
    print()
    print("Classification Report:")
    print(classification_report(y_test, y_pred_rf, target_names=['Not Effective', 'Effective']))
    
    print("Confusion Matrix:")
    cm_rf = confusion_matrix(y_test, y_pred_rf)
    print(f"  [[TN={cm_rf[0,0]}, FP={cm_rf[0,1]}],")
    print(f"   [FN={cm_rf[1,0]}, TP={cm_rf[1,1]}]]")
    print()
    
    # Cross-validation
    cv_scores_rf = cross_val_score(rf_model, X_train_scaled, y_train, cv=5, scoring='accuracy')
    print(f"CV Accuracy: {cv_scores_rf.mean():.4f} ± {cv_scores_rf.std():.4f}")
    print()
    
    # === Select Best Model ===
    if test_score_rf > test_score:
        print("=" * 70)
        print("✅ Random Forest Selected (Higher Test Accuracy)")
        best_model = rf_model
        best_name = "RandomForest"
        best_test_score = test_score_rf
        best_auc = auc_rf if len(np.unique(y_test)) > 1 else None
    else:
        print("=" * 70)
        print("✅ Gradient Boosting Selected (Higher Test Accuracy)")
        best_model = gb_model
        best_name = "GradientBoosting"
        best_test_score = test_score
        best_auc = auc_gb if len(np.unique(y_test)) > 1 else None
    
    print(f"Test Accuracy: {best_test_score:.4f}")
    if best_auc:
        print(f"ROC AUC: {best_auc:.4f}")
    print("=" * 70)
    print()
    
    # Save model
    model_path = PROJECT_ROOT / 'models' / 'therapeutic_binary_classifier_v1.joblib'
    PROJECT_ROOT.joinpath('models').mkdir(exist_ok=True)
    
    joblib.dump({
        'model': best_model,
        'scaler': scaler,
        'feature_names': list(X.columns),
        'model_type': best_name,
        'test_accuracy': best_test_score,
        'roc_auc': best_auc,
        'class_distribution': {
            'positive': stats['positive_samples'],
            'negative': stats['negative_samples']
        }
    }, model_path)
    
    print(f"💾 Model saved: {model_path}")
    print()
    
    # Feature importance
    if hasattr(best_model, 'feature_importances_'):
        importances = best_model.feature_importances_
        feature_importance = sorted(zip(X.columns, importances), key=lambda x: x[1], reverse=True)
        
        print("Top 10 Feature Importances:")
        for i, (feature, importance) in enumerate(feature_importance[:10], 1):
            print(f"  {i}. {feature}: {importance:.4f}")
    
    print()
    print("=" * 70)
    print("✅ Binary Classification Training Complete")
    print("=" * 70)


if __name__ == '__main__':
    train_binary_classifier()
