#!/usr/bin/env python3
"""
Train NeuroBotanica ML Models
Uses fully enriched dataset with clinical studies, synergies, and patient data.

Trains:
1. TherapeuticPredictionModel - Cannabinoid efficacy for conditions
2. DimerPotentialModel - Dimer synergy prediction
3. PatientResponseModel - Patient-specific response prediction
"""
import json
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime
from typing import Dict, List
import time

# Add backend to path
SCRIPT_DIR = Path(__file__).parent
PROJECT_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(PROJECT_ROOT / 'backend'))

from services.ml_models import TherapeuticPredictionModel, DimerPotentialModel, PatientResponseModel
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import StandardScaler
import joblib


def load_enriched_dataset() -> Dict:
    """Load fully enriched dataset."""
    # Use the fixed dataset with merged therapeutic targets
    data_file = PROJECT_ROOT / 'data' / 'training' / 'neurobotanica_fully_enriched_fixed.json'
    
    if not data_file.exists():
        # Fall back to original if fix hasn't been run
        data_file = PROJECT_ROOT / 'data' / 'training' / 'neurobotanica_fully_enriched.json'
    
    with open(data_file, 'r', encoding='utf-8') as f:
        return json.load(f)


def prepare_therapeutic_data(dataset: Dict) -> tuple:
    """Prepare data for TherapeuticPredictionModel.
    
    Loads augmented CSV if available, otherwise falls back to on-the-fly extraction.
    
    Returns:
        (X_features, y_efficacy, sample_weights)
    """
    # Try to load pre-augmented data first
    augmented_csv = PROJECT_ROOT / 'data' / 'training' / 'therapeutic_augmented_data.csv'
    
    if augmented_csv.exists():
        print("Loading pre-augmented dataset...")
        df = pd.read_csv(augmented_csv)
        
        # Split features, target, and weights
        feature_cols = ['molecular_weight', 'logP', 'h_bond_donors', 'h_bond_acceptors',
                       'rotatable_bonds', 'tpsa', 'num_aromatic_rings', 'num_stereocenters',
                       'sp3_fraction', 'complexity', 'qed', 'bbb_penetration', 'has_phenol',
                       'num_terpene_synergies', 'avg_synergy_evidence_tier', 'target_encoded']
        
        X = df[feature_cols]
        y = df['efficacy']
        weights = df['confidence']
        
        print(f"   Total samples: {len(df)}")
        print(f"   Synthetic samples: {df['is_synthetic'].sum()}")
        print(f"   Original samples: {(~df['is_synthetic']).sum()}")
        
        return X, y, weights
    
    # Fall back to original extraction logic
    print("Augmented data not found, extracting from dataset...")
    compounds = dataset['compounds']
    
    X_list = []
    y_list = []
    weights_list = []
    
    for compound in compounds:
        # Skip if no RDKit descriptors
        descriptors = compound.get('rdkit_descriptors')
        if not descriptors or 'calculation_note' in descriptors:
            continue  # Skip compounds with estimated descriptors
        
        # Get clinical studies for confidence weighting
        clinical_studies = compound.get('clinical_studies', [])
        num_studies = len(clinical_studies)
        
        # Base confidence from study count
        base_confidence = min(1.0, 0.3 + (num_studies * 0.05))
        
        # Extract features from RDKit descriptors
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
        
        # Create targets for each therapeutic target
        therapeutic_targets = compound.get('therapeutic_targets', [])
        
        if not therapeutic_targets:
            continue
        
        for target in therapeutic_targets:
            # Estimate efficacy from clinical studies
            efficacy_score = _estimate_efficacy(compound, target, clinical_studies)
            
            # Confidence based on evidence
            confidence = base_confidence
            
            # Check if target matches any study condition (improved matching)
            target_normalized = target.upper().replace(' ', '_').replace('-', '_')
            has_direct_evidence = any(
                target_normalized in s.get('condition', '').upper().replace(' ', '_').replace('-', '_') or
                s.get('condition', '').upper().replace(' ', '_').replace('-', '_') in target_normalized
                for s in clinical_studies
            )
            
            if has_direct_evidence:
                confidence = min(1.0, confidence + 0.3)  # Boost if direct evidence
            
            features_with_target = features.copy()
            features_with_target['target_encoded'] = hash(target) % 100  # Simple target encoding
            
            X_list.append(features_with_target)
            y_list.append(efficacy_score)
            weights_list.append(confidence)
    
    X = pd.DataFrame(X_list)
    y = pd.Series(y_list)
    weights = pd.Series(weights_list)
    
    return X, y, weights


def _avg_synergy_tier(synergies: List[Dict]) -> float:
    """Average evidence tier of synergies."""
    if not synergies:
        return 0.0
    return np.mean([s.get('evidence_tier', 3) for s in synergies])


def _estimate_efficacy(compound: Dict, target: str, studies: List[Dict]) -> float:
    """Estimate efficacy score (0-1) from clinical evidence."""
    # Normalize target name for matching
    target_normalized = target.upper().replace(' ', '_').replace('-', '_')
    
    # Find relevant studies with fuzzy matching
    relevant_studies = []
    for s in studies:
        condition = s.get('condition', '').upper().replace(' ', '_').replace('-', '_')
        # Check if target matches condition or is a substring
        if (target_normalized in condition or 
            condition in target_normalized or
            target_normalized.replace('_', '') in condition.replace('_', '')):
            relevant_studies.append(s)
    
    if not relevant_studies:
        return 0.5  # Neutral baseline
    
    # Map effect sizes to scores
    effect_map = {
        'large': 0.85,
        'medium': 0.70,
        'small': 0.55,
        'minimal': 0.40,
        None: 0.50
    }
    
    efficacies = [effect_map.get(s.get('effect_size'), 0.5) for s in relevant_studies]
    return np.mean(efficacies)


def prepare_dimer_data(dataset: Dict) -> tuple:
    """Prepare data for DimerPotentialModel."""
    # Check if dimeric predictions exist
    dimeric_file = PROJECT_ROOT / 'data' / 'training' / 'neurobotanica_dimeric_predictions.json'
    
    if not dimeric_file.exists():
        print("⚠️  No dimeric data found, skipping DimerPotentialModel")
        return None, None, None
    
    with open(dimeric_file, 'r', encoding='utf-8') as f:
        dimeric_data = json.load(f)
    
    dimers = dimeric_data.get('dimeric_cannabinoids', [])
    
    if not dimers:
        return None, None, None
    
    X_list = []
    y_list = []
    weights_list = []
    
    for dimer in dimers:
        features = {
            'enhanced_lipophilicity': dimer.get('enhanced_lipophilicity', 0),
            'increased_receptor_affinity': dimer.get('increased_receptor_affinity', 0),
            'prolonged_duration': dimer.get('prolonged_duration', 0),
            'reduced_metabolism': dimer.get('reduced_metabolism', 0),
            'therapeutic_potential_score': dimer.get('therapeutic_potential_score', 0.5),
            'confidence_score': dimer.get('confidence_score', 0.5)
        }
        
        X_list.append(features)
        y_list.append(dimer.get('therapeutic_potential_score', 0.5))
        weights_list.append(dimer.get('confidence_score', 0.5))
    
    return pd.DataFrame(X_list), pd.Series(y_list), pd.Series(weights_list)


def prepare_patient_data(dataset: Dict) -> tuple:
    """Prepare data for PatientResponseModel."""
    patient_data = dataset.get('patient_response_data', {})
    profiles = patient_data.get('profiles', [])
    
    if not profiles:
        print("⚠️  No patient data found, skipping PatientResponseModel")
        return None, None, None
    
    X_list = []
    y_list = []
    weights_list = []
    
    for profile in profiles:
        demographics = profile.get('demographics', {})
        treatment = profile.get('treatment', {})
        response = profile.get('response', {})
        quality = profile.get('data_quality', {})
        
        # Encode demographics
        age_map = {"18-25": 0, "26-35": 1, "36-45": 2, "46-55": 3, "56-65": 4, "66+": 5}
        exp_map = {"naive": 0, "occasional": 1, "regular": 2, "daily": 3}
        tol_map = {"low": 0, "medium": 1, "high": 2}
        
        features = {
            'age_range_encoded': age_map.get(demographics.get('age_range'), 2),
            'experience_level': exp_map.get(demographics.get('experience_level'), 1),
            'thc_tolerance': tol_map.get(demographics.get('thc_tolerance'), 1),
            'treatment_duration_days': treatment.get('duration_days', 30),
            'adherence_score': response.get('adherence_score', 7) / 10.0,
            'num_adverse_effects': len(response.get('adverse_effects', [])),
            'study_quality_encoded': {'high': 3, 'good': 2, 'moderate': 1, 'low': 0}.get(
                quality.get('source_study_quality'), 1
            )
        }
        
        # Target: overall efficacy score
        efficacy = response.get('overall_efficacy_score', 5) / 10.0
        
        # Weight by study quality
        quality_weight = features['study_quality_encoded'] / 3.0
        quality_weight = max(0.3, quality_weight)  # Minimum weight
        
        X_list.append(features)
        y_list.append(efficacy)
        weights_list.append(quality_weight)
    
    return pd.DataFrame(X_list), pd.Series(y_list), pd.Series(weights_list)


def train_all_models():
    """Train all NeuroBotanica models."""
    print("=" * 70)
    print("NeuroBotanica ML Model Training")
    print("=" * 70)
    print()
    
    # Load dataset
    print("Loading fully enriched dataset...")
    dataset = load_enriched_dataset()
    print(f"✅ Loaded {len(dataset['compounds'])} compounds")
    
    norml_meta = dataset.get('enrichment_metadata', {}).get('norml_integration', {})
    total_studies = norml_meta.get('total_studies', 0) if norml_meta else 0
    
    patient_meta = dataset.get('patient_response_data', {}).get('metadata', {})
    total_patients = patient_meta.get('total_patients', 0) if patient_meta else 0
    
    print(f"   Clinical studies: {total_studies}")
    print(f"   Patient profiles: {total_patients}")
    print()
    
    results = {}
    
    # === 1. Therapeutic Prediction Model ===
    print("=" * 70)
    print("1. THERAPEUTIC PREDICTION MODEL")
    print("=" * 70)
    
    X_therapeutic, y_therapeutic, weights_therapeutic = prepare_therapeutic_data(dataset)
    print(f"Prepared {len(X_therapeutic)} training samples")
    print(f"Features: {list(X_therapeutic.columns)}")
    print()
    
    # Try both approaches for comparison
    print("Approach 1: Regularized Gradient Boosting")
    # Regularized hyperparameters to prevent overfitting with limited data
    model_therapeutic_gb = TherapeuticPredictionModel(
        n_estimators=100,        # Reduced from 200
        learning_rate=0.1,       # Increased from 0.05 for faster convergence
        max_depth=3,             # Reduced from 6 to prevent deep trees
        random_state=42
    )
    
    print("Training Gradient Boosting...")
    metrics_therapeutic_gb = model_therapeutic_gb.train(
        X_therapeutic, 
        y_therapeutic, 
        sample_weights=weights_therapeutic,
        test_size=0.2
    )
    
    print(f"   Train R²: {metrics_therapeutic_gb.train_score:.4f}")
    print(f"   Test R²: {metrics_therapeutic_gb.test_score:.4f}")
    print(f"   CV: {metrics_therapeutic_gb.cv_mean:.4f} ± {metrics_therapeutic_gb.cv_std:.4f}")
    print()
    
    print("Approach 2: Random Forest (better for small datasets)")
    # Random Forest is more robust with limited data
    from sklearn.model_selection import train_test_split, cross_val_score
    
    scaler_rf = StandardScaler()
    X_scaled = scaler_rf.fit_transform(X_therapeutic)
    
    X_train, X_test, y_train, y_test, w_train, w_test = train_test_split(
        X_scaled, y_therapeutic, weights_therapeutic,
        test_size=0.2, random_state=42
    )
    
    model_rf = RandomForestRegressor(
        n_estimators=100,
        max_depth=5,
        min_samples_split=10,
        min_samples_leaf=4,
        max_features='sqrt',
        random_state=42,
        n_jobs=-1
    )
    
    print("Training Random Forest...")
    start = time.time()
    model_rf.fit(X_train, y_train, sample_weight=w_train)
    
    train_score_rf = model_rf.score(X_train, y_train)
    test_score_rf = model_rf.score(X_test, y_test)
    cv_scores_rf = cross_val_score(model_rf, X_scaled, y_therapeutic, cv=5, scoring='r2')
    
    print(f"   Train R²: {train_score_rf:.4f}")
    print(f"   Test R²: {test_score_rf:.4f}")
    print(f"   CV: {cv_scores_rf.mean():.4f} ± {cv_scores_rf.std():.4f}")
    print(f"   Time: {time.time() - start:.2f}s")
    print()
    
    # Use the better model
    if test_score_rf > metrics_therapeutic_gb.test_score:
        print("✅ Random Forest performs better - using RF model")
        model_therapeutic = model_rf
        metrics_therapeutic = metrics_therapeutic_gb  # Keep GB metrics structure but update
        metrics_therapeutic.model_name = "TherapeuticPrediction_RandomForest"
        metrics_therapeutic.train_score = train_score_rf
        metrics_therapeutic.test_score = test_score_rf
        metrics_therapeutic.cv_scores = cv_scores_rf.tolist()
        metrics_therapeutic.cv_mean = cv_scores_rf.mean()
        metrics_therapeutic.cv_std = cv_scores_rf.std()
        
        # Save RF model separately
        model_path_rf = PROJECT_ROOT / 'models' / 'therapeutic_prediction_rf_v1.joblib'
        joblib.dump({
            'model': model_rf,
            'scaler': scaler_rf,
            'feature_names': list(X_therapeutic.columns)
        }, model_path_rf)
        print(f"💾 Saved Random Forest: {model_path_rf}")
    else:
        print("✅ Gradient Boosting performs better - using GB model")
        model_therapeutic = model_therapeutic_gb
        metrics_therapeutic = metrics_therapeutic_gb
        model_path = PROJECT_ROOT / 'models' / 'therapeutic_prediction_v1.joblib'
        model_therapeutic.save(str(model_path))
        print(f"💾 Saved Gradient Boosting: {model_path}")
    
    print()
    
    results['therapeutic'] = metrics_therapeutic.to_dict()
    
    # === 2. Dimer Potential Model ===
    print("=" * 70)
    print("2. DIMER POTENTIAL MODEL")
    print("=" * 70)
    
    X_dimer, y_dimer, weights_dimer = prepare_dimer_data(dataset)
    
    if X_dimer is not None and len(X_dimer) > 10:
        print(f"Prepared {len(X_dimer)} dimer samples")
        
        model_dimer = DimerPotentialModel(random_state=42)
        
        print("Training...")
        metrics_dimer = model_dimer.train(
            X_dimer, 
            y_dimer, 
            sample_weights=weights_dimer,
            test_size=0.2
        )
        
        print(f"✅ Training complete:")
        print(f"   Train R²: {metrics_dimer.train_score:.4f}")
        print(f"   Test R²: {metrics_dimer.test_score:.4f}")
        print(f"   CV: {metrics_dimer.cv_mean:.4f} ± {metrics_dimer.cv_std:.4f}")
        print()
        
        model_path = PROJECT_ROOT / 'models' / 'dimer_potential_v1.joblib'
        model_dimer.save(str(model_path))
        print(f"💾 Saved: {model_path}")
        print()
        
        results['dimer'] = metrics_dimer.to_dict()
    else:
        print("⏭️  Skipped (insufficient data)")
        print()
    
    # === 3. Patient Response Model ===
    print("=" * 70)
    print("3. PATIENT RESPONSE MODEL")
    print("=" * 70)
    
    X_patient, y_patient, weights_patient = prepare_patient_data(dataset)
    
    if X_patient is not None and len(X_patient) > 10:
        print(f"Prepared {len(X_patient)} patient samples")
        
        model_patient = PatientResponseModel(random_state=42)
        
        print("Training...")
        metrics_patient = model_patient.train(
            X_patient, 
            y_patient, 
            sample_weights=weights_patient,
            test_size=0.2
        )
        
        print(f"✅ Training complete:")
        print(f"   Train R²: {metrics_patient.train_score:.4f}")
        print(f"   Test R²: {metrics_patient.test_score:.4f}")
        print(f"   CV: {metrics_patient.cv_mean:.4f} ± {metrics_patient.cv_std:.4f}")
        print()
        
        model_path = PROJECT_ROOT / 'models' / 'patient_response_v1.joblib'
        model_patient.save(str(model_path))
        print(f"💾 Saved: {model_path}")
        print()
        
        results['patient'] = metrics_patient.to_dict()
    else:
        print("⏭️  Skipped (insufficient data)")
        print()
    
    # === Save Training Report ===
    norml_meta = dataset.get('enrichment_metadata', {}).get('norml_integration', {})
    synergy_meta = dataset.get('enrichment_metadata', {}).get('synergy_integration', {})
    patient_meta = dataset.get('patient_response_data', {}).get('metadata', {})
    
    report = {
        "training_date": datetime.now().isoformat(),
        "dataset": {
            "total_compounds": len(dataset['compounds']),
            "clinical_studies": norml_meta.get('total_studies', 0) if norml_meta else 0,
            "patient_profiles": patient_meta.get('total_patients', 0) if patient_meta else 0,
            "terpene_synergies": synergy_meta.get('statistics', {}).get('total_synergies_added', 0) if synergy_meta else 0
        },
        "models": results
    }
    
    report_file = PROJECT_ROOT / 'models' / 'training_report.json'
    PROJECT_ROOT.joinpath('models').mkdir(exist_ok=True)
    with open(report_file, 'w') as f:
        json.dump(report, f, indent=2)
    
    print("=" * 70)
    print("✅ ALL MODELS TRAINED SUCCESSFULLY")
    print(f"📊 Training report: {report_file}")
    print("=" * 70)
    
    return results


if __name__ == '__main__':
    train_all_models()
