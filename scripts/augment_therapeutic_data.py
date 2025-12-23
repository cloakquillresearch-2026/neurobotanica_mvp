#!/usr/bin/env python3
"""
Synthetic Data Augmentation for TherapeuticPredictionModel
Generate 10x samples from existing clinical evidence with realistic perturbations.

Strategy:
- For each compound-target-study combination, create 10 synthetic variants
- Add Gaussian noise to molecular features (σ=2-5% of std)
- Add noise to efficacy scores based on effect size confidence
- Maintain realistic bounds and correlations
"""
import json
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple
from collections import defaultdict

# Add backend to path
SCRIPT_DIR = Path(__file__).parent
PROJECT_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(PROJECT_ROOT / 'backend'))


def load_enriched_dataset() -> Dict:
    """Load fixed dataset."""
    data_file = PROJECT_ROOT / 'data' / 'training' / 'neurobotanica_fully_enriched_fixed.json'
    
    with open(data_file, 'r', encoding='utf-8') as f:
        return json.load(f)


def prepare_base_samples(dataset: Dict) -> List[Dict]:
    """Extract base samples (same logic as training script)."""
    compounds = dataset['compounds']
    
    samples = []
    
    for compound in compounds:
        descriptors = compound.get('rdkit_descriptors')
        if not descriptors or 'calculation_note' in descriptors:
            continue
        
        clinical_studies = compound.get('clinical_studies', [])
        num_studies = len(clinical_studies)
        
        base_confidence = min(1.0, 0.3 + (num_studies * 0.05))
        
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
            efficacy_score = _estimate_efficacy(compound, target, clinical_studies)
            
            target_normalized = target.upper().replace(' ', '_').replace('-', '_')
            has_direct_evidence = any(
                target_normalized in s.get('condition', '').upper().replace(' ', '_').replace('-', '_') or
                s.get('condition', '').upper().replace(' ', '_').replace('-', '_') in target_normalized
                for s in clinical_studies
            )
            
            confidence = base_confidence
            if has_direct_evidence:
                confidence = min(1.0, confidence + 0.3)
            
            features_with_target = features.copy()
            features_with_target['target_encoded'] = hash(target) % 100
            
            samples.append({
                'features': features_with_target,
                'efficacy': efficacy_score,
                'confidence': confidence,
                'compound': compound['compound_name'],
                'target': target,
                'has_evidence': has_direct_evidence
            })
    
    return samples


def _avg_synergy_tier(synergies: List[Dict]) -> float:
    if not synergies:
        return 0.0
    return np.mean([s.get('evidence_tier', 3) for s in synergies])


def _estimate_efficacy(compound: Dict, target: str, studies: List[Dict]) -> float:
    target_normalized = target.upper().replace(' ', '_').replace('-', '_')
    
    relevant_studies = []
    for s in studies:
        condition = s.get('condition', '').upper().replace(' ', '_').replace('-', '_')
        if (target_normalized in condition or 
            condition in target_normalized or
            target_normalized.replace('_', '') in condition.replace('_', '')):
            relevant_studies.append(s)
    
    if not relevant_studies:
        return 0.5
    
    effect_map = {
        'large': 0.85,
        'medium': 0.70,
        'small': 0.55,
        'minimal': 0.40,
        None: 0.50
    }
    
    efficacies = [effect_map.get(s.get('effect_size'), 0.5) for s in relevant_studies]
    return np.mean(efficacies)


def augment_sample(sample: Dict, augmentation_id: int, feature_stats: Dict) -> Dict:
    """Create synthetic variant with realistic perturbations."""
    np.random.seed(hash((sample['compound'], sample['target'], augmentation_id)) % (2**32))
    
    augmented = {
        'features': {},
        'efficacy': sample['efficacy'],
        'confidence': sample['confidence'],
        'is_synthetic': True,
        'parent_compound': sample['compound'],
        'parent_target': sample['target']
    }
    
    # Determine noise levels based on evidence quality
    if sample['has_evidence']:
        efficacy_noise_std = 0.05  # ±5% for studies with direct evidence
        feature_noise_factor = 0.02  # ±2% feature perturbation
    else:
        efficacy_noise_std = 0.15  # ±15% for inferred efficacy
        feature_noise_factor = 0.05  # ±5% feature perturbation
    
    # Perturb efficacy score
    efficacy_perturbed = np.random.normal(sample['efficacy'], efficacy_noise_std)
    augmented['efficacy'] = np.clip(efficacy_perturbed, 0.0, 1.0)
    
    # Perturb molecular features
    for feature, value in sample['features'].items():
        if feature == 'target_encoded':
            # Keep target encoding unchanged
            augmented['features'][feature] = value
        elif feature in ['bbb_penetration', 'has_phenol']:
            # Binary features: keep unchanged
            augmented['features'][feature] = value
        elif feature in ['h_bond_donors', 'h_bond_acceptors', 'rotatable_bonds', 
                         'num_aromatic_rings', 'num_stereocenters', 'num_terpene_synergies']:
            # Integer features: add small integer noise
            noise = np.random.choice([-1, 0, 0, 0, 1])  # Bias toward no change
            augmented['features'][feature] = max(0, value + noise)
        else:
            # Continuous features: add Gaussian noise proportional to feature std
            feature_std = feature_stats.get(feature, {}).get('std', 1.0)
            noise = np.random.normal(0, feature_std * feature_noise_factor)
            
            # Apply noise and clip to reasonable bounds
            perturbed = value + noise
            
            # Feature-specific bounds
            if feature == 'sp3_fraction' or feature == 'qed':
                perturbed = np.clip(perturbed, 0.0, 1.0)
            elif feature == 'logP':
                perturbed = np.clip(perturbed, -2.0, 10.0)
            elif feature == 'avg_synergy_evidence_tier':
                perturbed = np.clip(perturbed, 0.0, 5.0)
            else:
                perturbed = max(0, perturbed)  # Non-negative
            
            augmented['features'][feature] = perturbed
    
    return augmented


def calculate_feature_statistics(samples: List[Dict]) -> Dict:
    """Calculate mean/std for each feature across base samples."""
    feature_values = defaultdict(list)
    
    for sample in samples:
        for feature, value in sample['features'].items():
            if feature not in ['bbb_penetration', 'has_phenol', 'target_encoded']:
                feature_values[feature].append(value)
    
    stats = {}
    for feature, values in feature_values.items():
        stats[feature] = {
            'mean': np.mean(values),
            'std': np.std(values),
            'min': np.min(values),
            'max': np.max(values)
        }
    
    return stats


def augment_dataset(base_samples: List[Dict], augmentation_factor: int = 10) -> List[Dict]:
    """Generate synthetic samples."""
    print(f"Calculating feature statistics...")
    feature_stats = calculate_feature_statistics(base_samples)
    
    print(f"Generating {augmentation_factor}x synthetic samples...")
    all_samples = []
    
    # Keep original samples
    for sample in base_samples:
        sample['is_synthetic'] = False
        all_samples.append(sample)
    
    # Generate synthetic variants
    for sample in base_samples:
        for aug_id in range(augmentation_factor):
            synthetic = augment_sample(sample, aug_id, feature_stats)
            all_samples.append(synthetic)
    
    return all_samples


def save_augmented_data(samples: List[Dict], output_path: Path):
    """Save to JSON for inspection and CSV for training."""
    # Convert to DataFrame
    X_list = [s['features'] for s in samples]
    y_list = [s['efficacy'] for s in samples]
    w_list = [s['confidence'] for s in samples]
    
    df = pd.DataFrame(X_list)
    df['efficacy'] = y_list
    df['confidence'] = w_list
    df['is_synthetic'] = [s.get('is_synthetic', False) for s in samples]
    df['compound'] = [s.get('parent_compound', s.get('compound')) for s in samples]
    df['target'] = [s.get('parent_target', s.get('target')) for s in samples]
    
    # Save CSV
    csv_path = output_path.with_suffix('.csv')
    df.to_csv(csv_path, index=False)
    
    # Save JSON summary
    summary = {
        'total_samples': len(samples),
        'original_samples': sum(1 for s in samples if not s.get('is_synthetic', False)),
        'synthetic_samples': sum(1 for s in samples if s.get('is_synthetic', False)),
        'augmentation_factor': len(samples) // sum(1 for s in samples if not s.get('is_synthetic', False)) - 1,
        'features': list(df.columns[:16]),
        'efficacy_distribution': {
            'mean': float(df['efficacy'].mean()),
            'std': float(df['efficacy'].std()),
            'min': float(df['efficacy'].min()),
            'max': float(df['efficacy'].max())
        },
        'samples_by_evidence': {
            'with_direct_evidence': len([s for s in samples if s.get('has_evidence', False)]),
            'inferred': len([s for s in samples if not s.get('has_evidence', False)])
        }
    }
    
    json_path = output_path.with_suffix('.json')
    with open(json_path, 'w') as f:
        json.dump(summary, f, indent=2)
    
    return csv_path, json_path


def main():
    """Main augmentation workflow."""
    print("=" * 70)
    print("Synthetic Data Augmentation for TherapeuticPredictionModel")
    print("=" * 70)
    print()
    
    # Load dataset
    print("Loading enriched dataset...")
    dataset = load_enriched_dataset()
    print(f"✅ Loaded {len(dataset['compounds'])} compounds")
    print()
    
    # Extract base samples
    print("Extracting base training samples...")
    base_samples = prepare_base_samples(dataset)
    print(f"✅ {len(base_samples)} base samples")
    print(f"   With direct evidence: {sum(1 for s in base_samples if s['has_evidence'])}")
    print(f"   Inferred (no direct evidence): {sum(1 for s in base_samples if not s['has_evidence'])}")
    print()
    
    # Augment
    augmentation_factor = 10
    print(f"Augmenting with factor={augmentation_factor}...")
    augmented_samples = augment_dataset(base_samples, augmentation_factor)
    print(f"✅ Generated {len(augmented_samples)} total samples")
    print(f"   Original: {len(base_samples)}")
    print(f"   Synthetic: {len(augmented_samples) - len(base_samples)}")
    print()
    
    # Check if target met
    target_min = 800
    if len(augmented_samples) >= target_min:
        print(f"✅ Target met: {len(augmented_samples)} ≥ {target_min} samples")
    else:
        print(f"⚠️  Below target: {len(augmented_samples)} < {target_min} samples")
        print(f"   Try augmentation_factor={target_min // len(base_samples) + 1}")
    print()
    
    # Save
    output_path = PROJECT_ROOT / 'data' / 'training' / 'therapeutic_augmented_data'
    csv_path, json_path = save_augmented_data(augmented_samples, output_path)
    
    print("=" * 70)
    print("✅ Augmentation Complete!")
    print(f"📊 Data: {csv_path}")
    print(f"📋 Summary: {json_path}")
    print()
    print("Next steps:")
    print("  1. Review augmented data quality")
    print("  2. Update train_neurobotanica.py to load CSV")
    print("  3. Retrain TherapeuticPredictionModel")
    print("=" * 70)


if __name__ == '__main__':
    main()
