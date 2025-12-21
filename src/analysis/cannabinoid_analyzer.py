"""
Cannabinoid Analyzer - Main Analysis Script
============================================

This script loads and analyzes the NeuroBotanica cannabinoid dataset,
providing insights into molecular properties, receptor binding, and
therapeutic predictions.

Usage:
    python src/analysis/cannabinoid_analyzer.py
"""

import json
import sys
from pathlib import Path
from typing import Dict, List, Any
import pandas as pd
import numpy as np

# Add project root to path for imports
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))


def load_dataset(filepath: str = "data/training/neurobotanica_complete_dataset_63compounds.json") -> Dict[str, Any]:
    """
    Load the complete NeuroBotanica dataset.
    
    Args:
        filepath: Path to the dataset JSON file
        
    Returns:
        Dictionary containing all cannabinoid data
        
    Raises:
        FileNotFoundError: If dataset file doesn't exist
        json.JSONDecodeError: If file is not valid JSON
    """
    full_path = project_root / filepath
    
    if not full_path.exists():
        raise FileNotFoundError(f"Dataset not found at: {full_path}")
    
    with open(full_path, 'r') as f:
        data = json.load(f)
    
    print(f"✓ Loaded dataset from: {filepath}")
    print(f"  Total compounds: {len(data.get('compounds', {}))}")
    print(f"  Clinical studies: {len(data.get('clinical_studies', []))}")
    
    return data


def analyze_molecular_properties(data: Dict[str, Any]) -> pd.DataFrame:
    """
    Analyze molecular properties of cannabinoids.
    
    Args:
        data: Complete dataset dictionary
        
    Returns:
        DataFrame with molecular property analysis
    """
    print("\n" + "="*60)
    print("MOLECULAR PROPERTY ANALYSIS")
    print("="*60)
    
    compounds = []
    for cannabinoid in data['compounds']:
        # Skip dimeric compounds for basic analysis (analyze separately)
        if cannabinoid.get('compound_type') == 'dimeric':
            continue
        
        name = cannabinoid.get('compound_name', 'unknown')
            
        rdkit = cannabinoid.get('rdkit_descriptors', {})
        
        compounds.append({
            'name': name,
            'type': cannabinoid.get('compound_type', 'unknown'),
            'molecular_weight': rdkit.get('molecular_weight'),
            'logP': rdkit.get('logP'),
            'TPSA': rdkit.get('tpsa'),
            'h_bond_donors': rdkit.get('h_bond_donors'),
            'h_bond_acceptors': rdkit.get('h_bond_acceptors'),
            'bbb_penetration': rdkit.get('bbb_penetration_predicted'),
            'psychoactive': cannabinoid.get('psychoactive', False),
            'num_stereocenters': rdkit.get('num_stereocenters', 0)
        })
    
    df = pd.DataFrame(compounds)
    
    print(f"\nMonomeric Cannabinoids: {len(df)}")
    print("\nDescriptive Statistics:")
    print(df[['molecular_weight', 'logP', 'TPSA', 'h_bond_donors', 'h_bond_acceptors']].describe())
    
    print("\nBlood-Brain Barrier Penetration:")
    print(df.groupby('bbb_penetration')['name'].count())
    
    print("\nPsychoactivity by LogP (lipophilicity):")
    print(df.groupby('psychoactive')['logP'].agg(['mean', 'std', 'count']))
    
    return df


def analyze_receptor_binding(data: Dict[str, Any]) -> pd.DataFrame:
    """
    Analyze receptor binding profiles.
    
    Args:
        data: Complete dataset dictionary
        
    Returns:
        DataFrame with receptor binding analysis
    """
    print("\n" + "="*60)
    print("RECEPTOR BINDING ANALYSIS")
    print("="*60)
    
    binding_data = []
    for cannabinoid in data['compounds']:
        if cannabinoid.get('compound_type') == 'dimeric':
            continue
        
        name = cannabinoid.get('compound_name', 'unknown')
            
        receptor_binding = cannabinoid.get('receptor_binding', {})
        
        for receptor, binding in receptor_binding.items():
            affinity = binding.get('affinity_nM')
            if affinity and affinity != 'inactive':
                binding_data.append({
                    'cannabinoid': name,
                    'receptor': receptor,
                    'affinity_nM': float(affinity),
                    'activity': binding.get('activity', 'unknown')
                })
    
    df = pd.DataFrame(binding_data)
    
    if len(df) > 0:
        print(f"\nTotal binding measurements: {len(df)}")
        print("\nTop 10 Strongest CB1 Binders (lowest Ki):")
        cb1_binders = df[df['receptor'] == 'CB1'].sort_values('affinity_nM')
        print(cb1_binders[['cannabinoid', 'affinity_nM', 'activity']].head(10).to_string(index=False))
        
        print("\nTop 10 Strongest CB2 Binders:")
        cb2_binders = df[df['receptor'] == 'CB2'].sort_values('affinity_nM')
        print(cb2_binders[['cannabinoid', 'affinity_nM', 'activity']].head(10).to_string(index=False))
        
        print("\nReceptor Coverage:")
        print(df.groupby('receptor')['cannabinoid'].count().sort_values(ascending=False))
    else:
        print("\n⚠ No receptor binding data found in dataset")
    
    return df


def analyze_clinical_studies(data: Dict[str, Any]) -> pd.DataFrame:
    """
    Analyze clinical study patterns.
    
    Args:
        data: Complete dataset dictionary
        
    Returns:
        DataFrame with clinical study analysis
    """
    print("\n" + "="*60)
    print("CLINICAL STUDY ANALYSIS")
    print("="*60)
    
    studies = data.get('clinical_studies', [])
    
    if not studies:
        print("\n⚠ No clinical studies found in dataset")
        return pd.DataFrame()
    
    study_data = []
    for study in studies:
        study_data.append({
            'condition': study.get('condition'),
            'year': study.get('year'),
            'sample_size': study.get('sample_size', 0),
            'num_cannabinoids': len(study.get('cannabinoid_profile', [])),
            'inference_status': study.get('inference_status', 'unknown'),
            'confidence': study.get('confidence_weight', 1.0)
        })
    
    df = pd.DataFrame(study_data)
    
    print(f"\nTotal Studies: {len(df)}")
    print("\nStudies by Condition:")
    print(df.groupby('condition').size().sort_values(ascending=False))
    
    print("\nData Quality (Inference Status):")
    print(df.groupby('inference_status').size())
    
    print("\nAverage Confidence Weight by Condition:")
    print(df.groupby('condition')['confidence'].mean().sort_values(ascending=False))
    
    print("\nStudy Timeline:")
    if 'year' in df.columns and df['year'].notna().any():
        print(df.groupby('year').size().sort_index())
    
    return df


def find_therapeutic_candidates(data: Dict[str, Any], condition: str = "Epilepsy") -> pd.DataFrame:
    """
    Find cannabinoid candidates for a specific therapeutic condition.
    
    Args:
        data: Complete dataset dictionary
        condition: Therapeutic condition to search for
        
    Returns:
        DataFrame with ranked cannabinoid candidates
    """
    print("\n" + "="*60)
    print(f"THERAPEUTIC CANDIDATES FOR: {condition}")
    print("="*60)
    
    candidates = []
    for cannabinoid in data['compounds']:
        name = cannabinoid.get('compound_name', 'unknown')
        targets = cannabinoid.get('therapeutic_targets', [])
        
        # Check if condition is in therapeutic targets (case-insensitive)
        targets_lower = [t.lower() for t in targets]
        if condition.lower() in targets_lower:
            
            rdkit = cannabinoid.get('rdkit_descriptors', {})
            receptor_binding = cannabinoid.get('receptor_binding', {})
            
            # Calculate a simple efficacy score based on CB1/CB2 binding
            cb1_affinity = receptor_binding.get('CB1', {}).get('affinity_nM')
            cb2_affinity = receptor_binding.get('CB2', {}).get('affinity_nM')
            
            # Lower affinity = stronger binding = higher score
            efficacy_score = 0.0
            if cb1_affinity and cb1_affinity != 'inactive':
                efficacy_score += 1.0 / (float(cb1_affinity) + 1.0) * 1000
            if cb2_affinity and cb2_affinity != 'inactive':
                efficacy_score += 1.0 / (float(cb2_affinity) + 1.0) * 1000
            
            candidates.append({
                'cannabinoid': name,
                'type': cannabinoid.get('compound_type', 'unknown'),
                'psychoactive': cannabinoid.get('psychoactive', False),
                'bbb_penetration': rdkit.get('bbb_penetration_predicted', False),
                'cb1_affinity_nM': cb1_affinity if cb1_affinity != 'inactive' else None,
                'cb2_affinity_nM': cb2_affinity if cb2_affinity != 'inactive' else None,
                'efficacy_score': efficacy_score
            })
    
    df = pd.DataFrame(candidates)
    
    if len(df) > 0:
        df = df.sort_values('efficacy_score', ascending=False)
        print(f"\nFound {len(df)} candidates for {condition}")
        print("\nTop 10 Candidates (ranked by efficacy score):")
        print(df[['cannabinoid', 'type', 'psychoactive', 'bbb_penetration', 'efficacy_score']].head(10).to_string(index=False))
    else:
        print(f"\n⚠ No cannabinoid candidates found for {condition}")
    
    return df


def analyze_dimeric_predictions(data: Dict[str, Any], condition: str = "Epilepsy", min_synergy: float = 0.7) -> pd.DataFrame:
    """
    Analyze dimeric cannabinoid predictions.
    
    Args:
        data: Complete dataset dictionary
        condition: Therapeutic condition to filter by
        min_synergy: Minimum synergy score threshold
        
    Returns:
        DataFrame with dimeric predictions
    """
    print("\n" + "="*60)
    print(f"DIMERIC CANNABINOID PREDICTIONS")
    print(f"Condition: {condition} | Min Synergy: {min_synergy}")
    print("="*60)
    
    dimers = []
    for cannabinoid in data['compounds']:
        if cannabinoid.get('compound_type') != 'dimeric':
            continue
        
        name = cannabinoid.get('compound_name', 'unknown')
        
        synergy = cannabinoid.get('synergy_score', 0.0)
        formation_prob = cannabinoid.get('formation_probability', 0.0)
        therapeutic_potential = cannabinoid.get('therapeutic_potential', {})
        
        condition_potential = therapeutic_potential.get(condition, 0.0)
        
        if synergy >= min_synergy and condition_potential > 0:
            dimers.append({
                'name': name,
                'monomer_1': cannabinoid.get('monomer_1'),
                'monomer_2': cannabinoid.get('monomer_2'),
                'synergy_score': synergy,
                'formation_probability': formation_prob,
                'therapeutic_potential': condition_potential,
                'formation_mechanism': cannabinoid.get('formation_mechanism', 'unknown'),
                'composite_score': synergy * 0.4 + formation_prob * 0.3 + condition_potential * 0.3
            })
    
    df = pd.DataFrame(dimers)
    
    if len(df) > 0:
        df = df.sort_values('composite_score', ascending=False)
        print(f"\nFound {len(df)} dimeric predictions matching criteria")
        print("\nTop 10 Dimeric Candidates:")
        print(df[['name', 'monomer_1', 'monomer_2', 'synergy_score', 'formation_probability', 'therapeutic_potential', 'composite_score']].head(10).to_string(index=False))
        
        print("\nFormation Mechanism Distribution:")
        print(df.groupby('formation_mechanism').size())
    else:
        print(f"\n⚠ No dimeric predictions found matching criteria")
    
    return df


def main():
    """Main analysis workflow."""
    print("="*60)
    print("NEUROBOTANICA CANNABINOID ANALYZER")
    print("Phase 1 Dataset Analysis")
    print("="*60)
    
    # Load dataset
    try:
        data = load_dataset()
    except FileNotFoundError as e:
        print(f"\n❌ Error: {e}")
        print("\nMake sure you're running from the project root directory:")
        print("  cd neurobotanica_project")
        print("  python src/analysis/cannabinoid_analyzer.py")
        return
    
    # Run analyses
    molecular_df = analyze_molecular_properties(data)
    receptor_df = analyze_receptor_binding(data)
    clinical_df = analyze_clinical_studies(data)
    
    # Therapeutic candidate search
    epilepsy_candidates = find_therapeutic_candidates(data, condition="Epilepsy")
    ptsd_candidates = find_therapeutic_candidates(data, condition="PTSD")
    
    # Dimeric predictions
    dimeric_df = analyze_dimeric_predictions(data, condition="Epilepsy", min_synergy=0.7)
    
    print("\n" + "="*60)
    print("ANALYSIS COMPLETE")
    print("="*60)
    print("\nDatasets created:")
    print(f"  - molecular_df: {len(molecular_df)} monomeric cannabinoids")
    print(f"  - receptor_df: {len(receptor_df)} receptor binding measurements")
    print(f"  - clinical_df: {len(clinical_df)} clinical studies")
    print(f"  - epilepsy_candidates: {len(epilepsy_candidates)} candidates")
    print(f"  - ptsd_candidates: {len(ptsd_candidates)} candidates")
    print(f"  - dimeric_df: {len(dimeric_df)} dimeric predictions")
    
    print("\nNext steps:")
    print("  1. Train ML models: python src/ml_models/sar_predictor.py")
    print("  2. Run dimeric screening: python src/ml_models/dimeric_predictor.py")
    print("  3. View integration guide: docs/INTEGRATION_GUIDE.md")


if __name__ == "__main__":
    main()
