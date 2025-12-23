#!/usr/bin/env python3
"""
Integrate Terpene-Cannabinoid Synergy Data with Evidence Tiers
NeuroBotanica Data Enrichment Step 2/4

Adds structured synergy interactions from terpene_analyzer.py to dataset.
"""
import json
import sys
from pathlib import Path
from typing import Dict, List, Set
from collections import defaultdict

# Add backend to path
SCRIPT_DIR = Path(__file__).parent
PROJECT_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(PROJECT_ROOT / 'backend'))

from services.terpene_analyzer import TerpeneDatabase, EvidenceTier


def extract_synergy_data() -> Dict[str, List[Dict]]:
    """Extract synergy data from TerpeneDatabase."""
    db = TerpeneDatabase()
    synergies = defaultdict(list)
    
    for terpene_name, terpene_data in db.terpenes.items():
        for cannabinoid, synergy_info in terpene_data.cannabinoid_synergy.items():
            synergies[cannabinoid].append({
                'partner_compound': terpene_name,
                'partner_type': 'terpene',
                'mechanism': terpene_data.mechanism,
                'enhancement_factor': synergy_info['enhancement_factor'],
                'evidence_tier': int(synergy_info['evidence_tier']),
                'evidence_quality': _tier_to_label(synergy_info['evidence_tier']),
                'therapeutic_effects': terpene_data.therapeutic_effects,
                'pubmed_ids': terpene_data.pubmed_ids[:3],  # Top 3 refs
                'boiling_point_c': terpene_data.boiling_point_c,
                'aroma_profile': terpene_data.aroma_profile
            })
    
    return dict(synergies)


def _tier_to_label(tier: int) -> str:
    """Convert evidence tier to quality label."""
    labels = {
        1: "Very Low - Anecdotal/Theoretical",
        2: "Low - In Vitro/Animal Studies",
        3: "Moderate - Small Studies/Case Series",
        4: "Good - Well-Designed Observational",
        5: "High - Multiple RCTs/Reviews"
    }
    return labels.get(tier, "Unknown")


def find_cannabinoid_matches(compound: Dict, synergy_data: Dict) -> List[Dict]:
    """Match compound to synergy data."""
    matches = []
    compound_name = compound.get('name', '').upper()
    
    # Check exact matches and common variants
    for cannabinoid_key in synergy_data.keys():
        if (cannabinoid_key.upper() in compound_name or 
            compound_name in cannabinoid_key.upper() or
            _normalize_name(cannabinoid_key) == _normalize_name(compound_name)):
            matches.extend(synergy_data[cannabinoid_key])
    
    return matches


def _normalize_name(name: str) -> str:
    """Normalize compound name for matching."""
    # Remove common prefixes/suffixes
    name = name.upper()
    name = name.replace('Δ9-', '').replace('DELTA-9-', '')
    name = name.replace('-ACID', 'A').replace('ACID', 'A')
    name = name.replace('-', '').replace(' ', '')
    return name


def enrich_dataset_with_synergies(dataset: Dict, synergy_data: Dict) -> Dict:
    """Add synergy data to dataset compounds."""
    enriched = dataset.copy()
    
    stats = {
        'compounds_enhanced': 0,
        'total_synergies_added': 0,
        'tier_5_synergies': 0,
        'tier_4_synergies': 0,
        'tier_3_synergies': 0
    }
    
    for compound in enriched['compounds']:
        # Find synergies for this compound
        synergies = find_cannabinoid_matches(compound, synergy_data)
        
        if synergies:
            # Add synergies to compound
            compound['terpene_synergies'] = synergies
            stats['compounds_enhanced'] += 1
            stats['total_synergies_added'] += len(synergies)
            
            # Count by tier
            for synergy in synergies:
                tier = synergy['evidence_tier']
                if tier == 5:
                    stats['tier_5_synergies'] += 1
                elif tier == 4:
                    stats['tier_4_synergies'] += 1
                elif tier == 3:
                    stats['tier_3_synergies'] += 1
    
    enriched['enrichment_metadata'] = enriched.get('enrichment_metadata', {})
    enriched['enrichment_metadata']['synergy_integration'] = {
        'date': '2025-12-23',
        'source': 'backend/services/terpene_analyzer.py',
        'statistics': stats
    }
    
    return enriched


def main():
    """Main integration workflow."""
    print("=" * 70)
    print("Terpene-Cannabinoid Synergy Integration")
    print("=" * 70)
    print()
    
    # Input/output paths
    input_file = PROJECT_ROOT / 'data' / 'training' / 'neurobotanica_enriched_dataset.json'
    output_file = PROJECT_ROOT / 'data' / 'training' / 'neurobotanica_enriched_synergies.json'
    
    # Load enriched dataset from step 1
    print(f"Loading enriched dataset: {input_file}")
    with open(input_file, 'r', encoding='utf-8') as f:
        dataset = json.load(f)
    print(f"✅ Loaded {len(dataset['compounds'])} compounds")
    print()
    
    # Extract synergy data
    print("Extracting synergy data from TerpeneDatabase...")
    synergy_data = extract_synergy_data()
    print(f"✅ Extracted synergies for {len(synergy_data)} cannabinoids:")
    for cannabinoid, synergies in synergy_data.items():
        print(f"   {cannabinoid}: {len(synergies)} terpene synergies")
    print()
    
    # Enrich dataset
    print("Adding synergy data to compounds...")
    enriched = enrich_dataset_with_synergies(dataset, synergy_data)
    
    stats = enriched['enrichment_metadata']['synergy_integration']['statistics']
    print(f"   Compounds enhanced: {stats['compounds_enhanced']}")
    print(f"   Total synergies added: {stats['total_synergies_added']}")
    print(f"   High evidence (Tier 5): {stats['tier_5_synergies']}")
    print(f"   Good evidence (Tier 4): {stats['tier_4_synergies']}")
    print(f"   Moderate evidence (Tier 3): {stats['tier_3_synergies']}")
    print()
    
    # Save
    print(f"Saving: {output_file}")
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(enriched, f, indent=2)
    
    print("=" * 70)
    print("✅ Synergy Integration Complete!")
    print(f"✅ Output: {output_file}")
    print("=" * 70)


if __name__ == '__main__':
    main()
