#!/usr/bin/env python3
"""
Add Cannabinoid Compound Variations
NeuroBotanica Data Enrichment Step 4/4

Adds emerging and rare cannabinoids plus structural variations.
"""
import json
from pathlib import Path
from typing import Dict, List


# New cannabinoid variants to add (emerging research compounds)
NEW_CANNABINOIDS = [
    {
        "compound_name": "CBDP",
        "full_name": "Cannabidiphorol",
        "compound_type": "cannabinoid",
        "molecular_formula": "C23H34O2",
        "smiles": "CCCCCCCc1cc(c2c(c1)OC(=C3C2C=C(CC3)C)C)O",
        "therapeutic_targets": ["pain", "inflammation", "anxiety"],
        "discovery_note": "Discovered 2019 - Seven-carbon alkyl side chain homolog of CBD",
        "potency_estimate": "Enhanced CB1/CB2 binding vs CBD",
        "research_status": "preclinical"
    },
    {
        "compound_name": "THCP",
        "full_name": "Tetrahydrocannabiphorol",
        "compound_type": "cannabinoid",
        "molecular_formula": "C23H34O2",
        "smiles": "CCCCCCCc1cc(c2c(c1)OC(C3C2C=C(CC3)C)(C)C)O",
        "therapeutic_targets": ["pain", "appetite_stimulation", "sleep"],
        "discovery_note": "Discovered 2019 - 30x more potent at CB1 than THC",
        "potency_estimate": "30x CB1 affinity vs THC",
        "research_status": "preclinical"
    },
    {
        "compound_name": "CBDB",
        "full_name": "Cannabidibutol",
        "compound_type": "cannabinoid",
        "molecular_formula": "C19H26O2",
        "smiles": "CCCCc1cc(c2c(c1)OC(=C3C2C=C(CC3)C)C)O",
        "therapeutic_targets": ["inflammation", "anxiety"],
        "discovery_note": "Four-carbon side chain homolog of CBD",
        "potency_estimate": "Lower CB receptor affinity than CBD",
        "research_status": "preclinical"
    },
    {
        "compound_name": "THCB",
        "full_name": "Tetrahydrocannabutol",
        "compound_type": "cannabinoid",
        "molecular_formula": "C19H26O2",
        "smiles": "CCCCc1cc(c2c(c1)OC(C3C2C=C(CC3)C)(C)C)O",
        "therapeutic_targets": ["pain", "inflammation"],
        "discovery_note": "Four-carbon side chain homolog of THC",
        "potency_estimate": "Similar CB1/CB2 binding to THC",
        "research_status": "preclinical"
    },
    {
        "compound_name": "CBT",
        "full_name": "Cannabitriol",
        "compound_type": "cannabinoid",
        "molecular_formula": "C21H28O3",
        "smiles": "CCCCCc1cc(c2c(c1)OC(C3C2C(C(CC3)(C)O)O)(C)C)O",
        "therapeutic_targets": ["inflammation", "neuroprotection"],
        "discovery_note": "THC metabolite with unique hydroxy groups",
        "potency_estimate": "Moderate CB receptor activity",
        "research_status": "preclinical"
    },
    {
        "compound_name": "CBND",
        "full_name": "Cannabinodiol",
        "compound_type": "cannabinoid",
        "molecular_formula": "C21H26O2",
        "smiles": "CCCCCc1cc2c(c(c1)O)C1C=C(C)CCC1C(C)O2",
        "therapeutic_targets": ["sleep", "inflammation"],
        "discovery_note": "CBD metabolite with sedative properties",
        "potency_estimate": "Low CB1/CB2 activity, unique pharmacology",
        "research_status": "preclinical"
    },
    {
        "compound_name": "CBE",
        "full_name": "Cannabielsoin",
        "compound_type": "cannabinoid",
        "molecular_formula": "C21H30O3",
        "smiles": "CCCCCc1cc(c2c(c1)OC(C3C2C(CC(C3)(C)O)O)(C)C)O",
        "therapeutic_targets": ["inflammation", "analgesic"],
        "discovery_note": "CBD metabolite with additional hydroxyl group",
        "potency_estimate": "Weak CB receptor activity",
        "research_status": "preclinical"
    },
    {
        "compound_name": "11-OH-THC",
        "full_name": "11-Hydroxy-Δ9-tetrahydrocannabinol",
        "compound_type": "cannabinoid_metabolite",
        "molecular_formula": "C21H30O3",
        "smiles": "CCCCCc1cc(c2c(c1)OC(C3C2C(C(CC3)C)O)(C)C)O",
        "therapeutic_targets": ["pain", "psychoactive_effects"],
        "discovery_note": "Primary THC metabolite - more potent than THC",
        "potency_estimate": "Higher CB1 potency than THC, 1.5-2x psychoactive",
        "research_status": "well_characterized"
    },
    {
        "compound_name": "CBGM",
        "full_name": "Cannabigerol Monomethyl Ether",
        "compound_type": "cannabinoid",
        "molecular_formula": "C22H34O2",
        "smiles": "CCCCCC=Cc1cc(cc(c1OC)O)CCCCC",
        "therapeutic_targets": ["inflammation", "antibacterial"],
        "discovery_note": "Methylated CBG variant",
        "potency_estimate": "Similar to CBG with enhanced lipophilicity",
        "research_status": "preclinical"
    },
    {
        "compound_name": "CBGV",
        "full_name": "Cannabigerivarin",
        "compound_type": "cannabinoid",
        "molecular_formula": "C19H28O2",
        "smiles": "CCCC=Cc1cc(cc(c1O)O)CCCC",
        "therapeutic_targets": ["inflammation", "neuroprotection"],
        "discovery_note": "Propyl chain variant of CBG",
        "potency_estimate": "Weaker CB receptor binding than CBG",
        "research_status": "preclinical"
    }
]


def add_placeholder_descriptors(compound: Dict) -> Dict:
    """Add placeholder RDKit descriptors for new compounds."""
    # Note: In production, these would be calculated via RDKit
    # For now, use estimates based on molecular formula
    enriched = compound.copy()
    
    # Extract carbon, hydrogen counts from formula
    formula = compound['molecular_formula']
    c_count = int(formula.split('C')[1].split('H')[0])
    h_count = int(formula.split('H')[1].split('O')[0]) if 'O' in formula else int(formula.split('H')[1])
    
    enriched['rdkit_descriptors'] = {
        "molecular_weight": _estimate_mw(formula),
        "logP": _estimate_logp(c_count, h_count),
        "h_bond_donors": _count_oh_groups(compound['smiles']),
        "h_bond_acceptors": formula.count('O'),
        "rotatable_bonds": _estimate_rotatable(compound['smiles']),
        "lipinski_violations": 1,
        "tpsa": 29.46 + (10.9 * max(0, formula.count('O') - 2)),  # Estimate
        "num_atoms": c_count + h_count + formula.count('O'),
        "num_heavy_atoms": c_count + formula.count('O'),
        "num_aromatic_rings": 1,
        "num_aliphatic_rings": 2,
        "bbb_penetration_predicted": True,
        "veber_rule_pass": True,
        "has_phenol": True,
        "has_pentyl_chain": 'CCCCC' in compound['smiles'],
        "terpenoid_scaffold": True,
        "calculation_note": "Estimated descriptors - requires RDKit validation"
    }
    
    enriched['source'] = 'emerging_research_2024'
    enriched['validation_status'] = 'pending_rdkit_calculation'
    
    return enriched


def _estimate_mw(formula: str) -> float:
    """Estimate molecular weight from formula."""
    mw = 0
    mw += int(formula.split('C')[1].split('H')[0]) * 12.01  # Carbon
    mw += int(formula.split('H')[1].split('O')[0] if 'O' in formula else formula.split('H')[1]) * 1.008  # Hydrogen
    if 'O' in formula:
        mw += int(formula.split('O')[1] if len(formula.split('O')) > 1 and formula.split('O')[1] else '1') * 16.00
    return round(mw, 2)


def _estimate_logp(c_count: int, h_count: int) -> float:
    """Estimate logP based on hydrophobicity."""
    # Rough estimate: cannabinoids are lipophilic
    return round(3.5 + (c_count * 0.15) - (h_count * 0.02), 2)


def _count_oh_groups(smiles: str) -> int:
    """Count OH groups in SMILES."""
    return smiles.count('O)') + smiles.count('OH')


def _estimate_rotatable(smiles: str) -> int:
    """Estimate rotatable bonds."""
    # Count acyclic single bonds (rough)
    rotatable = smiles.count('CC') + smiles.count('c1') - smiles.count('=')
    return min(max(rotatable, 4), 10)


def integrate_new_cannabinoids(dataset: Dict) -> Dict:
    """Add new cannabinoid variants to dataset."""
    enriched = dataset.copy()
    
    # Add descriptors to new compounds
    enriched_cannabinoids = [add_placeholder_descriptors(c) for c in NEW_CANNABINOIDS]
    
    # Append to compounds list
    original_count = len(enriched['compounds'])
    enriched['compounds'].extend(enriched_cannabinoids)
    
    # Update metadata
    enriched['enrichment_metadata'] = enriched.get('enrichment_metadata', {})
    enriched['enrichment_metadata']['cannabinoid_expansion'] = {
        'date': '2025-12-23',
        'compounds_added': len(NEW_CANNABINOIDS),
        'new_total': len(enriched['compounds']),
        'original_total': original_count,
        'new_cannabinoids': [c['compound_name'] for c in NEW_CANNABINOIDS],
        'note': 'Emerging research cannabinoids - descriptors are estimates pending RDKit calculation'
    }
    
    # Update top-level metadata
    if 'metadata' in enriched:
        enriched['metadata']['total_compounds'] = len(enriched['compounds'])
        enriched['metadata']['last_updated'] = '2025-12-23'
    
    return enriched


def main():
    """Main integration workflow."""
    print("=" * 70)
    print("Cannabinoid Compound Expansion")
    print("=" * 70)
    print()
    
    # Input/output paths
    SCRIPT_DIR = Path(__file__).parent
    PROJECT_ROOT = SCRIPT_DIR.parent
    
    input_file = PROJECT_ROOT / 'data' / 'training' / 'neurobotanica_enriched_patients.json'
    output_file = PROJECT_ROOT / 'data' / 'training' / 'neurobotanica_fully_enriched.json'
    
    # Load enriched dataset from step 3
    print(f"Loading patient-enriched dataset: {input_file}")
    with open(input_file, 'r', encoding='utf-8') as f:
        dataset = json.load(f)
    print(f"✅ Loaded {len(dataset['compounds'])} compounds")
    print()
    
    # Add new cannabinoids
    print(f"Adding {len(NEW_CANNABINOIDS)} emerging cannabinoid variants...")
    enriched = integrate_new_cannabinoids(dataset)
    
    print(f"   New compounds:")
    for cannabinoid in NEW_CANNABINOIDS:
        print(f"   • {cannabinoid['compound_name']} - {cannabinoid['discovery_note']}")
    print()
    print(f"   Total compounds: {len(enriched['compounds'])}")
    print()
    
    # Save
    print(f"Saving: {output_file}")
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(enriched, f, indent=2)
    
    print("=" * 70)
    print("✅ Cannabinoid Expansion Complete!")
    print(f"✅ Output: {output_file}")
    print(f"✅ Final compound count: {len(enriched['compounds'])}")
    print()
    print("📊 Full Enrichment Summary:")
    print(f"   Step 1: 368 clinical studies integrated")
    print(f"   Step 2: 1,764 terpene synergies added (Tiers 3-5)")
    print(f"   Step 3: 1,045 synthetic patient profiles")
    print(f"   Step 4: {len(NEW_CANNABINOIDS)} emerging cannabinoids")
    print("=" * 70)


if __name__ == '__main__':
    main()
