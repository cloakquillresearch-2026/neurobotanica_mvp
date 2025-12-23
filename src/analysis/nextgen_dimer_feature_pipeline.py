import numpy as np
from typing import List, Dict, Any
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

# --- Input Module: Compound Embeddings ---
def smiles_to_ecfp4(smiles: str, n_bits: int = 2048) -> np.ndarray:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return np.zeros(n_bits)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=n_bits)
    arr = np.zeros((n_bits,), dtype=int)
    AllChem.DataStructs.ConvertToNumpyArray(fp, arr)
    return arr

def get_mol_descriptors(smiles: str) -> Dict[str, float]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {k: 0.0 for k in [
            'MolWt', 'MolLogP', 'NumHDonors', 'NumHAcceptors', 'TPSA', 'NumRotatableBonds', 'NumAromaticRings']}
    return {
        'MolWt': Descriptors.MolWt(mol),
        'MolLogP': Descriptors.MolLogP(mol),
        'NumHDonors': Descriptors.NumHDonors(mol),
        'NumHAcceptors': Descriptors.NumHAcceptors(mol),
        'TPSA': Descriptors.TPSA(mol),
        'NumRotatableBonds': Descriptors.NumRotatableBonds(mol),
        'NumAromaticRings': Descriptors.NumAromaticRings(mol)
    }

# --- Input Module: Contextual Features ---
def extract_context_features(entry: Dict[str, Any]) -> Dict[str, Any]:
    return {
        'synergy_score': entry.get('synergy_score', 0.0),
        'ai_predicted': float(entry.get('ai_predicted', False)),
        'benefit_risk_ratio': entry.get('benefit_risk_ratio', 1.0),
        'population_group': entry.get('population_group', ''),
        'regulatory_status': entry.get('regulatory_status', ''),
        'effect_size': entry.get('effect_size', ''),
        'confidence': entry.get('confidence', ''),
        'omics_signature': entry.get('omics_signature', ''),
        'pathway': entry.get('pathway', ''),
        'patient_stratification': entry.get('patient_stratification', ''),
    }

# --- Pipeline: Feature Extraction for Model ---
def extract_features_for_model(entries: List[Dict[str, Any]], smiles_lookup: Dict[str, str]) -> List[Dict[str, Any]]:
    features = []
    for entry in entries:
        c1 = entry['compound_1']
        c2 = entry['compound_2']
        smiles1 = smiles_lookup.get(c1, '')
        smiles2 = smiles_lookup.get(c2, '')
        fp1 = smiles_to_ecfp4(smiles1)
        fp2 = smiles_to_ecfp4(smiles2)
        desc1 = get_mol_descriptors(smiles1)
        desc2 = get_mol_descriptors(smiles2)
        context = extract_context_features(entry)
        features.append({
            'compound_1': c1,
            'compound_2': c2,
            'fp1': fp1,
            'fp2': fp2,
            'desc1': desc1,
            'desc2': desc2,
            **context
        })
    return features

# Example usage:
# from src.analysis.nextgen_dimer_ingest import load_nextgen_dimer_evidence
# entries = [d.__dict__ for d in load_nextgen_dimer_evidence('data/processed/nextgen_dimers.json')]
# smiles_lookup = {'THC': 'CCCCCc1cc(c2c(c1)OC(C3C2C=C(CC3)C)(C)C)O', ...}
# features = extract_features_for_model(entries, smiles_lookup)
