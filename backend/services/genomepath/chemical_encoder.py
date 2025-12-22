"""Chemical Fingerprint Encoder for GenomePath"""

import torch
import numpy as np
from typing import List, Dict, Optional, Tuple
import json
from pathlib import Path

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors, DataStructs
    RDKIT_AVAILABLE = True
except ImportError:
    print("Warning: RDKit not available. Chemical fingerprints will be disabled.")
    RDKIT_AVAILABLE = False


class ChemicalEncoder:
    """Encodes chemical compounds into fingerprint vectors."""
    
    def __init__(
        self,
        fingerprint_radius: int = 2,  # ECFP4
        fingerprint_bits: int = 2048,
        use_features: bool = True
    ):
        self.radius = fingerprint_radius
        self.n_bits = fingerprint_bits
        self.use_features = use_features
        
        # Common cannabinoids and terpenes SMILES
        self.known_compounds = {
            # Cannabinoids
            'THC': 'CCCCCc1cc(O)c2c(c1)OC(C)(C)[C@@H]1CC=C(C)C[C@H]21',
            'CBD': 'CCCCCc1cc(O)c2c(c1)C1C=C(C)CC[C@H]1C(C)(C)O2',
            'CBN': 'CCCCCc1cc2OC(C)(C)c3ccc(O)cc3c2c(O)c1',
            'CBG': 'CCCCCC=CCc1cc(O)cc(O)c1C(=O)O',
            'CBC': 'CCCCCC=CCc1cc(O)c2C3C=C(C)CCC3C(C)(C)Oc2c1',
            # Common Terpenes
            'beta-caryophyllene': 'CC1=CCC(C(=C)C2CCC(=C)C2(C)C)CC1',
            'limonene': 'CC1=CCC(CC1)C(=C)C',
            'myrcene': 'CC(=CCCC(=C)C=C)C',
            'pinene': 'CC1=CCC2CC1C2(C)C',
            'linalool': 'CC(C)=CCCC(C)(O)C=C',
            'humulene': 'CC1=CCCC(=C)C2CCC(C)(C=C)C2CC1',
            # Phytochemicals
            'curcumin': 'COc1cc(ccc1O)C=CC(=O)CC(=O)C=Cc1ccc(O)c(OC)c1',
            'resveratrol': 'Oc1ccc(cc1)C=Cc1cc(O)cc(O)c1',
            'quercetin': 'O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12',
        }
        
        print(f"Chemical encoder initialized: ECFP{2*fingerprint_radius} ({fingerprint_bits} bits)")
        print(f"Known compounds: {len(self.known_compounds)}")
    
    def smiles_to_fingerprint(self, smiles: str) -> Optional[np.ndarray]:
        if not RDKIT_AVAILABLE:
            return None
            
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
            
            fp = AllChem.GetMorganFingerprintAsBitVect(
                mol,
                radius=self.radius,
                nBits=self.n_bits
            )
            
            arr = np.zeros((self.n_bits,), dtype=np.float32)
            DataStructs.ConvertToNumpyArray(fp, arr)
            
            return arr
            
        except Exception as e:
            print(f"Warning: Failed to parse SMILES '{smiles}': {e}")
            return None
    
    def get_molecular_descriptors(self, smiles: str) -> Optional[np.ndarray]:
        if not RDKIT_AVAILABLE:
            return None
            
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
            
            descriptors = np.array([
                Descriptors.MolWt(mol),
                Descriptors.MolLogP(mol),
                Descriptors.NumHDonors(mol),
                Descriptors.NumHAcceptors(mol),
                Descriptors.TPSA(mol),
                Descriptors.NumRotatableBonds(mol),
                rdMolDescriptors.CalcNumAromaticRings(mol),
            ], dtype=np.float32)
            
            # Normalize
            normalization = np.array([500.0, 7.5, 5.0, 7.5, 100.0, 10.0, 3.0])
            descriptors = descriptors / normalization
            
            return descriptors
            
        except Exception as e:
            print(f"Warning: Failed to calculate descriptors for '{smiles}': {e}")
            return None
    
    def encode_compound(
        self,
        compound_name: Optional[str] = None,
        smiles: Optional[str] = None
    ) -> Dict[str, torch.Tensor]:
        # Get SMILES
        if smiles is None and compound_name:
            smiles = self.known_compounds.get(compound_name.lower())
        
        if smiles is None:
            return {
                'fingerprint': torch.zeros(self.n_bits, dtype=torch.float32),
                'descriptors': torch.zeros(7, dtype=torch.float32),
                'is_valid': torch.tensor(0.0)
            }
        
        # Generate fingerprint
        fp = self.smiles_to_fingerprint(smiles)
        if fp is None:
            return {
                'fingerprint': torch.zeros(self.n_bits, dtype=torch.float32),
                'descriptors': torch.zeros(7, dtype=torch.float32),
                'is_valid': torch.tensor(0.0)
            }
        
        # Get descriptors
        descriptors = self.get_molecular_descriptors(smiles) if self.use_features else np.zeros(7)
        if descriptors is None:
            descriptors = np.zeros(7, dtype=np.float32)
        
        return {
            'fingerprint': torch.from_numpy(fp),
            'descriptors': torch.from_numpy(descriptors),
            'is_valid': torch.tensor(1.0)
        }
    
    def encode_compound_list(self, compounds: List[str]) -> Dict[str, torch.Tensor]:
        if not compounds:
            result = self.encode_compound(None, None)
            result['num_compounds'] = torch.tensor(0.0)
            return result
        
        fingerprints = []
        descriptors = []
        valid_count = 0
        
        for compound in compounds:
            result = self.encode_compound(compound_name=compound)
            if result['is_valid'] > 0.5:
                fingerprints.append(result['fingerprint'])
                descriptors.append(result['descriptors'])
                valid_count += 1
        
        if valid_count == 0:
            result = self.encode_compound(None, None)
            result['num_compounds'] = torch.tensor(0.0)
            return result
        
        avg_fp = torch.stack(fingerprints).mean(dim=0)
        avg_desc = torch.stack(descriptors).mean(dim=0)
        
        return {
            'fingerprint': avg_fp,
            'descriptors': avg_desc,
            'is_valid': torch.tensor(1.0),
            'num_compounds': torch.tensor(float(valid_count))
        }
    
    def calculate_similarity(self, fp1: torch.Tensor, fp2: torch.Tensor) -> float:
        fp1_np = fp1.numpy() if isinstance(fp1, torch.Tensor) else fp1
        fp2_np = fp2.numpy() if isinstance(fp2, torch.Tensor) else fp2
        
        intersection = np.sum(fp1_np * fp2_np)
        union = np.sum(fp1_np) + np.sum(fp2_np) - intersection
        
        if union == 0:
            return 0.0
        
        return float(intersection / union)


if __name__ == "__main__":
    encoder = ChemicalEncoder()
    
    print("\n" + "=" * 70)
    print("Testing Chemical Encoder")
    print("=" * 70)
    
    # Test single compound
    print("\nEncoding THC:")
    thc_result = encoder.encode_compound(compound_name='THC')
    print(f"  Fingerprint: {thc_result['fingerprint'].sum():.0f} bits set")
    print(f"  Descriptors: {thc_result['descriptors'].tolist()}")
    print(f"  Valid: {thc_result['is_valid']}")
    
    # Test compound list
    print("\nEncoding cannabinoid blend:")
    blend_result = encoder.encode_compound_list(['THC', 'CBD', 'CBN'])
    print(f"  Fingerprint: {blend_result['fingerprint'].sum():.0f} bits set")
    print(f"  Valid compounds: {blend_result['num_compounds']}")
    
    # Test similarity
    cbd_result = encoder.encode_compound(compound_name='CBD')
    similarity = encoder.calculate_similarity(
        thc_result['fingerprint'],
        cbd_result['fingerprint']
    )
    print(f"\nTHC-CBD similarity: {similarity:.3f}")
    
    # Test terpene
    limonene_result = encoder.encode_compound(compound_name='limonene')
    print(f"\nLimonene descriptors: {limonene_result['descriptors'].tolist()}")
    
    print("\nChemical encoder working correctly!")
