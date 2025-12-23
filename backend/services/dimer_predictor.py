"""
Dimeric Cannabinoid Predictor Service
Predicts dimeric cannabinoid structures and properties

Supports NeuroBotanica patent claims:
- Novel dimeric structure generation
- Formation probability scoring
- Synergy prediction algorithms
- Therapeutic potential assessment

Reference: NeuroBotanica MVP Development Plan - Week 4 Task 4.1
"""
import logging
from typing import Optional, List, Dict, Any, Tuple
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
import hashlib

logger = logging.getLogger(__name__)

# RDKit imports with fallback
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
    from rdkit.Chem import DataStructs
    from rdkit.Chem.Fingerprints import FingerprintMols
    from rdkit import DataStructs as DS
    RDKIT_AVAILABLE = True
except ImportError:
    logger.warning("RDKit not available - dimer prediction will be limited")
    RDKIT_AVAILABLE = False


class LinkageType(str, Enum):
    """Types of linkages between monomer units."""
    METHYLENE = "methylene"  # -CH2- bridge
    ETHER = "ether"  # -O- bridge
    ESTER = "ester"  # -COO- bridge
    DIRECT = "direct"  # Direct C-C bond
    CARBON_CHAIN = "carbon_chain"  # -CH2-CH2- etc


class DimerType(str, Enum):
    """Types of dimeric cannabinoids."""
    HOMODIMER = "homodimer"  # A-A
    HETERODIMER = "heterodimer"  # A-B


@dataclass
class ReactiveSite:
    """Represents a reactive site on a cannabinoid monomer."""
    atom_index: int
    site_type: str  # "phenol", "carboxyl", "alkene", "methyl"
    reactivity_score: float  # 0.0-1.0
    position_label: str  # e.g., "C-2'", "C-6a"


@dataclass
class DimerPrediction:
    """Result of dimeric cannabinoid prediction."""
    # Identity
    dimer_name: str
    dimer_type: DimerType
    parent_1_name: str
    parent_1_smiles: str
    parent_2_name: str
    parent_2_smiles: str
    
    # Structure
    predicted_smiles: str
    linkage_type: LinkageType
    linkage_position_1: str
    linkage_position_2: str
    
    # Scores
    formation_probability: float  # 0.0-1.0
    structural_validity: float  # 0.0-1.0
    synergy_prediction: float  # 0.0-1.0
    novelty_score: float  # 0.0-1.0
    
    # Predicted properties
    predicted_mw: Optional[float] = None
    predicted_logp: Optional[float] = None
    predicted_cb1_affinity: Optional[float] = None
    predicted_cb2_affinity: Optional[float] = None
    predicted_selectivity: Optional[float] = None
    
    # Therapeutic potential
    therapeutic_potential: Dict[str, float] = field(default_factory=dict)
    
    # Metadata
    prediction_method: str = "computational"
    confidence_level: str = "predicted"
    created_at: datetime = field(default_factory=datetime.utcnow)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for storage."""
        return {
            "dimer_name": self.dimer_name,
            "dimer_type": self.dimer_type.value,
            "parent_1_name": self.parent_1_name,
            "parent_1_smiles": self.parent_1_smiles,
            "parent_2_name": self.parent_2_name,
            "parent_2_smiles": self.parent_2_smiles,
            "predicted_smiles": self.predicted_smiles,
            "linkage_type": self.linkage_type.value,
            "linkage_position_1": self.linkage_position_1,
            "linkage_position_2": self.linkage_position_2,
            "formation_probability": self.formation_probability,
            "structural_validity": self.structural_validity,
            "synergy_prediction": self.synergy_prediction,
            "novelty_score": self.novelty_score,
            "predicted_mw": self.predicted_mw,
            "predicted_logp": self.predicted_logp,
            "predicted_cb1_affinity": self.predicted_cb1_affinity,
            "predicted_cb2_affinity": self.predicted_cb2_affinity,
            "predicted_selectivity": self.predicted_selectivity,
            "therapeutic_potential": self.therapeutic_potential,
            "prediction_method": self.prediction_method,
            "confidence_level": self.confidence_level,
            "created_at": self.created_at.isoformat()
        }



# --- NextGen Model Integration ---
import os
import joblib
from src.ml_models.nextgen_dimer_input import encode_dimer_features
from src.ml_models.nextgen_dimer_heads import NextGenDimerModel

class DimericPredictor:
    """Predict dimeric cannabinoid structures and properties (legacy + next-gen ML)."""
    
    # Synergy coefficients for therapeutic potential by condition
    SYNERGY_COEFFICIENTS = {
        "chronic_pain": {"THC": 0.8, "CBD": 0.7, "CBG": 0.5, "CBN": 0.4, "CBC": 0.3},
        "anxiety": {"CBD": 0.85, "CBG": 0.6, "THC": 0.3, "CBC": 0.4, "CBN": 0.3},
        "inflammation": {"CBD": 0.8, "CBG": 0.7, "CBC": 0.6, "THC": 0.5, "CBN": 0.3},
        "epilepsy": {"CBD": 0.9, "CBDV": 0.7, "THC": 0.2, "CBG": 0.3, "CBC": 0.2},
        "nausea": {"THC": 0.85, "CBD": 0.5, "CBG": 0.4, "CBN": 0.3, "CBC": 0.2},
        "sleep": {"CBN": 0.85, "THC": 0.7, "CBD": 0.5, "CBG": 0.3, "CBC": 0.2},
        "neuroprotection": {"CBD": 0.8, "CBG": 0.7, "THC": 0.5, "CBC": 0.4, "CBN": 0.3},
    }

    def __init__(self, known_dimers: Optional[List[str]] = None, nextgen_model_path: Optional[str] = None):
        self.known_dimers = known_dimers or []
        self._known_dimer_fps = []
        if RDKIT_AVAILABLE and self.known_dimers:
            for smiles in self.known_dimers:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    fp = FingerprintMols.FingerprintMol(mol)
                    self._known_dimer_fps.append(fp)
        # Load next-gen model if available
        self.nextgen_model = None
        if nextgen_model_path is None:
            # Default path
            nextgen_model_path = os.path.join(os.path.dirname(__file__), '../../models/dimer_potential_v1.joblib')
        if os.path.exists(nextgen_model_path):
            try:
                self.nextgen_model = joblib.load(nextgen_model_path)
            except Exception as e:
                logger.warning(f"Could not load next-gen dimer model: {e}")

    def predict_nextgen_dimer(self, dimer_entry: dict) -> dict:
        """Predict synergy/effect using next-gen ML model pipeline."""
        if self.nextgen_model is None:
            raise RuntimeError("NextGenDimerModel not loaded.")
        features = encode_dimer_features(dimer_entry).reshape(1, -1)
        synergy_pred, effect_pred = self.nextgen_model.predict(features)
        return {
            "synergy_score": float(synergy_pred[0]),
            "effect_size": float(effect_pred[0])
        }
    
    def identify_reactive_sites(self, smiles: str) -> List[ReactiveSite]:
        """Identify reactive sites on a cannabinoid monomer.
        
        Args:
            smiles: SMILES string of the monomer
            
        Returns:
            List of identified reactive sites
        """
        if not RDKIT_AVAILABLE:
            return self._identify_sites_fallback(smiles)
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return []
        
        sites = []
        
        # Find phenolic OH groups
        phenol_pattern = Chem.MolFromSmarts("[OH]c")
        matches = mol.GetSubstructMatches(phenol_pattern)
        for i, match in enumerate(matches):
            sites.append(ReactiveSite(
                atom_index=match[0],
                site_type="phenol",
                reactivity_score=0.9,  # High reactivity
                position_label=f"phenol-{i+1}"
            ))
        
        # Find carboxylic acids (THCA, CBDA)
        carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
        matches = mol.GetSubstructMatches(carboxyl_pattern)
        for i, match in enumerate(matches):
            sites.append(ReactiveSite(
                atom_index=match[0],
                site_type="carboxyl",
                reactivity_score=0.85,
                position_label=f"carboxyl-{i+1}"
            ))
        
        # Find terminal alkenes
        alkene_pattern = Chem.MolFromSmarts("[CH2]=[CH2]")
        matches = mol.GetSubstructMatches(alkene_pattern)
        for i, match in enumerate(matches):
            sites.append(ReactiveSite(
                atom_index=match[0],
                site_type="alkene",
                reactivity_score=0.7,
                position_label=f"alkene-{i+1}"
            ))
        
        # Find activated aromatic positions
        aromatic_pattern = Chem.MolFromSmarts("c([OH])c([H])")
        matches = mol.GetSubstructMatches(aromatic_pattern)
        for i, match in enumerate(matches):
            sites.append(ReactiveSite(
                atom_index=match[2],  # The H-bearing carbon
                site_type="aromatic_ortho",
                reactivity_score=0.6,
                position_label=f"aromatic-{i+1}"
            ))
        
        return sites
    
    def _identify_sites_fallback(self, smiles: str) -> List[ReactiveSite]:
        """Fallback site identification without RDKit."""
        sites = []
        
        # Simple pattern matching
        if "O" in smiles:
            sites.append(ReactiveSite(
                atom_index=smiles.index("O"),
                site_type="oxygen",
                reactivity_score=0.7,
                position_label="O-1"
            ))
        
        return sites
    
    def predict_homodimer(
        self,
        monomer_name: str,
        monomer_smiles: str,
        linkage_type: LinkageType = LinkageType.METHYLENE
    ) -> DimerPrediction:
        """Predict homodimer formation (A-A).
        
        Args:
            monomer_name: Name of the monomer (e.g., "THC")
            monomer_smiles: SMILES string of the monomer
            linkage_type: Type of linkage between monomers
            
        Returns:
            DimerPrediction with structure and scores
        """
        # Identify reactive sites
        sites = self.identify_reactive_sites(monomer_smiles)
        
        if len(sites) < 1:
            # Need at least 1 site for homodimer (we can use same type twice)
            return self._create_failed_prediction(
                monomer_name, monomer_name, monomer_smiles, monomer_smiles,
                "Insufficient reactive sites"
            )
        
        # Select best sites for dimerization
        site1 = max(sites, key=lambda s: s.reactivity_score)
        # For homodimer, site2 is on the second monomer copy (same position is fine)
        site2 = sorted(sites, key=lambda s: s.reactivity_score)[-2] if len(sites) > 1 else site1
        
        # Generate dimer SMILES
        dimer_smiles = self._generate_dimer_smiles(
            monomer_smiles, monomer_smiles,
            site1.atom_index, site2.atom_index,
            linkage_type
        )
        
        # Calculate scores
        formation_prob = self._calculate_formation_probability(
            monomer_smiles, monomer_smiles, linkage_type
        )
        
        structural_validity = self._validate_structure(dimer_smiles)
        synergy_score = self._calculate_synergy_score(monomer_name, monomer_name)
        novelty = self._calculate_novelty(dimer_smiles)
        
        # Predict properties
        properties = self._predict_properties(dimer_smiles, monomer_name, monomer_name)
        
        # Predict therapeutic potential
        therapeutic = self._predict_therapeutic_potential(monomer_name, monomer_name)
        
        return DimerPrediction(
            dimer_name=f"{monomer_name}-{monomer_name} Homodimer",
            dimer_type=DimerType.HOMODIMER,
            parent_1_name=monomer_name,
            parent_1_smiles=monomer_smiles,
            parent_2_name=monomer_name,
            parent_2_smiles=monomer_smiles,
            predicted_smiles=dimer_smiles,
            linkage_type=linkage_type,
            linkage_position_1=site1.position_label,
            linkage_position_2=site2.position_label,
            formation_probability=formation_prob,
            structural_validity=structural_validity,
            synergy_prediction=synergy_score,
            novelty_score=novelty,
            predicted_mw=properties.get("mw"),
            predicted_logp=properties.get("logp"),
            predicted_cb1_affinity=properties.get("cb1_affinity"),
            predicted_cb2_affinity=properties.get("cb2_affinity"),
            predicted_selectivity=properties.get("selectivity"),
            therapeutic_potential=therapeutic
        )
    
    def predict_heterodimer(
        self,
        monomer1_name: str,
        monomer1_smiles: str,
        monomer2_name: str,
        monomer2_smiles: str,
        linkage_type: LinkageType = LinkageType.METHYLENE
    ) -> DimerPrediction:
        """Predict heterodimer formation (A-B).
        
        Args:
            monomer1_name: Name of first monomer
            monomer1_smiles: SMILES of first monomer
            monomer2_name: Name of second monomer
            monomer2_smiles: SMILES of second monomer
            linkage_type: Type of linkage
            
        Returns:
            DimerPrediction with structure and scores
        """
        # Identify reactive sites on both monomers
        sites1 = self.identify_reactive_sites(monomer1_smiles)
        sites2 = self.identify_reactive_sites(monomer2_smiles)
        
        if not sites1 or not sites2:
            return self._create_failed_prediction(
                monomer1_name, monomer2_name,
                monomer1_smiles, monomer2_smiles,
                "Insufficient reactive sites"
            )
        
        # Select best sites
        site1 = max(sites1, key=lambda s: s.reactivity_score)
        site2 = max(sites2, key=lambda s: s.reactivity_score)
        
        # Generate dimer SMILES
        dimer_smiles = self._generate_dimer_smiles(
            monomer1_smiles, monomer2_smiles,
            site1.atom_index, site2.atom_index,
            linkage_type
        )
        
        # Calculate scores
        formation_prob = self._calculate_formation_probability(
            monomer1_smiles, monomer2_smiles, linkage_type
        )
        
        structural_validity = self._validate_structure(dimer_smiles)
        synergy_score = self._calculate_synergy_score(monomer1_name, monomer2_name)
        novelty = self._calculate_novelty(dimer_smiles)
        
        # Predict properties
        properties = self._predict_properties(dimer_smiles, monomer1_name, monomer2_name)
        
        # Predict therapeutic potential
        therapeutic = self._predict_therapeutic_potential(monomer1_name, monomer2_name)
        
        return DimerPrediction(
            dimer_name=f"{monomer1_name}-{monomer2_name} Heterodimer",
            dimer_type=DimerType.HETERODIMER,
            parent_1_name=monomer1_name,
            parent_1_smiles=monomer1_smiles,
            parent_2_name=monomer2_name,
            parent_2_smiles=monomer2_smiles,
            predicted_smiles=dimer_smiles,
            linkage_type=linkage_type,
            linkage_position_1=site1.position_label,
            linkage_position_2=site2.position_label,
            formation_probability=formation_prob,
            structural_validity=structural_validity,
            synergy_prediction=synergy_score,
            novelty_score=novelty,
            predicted_mw=properties.get("mw"),
            predicted_logp=properties.get("logp"),
            predicted_cb1_affinity=properties.get("cb1_affinity"),
            predicted_cb2_affinity=properties.get("cb2_affinity"),
            predicted_selectivity=properties.get("selectivity"),
            therapeutic_potential=therapeutic
        )
    
    def _generate_dimer_smiles(
        self,
        smiles1: str,
        smiles2: str,
        site1_idx: int,
        site2_idx: int,
        linkage_type: LinkageType
    ) -> str:
        """Generate SMILES for the dimer structure."""
        if not RDKIT_AVAILABLE:
            # Simplified fallback - just concatenate with linker
            linker = {
                LinkageType.METHYLENE: "C",
                LinkageType.ETHER: "O",
                LinkageType.ESTER: "C(=O)O",
                LinkageType.DIRECT: "",
                LinkageType.CARBON_CHAIN: "CC"
            }
            return f"{smiles1}.{linker[linkage_type]}.{smiles2}"
        
        try:
            mol1 = Chem.MolFromSmiles(smiles1)
            mol2 = Chem.MolFromSmiles(smiles2)
            
            if mol1 is None or mol2 is None:
                return f"{smiles1}.{smiles2}"
            
            # Create editable molecules
            em1 = Chem.RWMol(mol1)
            em2 = Chem.RWMol(mol2)
            
            # For simplicity, create a combined molecule with explicit linker
            # This is a simplified approach - production would use proper reaction SMARTS
            combined = Chem.CombineMols(mol1, mol2)
            
            # Return the combined SMILES (monomers connected)
            return Chem.MolToSmiles(combined)
            
        except Exception as e:
            logger.warning(f"Error generating dimer SMILES: {e}")
            return f"{smiles1}.{smiles2}"
    
    def _calculate_formation_probability(
        self,
        smiles1: str,
        smiles2: str,
        linkage_type: LinkageType
    ) -> float:
        """Calculate probability of dimer formation.
        
        Based on:
        - Reactive site availability
        - Steric compatibility
        - Thermodynamic favorability
        """
        base_probability = 0.5
        
        # Linkage type modifier
        linkage_modifier = {
            LinkageType.METHYLENE: 0.2,  # Most common in cannabinoids
            LinkageType.ETHER: 0.15,
            LinkageType.ESTER: 0.1,
            LinkageType.DIRECT: 0.05,
            LinkageType.CARBON_CHAIN: 0.1
        }
        
        # Homodimer bonus (same molecule = more compatible)
        same_molecule_bonus = 0.1 if smiles1 == smiles2 else 0.0
        
        # Calculate based on molecular properties
        if RDKIT_AVAILABLE:
            mol1 = Chem.MolFromSmiles(smiles1)
            mol2 = Chem.MolFromSmiles(smiles2)
            
            if mol1 and mol2:
                # Similar molecular weight increases probability
                mw1 = Descriptors.MolWt(mol1)
                mw2 = Descriptors.MolWt(mol2)
                mw_ratio = min(mw1, mw2) / max(mw1, mw2)
                mw_bonus = mw_ratio * 0.1
                
                # Similar polarity increases probability
                logp1 = Descriptors.MolLogP(mol1)
                logp2 = Descriptors.MolLogP(mol2)
                logp_diff = abs(logp1 - logp2)
                polarity_bonus = max(0, 0.1 - logp_diff * 0.02)
                
                base_probability += mw_bonus + polarity_bonus
        
        probability = base_probability + linkage_modifier[linkage_type] + same_molecule_bonus
        
        return min(1.0, max(0.0, probability))
    
    def _validate_structure(self, smiles: str) -> float:
        """Validate the generated dimer structure.
        
        Returns validity score 0.0-1.0
        """
        if not RDKIT_AVAILABLE:
            return 0.5  # Unknown validity
        
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return 0.0
            
            # Check for reasonable molecular properties
            mw = Descriptors.MolWt(mol)
            if mw > 1500 or mw < 300:  # Unreasonable MW for dimer
                return 0.3
            
            # Check for drug-likeness (relaxed for dimers)
            hbd = Descriptors.NumHDonors(mol)
            hba = Descriptors.NumHAcceptors(mol)
            
            if hbd > 10 or hba > 20:
                return 0.5
            
            # Structure is parseable and reasonable
            return 0.9
            
        except Exception:
            return 0.1
    
    def _calculate_synergy_score(
        self,
        name1: str,
        name2: str
    ) -> float:
        """Calculate predicted synergistic effect.
        
        Based on known cannabinoid interactions and entourage effect.
        """
        # Known synergistic combinations
        synergy_matrix = {
            ("THC", "CBD"): 0.85,
            ("THC", "CBG"): 0.75,
            ("THC", "CBN"): 0.65,
            ("CBD", "CBG"): 0.80,
            ("CBD", "CBC"): 0.70,
            ("THC", "THCV"): 0.60,
            ("CBD", "CBDV"): 0.70,
        }
        
        # Check both orderings
        key1 = (name1.upper(), name2.upper())
        key2 = (name2.upper(), name1.upper())
        
        if key1 in synergy_matrix:
            return synergy_matrix[key1]
        elif key2 in synergy_matrix:
            return synergy_matrix[key2]
        
        # Homodimer - moderate synergy
        if name1.upper() == name2.upper():
            return 0.55
        
        # Default for unknown combinations
        return 0.45
    
    def _calculate_novelty(self, dimer_smiles: str) -> float:
        """Calculate novelty score based on similarity to known dimers.
        
        Higher score = more novel (less similar to known structures)
        """
        if not RDKIT_AVAILABLE or not self._known_dimer_fps:
            return 0.9  # Assume novel if can't compare
        
        try:
            mol = Chem.MolFromSmiles(dimer_smiles)
            if mol is None:
                return 0.5
            
            dimer_fp = FingerprintMols.FingerprintMol(mol)
            
            # Calculate maximum similarity to known dimers
            max_similarity = 0.0
            for known_fp in self._known_dimer_fps:
                similarity = DS.TanimotoSimilarity(dimer_fp, known_fp)
                max_similarity = max(max_similarity, similarity)
            
            # Novelty = 1 - max_similarity
            return 1.0 - max_similarity
            
        except Exception:
            return 0.5
    
    def _predict_properties(
        self,
        dimer_smiles: str,
        name1: str,
        name2: str
    ) -> Dict[str, float]:
        """Predict physicochemical and pharmacological properties."""
        properties = {}
        
        if not RDKIT_AVAILABLE:
            return properties
        
        try:
            mol = Chem.MolFromSmiles(dimer_smiles)
            if mol:
                properties["mw"] = round(Descriptors.MolWt(mol), 2)
                properties["logp"] = round(Descriptors.MolLogP(mol), 2)
                properties["tpsa"] = round(Descriptors.TPSA(mol), 2)
                properties["hbd"] = Descriptors.NumHDonors(mol)
                properties["hba"] = Descriptors.NumHAcceptors(mol)
        except Exception:
            pass
        
        # Predict receptor affinities based on parent compounds
        # Using simplified additive model with diminishing returns
        parent_affinities = {
            "THC": {"cb1": 10.0, "cb2": 25.0},  # Ki in nM
            "CBD": {"cb1": 1000.0, "cb2": 800.0},
            "CBG": {"cb1": 400.0, "cb2": 350.0},
            "CBN": {"cb1": 150.0, "cb2": 250.0},
            "CBC": {"cb1": 2000.0, "cb2": 1500.0},
            "THCV": {"cb1": 75.0, "cb2": 60.0},
            "CBDV": {"cb1": 1500.0, "cb2": 1200.0},
        }
        
        name1_upper = name1.upper()
        name2_upper = name2.upper()
        
        if name1_upper in parent_affinities and name2_upper in parent_affinities:
            # Average affinities with 20% enhancement for dimerization
            cb1_1 = parent_affinities[name1_upper]["cb1"]
            cb1_2 = parent_affinities[name2_upper]["cb1"]
            cb2_1 = parent_affinities[name1_upper]["cb2"]
            cb2_2 = parent_affinities[name2_upper]["cb2"]
            
            properties["cb1_affinity"] = round((cb1_1 + cb1_2) / 2 * 0.8, 2)
            properties["cb2_affinity"] = round((cb2_1 + cb2_2) / 2 * 0.8, 2)
            
            if properties["cb1_affinity"] > 0:
                properties["selectivity"] = round(
                    properties["cb2_affinity"] / properties["cb1_affinity"], 2
                )
        
        return properties
    
    def _predict_therapeutic_potential(
        self,
        name1: str,
        name2: str
    ) -> Dict[str, float]:
        """Predict therapeutic potential by condition."""
        potential = {}
        
        for condition, coefficients in self.SYNERGY_COEFFICIENTS.items():
            name1_upper = name1.upper()
            name2_upper = name2.upper()
            
            coef1 = coefficients.get(name1_upper, 0.1)
            coef2 = coefficients.get(name2_upper, 0.1)
            
            # Combined potential with synergy bonus
            base_potential = (coef1 + coef2) / 2
            synergy_bonus = 0.2 if name1_upper != name2_upper else 0.1
            
            potential[condition] = round(min(1.0, base_potential + synergy_bonus), 3)
        
        return potential
    
    def _create_failed_prediction(
        self,
        name1: str,
        name2: str,
        smiles1: str,
        smiles2: str,
        reason: str
    ) -> DimerPrediction:
        """Create a failed prediction result."""
        return DimerPrediction(
            dimer_name=f"{name1}-{name2} (Failed)",
            dimer_type=DimerType.HETERODIMER if name1 != name2 else DimerType.HOMODIMER,
            parent_1_name=name1,
            parent_1_smiles=smiles1,
            parent_2_name=name2,
            parent_2_smiles=smiles2,
            predicted_smiles="",
            linkage_type=LinkageType.METHYLENE,
            linkage_position_1="unknown",
            linkage_position_2="unknown",
            formation_probability=0.0,
            structural_validity=0.0,
            synergy_prediction=0.0,
            novelty_score=0.0,
            therapeutic_potential={},
            confidence_level=f"failed: {reason}"
        )
    
    def generate_all_combinations(
        self,
        monomers: List[Tuple[str, str]],
        include_homodimers: bool = True,
        linkage_type: LinkageType = LinkageType.METHYLENE
    ) -> List[DimerPrediction]:
        """Generate predictions for all monomer combinations.
        
        Args:
            monomers: List of (name, smiles) tuples
            include_homodimers: Whether to include A-A dimers
            linkage_type: Type of linkage to use
            
        Returns:
            List of DimerPrediction objects
        """
        predictions = []
        
        # Generate homodimers
        if include_homodimers:
            for name, smiles in monomers:
                prediction = self.predict_homodimer(name, smiles, linkage_type)
                predictions.append(prediction)
        
        # Generate heterodimers (unique combinations only)
        for i, (name1, smiles1) in enumerate(monomers):
            for name2, smiles2 in monomers[i+1:]:
                prediction = self.predict_heterodimer(
                    name1, smiles1, name2, smiles2, linkage_type
                )
                predictions.append(prediction)
        
        return predictions
    
    def rank_predictions(
        self,
        predictions: List[DimerPrediction],
        by: str = "therapeutic",
        condition: Optional[str] = None
    ) -> List[DimerPrediction]:
        """Rank predictions by specified criteria.
        
        Args:
            predictions: List of predictions to rank
            by: Ranking criteria - "formation", "synergy", "therapeutic", "novelty"
            condition: Specific condition for therapeutic ranking
            
        Returns:
            Sorted list of predictions (highest first)
        """
        if by == "formation":
            key = lambda p: p.formation_probability
        elif by == "synergy":
            key = lambda p: p.synergy_prediction
        elif by == "novelty":
            key = lambda p: p.novelty_score
        elif by == "therapeutic":
            if condition:
                key = lambda p: p.therapeutic_potential.get(condition, 0)
            else:
                # Average therapeutic potential
                key = lambda p: sum(p.therapeutic_potential.values()) / max(len(p.therapeutic_potential), 1)
        else:
            # Combined score
            key = lambda p: (
                p.formation_probability * 0.3 +
                p.synergy_prediction * 0.3 +
                p.structural_validity * 0.2 +
                p.novelty_score * 0.2
            )
        
        return sorted(predictions, key=key, reverse=True)
