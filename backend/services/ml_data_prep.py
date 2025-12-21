"""
ML Data Preparation - NeuroBotanica Week 6
Prepare training data with 2D and 3D molecular features for ML models.

Supports:
- 2D RDKit descriptors (40+)
- 3D descriptors (PMI, asphericity, etc.)
- Receptor affinity features with provenance weighting
- Clinical evidence integration
- Confidence-weighted training samples
"""
import pandas as pd
import numpy as np
from typing import List, Dict, Optional, Tuple, Any
from dataclasses import dataclass, field
from datetime import datetime
import logging
import json

logger = logging.getLogger(__name__)


@dataclass
class FeatureSet:
    """Complete feature set for ML training."""
    features: pd.DataFrame
    feature_names: List[str]
    target_column: str
    weight_column: Optional[str]
    metadata: Dict[str, Any] = field(default_factory=dict)
    
    def to_dict(self) -> Dict:
        return {
            "num_samples": len(self.features),
            "num_features": len(self.feature_names),
            "target_column": self.target_column,
            "has_weights": self.weight_column is not None,
            "feature_names": self.feature_names[:20],  # First 20
            "metadata": self.metadata
        }


@dataclass
class DatasetStats:
    """Statistics about prepared dataset."""
    total_samples: int
    num_features: int
    num_2d_features: int
    num_3d_features: int
    num_receptor_features: int
    num_conditions: int
    condition_distribution: Dict[str, int]
    missing_3d_count: int
    avg_confidence_weight: float
    
    def to_dict(self) -> Dict:
        return {
            "total_samples": self.total_samples,
            "num_features": self.num_features,
            "num_2d_features": self.num_2d_features,
            "num_3d_features": self.num_3d_features,
            "num_receptor_features": self.num_receptor_features,
            "num_conditions": self.num_conditions,
            "condition_distribution": self.condition_distribution,
            "missing_3d_count": self.missing_3d_count,
            "avg_confidence_weight": round(self.avg_confidence_weight, 3)
        }


# Standard 2D RDKit descriptors to use
STANDARD_2D_DESCRIPTORS = [
    "molecular_weight", "logp", "tpsa", "h_bond_donors", "h_bond_acceptors",
    "rotatable_bonds", "aromatic_rings", "heavy_atoms", "fraction_csp3",
    "num_aliphatic_rings", "num_saturated_rings", "num_heteroatoms",
    "num_amide_bonds", "num_bridgehead_atoms", "molar_refractivity",
    "exact_mass", "formal_charge", "num_radical_electrons", "num_valence_electrons"
]

# Standard 3D descriptors from conformer generation
STANDARD_3D_DESCRIPTORS = [
    "pmi1", "pmi2", "pmi3",  # Principal moments of inertia
    "npr1", "npr2",  # Normalized PMI ratios
    "asphericity", "eccentricity", "inertial_shape_factor",
    "radius_of_gyration", "spherocity_index",
    "lowest_energy", "energy_range", "num_conformers"
]

# Receptor targets for affinity features
RECEPTOR_TARGETS = [
    "CB1", "CB2", "GPR55", "TRPV1", "TRPA1",
    "5-HT1A", "PPARα", "PPARγ", "FAAH", "MAGL"
]


class MLDataPreparator:
    """Prepare training data with 2D and 3D features for ML models.
    
    Handles:
    - Feature extraction from cannabinoid compounds
    - Clinical study outcome parsing
    - Evidence-weighted sample generation
    - Missing data imputation
    """
    
    def __init__(
        self,
        include_3d: bool = True,
        min_evidence_tier: int = 3,
        impute_missing: bool = True
    ):
        """Initialize data preparator.
        
        Args:
            include_3d: Whether to include 3D conformer descriptors
            min_evidence_tier: Minimum evidence tier for clinical studies (1-5, 1=highest)
            impute_missing: Whether to impute missing values
        """
        self.include_3d = include_3d
        self.min_evidence_tier = min_evidence_tier
        self.impute_missing = impute_missing
        
        self._2d_features = STANDARD_2D_DESCRIPTORS.copy()
        self._3d_features = STANDARD_3D_DESCRIPTORS.copy() if include_3d else []
        self._receptor_features = [f"{r}_affinity" for r in RECEPTOR_TARGETS]
        
        logger.info(f"MLDataPreparator initialized (include_3d={include_3d})")
    
    def prepare_therapeutic_dataset(
        self,
        cannabinoids: List[Dict],
        studies: List[Dict]
    ) -> Tuple[FeatureSet, DatasetStats]:
        """Prepare dataset for therapeutic prediction model.
        
        Args:
            cannabinoids: List of cannabinoid compound dicts
            studies: List of clinical study dicts
            
        Returns:
            Tuple of (FeatureSet, DatasetStats)
        """
        logger.info(f"Preparing therapeutic dataset: {len(cannabinoids)} compounds, {len(studies)} studies")
        
        data_rows = []
        condition_counts = {}
        missing_3d = 0
        
        # Create lookup for cannabinoids
        cb_lookup = {cb.get("name", "").upper(): cb for cb in cannabinoids}
        cb_lookup.update({cb.get("abbreviation", "").upper(): cb for cb in cannabinoids if cb.get("abbreviation")})
        
        for study in studies:
            # Skip low-evidence studies
            evidence_tier = study.get("evidence_tier", 3)
            if evidence_tier > self.min_evidence_tier:
                continue
            
            # Get cannabinoid(s) studied
            cannabinoid_name = study.get("cannabinoid", "")
            compound = cb_lookup.get(cannabinoid_name.upper())
            
            if not compound:
                continue
            
            # Extract features
            features = self._extract_compound_features(compound)
            
            # Track 3D availability
            if not compound.get("has_conformers", False):
                missing_3d += 1
            
            # Parse study outcomes
            efficacy = self._parse_efficacy(study)
            condition = study.get("condition", "unknown")
            
            # Track condition distribution
            condition_counts[condition] = condition_counts.get(condition, 0) + 1
            
            # Calculate confidence weight
            confidence_weight = self._calculate_confidence_weight(study)
            
            # Build row
            row = features.copy()
            row["condition"] = condition
            row["efficacy"] = efficacy
            row["confidence_weight"] = confidence_weight
            row["study_id"] = study.get("study_id", "")
            row["compound_name"] = compound.get("name", "")
            
            data_rows.append(row)
        
        # Create DataFrame
        df = pd.DataFrame(data_rows)
        
        # Impute missing values if enabled
        if self.impute_missing and len(df) > 0:
            df = self._impute_missing_values(df)
        
        # Get feature column names
        feature_cols = self._get_feature_columns(df)
        
        # Create stats
        stats = DatasetStats(
            total_samples=len(df),
            num_features=len(feature_cols),
            num_2d_features=len([f for f in feature_cols if f in self._2d_features]),
            num_3d_features=len([f for f in feature_cols if f in self._3d_features]),
            num_receptor_features=len([f for f in feature_cols if f in self._receptor_features]),
            num_conditions=len(condition_counts),
            condition_distribution=condition_counts,
            missing_3d_count=missing_3d,
            avg_confidence_weight=df["confidence_weight"].mean() if len(df) > 0 else 0.0
        )
        
        feature_set = FeatureSet(
            features=df,
            feature_names=feature_cols,
            target_column="efficacy",
            weight_column="confidence_weight",
            metadata={
                "created_at": datetime.now().isoformat(),
                "include_3d": self.include_3d,
                "min_evidence_tier": self.min_evidence_tier
            }
        )
        
        logger.info(f"Dataset prepared: {stats.total_samples} samples, {stats.num_features} features")
        return feature_set, stats
    
    def prepare_dimer_dataset(
        self,
        dimer_predictions: List[Dict],
        cannabinoids: List[Dict]
    ) -> Tuple[FeatureSet, DatasetStats]:
        """Prepare dataset for dimer therapeutic potential model.
        
        Args:
            dimer_predictions: List of dimer prediction dicts
            cannabinoids: List of parent cannabinoid compounds
            
        Returns:
            Tuple of (FeatureSet, DatasetStats)
        """
        logger.info(f"Preparing dimer dataset: {len(dimer_predictions)} dimers")
        
        data_rows = []
        cb_lookup = {cb.get("name", "").upper(): cb for cb in cannabinoids}
        
        for dimer in dimer_predictions:
            # Get parent compound features
            parent1_name = dimer.get("monomer_a", dimer.get("parent_1_name", ""))
            parent2_name = dimer.get("monomer_b", dimer.get("parent_2_name", ""))
            
            parent1 = cb_lookup.get(parent1_name.upper(), {})
            parent2 = cb_lookup.get(parent2_name.upper(), {})
            
            # Extract combined features
            features = self._extract_dimer_features(dimer, parent1, parent2)
            
            # Use triangulation score as target
            therapeutic_potential = dimer.get("triangulation_score", 0.5)
            
            row = features.copy()
            row["therapeutic_potential"] = therapeutic_potential
            row["dimer_name"] = dimer.get("dimer_name", f"{parent1_name}-{parent2_name}")
            row["is_homodimer"] = 1 if parent1_name == parent2_name else 0
            
            data_rows.append(row)
        
        df = pd.DataFrame(data_rows)
        
        if self.impute_missing and len(df) > 0:
            df = self._impute_missing_values(df)
        
        feature_cols = self._get_feature_columns(df)
        
        stats = DatasetStats(
            total_samples=len(df),
            num_features=len(feature_cols),
            num_2d_features=len([f for f in feature_cols if "_mw" in f or "_logp" in f]),
            num_3d_features=len([f for f in feature_cols if "pmi" in f or "energy" in f]),
            num_receptor_features=len([f for f in feature_cols if "_affinity" in f]),
            num_conditions=0,
            condition_distribution={},
            missing_3d_count=0,
            avg_confidence_weight=1.0
        )
        
        feature_set = FeatureSet(
            features=df,
            feature_names=feature_cols,
            target_column="therapeutic_potential",
            weight_column=None,
            metadata={
                "created_at": datetime.now().isoformat(),
                "dataset_type": "dimer_therapeutic_potential"
            }
        )
        
        return feature_set, stats
    
    def prepare_patient_response_dataset(
        self,
        cannabinoids: List[Dict],
        studies: List[Dict]
    ) -> Tuple[FeatureSet, DatasetStats]:
        """Prepare dataset for patient treatment response model.
        
        Focuses on response probability prediction.
        """
        logger.info("Preparing patient response dataset")
        
        # Similar to therapeutic dataset but with response as binary outcome
        feature_set, stats = self.prepare_therapeutic_dataset(cannabinoids, studies)
        
        # Convert efficacy to binary response
        if len(feature_set.features) > 0:
            feature_set.features["response"] = (
                feature_set.features["efficacy"] > 0.5
            ).astype(int)
            feature_set.target_column = "response"
        
        return feature_set, stats
    
    def _extract_compound_features(self, compound: Dict) -> Dict[str, float]:
        """Extract all features from a compound."""
        features = {}
        
        # 2D descriptors
        rdkit_2d = compound.get("rdkit_descriptors", {})
        for feat in self._2d_features:
            value = rdkit_2d.get(feat)
            if value is None:
                # Try direct compound attributes
                value = compound.get(feat)
            features[feat] = value if value is not None else np.nan
        
        # 3D descriptors (if available)
        if self.include_3d:
            rdkit_3d = compound.get("rdkit_descriptors_3d", {})
            for feat in self._3d_features:
                features[feat] = rdkit_3d.get(feat, np.nan)
            features["has_3d"] = 1.0 if compound.get("has_conformers", False) else 0.0
        
        # Receptor affinities
        affinities = compound.get("receptor_affinities", {})
        for receptor in RECEPTOR_TARGETS:
            receptor_data = affinities.get(receptor, [])
            if isinstance(receptor_data, list) and len(receptor_data) > 0:
                # Average all measurements
                values = [m.get("affinity_value", 0) for m in receptor_data if m.get("affinity_value")]
                features[f"{receptor}_affinity"] = np.mean(values) if values else np.nan
            elif isinstance(receptor_data, (int, float)):
                features[f"{receptor}_affinity"] = float(receptor_data)
            else:
                features[f"{receptor}_affinity"] = np.nan
        
        return features
    
    def _extract_dimer_features(
        self,
        dimer: Dict,
        parent1: Dict,
        parent2: Dict
    ) -> Dict[str, float]:
        """Extract features from dimer and parent compounds."""
        features = {}
        
        # Dimer-specific features
        features["synergy_score"] = dimer.get("synergy_score", 0.5)
        features["formation_probability"] = dimer.get("formation_probability", 0.5)
        features["novelty_score"] = dimer.get("novelty_score", 0.5)
        features["structural_validity"] = dimer.get("structural_validity", 0.5)
        
        # Predicted properties
        pred_props = dimer.get("predicted_properties", {})
        features["dimer_mw"] = pred_props.get("mw", np.nan)
        features["dimer_logp"] = pred_props.get("logp", np.nan)
        features["dimer_tpsa"] = pred_props.get("tpsa", np.nan)
        
        # Parent 1 features (prefixed)
        p1_features = self._extract_compound_features(parent1)
        for key, value in p1_features.items():
            features[f"p1_{key}"] = value
        
        # Parent 2 features (prefixed)
        p2_features = self._extract_compound_features(parent2)
        for key, value in p2_features.items():
            features[f"p2_{key}"] = value
        
        # Combined features
        if not np.isnan(features.get("p1_molecular_weight", np.nan)) and not np.isnan(features.get("p2_molecular_weight", np.nan)):
            features["mw_ratio"] = features["p1_molecular_weight"] / max(features["p2_molecular_weight"], 0.001)
            features["mw_sum"] = features["p1_molecular_weight"] + features["p2_molecular_weight"]
        
        return features
    
    def _parse_efficacy(self, study: Dict) -> float:
        """Parse efficacy score from study data."""
        # Try numeric efficacy first
        efficacy = study.get("efficacy", study.get("effect_size_numeric"))
        if isinstance(efficacy, (int, float)):
            return float(np.clip(efficacy, 0.0, 1.0))
        
        # Parse from text summary
        summary = study.get("results_summary", study.get("efficacy_summary", "")).lower()
        
        if "significant improvement" in summary or "highly effective" in summary:
            return 0.85
        elif "significant" in summary:
            return 0.75
        elif "moderate" in summary:
            return 0.55
        elif "mild" in summary or "slight" in summary:
            return 0.35
        elif "no effect" in summary or "ineffective" in summary:
            return 0.15
        else:
            return 0.5  # Default neutral
    
    def _calculate_confidence_weight(self, study: Dict) -> float:
        """Calculate confidence weight for a study."""
        weight = 1.0
        
        # Evidence tier adjustment (1=best, 5=weakest)
        tier = study.get("evidence_tier", 3)
        tier_weights = {1: 1.5, 2: 1.2, 3: 1.0, 4: 0.7, 5: 0.4}
        weight *= tier_weights.get(tier, 1.0)
        
        # Sample size adjustment
        sample_size = study.get("sample_size", 50)
        if sample_size > 200:
            weight *= 1.3
        elif sample_size > 100:
            weight *= 1.1
        elif sample_size < 30:
            weight *= 0.7
        
        # Confidence score from study
        confidence = study.get("confidence_weight", study.get("confidence_score", 1.0))
        if isinstance(confidence, (int, float)):
            weight *= confidence
        
        return np.clip(weight, 0.1, 2.0)
    
    def _impute_missing_values(self, df: pd.DataFrame) -> pd.DataFrame:
        """Impute missing values in dataframe."""
        numeric_cols = df.select_dtypes(include=[np.number]).columns
        
        for col in numeric_cols:
            if df[col].isna().any():
                # Use median for imputation
                median_val = df[col].median()
                if pd.isna(median_val):
                    median_val = 0.0
                df[col] = df[col].fillna(median_val)
        
        return df
    
    def _get_feature_columns(self, df: pd.DataFrame) -> List[str]:
        """Get list of feature columns (excluding target and metadata)."""
        exclude_cols = {
            "efficacy", "therapeutic_potential", "response",
            "confidence_weight", "study_id", "compound_name",
            "condition", "dimer_name", "is_homodimer"
        }
        
        return [
            col for col in df.columns
            if col not in exclude_cols and df[col].dtype in [np.float64, np.int64, np.float32, np.int32]
        ]
    
    def get_feature_importance_groups(self) -> Dict[str, List[str]]:
        """Get feature groups for interpretability."""
        return {
            "2d_structural": self._2d_features,
            "3d_conformational": self._3d_features,
            "receptor_affinity": self._receptor_features,
            "experimental": ["has_3d", "confidence_weight"]
        }
