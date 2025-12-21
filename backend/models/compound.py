"""
Cannabinoid Compound Model
Schema for 63-compound dataset with dimeric predictions

Supports NeuroBotanica patent claims:
- 3D conformer data for molecular analysis
- Receptor affinity with provenance
- Dimeric prediction scoring
"""
from sqlalchemy import Column, Integer, String, Float, Boolean, JSON, Text, DateTime
from sqlalchemy.sql import func
from .database import Base


class Cannabinoid(Base):
    """Cannabinoid compound model with 2D/3D molecular descriptors.
    
    Supports:
    - Dimeric triangulation scoring
    - ChEMBL/PubChem integration
    - Receptor affinity prediction
    """
    __tablename__ = "cannabinoids"
    
    id = Column(Integer, primary_key=True, index=True)
    
    # Identity
    name = Column(String(100), nullable=False, unique=True, index=True)
    abbreviation = Column(String(20), index=True)
    smiles = Column(String(500), nullable=False)
    inchi = Column(String(500))
    inchi_key = Column(String(100), index=True)
    
    # Classification
    compound_class = Column(String(50))  # Phytocannabinoid, Endocannabinoid, Synthetic
    is_natural = Column(Boolean, default=True)
    is_metabolite = Column(Boolean, default=False)
    is_dimeric = Column(Boolean, default=False)
    parent_compounds = Column(JSON)  # For dimeric compounds
    
    # 2D RDKit Descriptors
    molecular_weight = Column(Float)
    exact_mass = Column(Float)
    logp = Column(Float)
    tpsa = Column(Float)  # Topological polar surface area
    h_bond_donors = Column(Integer)
    h_bond_acceptors = Column(Integer)
    rotatable_bonds = Column(Integer)
    aromatic_rings = Column(Integer)
    heavy_atoms = Column(Integer)
    fraction_csp3 = Column(Float)
    rdkit_descriptors = Column(JSON)  # Full 2D descriptor set (40+ features)
    
    # 3D Conformer Data (ETKDG method)
    has_conformers = Column(Boolean, default=False)
    conformer_generation_method = Column(String(50))  # "ETKDG"
    num_conformers_generated = Column(Integer)
    conformer_ensemble_id = Column(String(100))
    conformer_metadata = Column(JSON)  # Energy ranges, RMSD clusters
    rdkit_descriptors_3d = Column(JSON)  # PMI, asphericity, eccentricity, etc.
    conformers_generated_at = Column(DateTime)
    
    # Receptor Affinities with Provenance
    cb1_affinity_ki = Column(Float)  # Ki in nM
    cb1_source = Column(String(255))
    cb2_affinity_ki = Column(Float)
    cb2_source = Column(String(255))
    trpv1_affinity = Column(Float)
    gpr55_affinity = Column(Float)
    receptor_affinities = Column(JSON)  # Full provenance-rich structure
    
    # Pharmacology
    mechanism_of_action = Column(Text)
    therapeutic_categories = Column(JSON)  # ["analgesic", "anti-inflammatory", etc.]
    conditions_treated = Column(JSON)  # From clinical evidence
    
    # FDA Drug Status
    fda_approved = Column(Boolean, default=False)
    fda_drug_name = Column(String(100))  # Epidiolex, Marinol, Cesamet, Sativex
    fda_indications = Column(JSON)
    schedule_classification = Column(String(20))  # I, II, III, IV, V, unscheduled
    
    # External Database IDs (for ChEMBL/PubChem sync)
    chembl_id = Column(String(50))
    pubchem_cid = Column(String(50))
    drugbank_id = Column(String(50))
    last_external_sync = Column(DateTime)
    
    # Dimeric Prediction Scores
    dimeric_potential_score = Column(Float)
    predicted_synergy_score = Column(Float)
    triangulation_confidence = Column(Float)
    dimeric_predictions = Column(JSON)  # Full prediction data
    
    # OmniPath Integration
    omnipath_manifest_id = Column(String(100))
    tk_attribution_required = Column(Boolean, default=False)
    benefit_sharing_percentage = Column(Float)
    
    # Timestamps
    created_at = Column(DateTime(timezone=True), server_default=func.now())
    updated_at = Column(DateTime(timezone=True), onupdate=func.now())
    
    def to_cmc_profile(self) -> dict:
        """Generate CMC (Chemistry Manufacturing Controls) profile.
        
        Patent Claim 1(j): CMC templates for FDA submission
        """
        return {
            "compound_identity": {
                "name": self.name,
                "abbreviation": self.abbreviation,
                "smiles": self.smiles,
                "inchi": self.inchi,
                "inchi_key": self.inchi_key
            },
            "physicochemical_properties": {
                "molecular_weight": self.molecular_weight,
                "exact_mass": self.exact_mass,
                "logP": self.logp,
                "tpsa": self.tpsa,
                "h_bond_donors": self.h_bond_donors,
                "h_bond_acceptors": self.h_bond_acceptors,
                "rotatable_bonds": self.rotatable_bonds
            },
            "classification": {
                "compound_class": self.compound_class,
                "is_natural": self.is_natural,
                "fda_approved": self.fda_approved,
                "schedule": self.schedule_classification
            },
            "3d_structure": {
                "has_conformers": self.has_conformers,
                "generation_method": self.conformer_generation_method,
                "num_conformers": self.num_conformers_generated
            }
        }
    
    def to_pharmacology_profile(self) -> dict:
        """Generate pharmacology profile for FDA drug submission."""
        return {
            "compound": self.name,
            "receptor_binding": {
                "CB1": {"ki_nm": self.cb1_affinity_ki, "source": self.cb1_source},
                "CB2": {"ki_nm": self.cb2_affinity_ki, "source": self.cb2_source},
                "TRPV1": {"affinity": self.trpv1_affinity},
                "GPR55": {"affinity": self.gpr55_affinity}
            },
            "mechanism_of_action": self.mechanism_of_action,
            "therapeutic_applications": self.therapeutic_categories,
            "clinical_evidence": self.conditions_treated
        }


class DimericPrediction(Base):
    """Dimeric cannabinoid prediction results.
    
    Stores triangulation scoring for predicted novel dimers.
    """
    __tablename__ = "dimeric_predictions"
    
    id = Column(Integer, primary_key=True, index=True)
    
    # Dimer Identity
    dimer_name = Column(String(200), nullable=False, index=True)
    parent_1_name = Column(String(100), nullable=False)
    parent_2_name = Column(String(100), nullable=False)
    predicted_smiles = Column(String(1000))
    
    # Triangulation Scores
    triangulation_score = Column(Float)
    synergy_prediction = Column(Float)
    novelty_score = Column(Float)
    therapeutic_potential = Column(JSON)  # By condition
    
    # Predicted Properties
    predicted_cb1_affinity = Column(Float)
    predicted_cb2_affinity = Column(Float)
    predicted_selectivity = Column(Float)  # CB2/CB1 ratio
    
    # Evidence Support
    supporting_studies = Column(JSON)  # List of study_ids
    confidence_score = Column(Float)
    
    # PatentPath Lite Integration
    prior_art_score = Column(Float)
    novelty_assessment = Column(Text)
    fto_status = Column(String(50))  # Freedom to Operate
    
    # Timestamps
    created_at = Column(DateTime(timezone=True), server_default=func.now())
