"""
Receptor Affinity Provenance Schema
Enhanced data model for receptor binding data with full provenance tracking

Supports NeuroBotanica patent claims:
- Assay-specific provenance for receptor affinities
- Multi-source evidence triangulation
- Data quality and confidence assessment

Reference: NeuroBotanica MVP Development Plan - Week 3 Task 3.2
"""
from typing import Optional, List, Dict, Any
from datetime import datetime
from enum import Enum
from dataclasses import dataclass, field, asdict
from pydantic import BaseModel, Field
from sqlalchemy import Column, Integer, String, Float, Boolean, JSON, Text, DateTime, ForeignKey, Enum as SQLEnum
from sqlalchemy.orm import relationship
from sqlalchemy.sql import func

from .database import Base


# ==================== Enumerations ====================

class AffinityUnit(str, Enum):
    """Standard units for receptor binding affinity."""
    NANOMOLAR = "nM"
    MICROMOLAR = "µM"
    PICOMOLAR = "pM"
    MILLIMOLAR = "mM"
    PERCENT_INHIBITION = "% inhibition"
    EC50 = "EC50"
    IC50 = "IC50"


class AssayType(str, Enum):
    """Standard assay types for receptor binding studies."""
    RADIOLIGAND_BINDING = "radioligand_binding"
    FUNCTIONAL_CAMP = "functional_cAMP"
    FUNCTIONAL_CALCIUM = "functional_calcium"
    BETA_ARRESTIN = "beta_arrestin"
    GTPGAMMA = "GTPγS"
    BRET = "BRET"
    FRET = "FRET"
    SPR = "SPR"  # Surface Plasmon Resonance
    ITC = "ITC"  # Isothermal Titration Calorimetry
    OTHER = "other"


class ReceptorType(str, Enum):
    """Cannabinoid-related receptors."""
    CB1 = "CB1"
    CB2 = "CB2"
    GPR55 = "GPR55"
    GPR18 = "GPR18"
    GPR119 = "GPR119"
    TRPV1 = "TRPV1"
    TRPV2 = "TRPV2"
    TRPA1 = "TRPA1"
    TRPM8 = "TRPM8"
    PPAR_ALPHA = "PPAR-α"
    PPAR_GAMMA = "PPAR-γ"
    ADENOSINE_A1 = "A1"
    ADENOSINE_A2A = "A2A"
    SEROTONIN_5HT1A = "5-HT1A"
    DOPAMINE_D2 = "D2"


class SourceQuality(str, Enum):
    """Data source quality classification."""
    PEER_REVIEWED = "peer_reviewed"
    PREPRINT = "preprint"
    DATABASE_CURATED = "database_curated"
    VENDOR_DATA = "vendor_data"
    PREDICTED = "predicted"
    INTERNAL = "internal"


# ==================== SQLAlchemy Models ====================

class ReceptorAffinity(Base):
    """Receptor affinity measurement with full provenance.
    
    Stores individual binding measurements with:
    - Assay details and conditions
    - Source attribution
    - Confidence scoring
    - Data quality metrics
    """
    __tablename__ = "receptor_affinities"
    
    id = Column(Integer, primary_key=True, index=True)
    
    # Compound Reference
    cannabinoid_id = Column(Integer, ForeignKey("cannabinoids.id"), index=True)
    compound_name = Column(String(100), index=True)  # Denormalized for query efficiency
    
    # Receptor Target
    receptor = Column(String(50), nullable=False, index=True)
    receptor_subtype = Column(String(50))  # e.g., CB1_human, CB1_rat
    target_organism = Column(String(50), default="Homo sapiens")
    target_tissue = Column(String(100))  # e.g., brain membrane
    target_cell_type = Column(String(100))  # e.g., CHO cells, HEK293
    
    # Affinity Measurement
    affinity_value = Column(Float, nullable=False)
    affinity_unit = Column(String(20), nullable=False, default="nM")
    affinity_type = Column(String(20), default="Ki")  # Ki, Kd, EC50, IC50
    affinity_modifier = Column(String(10))  # '<', '>', '≈', '±'
    affinity_error = Column(Float)  # Standard error or deviation
    error_type = Column(String(20))  # "SEM", "SD", "95% CI"
    
    # Functional Data
    efficacy = Column(Float)  # % efficacy (for functional assays)
    efficacy_type = Column(String(50))  # "Emax", "intrinsic activity"
    functional_response = Column(String(50))  # "agonist", "antagonist", "inverse_agonist", "partial_agonist"
    
    # Assay Information
    assay_type = Column(String(50), nullable=False)
    assay_description = Column(Text)
    radioligand = Column(String(100))  # For radioligand binding assays
    radioligand_concentration = Column(String(50))
    incubation_time = Column(String(50))
    temperature = Column(String(20))
    buffer_conditions = Column(Text)
    
    # Source Attribution
    source_type = Column(String(50), nullable=False)  # "pubmed", "chembl", "vendor", etc.
    pubmed_id = Column(String(20), index=True)
    doi = Column(String(200))
    chembl_assay_id = Column(String(50), index=True)
    chembl_document_id = Column(String(50))
    source_database = Column(String(50))  # "ChEMBL", "BindingDB", "PDSP", "IUPHAR"
    source_url = Column(String(500))
    
    # Publication Details
    authors = Column(Text)
    journal = Column(String(200))
    publication_year = Column(Integer)
    study_title = Column(Text)
    
    # Confidence and Quality
    confidence_score = Column(Float, default=0.5)  # 0.0-1.0
    data_quality = Column(String(50), default="peer_reviewed")
    curation_level = Column(String(50))  # "expert", "automated", "uncurated"
    assay_heterogeneity_flag = Column(Boolean, default=False)
    heterogeneity_notes = Column(Text)
    
    # OmniPath Provenance
    provenance_manifest_id = Column(String(100))
    tk_attribution_id = Column(String(100))
    extraction_timestamp = Column(DateTime)
    
    # Timestamps
    created_at = Column(DateTime(timezone=True), server_default=func.now())
    updated_at = Column(DateTime(timezone=True), onupdate=func.now())
    
    def to_provenance_dict(self) -> Dict[str, Any]:
        """Generate provenance-rich dictionary representation."""
        return {
            "compound": self.compound_name,
            "receptor": self.receptor,
            "affinity": {
                "value": self.affinity_value,
                "unit": self.affinity_unit,
                "type": self.affinity_type,
                "modifier": self.affinity_modifier,
                "error": self.affinity_error,
                "error_type": self.error_type
            },
            "assay": {
                "type": self.assay_type,
                "description": self.assay_description,
                "organism": self.target_organism,
                "cell_type": self.target_cell_type,
                "radioligand": self.radioligand
            },
            "source": {
                "type": self.source_type,
                "pubmed_id": self.pubmed_id,
                "doi": self.doi,
                "database": self.source_database,
                "chembl_assay_id": self.chembl_assay_id
            },
            "quality": {
                "confidence_score": self.confidence_score,
                "data_quality": self.data_quality,
                "curation_level": self.curation_level,
                "heterogeneity_flag": self.assay_heterogeneity_flag
            },
            "provenance": {
                "manifest_id": self.provenance_manifest_id,
                "extraction_timestamp": self.extraction_timestamp.isoformat() if self.extraction_timestamp else None
            }
        }


class ReceptorAffinityAggregate(Base):
    """Aggregated receptor affinity across multiple measurements.
    
    Provides consensus values with:
    - Weighted mean calculation
    - Variance tracking
    - Source triangulation
    """
    __tablename__ = "receptor_affinity_aggregates"
    
    id = Column(Integer, primary_key=True, index=True)
    
    # Compound Reference
    cannabinoid_id = Column(Integer, ForeignKey("cannabinoids.id"), index=True)
    compound_name = Column(String(100), index=True)
    
    # Receptor Target
    receptor = Column(String(50), nullable=False, index=True)
    
    # Consensus Values
    consensus_ki = Column(Float)  # Weighted mean Ki
    consensus_unit = Column(String(20), default="nM")
    measurement_count = Column(Integer, default=0)
    
    # Statistical Measures
    geometric_mean = Column(Float)
    arithmetic_mean = Column(Float)
    median_value = Column(Float)
    min_value = Column(Float)
    max_value = Column(Float)
    standard_deviation = Column(Float)
    coefficient_of_variation = Column(Float)
    interquartile_range = Column(Float)
    
    # Source Analysis
    source_count = Column(Integer)
    pubmed_sources = Column(Integer)
    chembl_sources = Column(Integer)
    unique_assay_types = Column(Integer)
    source_agreement_score = Column(Float)  # 0-1 how well sources agree
    
    # Quality Metrics
    overall_confidence = Column(Float)  # Weighted confidence
    heterogeneity_detected = Column(Boolean, default=False)
    heterogeneity_type = Column(String(100))  # "unit_mismatch", "assay_variance", etc.
    heterogeneity_score = Column(Float)  # Measure of data heterogeneity
    
    # Reference to individual measurements
    measurement_ids = Column(JSON)  # List of receptor_affinity IDs
    
    # Timestamps
    calculated_at = Column(DateTime(timezone=True), server_default=func.now())
    updated_at = Column(DateTime(timezone=True), onupdate=func.now())


# ==================== Pydantic Models for API ====================

class AffinityMeasurement(BaseModel):
    """Single affinity measurement with provenance."""
    value: float
    unit: str = "nM"
    type: str = "Ki"  # Ki, Kd, EC50, IC50
    modifier: Optional[str] = None
    error: Optional[float] = None
    error_type: Optional[str] = None


class AssayDetails(BaseModel):
    """Assay methodology details."""
    type: str
    description: Optional[str] = None
    organism: str = "Homo sapiens"
    cell_type: Optional[str] = None
    tissue: Optional[str] = None
    radioligand: Optional[str] = None
    temperature: Optional[str] = None


class SourceAttribution(BaseModel):
    """Source attribution for data."""
    type: str  # "pubmed", "chembl", "vendor"
    pubmed_id: Optional[str] = None
    doi: Optional[str] = None
    database: Optional[str] = None
    chembl_assay_id: Optional[str] = None
    url: Optional[str] = None


class QualityMetrics(BaseModel):
    """Data quality assessment."""
    confidence_score: float = Field(ge=0.0, le=1.0)
    data_quality: str = "peer_reviewed"
    curation_level: Optional[str] = None
    heterogeneity_flag: bool = False
    heterogeneity_notes: Optional[str] = None


class ReceptorAffinityCreate(BaseModel):
    """Create new receptor affinity record."""
    cannabinoid_id: Optional[int] = None
    compound_name: str
    receptor: str
    affinity: AffinityMeasurement
    assay: AssayDetails
    source: SourceAttribution
    quality: Optional[QualityMetrics] = None


class ReceptorAffinityResponse(BaseModel):
    """Response model for receptor affinity queries."""
    id: int
    compound_name: str
    receptor: str
    affinity: AffinityMeasurement
    assay: AssayDetails
    source: SourceAttribution
    quality: QualityMetrics
    created_at: datetime


class AggregatedAffinityResponse(BaseModel):
    """Aggregated affinity with statistics."""
    compound_name: str
    receptor: str
    consensus_ki: float
    consensus_unit: str
    measurement_count: int
    statistics: Dict[str, float]
    source_analysis: Dict[str, Any]
    quality_metrics: Dict[str, Any]
    individual_measurements: List[ReceptorAffinityResponse]


# ==================== Helper Functions ====================

def calculate_confidence_score(
    source_type: str,
    has_pubmed: bool,
    has_chembl: bool,
    data_quality: str,
    measurement_count: int = 1
) -> float:
    """Calculate confidence score for receptor affinity data.
    
    Factors:
    - Source type (peer-reviewed highest)
    - Database presence (ChEMBL, PubMed)
    - Data quality classification
    - Number of corroborating measurements
    """
    score = 0.0
    
    # Source type contribution (0-0.4)
    source_weights = {
        "peer_reviewed": 0.4,
        "database_curated": 0.35,
        "preprint": 0.25,
        "vendor_data": 0.2,
        "predicted": 0.1,
        "internal": 0.15
    }
    score += source_weights.get(source_type, 0.1)
    
    # Database presence (0-0.3)
    if has_pubmed:
        score += 0.15
    if has_chembl:
        score += 0.15
    
    # Quality classification (0-0.2)
    quality_weights = {
        "peer_reviewed": 0.2,
        "expert_curated": 0.18,
        "automated_curated": 0.1,
        "uncurated": 0.05
    }
    score += quality_weights.get(data_quality, 0.05)
    
    # Measurement count bonus (0-0.1)
    count_bonus = min(measurement_count / 10, 1.0) * 0.1
    score += count_bonus
    
    return min(score, 1.0)


def normalize_affinity_unit(value: float, from_unit: str, to_unit: str = "nM") -> float:
    """Convert affinity values between units.
    
    Supports: pM, nM, µM, mM
    """
    # Convert to nM first
    to_nm = {
        "pM": 0.001,
        "nM": 1.0,
        "µM": 1000.0,
        "uM": 1000.0,
        "mM": 1000000.0
    }
    
    nm_value = value * to_nm.get(from_unit, 1.0)
    
    # Convert from nM to target
    from_nm = {
        "pM": 1000.0,
        "nM": 1.0,
        "µM": 0.001,
        "uM": 0.001,
        "mM": 0.000001
    }
    
    return nm_value * from_nm.get(to_unit, 1.0)


def detect_assay_heterogeneity(measurements: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Detect heterogeneity in a set of affinity measurements.
    
    Checks for:
    - Unit mismatches
    - Assay type variance
    - Value variance (coefficient of variation)
    - Source disagreement
    """
    if not measurements:
        return {"detected": False, "issues": []}
    
    issues = []
    
    # Check unit consistency
    units = set(m.get("affinity", {}).get("unit", "nM") for m in measurements)
    if len(units) > 1:
        issues.append({
            "type": "unit_mismatch",
            "description": f"Multiple units detected: {units}",
            "severity": "high"
        })
    
    # Check assay type consistency
    assay_types = set(m.get("assay", {}).get("type", "unknown") for m in measurements)
    if len(assay_types) > 1:
        issues.append({
            "type": "assay_variance",
            "description": f"Multiple assay types: {assay_types}",
            "severity": "medium"
        })
    
    # Check value variance (coefficient of variation)
    values = [m.get("affinity", {}).get("value", 0) for m in measurements if m.get("affinity", {}).get("value")]
    if len(values) > 1:
        mean_val = sum(values) / len(values)
        if mean_val > 0:
            variance = sum((v - mean_val) ** 2 for v in values) / len(values)
            cv = (variance ** 0.5) / mean_val
            
            if cv > 1.0:  # CV > 100%
                issues.append({
                    "type": "high_variance",
                    "description": f"High coefficient of variation: {cv:.2f}",
                    "severity": "high"
                })
            elif cv > 0.5:  # CV > 50%
                issues.append({
                    "type": "moderate_variance",
                    "description": f"Moderate coefficient of variation: {cv:.2f}",
                    "severity": "medium"
                })
    
    return {
        "detected": len(issues) > 0,
        "issues": issues,
        "heterogeneity_score": len(issues) / 3.0  # 3 possible issue types
    }
