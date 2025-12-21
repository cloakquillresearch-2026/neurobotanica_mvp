"""
Clinical Study Model
Schema for 320-study database with FDA regulatory support

Supports Patent Claim 1(j):
- FDA Schedule III compliance documentation
- Pharmacology data packages
- Comparative efficacy analysis
"""
from sqlalchemy import Column, Integer, String, Float, Boolean, JSON, Text, DateTime, Enum
from sqlalchemy.sql import func
from .database import Base
import enum


class StudyType(enum.Enum):
    """Clinical study design types."""
    RCT = "RCT"
    OBSERVATIONAL = "Observational"
    SYSTEMATIC_REVIEW = "Systematic Review"
    META_ANALYSIS = "Meta-Analysis"
    CASE_SERIES = "Case Series"
    GUIDELINE = "Guideline"
    MECHANISTIC = "Mechanistic"
    PRECLINICAL = "Preclinical"


class EvidenceGrade(enum.Enum):
    """Evidence quality grades for FDA documentation."""
    LEVEL_1 = "Level 1"  # High-quality RCT
    LEVEL_2 = "Level 2"  # Lower-quality RCT
    LEVEL_3 = "Level 3"  # Observational
    LEVEL_4 = "Level 4"  # Case series/reports
    LEVEL_5 = "Level 5"  # Expert opinion


class ClinicalStudy(Base):
    """Clinical study model for 320-study evidence database.
    
    Supports FDA Schedule III documentation:
    - Pharmacology data packages (mechanism of action)
    - Comparative efficacy analysis (RCT data)
    - CMC support (compound characterization)
    """
    __tablename__ = "clinical_studies"
    
    id = Column(Integer, primary_key=True, index=True)
    
    # Study Identification
    study_id = Column(String(100), unique=True, nullable=False, index=True)
    condition = Column(String(100), nullable=False, index=True)
    study_type = Column(String(50))  # Made nullable for data flexibility
    study_title = Column(Text)  # Made nullable for data flexibility
    
    # Citation Data
    citation = Column(Text)
    authors = Column(Text)
    year = Column(Integer, index=True)
    journal = Column(String(255))
    doi = Column(String(100))
    pubmed_id = Column(String(50))
    
    # Population
    sample_size = Column(String(100))
    population_description = Column(Text)
    
    # Intervention Details
    cannabinoid = Column(String(100), index=True)  # Primary cannabinoid studied
    cannabinoid_type = Column(String(50))  # THC, CBD, Mixed, Synthetic
    dosage = Column(String(255))
    duration = Column(String(100))
    delivery_method = Column(String(100))
    intervention_details = Column(JSON)  # Full intervention object
    
    # Outcomes
    primary_measure = Column(Text)
    results_summary = Column(Text)
    effect_size = Column(String(100))
    effect_size_numeric = Column(Float)  # For meta-analysis
    secondary_outcomes = Column(Text)
    outcomes_json = Column(JSON)  # Full outcomes object
    
    # Safety Data (critical for FDA)
    adverse_events = Column(Text)
    serious_adverse_events = Column(Text)
    dropout_rate = Column(String(50))
    safety_json = Column(JSON)  # Full safety object
    
    # Quality Metrics
    randomization_method = Column(String(100))
    blinding = Column(String(50))
    funding_source = Column(String(255))
    conflicts_of_interest = Column(Text)
    quality_metrics_json = Column(JSON)
    
    # FDA Regulatory Relevance
    regulatory_relevance = Column(Text)
    priority_rating = Column(String(20))  # Star rating
    fda_drug_reference = Column(String(100))  # Epidiolex, Marinol, etc.
    is_pivotal_trial = Column(Boolean, default=False)
    
    # Evidence Grading (for FDA documentation)
    evidence_grade = Column(String(20))
    confidence_weight = Column(Float, default=0.7)  # 0.5-0.9 for ML training
    
    # Mechanism of Action (for pharmacology packages)
    mechanism_of_action = Column(Text)
    receptor_targets = Column(JSON)  # CB1, CB2, TRPV1, etc.
    pharmacodynamics = Column(Text)
    
    # Timestamps
    created_at = Column(DateTime(timezone=True), server_default=func.now())
    updated_at = Column(DateTime(timezone=True), onupdate=func.now())
    
    def to_pharmacology_package(self) -> dict:
        """Generate pharmacology data package format for FDA CMC support.
        
        Patent Claim 1(j): pharmacology data packages demonstrating mechanism of action
        """
        return {
            "study_id": self.study_id,
            "cannabinoid": self.cannabinoid,
            "mechanism_of_action": self.mechanism_of_action,
            "receptor_targets": self.receptor_targets or [],
            "pharmacodynamics": self.pharmacodynamics,
            "clinical_evidence": {
                "condition": self.condition,
                "effect_size": self.effect_size,
                "sample_size": self.sample_size,
                "study_type": self.study_type,
                "evidence_grade": self.evidence_grade
            },
            "safety_profile": {
                "adverse_events": self.adverse_events,
                "serious_adverse_events": self.serious_adverse_events,
                "dropout_rate": self.dropout_rate
            }
        }
    
    def to_efficacy_comparison(self) -> dict:
        """Generate comparative efficacy format for FDA analysis.
        
        Patent Claim 1(j): comparative efficacy analysis
        """
        return {
            "study_id": self.study_id,
            "condition": self.condition,
            "intervention": {
                "cannabinoid": self.cannabinoid,
                "dosage": self.dosage,
                "duration": self.duration
            },
            "efficacy": {
                "primary_outcome": self.primary_measure,
                "effect_size": self.effect_size,
                "effect_size_numeric": self.effect_size_numeric,
                "results": self.results_summary
            },
            "comparator": "placebo" if "placebo" in (self.blinding or "").lower() else "active",
            "evidence_quality": {
                "grade": self.evidence_grade,
                "study_type": self.study_type,
                "randomization": self.randomization_method,
                "blinding": self.blinding
            }
        }
