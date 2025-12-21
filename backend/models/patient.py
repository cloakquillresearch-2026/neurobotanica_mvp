"""
Patient Model
Schema for patient profiles and clinical data

Supports NeuroBotanica patent claims:
- Personalized treatment recommendations
- OmniPath consent and data provenance
- HIPAA-compliant anonymized profiles
"""
from sqlalchemy import Column, Integer, String, Float, Boolean, JSON, Text, DateTime, ForeignKey, Enum
from sqlalchemy.orm import relationship
from sqlalchemy.sql import func
from .database import Base
import enum
from datetime import datetime
from typing import Optional, List, Dict


class ConsentStatus(enum.Enum):
    """OmniPath consent status for data sharing."""
    FULL_CONSENT = "full_consent"  # All data can be used
    LIMITED_CONSENT = "limited_consent"  # Only anonymized aggregated
    RESEARCH_ONLY = "research_only"  # Research use only
    REVOKED = "revoked"  # Consent withdrawn


class PatientStatus(enum.Enum):
    """Patient profile status."""
    ACTIVE = "active"
    INACTIVE = "inactive"
    ARCHIVED = "archived"


class Patient(Base):
    """Patient profile model for personalized treatment recommendations.
    
    Supports:
    - Anonymized demographic data
    - Medical history and conditions
    - Cannabinoid sensitivity profiles
    - OmniPath consent tracking
    - HIPAA-compliant data handling
    """
    __tablename__ = "patients"
    
    id = Column(Integer, primary_key=True, index=True)
    
    # Anonymized Identity
    patient_hash = Column(String(64), unique=True, nullable=False, index=True)  # SHA-256 hash
    profile_code = Column(String(20), unique=True, index=True)  # NB-XXXXX format
    
    # Demographics (anonymized ranges)
    age_range = Column(String(20))  # "18-25", "26-35", etc.
    biological_sex = Column(String(20))  # "male", "female", "other", "unspecified"
    weight_category = Column(String(20))  # "underweight", "normal", "overweight", "obese"
    bmi_range = Column(String(20))  # "<18.5", "18.5-24.9", "25-29.9", "30+"
    
    # Geographic (for regulatory compliance)
    state_code = Column(String(2), index=True)  # US state code
    country_code = Column(String(3), default="USA")
    jurisdiction_type = Column(String(20))  # "recreational", "medical", "illegal"
    
    # Medical Profile
    primary_conditions = Column(JSON)  # List of primary conditions being treated
    secondary_conditions = Column(JSON)  # Comorbidities
    condition_severity = Column(JSON)  # Condition -> severity mapping
    medical_history_summary = Column(Text)  # Anonymized summary
    
    # Current Medications (for interaction checking)
    current_medications = Column(JSON)  # List of medication classes
    medication_classes = Column(JSON)  # ["SSRI", "NSAID", "Opioid", etc.]
    contraindications = Column(JSON)  # Known drug interactions
    
    # Cannabinoid Experience Profile
    cannabinoid_experience_level = Column(String(20))  # "naive", "occasional", "regular", "daily"
    thc_tolerance = Column(String(20))  # "low", "medium", "high"
    cbd_preference = Column(Float)  # Preferred CBD:THC ratio
    preferred_delivery_methods = Column(JSON)  # ["inhalation", "sublingual", "edible", etc.]
    adverse_reactions_history = Column(JSON)  # Past negative experiences
    
    # Genetic Markers (anonymized)
    cyp450_profile = Column(JSON)  # CYP2C9, CYP2C19, CYP3A4 variants
    cb1_receptor_variant = Column(String(50))
    has_genetic_data = Column(Boolean, default=False)
    
    # Treatment Goals
    primary_treatment_goal = Column(String(100))  # "pain relief", "anxiety reduction", etc.
    secondary_goals = Column(JSON)
    outcome_priorities = Column(JSON)  # Ranked priorities
    
    # OmniPath Consent & Attribution
    omnipath_consent_status = Column(String(20), default="limited_consent")
    consent_granted_at = Column(DateTime(timezone=True))
    consent_revoked_at = Column(DateTime(timezone=True))
    data_sharing_preferences = Column(JSON)  # Granular consent settings
    attribution_required = Column(Boolean, default=False)
    benefit_sharing_eligible = Column(Boolean, default=False)
    
    # Status
    status = Column(String(20), default="active")
    profile_completeness = Column(Float, default=0.0)  # 0.0 to 1.0
    last_assessment_date = Column(DateTime(timezone=True))
    
    # Timestamps
    created_at = Column(DateTime(timezone=True), server_default=func.now())
    updated_at = Column(DateTime(timezone=True), onupdate=func.now())
    
    # Relationships
    treatments = relationship("Treatment", back_populates="patient")
    treatment_responses = relationship("TreatmentResponse", back_populates="patient")
    
    def calculate_completeness(self) -> float:
        """Calculate profile completeness score."""
        required_fields = [
            self.age_range, self.biological_sex, self.primary_conditions,
            self.cannabinoid_experience_level, self.primary_treatment_goal,
            self.state_code
        ]
        optional_fields = [
            self.weight_category, self.current_medications, self.thc_tolerance,
            self.preferred_delivery_methods, self.medical_history_summary
        ]
        
        required_complete = sum(1 for f in required_fields if f is not None) / len(required_fields)
        optional_complete = sum(1 for f in optional_fields if f is not None) / len(optional_fields)
        
        # Required fields weighted 70%, optional 30%
        return round(required_complete * 0.7 + optional_complete * 0.3, 2)
    
    def to_recommendation_profile(self) -> dict:
        """Generate profile data for treatment recommendation engine."""
        return {
            "patient_code": self.profile_code,
            "demographics": {
                "age_range": self.age_range,
                "biological_sex": self.biological_sex,
                "weight_category": self.weight_category,
                "jurisdiction": self.jurisdiction_type
            },
            "conditions": {
                "primary": self.primary_conditions or [],
                "secondary": self.secondary_conditions or [],
                "severity": self.condition_severity or {}
            },
            "cannabinoid_profile": {
                "experience_level": self.cannabinoid_experience_level,
                "thc_tolerance": self.thc_tolerance,
                "cbd_preference_ratio": self.cbd_preference,
                "delivery_preferences": self.preferred_delivery_methods or []
            },
            "contraindications": {
                "medications": self.medication_classes or [],
                "adverse_reactions": self.adverse_reactions_history or [],
                "known_interactions": self.contraindications or []
            },
            "treatment_goals": {
                "primary": self.primary_treatment_goal,
                "secondary": self.secondary_goals or [],
                "priorities": self.outcome_priorities or []
            },
            "genetic_factors": {
                "has_data": self.has_genetic_data,
                "cyp450": self.cyp450_profile if self.has_genetic_data else None
            }
        }
    
    def to_anonymized_research_profile(self) -> dict:
        """Generate fully anonymized profile for research aggregation."""
        if self.omnipath_consent_status == "revoked":
            return {"error": "Consent revoked - data cannot be used"}
        
        return {
            "demographics": {
                "age_range": self.age_range,
                "biological_sex": self.biological_sex,
                "region": self.state_code[:1] + "*" if self.state_code else None  # Partial
            },
            "conditions": self.primary_conditions or [],
            "experience_level": self.cannabinoid_experience_level,
            "consent_level": self.omnipath_consent_status
        }


class TreatmentResponse(Base):
    """Track patient responses to cannabinoid treatments.
    
    Supports:
    - Outcome tracking for efficacy analysis
    - Safety monitoring
    - Personalization model training
    """
    __tablename__ = "treatment_responses"
    
    id = Column(Integer, primary_key=True, index=True)
    
    # References
    patient_id = Column(Integer, ForeignKey("patients.id"), nullable=False, index=True)
    treatment_id = Column(Integer, ForeignKey("treatments.id"), nullable=False, index=True)
    
    # Response Assessment
    assessment_date = Column(DateTime(timezone=True), nullable=False)
    assessment_type = Column(String(50))  # "initial", "follow_up", "final"
    days_on_treatment = Column(Integer)
    
    # Efficacy Outcomes (0-10 scales)
    primary_symptom_baseline = Column(Float)  # Pre-treatment score
    primary_symptom_current = Column(Float)  # Current score
    primary_symptom_change = Column(Float)  # Calculated change
    symptom_relief_score = Column(Float)  # 0-10 subjective relief
    
    # Specific Outcome Measures
    pain_vas_score = Column(Float)  # Visual Analog Scale 0-10
    anxiety_gad7_score = Column(Float)  # GAD-7 score
    sleep_quality_score = Column(Float)  # 0-10
    functional_improvement = Column(Float)  # 0-10
    quality_of_life_change = Column(Float)  # -5 to +5
    
    # Safety & Tolerability
    adverse_events = Column(JSON)  # List of AEs experienced
    severity_of_aes = Column(String(20))  # "none", "mild", "moderate", "severe"
    treatment_discontinued = Column(Boolean, default=False)
    discontinuation_reason = Column(String(255))
    
    # Patient Satisfaction
    satisfaction_score = Column(Float)  # 0-10
    would_recommend = Column(Boolean)
    patient_comments = Column(Text)
    
    # Clinical Notes (anonymized)
    clinician_assessment = Column(Text)
    dosing_adjustments = Column(JSON)
    
    # Data Quality
    response_completeness = Column(Float)  # 0.0 to 1.0
    verified = Column(Boolean, default=False)
    
    # Timestamps
    created_at = Column(DateTime(timezone=True), server_default=func.now())
    
    # Relationships
    patient = relationship("Patient", back_populates="treatment_responses")
    treatment = relationship("Treatment", back_populates="responses")
    
    def calculate_efficacy_score(self) -> float:
        """Calculate overall efficacy score from multiple measures."""
        scores = []
        
        # Primary symptom change (weighted 40%)
        if self.primary_symptom_change is not None:
            # Normalize to 0-10 where higher = better improvement
            normalized = min(10, max(0, -self.primary_symptom_change))
            scores.append(("primary_change", normalized, 0.4))
        
        # Symptom relief (weighted 30%)
        if self.symptom_relief_score is not None:
            scores.append(("relief", self.symptom_relief_score, 0.3))
        
        # Functional improvement (weighted 20%)
        if self.functional_improvement is not None:
            scores.append(("functional", self.functional_improvement, 0.2))
        
        # Quality of life (weighted 10%)
        if self.quality_of_life_change is not None:
            normalized = min(10, max(0, (self.quality_of_life_change + 5)))  # -5 to +5 -> 0 to 10
            scores.append(("qol", normalized, 0.1))
        
        if not scores:
            return None
        
        # Normalize weights
        total_weight = sum(s[2] for s in scores)
        weighted_sum = sum(s[1] * s[2] for s in scores)
        
        return round(weighted_sum / total_weight, 2)
