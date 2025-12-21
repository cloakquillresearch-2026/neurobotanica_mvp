"""
Treatment Model
Schema for cannabinoid treatment recommendations and protocols

Supports NeuroBotanica patent claims:
- Personalized dosing recommendations
- Dimeric cannabinoid optimization
- FDA Schedule III compliance tracking
"""
from sqlalchemy import Column, Integer, String, Float, Boolean, JSON, Text, DateTime, ForeignKey, Enum
from sqlalchemy.orm import relationship
from sqlalchemy.sql import func
from .database import Base
import enum
from datetime import datetime
from typing import Optional, List, Dict


class TreatmentStatus(enum.Enum):
    """Treatment protocol status."""
    DRAFT = "draft"  # Being developed
    PROPOSED = "proposed"  # Ready for patient review
    ACTIVE = "active"  # Currently being followed
    COMPLETED = "completed"  # Treatment course finished
    DISCONTINUED = "discontinued"  # Stopped early
    MODIFIED = "modified"  # Changed from original


class DeliveryMethod(enum.Enum):
    """Cannabinoid delivery methods."""
    INHALATION = "inhalation"  # Smoking/vaporization
    SUBLINGUAL = "sublingual"  # Under tongue
    ORAL = "oral"  # Edibles/capsules
    TOPICAL = "topical"  # Creams/patches
    TRANSDERMAL = "transdermal"  # Patches
    SUPPOSITORY = "suppository"


class Treatment(Base):
    """Treatment recommendation model for personalized cannabinoid therapy.
    
    Supports:
    - Evidence-based dosing protocols
    - Dimeric compound recommendations
    - Entourage effect optimization
    - Schedule III compliance documentation
    """
    __tablename__ = "treatments"
    
    id = Column(Integer, primary_key=True, index=True)
    
    # Treatment Identity
    treatment_code = Column(String(30), unique=True, nullable=False, index=True)  # NB-TRT-XXXXX
    version = Column(Integer, default=1)
    parent_treatment_id = Column(Integer, ForeignKey("treatments.id"))  # For modifications
    
    # Patient Reference
    patient_id = Column(Integer, ForeignKey("patients.id"), nullable=False, index=True)
    
    # Target Condition
    primary_condition = Column(String(100), nullable=False, index=True)
    condition_severity = Column(String(20))  # "mild", "moderate", "severe"
    secondary_conditions = Column(JSON)  # Additional conditions being addressed
    treatment_goals = Column(JSON)  # Specific outcome targets
    
    # Primary Cannabinoid Recommendation
    primary_cannabinoid_id = Column(Integer, ForeignKey("cannabinoids.id"))
    primary_cannabinoid_name = Column(String(100), nullable=False)
    primary_cannabinoid_dosage = Column(String(100))  # "10mg", "2.5-5mg", etc.
    
    # Secondary/Adjunct Cannabinoids
    secondary_cannabinoids = Column(JSON)  # List of {name, dosage, rationale}
    cannabinoid_ratio = Column(String(50))  # "THC:CBD 1:20", "CBD only", etc.
    
    # Dimeric Recommendations (NeuroBotanica specialty)
    dimeric_compound_id = Column(Integer, ForeignKey("cannabinoids.id"))
    dimeric_compound_name = Column(String(100))
    dimeric_prediction_score = Column(Float)
    dimeric_rationale = Column(Text)  # Why this dimeric compound was selected
    
    # Dosing Protocol
    delivery_method = Column(String(50))
    initial_dose = Column(String(100))
    maintenance_dose = Column(String(100))
    max_daily_dose = Column(String(100))
    dosing_frequency = Column(String(100))  # "BID", "TID", "PRN", etc.
    titration_schedule = Column(JSON)  # Week-by-week titration plan
    duration_weeks = Column(Integer)
    
    # Timing & Administration
    administration_timing = Column(JSON)  # ["morning", "evening", "before_bed"]
    food_requirements = Column(String(100))  # "with food", "empty stomach", etc.
    special_instructions = Column(Text)
    
    # Clinical Evidence Basis
    evidence_citations = Column(JSON)  # Study IDs supporting this treatment
    evidence_strength = Column(String(20))  # "strong", "moderate", "limited"
    evidence_summary = Column(Text)
    
    # Safety Considerations
    contraindications_checked = Column(JSON)  # Drug interactions verified
    warnings = Column(JSON)  # Patient-specific warnings
    monitoring_requirements = Column(JSON)  # Required follow-up tests/visits
    adverse_event_risk_profile = Column(JSON)  # Expected AE probabilities
    
    # Expected Outcomes
    expected_efficacy_score = Column(Float)  # Predicted 0-10 based on evidence
    time_to_effect_days = Column(Integer)  # Expected days until benefit
    response_probability = Column(Float)  # Probability of meaningful response
    
    # Regulatory Compliance
    schedule_classification = Column(String(20))  # "III", "IV", "unscheduled"
    requires_prescription = Column(Boolean, default=True)
    state_compliance_notes = Column(Text)  # State-specific requirements
    fda_drug_reference = Column(String(100))  # Comparable FDA-approved drug
    
    # Alternative Recommendations
    alternative_treatments = Column(JSON)  # Ranked alternatives if primary fails
    escalation_plan = Column(JSON)  # What to do if insufficient response
    de_escalation_plan = Column(JSON)  # What to do if over-response/AEs
    
    # Status & Workflow
    status = Column(String(20), default="draft")
    approved_by = Column(String(100))  # Clinician approval
    approved_at = Column(DateTime(timezone=True))
    start_date = Column(DateTime(timezone=True))
    end_date = Column(DateTime(timezone=True))
    
    # OmniPath Attribution
    traditional_knowledge_attribution = Column(JSON)  # TK sources if applicable
    benefit_sharing_required = Column(Boolean, default=False)
    
    # Timestamps
    created_at = Column(DateTime(timezone=True), server_default=func.now())
    updated_at = Column(DateTime(timezone=True), onupdate=func.now())
    
    # Relationships
    patient = relationship("Patient", back_populates="treatments")
    responses = relationship("TreatmentResponse", back_populates="treatment")
    
    def to_patient_summary(self) -> dict:
        """Generate patient-friendly treatment summary."""
        return {
            "treatment_id": self.treatment_code,
            "condition": self.primary_condition,
            "recommendation": {
                "primary_compound": self.primary_cannabinoid_name,
                "dose": self.primary_cannabinoid_dosage,
                "ratio": self.cannabinoid_ratio,
                "delivery": self.delivery_method,
                "frequency": self.dosing_frequency
            },
            "schedule": {
                "duration_weeks": self.duration_weeks,
                "timing": self.administration_timing,
                "special_instructions": self.special_instructions
            },
            "expected_outcomes": {
                "time_to_effect_days": self.time_to_effect_days,
                "response_probability": f"{(self.response_probability or 0) * 100:.0f}%"
            },
            "important_warnings": self.warnings,
            "follow_up_required": self.monitoring_requirements
        }
    
    def to_clinician_protocol(self) -> dict:
        """Generate clinical protocol for healthcare providers."""
        return {
            "treatment_code": self.treatment_code,
            "version": self.version,
            "patient_id": self.patient_id,
            "diagnosis": {
                "primary_condition": self.primary_condition,
                "severity": self.condition_severity,
                "comorbidities": self.secondary_conditions
            },
            "therapeutic_regimen": {
                "primary_agent": {
                    "compound": self.primary_cannabinoid_name,
                    "initial_dose": self.initial_dose,
                    "maintenance_dose": self.maintenance_dose,
                    "max_daily": self.max_daily_dose,
                    "route": self.delivery_method,
                    "frequency": self.dosing_frequency
                },
                "adjunct_agents": self.secondary_cannabinoids,
                "cannabinoid_ratio": self.cannabinoid_ratio,
                "titration": self.titration_schedule
            },
            "dimeric_optimization": {
                "recommended_compound": self.dimeric_compound_name,
                "prediction_score": self.dimeric_prediction_score,
                "rationale": self.dimeric_rationale
            } if self.dimeric_compound_name else None,
            "clinical_evidence": {
                "citations": self.evidence_citations,
                "strength": self.evidence_strength,
                "summary": self.evidence_summary
            },
            "safety": {
                "contraindications_verified": self.contraindications_checked,
                "warnings": self.warnings,
                "expected_adverse_events": self.adverse_event_risk_profile,
                "monitoring": self.monitoring_requirements
            },
            "regulatory": {
                "schedule": self.schedule_classification,
                "prescription_required": self.requires_prescription,
                "fda_reference": self.fda_drug_reference,
                "state_notes": self.state_compliance_notes
            },
            "outcome_prediction": {
                "expected_efficacy": self.expected_efficacy_score,
                "time_to_effect_days": self.time_to_effect_days,
                "response_probability": self.response_probability
            },
            "contingency_plans": {
                "alternatives": self.alternative_treatments,
                "escalation": self.escalation_plan,
                "de_escalation": self.de_escalation_plan
            }
        }
    
    def to_fda_documentation(self) -> dict:
        """Generate FDA Schedule III compliance documentation.
        
        Patent Claim 1(j): FDA Schedule III compliance documentation
        """
        return {
            "treatment_protocol_id": self.treatment_code,
            "regulatory_classification": {
                "schedule": self.schedule_classification,
                "fda_reference_drug": self.fda_drug_reference,
                "prescription_required": self.requires_prescription
            },
            "indication": {
                "primary_condition": self.primary_condition,
                "severity": self.condition_severity,
                "icd_codes": None  # To be populated from condition mapping
            },
            "dosing_regimen": {
                "initial_dose": self.initial_dose,
                "maintenance_dose": self.maintenance_dose,
                "maximum_dose": self.max_daily_dose,
                "route_of_administration": self.delivery_method,
                "frequency": self.dosing_frequency,
                "duration": f"{self.duration_weeks} weeks"
            },
            "clinical_evidence_package": {
                "supporting_studies": self.evidence_citations,
                "evidence_grade": self.evidence_strength,
                "efficacy_summary": self.evidence_summary
            },
            "safety_profile": {
                "contraindications": self.contraindications_checked,
                "warnings_and_precautions": self.warnings,
                "adverse_events": self.adverse_event_risk_profile,
                "monitoring_plan": self.monitoring_requirements
            },
            "comparator_analysis": {
                "reference_drug": self.fda_drug_reference,
                "expected_efficacy": self.expected_efficacy_score,
                "safety_comparison": None  # To be populated from comparative analysis
            },
            "state_compliance": {
                "notes": self.state_compliance_notes
            },
            "document_version": self.version,
            "approval_status": self.status,
            "approved_by": self.approved_by,
            "approval_date": self.approved_at.isoformat() if self.approved_at else None
        }


class TreatmentTemplate(Base):
    """Pre-defined treatment templates for common conditions.
    
    Supports rapid treatment generation based on evidence-based protocols.
    """
    __tablename__ = "treatment_templates"
    
    id = Column(Integer, primary_key=True, index=True)
    
    # Template Identity
    template_code = Column(String(30), unique=True, nullable=False, index=True)
    name = Column(String(200), nullable=False)
    description = Column(Text)
    
    # Target Condition
    condition = Column(String(100), nullable=False, index=True)
    condition_subtypes = Column(JSON)  # Applicable subtypes
    severity_range = Column(JSON)  # ["mild", "moderate"] etc.
    
    # Recommended Protocol
    primary_cannabinoid = Column(String(100), nullable=False)
    cannabinoid_ratio = Column(String(50))
    delivery_method = Column(String(50))
    
    # Dosing
    starting_dose = Column(String(100))
    target_dose = Column(String(100))
    max_dose = Column(String(100))
    titration_protocol = Column(JSON)
    
    # Evidence
    evidence_basis = Column(JSON)  # Study IDs
    evidence_grade = Column(String(20))
    
    # Applicability
    patient_criteria = Column(JSON)  # Requirements for template use
    exclusion_criteria = Column(JSON)  # When NOT to use
    
    # Status
    is_active = Column(Boolean, default=True)
    version = Column(Integer, default=1)
    
    # Timestamps
    created_at = Column(DateTime(timezone=True), server_default=func.now())
    updated_at = Column(DateTime(timezone=True), onupdate=func.now())
    
    def generate_treatment(self, patient_id: int, customizations: dict = None) -> dict:
        """Generate a Treatment from this template."""
        treatment_data = {
            "patient_id": patient_id,
            "primary_condition": self.condition,
            "primary_cannabinoid_name": self.primary_cannabinoid,
            "cannabinoid_ratio": self.cannabinoid_ratio,
            "delivery_method": self.delivery_method,
            "initial_dose": self.starting_dose,
            "maintenance_dose": self.target_dose,
            "max_daily_dose": self.max_dose,
            "titration_schedule": self.titration_protocol,
            "evidence_citations": self.evidence_basis,
            "evidence_strength": self.evidence_grade,
            "status": "draft"
        }
        
        # Apply customizations if provided
        if customizations:
            treatment_data.update(customizations)
        
        return treatment_data
