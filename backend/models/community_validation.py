"""
Community Healer Validation Protocol Schema

Models for community-based validation of therapeutic claims,
critical for BioPath's bias-aware validation engine.

Trade Secret Component: TS-BIO-001 ($2.0B)
- Community validation provides 1.8-2.5x evidence weighting
- Minimum 5 healers for quorum
- Bias correction factor applied based on representation gaps

Reference: TRADE_SECRET_IMPLEMENTATION_REPORT.md
"""

from sqlalchemy import Column, Integer, String, Float, Boolean, JSON, Text, DateTime, ForeignKey, Enum
from sqlalchemy.orm import relationship
from datetime import datetime
import enum

from backend.models.database import Base


class HealerCredentialType(str, enum.Enum):
    """Types of healer credentials recognized."""
    TRADITIONAL_MEDICINE_PRACTITIONER = "traditional_medicine_practitioner"
    COMMUNITY_HEALER = "community_healer"
    INDIGENOUS_KNOWLEDGE_KEEPER = "indigenous_knowledge_keeper"
    HERBALIST = "herbalist"
    NATUROPATHIC_DOCTOR = "naturopathic_doctor"
    INTEGRATIVE_MEDICINE = "integrative_medicine"
    AYURVEDIC_PRACTITIONER = "ayurvedic_practitioner"
    TRADITIONAL_CHINESE_MEDICINE = "traditional_chinese_medicine"


class ValidationQuorumStatus(str, enum.Enum):
    """Status of validation quorum."""
    PENDING = "pending"
    QUORUM_MET = "quorum_met"
    QUORUM_NOT_MET = "quorum_not_met"
    VALIDATED = "validated"
    REJECTED = "rejected"


class CommunityRegion(str, enum.Enum):
    """Geographic regions for community networks."""
    NORTH_AMERICA = "north_america"
    SOUTH_AMERICA = "south_america"
    EUROPE = "europe"
    AFRICA = "africa"
    ASIA_PACIFIC = "asia_pacific"
    MIDDLE_EAST = "middle_east"
    CARIBBEAN = "caribbean"


class CommunityNetwork(Base):
    """Federated community network for validation.
    
    Each network represents a geographic/cultural region with
    its own healer pool and governance structure.
    """
    __tablename__ = "community_networks"
    
    id = Column(Integer, primary_key=True, index=True)
    
    # Network Identity
    network_code = Column(String(50), unique=True, nullable=False)
    network_name = Column(String(200), nullable=False)
    region = Column(String(50), default=CommunityRegion.NORTH_AMERICA.value)
    
    # Governance
    governance_type = Column(String(100))  # dao, council, federation
    voting_mechanism = Column(String(100))  # quadratic, simple_majority, consensus
    quorum_requirement = Column(Integer, default=5)
    
    # Status
    is_active = Column(Boolean, default=True)
    established_date = Column(DateTime)
    
    # Metrics
    total_healers = Column(Integer, default=0)
    active_healers = Column(Integer, default=0)
    validations_completed = Column(Integer, default=0)
    average_response_time_hours = Column(Float)
    
    # OmniPath Integration
    omnipath_node_id = Column(String(100))
    
    # Timestamps
    created_at = Column(DateTime, default=datetime.utcnow)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow)
    
    # Relationships
    healers = relationship("CommunityHealer", back_populates="network")
    validations = relationship("CommunityValidation", back_populates="network")


class CommunityHealer(Base):
    """Registered community healer for validation."""
    __tablename__ = "community_healers"
    
    id = Column(Integer, primary_key=True, index=True)
    network_id = Column(Integer, ForeignKey("community_networks.id"))
    
    # Identity (anonymized for privacy)
    healer_code = Column(String(50), unique=True, nullable=False)
    display_name = Column(String(200))
    
    # Credentials
    credential_type = Column(String(100))
    credential_verified = Column(Boolean, default=False)
    verification_date = Column(DateTime)
    verified_by = Column(String(100))  # Admin or network lead
    
    # Expertise
    specializations = Column(JSON)  # ["chronic_pain", "anxiety", "sleep"]
    plant_expertise = Column(JSON)  # ["cannabis", "kratom", "kava"]
    years_experience = Column(Integer)
    community_affiliation = Column(String(200))
    
    # Participation
    is_active = Column(Boolean, default=True)
    validations_completed = Column(Integer, default=0)
    average_response_time_hours = Column(Float)
    validation_accuracy = Column(Float)  # Correlation with final outcomes
    
    # Weighting (Trade Secret: contributes to bias correction)
    expertise_weight = Column(Float, default=1.0)  # 0.5-2.0 based on track record
    community_weight = Column(Float, default=1.0)  # Representation factor
    
    # Consent and Privacy
    consent_to_participate = Column(Boolean, default=True)
    consent_date = Column(DateTime)
    privacy_level = Column(String(50), default="anonymous")  # anonymous, pseudonymous, public
    
    # Timestamps
    created_at = Column(DateTime, default=datetime.utcnow)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow)
    
    # Relationships
    network = relationship("CommunityNetwork", back_populates="healers")
    responses = relationship("HealerValidationResponse", back_populates="healer")


class CommunityValidation(Base):
    """Community validation request for a therapeutic claim."""
    __tablename__ = "community_validations"
    
    id = Column(Integer, primary_key=True, index=True)
    network_id = Column(Integer, ForeignKey("community_networks.id"))
    
    # Validation Request
    validation_code = Column(String(50), unique=True, nullable=False)
    compound_name = Column(String(100), nullable=False)
    condition = Column(String(200), nullable=False)
    claim_description = Column(Text)
    
    # Evidence Context
    clinical_evidence_score = Column(Float)  # From clinical studies
    rwe_evidence_score = Column(Float)  # From real-world evidence
    evidence_summary = Column(Text)
    
    # Quorum Status
    quorum_status = Column(String(50), default=ValidationQuorumStatus.PENDING.value)
    required_responses = Column(Integer, default=5)
    received_responses = Column(Integer, default=0)
    
    # Results (populated after quorum met)
    community_score = Column(Float)  # Raw community score
    bias_corrected_score = Column(Float)  # After bias correction
    bias_correction_factor = Column(Float)  # Applied factor
    confidence_interval = Column(JSON)  # [lower, upper]
    
    # Final Status
    validation_result = Column(String(50))  # validated, rejected, conditional
    validation_notes = Column(Text)
    
    # Timing
    requested_at = Column(DateTime, default=datetime.utcnow)
    quorum_reached_at = Column(DateTime)
    completed_at = Column(DateTime)
    
    # Timestamps
    created_at = Column(DateTime, default=datetime.utcnow)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow)
    
    # Relationships
    network = relationship("CommunityNetwork", back_populates="validations")
    responses = relationship("HealerValidationResponse", back_populates="validation")


class HealerValidationResponse(Base):
    """Individual healer's response to a validation request."""
    __tablename__ = "healer_validation_responses"
    
    id = Column(Integer, primary_key=True, index=True)
    validation_id = Column(Integer, ForeignKey("community_validations.id"))
    healer_id = Column(Integer, ForeignKey("community_healers.id"))
    
    # Response
    validates_claim = Column(Boolean)  # True = supports, False = rejects
    confidence_level = Column(Float)  # 0-1 how confident they are
    
    # Detailed Assessment
    efficacy_rating = Column(Float)  # 0-10
    safety_rating = Column(Float)  # 0-10
    traditional_use_confirmed = Column(Boolean)
    personal_experience = Column(Boolean)  # Has direct experience
    
    # Supporting Information
    reasoning = Column(Text)
    traditional_context = Column(Text)  # Traditional use context
    preparation_notes = Column(Text)  # Recommended preparation
    dosing_guidance = Column(Text)
    contraindications = Column(JSON)  # Known contraindications
    
    # Quality Metrics
    response_time_hours = Column(Float)
    completeness_score = Column(Float)  # How complete the response is
    
    # Weighting (calculated)
    weighted_contribution = Column(Float)  # Final weighted score contribution
    
    # Timestamps
    submitted_at = Column(DateTime, default=datetime.utcnow)
    
    # Relationships
    validation = relationship("CommunityValidation", back_populates="responses")
    healer = relationship("CommunityHealer", back_populates="responses")


# Community validation configuration
COMMUNITY_VALIDATION_CONFIG = {
    # Quorum requirements
    "min_healers_for_quorum": 5,
    "min_response_rate": 0.6,  # 60% of invited healers must respond
    
    # Weighting (Trade Secret values)
    "base_community_weight": 1.8,  # Community evidence weight multiplier
    "max_community_weight": 2.5,  # Maximum after bonuses
    "experience_bonus_per_year": 0.02,  # Per year of experience
    "accuracy_bonus_max": 0.3,  # Based on validation accuracy
    
    # Bias correction thresholds (Trade Secret)
    "representation_bias_threshold": 0.3,
    "correction_factor_min": 1.0,
    "correction_factor_max": 2.5,
    
    # Timing
    "validation_timeout_days": 7,
    "emergency_timeout_hours": 24,
    
    # Confidence adjustment
    "confidence_boost_for_consensus": 0.1,  # If all healers agree
    "confidence_penalty_for_disagreement": 0.05,  # Per dissenting healer
}

# Supported conditions for community validation
VALIDATED_CONDITIONS = [
    "chronic_pain",
    "anxiety",
    "depression",
    "ptsd",
    "insomnia",
    "nausea",
    "inflammation",
    "epilepsy",
    "neuroprotection",
    "appetite_stimulation",
    "muscle_spasticity",
    "glaucoma",
]

# Geographic network seed data
COMMUNITY_NETWORKS_SEED = [
    {
        "network_code": "NA-001",
        "network_name": "North American Traditional Medicine Network",
        "region": CommunityRegion.NORTH_AMERICA.value,
        "governance_type": "federation",
        "voting_mechanism": "quadratic"
    },
    {
        "network_code": "SA-001",
        "network_name": "Amazonian Traditional Knowledge Network",
        "region": CommunityRegion.SOUTH_AMERICA.value,
        "governance_type": "council",
        "voting_mechanism": "consensus"
    },
    {
        "network_code": "AP-001",
        "network_name": "Asia-Pacific Traditional Medicine Network",
        "region": CommunityRegion.ASIA_PACIFIC.value,
        "governance_type": "federation",
        "voting_mechanism": "simple_majority"
    },
    {
        "network_code": "AF-001",
        "network_name": "African Traditional Healing Network",
        "region": CommunityRegion.AFRICA.value,
        "governance_type": "council",
        "voting_mechanism": "consensus"
    }
]
