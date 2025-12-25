"""
Nevada Regulatory Data Models

Schema for Nevada-specific cannabis regulatory compliance data
required for the Nevada pilot launch (March 2026).

Reference: Nevada Administrative Code (NAC) Chapter 453D
"""

from sqlalchemy import Column, Integer, String, Float, Boolean, JSON, Text, DateTime, ForeignKey, Enum
from sqlalchemy.orm import relationship
from datetime import datetime
import enum

from backend.models.database import Base


class NevadaLicenseType(str, enum.Enum):
    """Nevada cannabis license types."""
    CULTIVATION = "cultivation"
    PRODUCTION = "production"
    DISPENSARY = "dispensary"
    DISTRIBUTION = "distribution"
    TESTING_LAB = "testing_lab"


class ComplianceStatus(str, enum.Enum):
    """Compliance status for products/dispensaries."""
    COMPLIANT = "compliant"
    PENDING_REVIEW = "pending_review"
    NON_COMPLIANT = "non_compliant"
    CONDITIONAL = "conditional"


class NevadaDispensary(Base):
    """Nevada dispensary registration for pilot program."""
    __tablename__ = "nevada_dispensaries"
    
    id = Column(Integer, primary_key=True, index=True)
    
    # License Information
    license_number = Column(String(50), unique=True, nullable=False)
    license_type = Column(String(50), default=NevadaLicenseType.DISPENSARY.value)
    license_expiry = Column(DateTime)
    
    # Business Information
    business_name = Column(String(200), nullable=False)
    dba_name = Column(String(200))
    address = Column(String(500))
    city = Column(String(100))
    county = Column(String(100))
    zip_code = Column(String(10))
    
    # Contact
    contact_name = Column(String(200))
    contact_email = Column(String(200))
    contact_phone = Column(String(20))
    
    # Pilot Program Status
    pilot_enrolled = Column(Boolean, default=False)
    pilot_enrollment_date = Column(DateTime)
    pilot_status = Column(String(50))  # active, paused, completed
    
    # Compliance
    compliance_status = Column(String(50), default=ComplianceStatus.PENDING_REVIEW.value)
    last_inspection_date = Column(DateTime)
    inspection_notes = Column(Text)
    
    # Integration
    omnipath_registered = Column(Boolean, default=False)
    api_key_hash = Column(String(256))
    
    # Timestamps
    created_at = Column(DateTime, default=datetime.utcnow)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow)
    
    # Relationships
    products = relationship("NevadaProduct", back_populates="dispensary")
    rwe_entries = relationship("RealWorldEvidence", back_populates="dispensary")


class NevadaProduct(Base):
    """Cannabis products registered in Nevada."""
    __tablename__ = "nevada_products"
    
    id = Column(Integer, primary_key=True, index=True)
    dispensary_id = Column(Integer, ForeignKey("nevada_dispensaries.id"))
    
    # Product Information
    product_name = Column(String(200), nullable=False)
    product_type = Column(String(100))  # flower, concentrate, edible, topical, tincture
    brand_name = Column(String(200))
    batch_number = Column(String(100))
    
    # Cannabinoid Profile
    thc_percent = Column(Float)
    cbd_percent = Column(Float)
    cbn_percent = Column(Float)
    cbg_percent = Column(Float)
    total_cannabinoids = Column(Float)
    
    # Terpene Profile
    terpene_profile = Column(JSON)  # {terpene_name: percent}
    dominant_terpenes = Column(JSON)  # ["myrcene", "limonene"]
    
    # Testing
    testing_lab = Column(String(200))
    coa_number = Column(String(100))
    test_date = Column(DateTime)
    
    # Compliance
    compliance_status = Column(String(50), default=ComplianceStatus.PENDING_REVIEW.value)
    meets_nevada_limits = Column(Boolean, default=True)
    
    # NeuroBotanica Integration
    neurobotanica_compound_id = Column(Integer)  # Links to cannabinoids table
    predicted_effects = Column(JSON)  # From therapeutic model
    predicted_conditions = Column(JSON)  # Conditions this may help
    
    # Timestamps
    created_at = Column(DateTime, default=datetime.utcnow)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow)
    
    # Relationships
    dispensary = relationship("NevadaDispensary", back_populates="products")


class NevadaComplianceRule(Base):
    """Nevada regulatory compliance rules."""
    __tablename__ = "nevada_compliance_rules"
    
    id = Column(Integer, primary_key=True, index=True)
    
    # Rule Information
    rule_code = Column(String(50), unique=True, nullable=False)
    nac_reference = Column(String(100))  # e.g., "NAC 453D.810"
    title = Column(String(500), nullable=False)
    description = Column(Text)
    
    # Rule Parameters
    rule_type = Column(String(100))  # limit, requirement, prohibition, labeling
    applies_to = Column(JSON)  # ["flower", "edible", "concentrate"]
    
    # Limits (if applicable)
    limit_type = Column(String(50))  # max_thc, min_cbd, max_serving
    limit_value = Column(Float)
    limit_unit = Column(String(20))  # percent, mg, count
    
    # Enforcement
    penalty_type = Column(String(100))  # warning, fine, suspension, revocation
    penalty_amount = Column(Float)
    
    # Status
    effective_date = Column(DateTime)
    expiry_date = Column(DateTime)
    is_active = Column(Boolean, default=True)
    
    # Timestamps
    created_at = Column(DateTime, default=datetime.utcnow)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow)


class RealWorldEvidence(Base):
    """Real-world evidence from Nevada dispensary data.
    
    This is critical for BioPath validation and ClinPath optimization.
    """
    __tablename__ = "real_world_evidence"
    
    id = Column(Integer, primary_key=True, index=True)
    dispensary_id = Column(Integer, ForeignKey("nevada_dispensaries.id"))
    product_id = Column(Integer, ForeignKey("nevada_products.id"))
    
    # Evidence Type
    evidence_type = Column(String(100))  # sale, feedback, adverse_event, outcome
    
    # Anonymized Patient Info (no PII)
    patient_hash = Column(String(256))  # Hashed patient ID for tracking
    age_range = Column(String(20))  # "18-25", "26-35", etc.
    gender = Column(String(20))
    condition_category = Column(String(100))  # What they're treating
    
    # Product Usage
    cannabinoid_used = Column(String(50))
    dose_mg = Column(Float)
    frequency = Column(String(50))  # daily, weekly, as_needed
    duration_days = Column(Integer)
    
    # Outcome
    reported_outcome = Column(String(100))  # improved, no_change, worsened
    outcome_score = Column(Float)  # 0-10 patient-reported
    side_effects = Column(JSON)  # List of reported side effects
    
    # Quality Metrics
    confidence_score = Column(Float)  # Data quality score
    is_validated = Column(Boolean, default=False)
    validated_by = Column(String(100))  # Healer/clinician who validated
    
    # Provenance
    consent_verified = Column(Boolean, default=True)
    omnipath_manifest_id = Column(String(100))
    
    # Timestamps
    recorded_at = Column(DateTime, default=datetime.utcnow)
    created_at = Column(DateTime, default=datetime.utcnow)
    
    # Relationships
    dispensary = relationship("NevadaDispensary", back_populates="rwe_entries")


# Nevada-specific THC limits by product type (NAC 453D)
NEVADA_THC_LIMITS = {
    "flower": {"max_thc_percent": None},  # No limit for flower
    "concentrate": {"max_thc_percent": None},  # No limit for concentrates
    "edible": {"max_thc_mg_per_serving": 10, "max_thc_mg_per_package": 100},
    "tincture": {"max_thc_mg_per_serving": 10, "max_thc_mg_per_package": 100},
    "topical": {"max_thc_percent": None},  # No limit for topicals
}

# Nevada testing requirements
NEVADA_TESTING_REQUIREMENTS = [
    "potency",
    "terpenes",
    "residual_solvents",
    "pesticides",
    "heavy_metals",
    "microbials",
    "mycotoxins",
    "moisture_content",
    "water_activity",
]

# Nevada pilot dispensary quota
NEVADA_PILOT_CONFIG = {
    "target_dispensaries": 5,
    "pilot_start_date": "2026-03-26",
    "pilot_end_date": "2026-09-30",
    "min_rwe_entries_per_dispensary": 100,
    "required_conditions": [
        "chronic_pain",
        "anxiety",
        "ptsd",
        "insomnia",
        "nausea"
    ]
}
