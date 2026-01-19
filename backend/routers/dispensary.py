"""
Dispensary API Router - NeuroBotanica Nevada Pilot

Use Case 4: Dispensary Personalized Recommendation System

Endpoints:
- POST /api/dispensary/recommend: Generate personalized product recommendations
- POST /api/dispensary/profile: Create/update customer profile
- POST /api/dispensary/profile/inflammatory: Create profile with inflammatory biomarkers
- GET /api/dispensary/profile/{profile_id}: Get customer profile
- POST /api/dispensary/feedback: Submit recommendation feedback
- GET /api/dispensary/statistics: Get dispensary analytics
- POST /api/dispensary/inflammatory-synergy: TS-PS-001 cross-kingdom synergy prediction
- POST /api/dispensary/adjuvants/optimize: Optimize adjuvant selection

Reference: NeuroBotanica Patent - Section 6: Use Case 4
TS-PS-001: Cross-Kingdom Inflammatory Synergy Engine
"""

from typing import Dict, List, Optional, Any, Tuple
from fastapi import APIRouter, HTTPException, Depends
from pydantic import BaseModel, Field
from sqlalchemy.orm import Session
from datetime import datetime
import uuid
import logging

from backend.models.patient import Patient
from backend.models.database import get_db
from backend.services.engine_loader import load_inflammatory_engine
from backend.services.adjuvant_optimizer import get_adjuvant_optimizer, PatientProfile as AdjuvantPatientProfile
from backend.dependencies.auth import get_current_user, User
from backend.services.mapping import get_cannabinoid_scores_for_condition
from backend.models.database import SessionLocal
from backend.models.study import ClinicalStudy
from backend.services.confidence import compute_confidence_for_study, compute_overall_confidence

logger = logging.getLogger(__name__)
router = APIRouter(tags=["Dispensary"])


# =============================================================================
# Request/Response Models
# =============================================================================

class ConditionInput(BaseModel):
    """Patient condition with severity."""
    name: str = Field(..., description="Condition name (e.g., chronic_pain, anxiety)")
    severity: int = Field(..., ge=1, le=10, description="Severity score 1-10")
    is_primary: bool = Field(default=False, description="Is this the primary condition?")


class MedicationInput(BaseModel):
    """Current medication information."""
    name: str
    drug_class: Optional[str] = None


class NegativeExperience(BaseModel):
    """Previous negative cannabis experience."""
    strain: Optional[str] = None
    product_type: Optional[str] = None
    issue: str = Field(..., description="What went wrong (e.g., increased_anxiety, paranoia)")


class PositiveExperience(BaseModel):
    """Previous positive cannabis experience."""
    strain: Optional[str] = None
    product_type: Optional[str] = None
    effect: str = Field(..., description="What worked well")


class InflammatoryBiomarkers(BaseModel):
    """Inflammatory biomarkers for TS-PS-001 synergy prediction."""
    tnf_alpha: Optional[float] = Field(None, ge=0, le=100, description="TNF-α level (pg/mL)")
    il6: Optional[float] = Field(None, ge=0, le=50, description="IL-6 level (pg/mL)")
    crp: Optional[float] = Field(None, ge=0, le=20, description="CRP level (mg/L)")
    il1b: Optional[float] = Field(None, ge=0, le=10, description="IL-1β level (pg/mL)")
    esr: Optional[float] = Field(None, ge=0, le=100, description="ESR (mm/hr)")
    fibrinogen: Optional[float] = Field(None, ge=0, le=700, description="Fibrinogen (mg/dL)")
    homocysteine: Optional[float] = Field(None, ge=0, le=50, description="Homocysteine (µmol/L)")


class InflammatoryProfileInput(BaseModel):
    """Enhanced profile input with inflammatory biomarkers for TS-PS-001."""
    # Basic demographics (from existing TabletProfileInput)
    age: Optional[int] = Field(None, ge=18, le=120)
    biological_sex: Optional[str] = Field(None, description="male, female, other, unspecified")

    # Medical conditions
    conditions: List[ConditionInput] = Field(default_factory=list)

    # Inflammatory biomarkers (TS-PS-001 addition)
    biomarkers: Optional[InflammatoryBiomarkers] = None

    # Cannabis experience
    experience_level: str = Field(
        default="beginner",
        description="naive, beginner, intermediate, regular, experienced"
    )

    # Preferences
    administration_preferences: List[str] = Field(
        default_factory=lambda: ["inhalation"],
        description="inhalation, oral, topical"
    )

    # Goals
    primary_goal: Optional[str] = Field(None, description="pain_relief, anxiety, sleep, inflammation, etc.")


class CustomerProfileInput(BaseModel):
    """Full customer profile for recommendation engine."""
    # Demographics
    age: int = Field(..., ge=18, le=120)
    weight_kg: float = Field(..., gt=0)
    sex: Optional[str] = Field(None, description="male, female, other")
    
    # Medical
    conditions: List[ConditionInput] = Field(default_factory=list)
    current_medications: List[str] = Field(default_factory=list)
    supplements: List[str] = Field(default_factory=list)
    allergies: List[str] = Field(default_factory=list)
    
    # Cannabis Experience
    experience_level: str = Field(
        default="occasional",
        description="naive, occasional, regular, experienced"
    )
    tolerance_level: str = Field(default="medium", description="low, medium, high")
    
    # Preferences
    time_of_use: str = Field(default="evening", description="morning, afternoon, evening, night")
    must_remain_functional: bool = Field(default=False)
    administration_preferences: List[str] = Field(
        default_factory=lambda: ["vape", "tincture"],
        description="vape, flower, edible, tincture, topical"
    )
    
    # Previous Experiences
    negative_experiences: List[NegativeExperience] = Field(default_factory=list)
    positive_experiences: List[PositiveExperience] = Field(default_factory=list)

    # Goals
    primary_goal: Optional[str] = Field(None, description="pain_relief, anxiety, sleep, inflammation, etc.")


class ProductInput(BaseModel):
    """Dispensary product from inventory."""
    product_id: str
    product_name: str
    product_type: str = Field(..., description="flower, vape, edible, tincture, topical")
    strain_name: Optional[str] = None
    cbd_percent: float = Field(default=0, ge=0)
    thc_percent: float = Field(default=0, ge=0)
    cbg_percent: float = Field(default=0, ge=0)
    cbn_percent: float = Field(default=0, ge=0)
    terpenes: Optional[Dict[str, float]] = Field(default=None)
    price: Optional[float] = None


class RecommendationRequest(BaseModel):
    """Request for personalized recommendations."""
    customer_profile: CustomerProfileInput
    available_inventory: List[ProductInput] = Field(
        default_factory=list,
        description="Dispensary's current product catalog"
    )
    max_recommendations: int = Field(default=3, ge=1, le=10)


class DosageGuidance(BaseModel):
    """Dosage guidance for a product."""
    starting_dose: str
    wait_time_minutes: int
    max_dose_per_session: str
    frequency: str


class AdjuvantOptimization(BaseModel):
    """Adjuvant enhancement recommendation."""
    note: str
    timing: Optional[str] = None
    expected_synergy: Optional[str] = None
    recommended_addition: Optional[str] = None


class ProductRecommendation(BaseModel):
    """Single product recommendation with full details."""
    rank: int
    time_of_day: str
    product_id: str
    product_name: str
    match_score: float = Field(..., ge=0, le=1)
    cannabinoid_profile: Dict[str, Any]
    key_terpenes: List[Dict[str, Any]] = Field(default_factory=list)
    why_recommended: str
    expected_benefits: List[str]
    dosage_guidance: DosageGuidance
    confidence: float = Field(default=0.0, ge=0.0, le=1.0, description="Confidence score based on clinical evidence")
    adjuvant_optimization: Optional[AdjuvantOptimization] = None
    contraindications: Optional[str] = None


class ProductToAvoid(BaseModel):
    """Product to avoid with reason."""
    product_id: str
    product_name: str
    reason: str


class RecommendationResponse(BaseModel):
    """Full recommendation response."""
    recommendation_id: str
    customer_profile_id: Optional[str] = None
    recommendations: List[ProductRecommendation]
    products_to_avoid: List[ProductToAvoid] = Field(default_factory=list)
    education_notes: List[str] = Field(default_factory=list)
    follow_up_recommended: Optional[Dict[str, Any]] = None
    created_at: str


class FeedbackInput(BaseModel):
    """Feedback on a recommendation."""
    recommendation_id: str
    product_id: str
    purchased: bool
    effectiveness_score: Optional[int] = Field(None, ge=1, le=10)
    usage_frequency: Optional[str] = None
    side_effects: List[str] = Field(default_factory=list)
    customer_feedback: Optional[str] = None


class ProfileResponse(BaseModel):
    """Customer profile response."""
    profile_id: str
    created_at: str
    completeness_score: float
    primary_condition: Optional[str] = None
    recommendation_count: int = 0


# =============================================================================
# In-Memory Storage (for MVP - replace with database in production)
# =============================================================================

# Temporary storage for profiles and recommendations
_profiles: Dict[str, Dict[str, Any]] = {}
_recommendations: Dict[str, Dict[str, Any]] = {}
_feedback: Dict[str, List[Dict[str, Any]]] = {}


# =============================================================================
# Recommendation Engine
# =============================================================================

class DispensaryRecommendationEngine:
    """
    Core recommendation engine for dispensary use case.
    
    Matches customer profiles to products based on:
    - Condition-cannabinoid efficacy data
    - Terpene-condition correlations
    - Previous experience patterns
    - Adjuvant optimization
    """
    
    # Cannabinoid-condition efficacy scores (from 398 clinical studies)
    CONDITION_EFFICACY: Dict[str, Dict[str, float]] = {
        "chronic_pain": {"THC": 0.8, "CBD": 0.7, "CBG": 0.6, "CBN": 0.5},
        "anxiety": {"CBD": 0.85, "CBG": 0.6, "THC": 0.3, "CBN": 0.4},
        "insomnia": {"CBN": 0.85, "THC": 0.7, "CBD": 0.5, "CBG": 0.4},
        "inflammation": {"CBD": 0.8, "CBG": 0.75, "THC": 0.6, "CBC": 0.7},
        "ptsd": {"THC": 0.7, "CBD": 0.75, "CBG": 0.5},
        "nausea": {"THC": 0.9, "CBD": 0.6, "CBG": 0.5},
        "epilepsy": {"CBD": 0.95, "CBDV": 0.8, "THC": 0.3},
        "migraine": {"THC": 0.7, "CBD": 0.65, "CBG": 0.5},
        "arthritis": {"CBD": 0.75, "THC": 0.7, "CBG": 0.6},
        "appetite": {"THC": 0.85, "CBG": 0.6, "CBD": 0.3},
        "muscle_spasms": {"THC": 0.8, "CBD": 0.7, "CBG": 0.6, "CBN": 0.5},
        "seizures": {"CBD": 0.95, "CBDV": 0.8, "THC": 0.3},
        "weight_management": {"THC": 0.7, "CBG": 0.6, "CBD": 0.5},
    }
    
    # Terpene-condition correlations
    TERPENE_EFFECTS: Dict[str, Dict[str, Any]] = {
        "myrcene": {
            "effects": ["sedative", "muscle_relaxant", "anti_inflammatory"],
            "conditions": ["insomnia", "chronic_pain", "inflammation", "muscle_spasms"],
            "anxiety_risk": 0.1,
        },
        "limonene": {
            "effects": ["mood_elevation", "stress_relief", "energizing"],
            "conditions": ["anxiety", "depression"],
            "anxiety_risk": 0.2,
        },
        "beta_caryophyllene": {
            "effects": ["anti_inflammatory", "pain_relief", "cb2_activation"],
            "conditions": ["chronic_pain", "inflammation", "arthritis"],
            "anxiety_risk": 0.0,
        },
        "linalool": {
            "effects": ["anxiolytic", "sedative", "anticonvulsant"],
            "conditions": ["anxiety", "insomnia", "epilepsy", "seizures"],
            "anxiety_risk": 0.0,
        },
        "pinene": {
            "effects": ["alertness", "memory_retention", "bronchodilator"],
            "conditions": ["asthma", "focus"],
            "anxiety_risk": 0.3,
        },
        "terpinolene": {
            "effects": ["sedative", "antioxidant"],
            "conditions": ["insomnia"],
            "anxiety_risk": 0.15,
        },
        "humulene": {
            "effects": ["appetite_suppressant", "anti_inflammatory"],
            "conditions": ["inflammation"],
            "anxiety_risk": 0.05,
        },
    }
    
    # Anxiety risk factors
    HIGH_ANXIETY_RISK_PATTERNS = {
        "high_thc_low_cbd": lambda p: p.thc_percent > 20 and p.cbd_percent < 2,
        "high_limonene": lambda p: p.terpenes and p.terpenes.get("limonene", 0) > 1.5,
        "high_pinene": lambda p: p.terpenes and p.terpenes.get("pinene", 0) > 1.0,
    }
    
    def __init__(self):
        self.adjuvant_optimizer = get_adjuvant_optimizer()
    
    def generate_recommendations(
        self,
        profile: CustomerProfileInput,
        inventory: List[ProductInput],
        max_recommendations: int = 3,
    ) -> Tuple[List[ProductRecommendation], List[ProductToAvoid], List[str]]:
        """Generate personalized product recommendations."""
        
        # Identify primary condition
        primary_condition = None
        for cond in profile.conditions:
            if cond.is_primary:
                primary_condition = cond.name.lower().replace(" ", "_")
                break
        if not primary_condition and profile.conditions:
            primary_condition = profile.conditions[0].name.lower().replace(" ", "_")
        
        # Score all products
        scored_products = []
        avoid_products = []
        
        for product in inventory:
            score, reasons, warnings = self._score_product(product, profile, primary_condition)
            
            # Check for avoid reasons
            avoid_reason = self._should_avoid(product, profile)
            if avoid_reason:
                avoid_products.append(ProductToAvoid(
                    product_id=product.product_id,
                    product_name=product.product_name,
                    reason=avoid_reason,
                ))
                continue
            
            scored_products.append((product, score, reasons, warnings))
        
        # Sort by score
        scored_products.sort(key=lambda x: x[1], reverse=True)
        
        # Take top recommendations
        recommendations = []
        for rank, (product, score, reasons, warnings) in enumerate(scored_products[:max_recommendations], 1):
            rec = self._build_recommendation(
                rank, product, score, reasons, warnings, profile, primary_condition
            )
            recommendations.append(rec)
        
        # Generate education notes
        education_notes = self._generate_education_notes(profile, primary_condition)
        
        return recommendations, avoid_products, education_notes
    
    def _score_product(
        self,
        product: ProductInput,
        profile: CustomerProfileInput,
        primary_condition: Optional[str],
    ) -> Tuple[float, List[str], List[str]]:
        """Score a product for a customer profile."""
        score = 0.0
        reasons = []
        warnings = []
        
        # 1. Cannabinoid-condition match (0-40 points)
        if primary_condition:
            # Prefer the curated mapping file; fall back to built-in table
            efficacy = get_cannabinoid_scores_for_condition(primary_condition) or self.CONDITION_EFFICACY.get(primary_condition, {})

            # THC contribution
            if product.thc_percent > 0:
                thc_score = efficacy.get("THC", 0.5) * min(product.thc_percent / 25, 1) * 20
                score += thc_score

            # CBD contribution
            if product.cbd_percent > 0:
                cbd_score = efficacy.get("CBD", 0.5) * min(product.cbd_percent / 25, 1) * 20
                score += cbd_score
                reasons.append(f"CBD beneficial for {primary_condition}")

            # CBG/CBN contribution
            if product.cbg_percent > 0:
                score += efficacy.get("CBG", 0.4) * 5
            if product.cbn_percent > 0:
                score += efficacy.get("CBN", 0.4) * 5
        
        # 2. Terpene match (0-20 points)
        if product.terpenes:
            for terp, percent in product.terpenes.items():
                terp_lower = terp.lower().replace("-", "_")
                if terp_lower in self.TERPENE_EFFECTS:
                    terp_data = self.TERPENE_EFFECTS[terp_lower]
                    if primary_condition in terp_data.get("conditions", []):
                        score += percent * 5
                        reasons.append(f"{terp} aids {primary_condition}")
        
        # 3. CBD:THC ratio for anxiety-prone profiles (0-15 points)
        has_anxiety = any(c.name.lower() == "anxiety" for c in profile.conditions)
        if has_anxiety:
            if product.cbd_percent > 0 and product.thc_percent > 0:
                ratio = product.cbd_percent / max(product.thc_percent, 0.1)
                if ratio >= 2:
                    score += 15
                    reasons.append("High CBD:THC ratio reduces anxiety risk")
                elif ratio >= 1:
                    score += 10
        
        # 4. Administration preference match (0-10 points)
        if product.product_type.lower() in [p.lower() for p in profile.administration_preferences]:
            score += 10
        
        # 5. Functional requirement check (0-10 points)
        if profile.must_remain_functional:
            if product.thc_percent < 10:
                score += 10
                reasons.append("Low THC maintains functionality")
            elif product.thc_percent > 20:
                score -= 15
                warnings.append("High THC may impair functionality")
        
        # 6. Time of day optimization (0-10 points)
        if profile.time_of_use in ["evening", "night"]:
            if product.terpenes and product.terpenes.get("myrcene", 0) > 0.5:
                score += 10
                reasons.append("Sedative terpenes for evening use")
        elif profile.time_of_use in ["morning", "afternoon"]:
            if product.thc_percent < 15 and product.cbd_percent > 5:
                score += 10
        
        # 7. Experience level adjustment
        if profile.experience_level == "naive":
            if product.thc_percent > 15:
                score -= 10
                warnings.append("Consider lower THC for new users")
        
        # Normalize to 0-1
        score = max(0, min(100, score)) / 100
        
        return score, reasons, warnings
    
    def _should_avoid(
        self, product: ProductInput, profile: CustomerProfileInput
    ) -> Optional[str]:
        """Check if product should be avoided based on profile."""
        
        # Check negative experiences
        for neg in profile.negative_experiences:
            if neg.strain and product.strain_name:
                if neg.strain.lower() in product.strain_name.lower():
                    return f"Customer reported {neg.issue} with similar strain"
        
        # Check anxiety-inducing patterns for anxiety-prone customers
        has_anxiety_concern = any(
            c.name.lower() == "anxiety" or 
            any(n.issue == "increased_anxiety" for n in profile.negative_experiences)
            for c in profile.conditions
        )
        
        if has_anxiety_concern:
            if product.thc_percent > 25 and product.cbd_percent < 2:
                return "High THC without CBD may trigger anxiety"
        
        return None
    
    def _build_recommendation(
        self,
        rank: int,
        product: ProductInput,
        score: float,
        reasons: List[str],
        warnings: List[str],
        profile: CustomerProfileInput,
        primary_condition: Optional[str],
    ) -> ProductRecommendation:
        """Build a full product recommendation."""
        
        # Determine time of day
        if profile.time_of_use == "morning" or (profile.must_remain_functional and rank == 1):
            time_of_day = "daytime"
        elif profile.time_of_use in ["evening", "night"]:
            time_of_day = "evening"
        else:
            time_of_day = "anytime"
        
        # Build cannabinoid profile
        cannabinoid_profile = {
            "thc": product.thc_percent,
            "cbd": product.cbd_percent,
        }
        if product.cbg_percent:
            cannabinoid_profile["cbg"] = product.cbg_percent
        if product.cbn_percent:
            cannabinoid_profile["cbn"] = product.cbn_percent
        if product.thc_percent > 0 and product.cbd_percent > 0:
            cannabinoid_profile["ratio"] = f"{round(product.cbd_percent/product.thc_percent, 1)}:1 CBD:THC"
        
        # Build terpene info
        key_terpenes = []
        if product.terpenes:
            for terp, percent in sorted(product.terpenes.items(), key=lambda x: x[1], reverse=True)[:3]:
                terp_lower = terp.lower().replace("-", "_")
                effect = ""
                if terp_lower in self.TERPENE_EFFECTS:
                    effects = self.TERPENE_EFFECTS[terp_lower]["effects"]
                    effect = ", ".join(effects[:2])
                key_terpenes.append({
                    "name": terp,
                    "percent": percent,
                    "effect": effect,
                })
        
        # Why recommended
        why = ". ".join(reasons[:3]) if reasons else f"Matches {primary_condition} profile"
        
        # Expected benefits
        benefits = []
        if primary_condition:
            benefits.append(f"Targets {primary_condition.replace('_', ' ')}")
        if product.cbd_percent > 10:
            benefits.append("Anti-inflammatory and anxiolytic effects")
        if product.terpenes and product.terpenes.get("myrcene", 0) > 0.5:
            benefits.append("Muscle relaxation and sedation")
        
        # Dosage guidance based on product type
        dosage = self._get_dosage_guidance(product, profile)
        
        # Adjuvant optimization
        adjuvant_opt = None
        if primary_condition:
            adj_result = self.adjuvant_optimizer.optimize(
                primary_compound="CBD" if product.cbd_percent > product.thc_percent else "THC",
                therapeutic_target=primary_condition,
                patient_profile=AdjuvantPatientProfile(
                    age=profile.age,
                    weight_kg=profile.weight_kg,
                    conditions=[c.name for c in profile.conditions],
                    current_medications=profile.current_medications,
                    supplements=profile.supplements,
                    allergies=profile.allergies,
                ),
                max_adjuvants=1,
            )
            
            if adj_result.recommendations:
                top_adj = adj_result.recommendations[0]
                adjuvant_opt = AdjuvantOptimization(
                    note=f"Consider {top_adj.adjuvant_name} for enhanced effect",
                    timing=f"{abs(top_adj.timing_offset_minutes)} minutes before",
                    expected_synergy=f"{top_adj.expected_enhancement_percent}% enhancement",
                    recommended_addition=f"{top_adj.adjuvant_name} {top_adj.dosage_mg}mg",
                )
        
        # Contraindications
        contraindications = None
        if warnings:
            contraindications = ". ".join(warnings)
        
        # Compute per-recommendation confidence using clinical studies
        conf = self._compute_product_confidence(product, primary_condition)

        return ProductRecommendation(
            rank=rank,
            time_of_day=time_of_day,
            product_id=product.product_id,
            product_name=product.product_name,
            match_score=round(score, 2),
            cannabinoid_profile=cannabinoid_profile,
            key_terpenes=key_terpenes,
            why_recommended=why,
            expected_benefits=benefits,
            dosage_guidance=dosage,
            confidence=round(conf, 2),
            adjuvant_optimization=adjuvant_opt,
            contraindications=contraindications,
        )

    def _compute_product_confidence(self, product: ProductInput, primary_condition: Optional[str]) -> float:
        """Compute a confidence score for a product recommendation based on matching clinical studies.

        Strategy:
        - Identify the prominent cannabinoids in the product (THC, CBD, CBG, CBN, etc.)
        - Query the `clinical_studies` table for studies matching the `primary_condition` and those cannabinoids
        - Compute per-study confidence via `compute_confidence_for_study` and average
        - Fallback to `compute_overall_confidence(primary_condition)` when no specific studies found
        """
        try:
            if not primary_condition:
                return compute_overall_confidence()

            # Determine present cannabinoids
            present = []
            if product.thc_percent and product.thc_percent > 0:
                present.append("THC")
            if product.cbd_percent and product.cbd_percent > 0:
                present.append("CBD")
            if product.cbg_percent and product.cbg_percent > 0:
                present.append("CBG")
            if product.cbn_percent and product.cbn_percent > 0:
                present.append("CBN")

            db = SessionLocal()
            try:
                query = db.query(ClinicalStudy).filter(ClinicalStudy.condition == primary_condition.upper())
                if present:
                    query = query.filter(ClinicalStudy.cannabinoid.in_(present))
                studies = query.limit(200).all()
                if not studies:
                    return compute_overall_confidence(primary_condition)

                weights = [compute_confidence_for_study(s) for s in studies]
                return sum(weights) / len(weights)
            finally:
                db.close()
        except Exception:
            return compute_overall_confidence(primary_condition)
    
    def _get_dosage_guidance(
        self, product: ProductInput, profile: CustomerProfileInput
    ) -> DosageGuidance:
        """Get dosage guidance based on product type and profile."""
        
        product_type = product.product_type.lower()
        experience = profile.experience_level
        
        if product_type in ["vape", "flower"]:
            if experience == "naive":
                return DosageGuidance(
                    starting_dose="1 small puff",
                    wait_time_minutes=15,
                    max_dose_per_session="2-3 puffs",
                    frequency="As needed",
                )
            else:
                return DosageGuidance(
                    starting_dose="1-2 puffs",
                    wait_time_minutes=10,
                    max_dose_per_session="4-5 puffs",
                    frequency="As needed",
                )
        
        elif product_type == "edible":
            if experience == "naive":
                return DosageGuidance(
                    starting_dose="2.5mg THC",
                    wait_time_minutes=120,
                    max_dose_per_session="5mg THC",
                    frequency="Once daily",
                )
            else:
                return DosageGuidance(
                    starting_dose="5-10mg THC",
                    wait_time_minutes=90,
                    max_dose_per_session="20mg THC",
                    frequency="As needed, max twice daily",
                )
        
        elif product_type == "tincture":
            return DosageGuidance(
                starting_dose="0.25 mL sublingual",
                wait_time_minutes=30,
                max_dose_per_session="1 mL",
                frequency="2-3 times daily as needed",
            )
        
        else:  # topical or other
            return DosageGuidance(
                starting_dose="Apply thin layer to affected area",
                wait_time_minutes=30,
                max_dose_per_session="Reapply as needed",
                frequency="Up to 4 times daily",
            )
    
    def _generate_education_notes(
        self, profile: CustomerProfileInput, primary_condition: Optional[str]
    ) -> List[str]:
        """Generate educational notes for the customer."""
        notes = []
        
        if primary_condition:
            notes.append(f"CBD:THC ratio is more important than total THC for {primary_condition.replace('_', ' ')} management")
        
        # Supplement synergy
        if "magnesium" in [s.lower() for s in profile.supplements]:
            notes.append("Your magnesium supplementation enhances cannabinoid effectiveness by 25-35%")
        
        # Terpene education
        notes.append("Terpenes like beta-caryophyllene and myrcene specifically target pain pathways")
        
        # Consistency note
        if primary_condition in ["chronic_pain", "anxiety", "insomnia"]:
            notes.append("Consistent daily use often works better than as-needed dosing for chronic conditions")
        
        return notes


# Singleton engine
_engine: Optional[DispensaryRecommendationEngine] = None


def get_recommendation_engine() -> DispensaryRecommendationEngine:
    """Get or create recommendation engine singleton."""
    global _engine
    if _engine is None:
        _engine = DispensaryRecommendationEngine()
    return _engine


# =============================================================================
# API Endpoints
# =============================================================================

@router.post("/recommend", response_model=RecommendationResponse)
async def generate_recommendations(
    request: RecommendationRequest,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
) -> RecommendationResponse:
    """
    Generate personalized product recommendations.
    
    This is the core dispensary recommendation endpoint that matches
    customer profiles to available products using clinical evidence
    from 398 studies.
    
    Pricing: Included in dispensary subscription ($500-1,500/month)
    """
    try:
        engine = get_recommendation_engine()
        
        recommendations, avoid_products, education_notes = engine.generate_recommendations(
            profile=request.customer_profile,
            inventory=request.available_inventory,
            max_recommendations=request.max_recommendations,
        )
        
        recommendation_id = f"rec_{uuid.uuid4().hex[:12]}"
        
        # Store recommendation for feedback tracking
        _recommendations[recommendation_id] = {
            "profile": request.customer_profile.model_dump(),
            "recommendations": [r.model_dump() for r in recommendations],
            "created_at": datetime.utcnow().isoformat(),
        }
        
        return RecommendationResponse(
            recommendation_id=recommendation_id,
            recommendations=recommendations,
            products_to_avoid=avoid_products,
            education_notes=education_notes,
            follow_up_recommended={
                "days": 14,
                "reason": "Assess effectiveness and optimize dosing",
            },
            created_at=datetime.utcnow().isoformat(),
        )
        
    except Exception as e:
        logger.error(f"Recommendation generation failed: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/profile", response_model=ProfileResponse)
async def create_customer_profile(
    profile: CustomerProfileInput,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
) -> ProfileResponse:
    """
    Create or update a customer profile.

    Profiles are stored in the database and used for personalized recommendations.
    """
    try:
        # Generate unique profile code
        profile_code = f"NB-{uuid.uuid4().hex[:5].upper()}"

        # Create new patient record
        patient = Patient(
            patient_hash=f"tablet_{uuid.uuid4().hex}",  # Temporary hash for tablet-created profiles
            profile_code=profile_code,
            age_range=profile.age and f"{(profile.age//10)*10}-{((profile.age//10)*10)+9}" or None,  # Convert to range
            biological_sex=(getattr(profile, "biological_sex", None) or getattr(profile, "sex", None) or "unspecified"),
            state_code="NV",  # Nevada pilot
            jurisdiction_type="recreational",
            primary_conditions=[cond.name for cond in profile.conditions if cond.is_primary],
            secondary_conditions=[cond.name for cond in profile.conditions if not cond.is_primary],
            condition_severity={cond.name: cond.severity for cond in profile.conditions},
            cannabinoid_experience_level=profile.experience_level,
            preferred_delivery_methods=profile.administration_preferences,
            primary_treatment_goal=profile.primary_goal,
            omnipath_consent_status="limited_consent",  # Default for dispensary
            status="active"
        )

        # Calculate completeness
        patient.profile_completeness = patient.calculate_completeness()

        # Save to database
        db.add(patient)
        db.commit()
        db.refresh(patient)

        primary_condition = patient.primary_conditions[0] if patient.primary_conditions else None

        return ProfileResponse(
            profile_id=f"prof_{uuid.uuid4().hex[:6].upper()}",
            created_at=patient.created_at.isoformat() if patient.created_at else datetime.utcnow().isoformat(),
            completeness_score=patient.profile_completeness,
            primary_condition=primary_condition,
        )
    except Exception as e:
        logger.error(f"Profile creation failed: {e}")
        db.rollback()  # Rollback on error
        raise HTTPException(status_code=500, detail=f"Failed to create customer profile: {str(e)}")

    return {
        "profile": patient.to_recommendation_profile(),
        "created_at": patient.created_at.isoformat() if patient.created_at else None,
        "completeness": patient.profile_completeness,
        "status": patient.status
    }


@router.post("/feedback")
async def submit_feedback(
    feedback: FeedbackInput,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
) -> Dict[str, str]:
    """
    Submit feedback on a recommendation.
    
    Feedback is used to improve model accuracy and personalization.
    """
    try:
        if feedback.recommendation_id not in _recommendations:
            raise HTTPException(status_code=404, detail="Recommendation not found")
        
        # Store feedback
        if feedback.recommendation_id not in _feedback:
            _feedback[feedback.recommendation_id] = []
        
        _feedback[feedback.recommendation_id].append(feedback.model_dump())
        
        return {"status": "Feedback recorded", "recommendation_id": feedback.recommendation_id}
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Feedback submission failed: {e}")
        raise HTTPException(status_code=500, detail=f"Failed to submit feedback: {str(e)}")


@router.get("/statistics")
async def get_dispensary_statistics(
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
) -> Dict[str, Any]:
    """Get dispensary recommendation statistics."""
    try:
        total_profiles = db.query(Patient).count()

        return {
            "total_profiles": total_profiles,
            "total_recommendations": len(_recommendations),  # Keep in-memory for now
            "total_feedback_entries": sum(len(f) for f in _feedback.values()),  # Keep in-memory for now
            "clinical_studies_used": 398,
            "conditions_covered": 22,
            "trade_secret_engines": 6,
            "pricing": {
                "single_location": "$1,200/month",
                "multi_location": "$900/month per location",
                "enterprise": "$700/month per location",
            },
        }
    except Exception as e:
        logger.error(f"Statistics retrieval failed: {e}")
        raise HTTPException(status_code=500, detail=f"Failed to retrieve statistics: {str(e)}")


# =============================================================================
# Adjuvant Optimization Endpoint
# =============================================================================

class AdjuvantRequest(BaseModel):
    """Request for adjuvant optimization."""
    primary_compound: str = Field(..., description="Primary cannabinoid (CBD, THC, CBN, etc.)")
    therapeutic_target: str = Field(..., description="Target condition (insomnia, chronic_pain, etc.)")
    patient_profile: Optional[CustomerProfileInput] = None


class AdjuvantResponse(BaseModel):
    """Response with adjuvant recommendations."""
    recommended_adjuvant: str
    dosage_mg: float
    timing_offset_minutes: int
    expected_enhancement_percent: float
    confidence_interval: List[float]
    mechanism: str
    evidence_tier: int
    citations: List[str]
    protocol_summary: str


@router.post("/adjuvants/optimize", response_model=AdjuvantResponse)
async def optimize_adjuvants(
    request: AdjuvantRequest,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
) -> AdjuvantResponse:
    """
    Optimize adjuvant selection for a cannabinoid therapy.
    
    Implements patent claims 14-17: Adjuvant Enhancement.
    
    Example: Magnesium glycinate 400mg taken 45 minutes before
    a CBD product can enhance insomnia relief by 50%.
    """
    try:
        optimizer = get_adjuvant_optimizer()
        
        patient = None
        if request.patient_profile:
            patient = AdjuvantPatientProfile(
                age=request.patient_profile.age,
                weight_kg=request.patient_profile.weight_kg,
                conditions=[c.name for c in request.patient_profile.conditions],
                current_medications=request.patient_profile.current_medications,
                supplements=request.patient_profile.supplements,
                allergies=request.patient_profile.allergies,
            )
        
        result = optimizer.optimize(
            primary_compound=request.primary_compound,
            therapeutic_target=request.therapeutic_target,
            patient_profile=patient,
            max_adjuvants=1,
        )
        
        if not result.recommendations:
            raise HTTPException(
                status_code=404,
                detail=f"No adjuvants found for {request.therapeutic_target}"
            )
        
        top = result.recommendations[0]
        
        return AdjuvantResponse(
            recommended_adjuvant=top.adjuvant_name,
            dosage_mg=top.dosage_mg,
            timing_offset_minutes=top.timing_offset_minutes,
            expected_enhancement_percent=top.expected_enhancement_percent,
            confidence_interval=list(top.confidence_interval),
            mechanism=top.mechanism,
            evidence_tier=top.evidence_tier,
            citations=top.citations,
            protocol_summary=result.protocol_summary,
        )
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Adjuvant optimization failed: {e}")
        raise HTTPException(status_code=500, detail=str(e))


# =============================================================================
# TS-PS-001 Cross-Kingdom Inflammatory Synergy Engine Integration
# =============================================================================

class InflammatorySynergyRequest(BaseModel):
    """Request for TS-PS-001 inflammatory synergy prediction."""
    biomarkers: InflammatoryBiomarkers = Field(..., description="Inflammatory biomarker panel")
    condition_profile: InflammatoryProfileInput = Field(..., description="Patient inflammatory profile")
    available_kingdoms: List[str] = Field(
        default_factory=lambda: ["cannabis", "fungal", "marine", "plant"],
        description="Kingdoms to consider for synergy prediction"
    )


class InflammatorySynergyResponse(BaseModel):
    """TS-PS-001 inflammatory synergy prediction response."""
    primary_kingdom: str
    secondary_kingdoms: List[str]
    synergy_score: float
    confidence_level: float
    recommended_compounds: List[str]
    dosing_guidance: Dict[str, Any]
    expected_reduction: Dict[str, float]
    warning: Optional[str] = None


@router.post("/inflammatory-synergy", response_model=InflammatorySynergyResponse)
async def predict_inflammatory_synergy(
    request: InflammatorySynergyRequest,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
) -> InflammatorySynergyResponse:
    """
    TS-PS-001 Cross-Kingdom Inflammatory Synergy Prediction

    Predicts optimal anti-inflammatory formulations using proprietary algorithms
    that analyze biomarker profiles across Cannabis, Fungal, Marine, and Plant kingdoms.

    **TRADE SECRET PROTECTION**: This endpoint implements TS-PS-001 algorithms
    protected under 18 U.S.C. § 1836 (DTSA). Access is logged and monitored.

    **Biomarker Requirements**: At minimum TNF-α, IL-6, CRP recommended for optimal predictions.

    **Reference**: TS-PS-001 Cross-Kingdom Inflammatory Synergy Engine Documentation
    """
    try:
        # Initialize TS-PS-001 engine with access verification (loader picks prod or stub)
        engine = load_inflammatory_engine()

        # Extract biomarkers for prediction
        biomarkers = {
            'tnf_alpha': request.biomarkers.tnf_alpha or 0,
            'il6': request.biomarkers.il6 or 0,
            'crp': request.biomarkers.crp or 0,
            'il1b': request.biomarkers.il1b or 0,
        }

        # Get inflammatory synergy prediction
        prediction = engine.predict_inflammatory_synergy(
            biomarkers=biomarkers,
            condition_profile=request.condition_profile.model_dump(),
            available_kingdoms=request.available_kingdoms
        )

        return InflammatorySynergyResponse(**prediction)

    except PermissionError as e:
        logger.warning(f"TS-PS-001 access denied: {e}")
        raise HTTPException(status_code=403, detail="Access denied: TS-PS-001 trade secret protection")
    except Exception as e:
        logger.error(f"TS-PS-001 prediction failed: {e}")
        raise HTTPException(status_code=500, detail=f"Inflammatory synergy prediction failed: {str(e)}")


@router.post("/profile/inflammatory", response_model=ProfileResponse)
async def create_inflammatory_profile(
    profile: InflammatoryProfileInput,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
) -> ProfileResponse:
    """
    Create customer profile with inflammatory biomarkers for TS-PS-001 integration.

    Enhanced profile creation that includes inflammatory biomarker data for
    advanced synergy predictions using TS-PS-001 algorithms.
    """
    try:
        # Generate unique profile code
        profile_code = f"NB-IF-{uuid.uuid4().hex[:6].upper()}"

        # Create patient record with inflammatory biomarkers
        patient = Patient(
            patient_hash=f"inflammatory_{uuid.uuid4().hex}",
            profile_code=profile_code,
            age_range=profile.age and f"{(profile.age//10)*10}-{((profile.age//10)*10)+9}" or None,
            biological_sex=(getattr(profile, "biological_sex", None) or getattr(profile, "sex", None) or "unspecified"),
            state_code="NV",
            jurisdiction_type="recreational",
            primary_conditions=[cond.name for cond in profile.conditions if cond.is_primary],
            secondary_conditions=[cond.name for cond in profile.conditions if not cond.is_primary],
            condition_severity={cond.name: cond.severity for cond in profile.conditions},
            cannabinoid_experience_level=profile.experience_level,
            preferred_delivery_methods=profile.administration_preferences,
            primary_treatment_goal=profile.primary_goal,
            omnipath_consent_status="limited_consent",
            status="active"
        )

        # Add inflammatory biomarker data if provided
        if profile.biomarkers:
            # Store biomarkers in JSON field (would be normalized in production)
            biomarker_data = {
                'tnf_alpha': profile.biomarkers.tnf_alpha,
                'il6': profile.biomarkers.il6,
                'crp': profile.biomarkers.crp,
                'il1b': profile.biomarkers.il1b,
                'esr': profile.biomarkers.esr,
                'fibrinogen': profile.biomarkers.fibrinogen,
                'homocysteine': profile.biomarkers.homocysteine,
            }
            patient.inflammatory_biomarkers = biomarker_data

        # Calculate completeness including biomarkers
        patient.profile_completeness = patient.calculate_completeness()

        # Save to database
        db.add(patient)
        db.commit()
        db.refresh(patient)

        primary_condition = patient.primary_conditions[0] if patient.primary_conditions else None

        return ProfileResponse(
            profile_id=patient.profile_code,
            created_at=patient.created_at.isoformat() if patient.created_at else datetime.utcnow().isoformat(),
            completeness_score=patient.profile_completeness,
            primary_condition=primary_condition,
        )

    except Exception as e:
        logger.error(f"Inflammatory profile creation failed: {e}")
        db.rollback()
        raise HTTPException(status_code=500, detail=f"Failed to create inflammatory profile: {str(e)}")
