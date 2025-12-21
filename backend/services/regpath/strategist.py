"""
RegPath Strategist - Regulatory Pathway Optimization Engine

MVP Implementation for regulatory pathway analysis and strategy generation.
Trade secret: Decision matrices, scoring weights, and timeline templates are proprietary.
"""

from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
from typing import Optional, List, Dict, Any, Tuple
import uuid


class ProductType(Enum):
    """Product type classification."""
    ISOLATED_CANNABINOID = "isolated_cannabinoid"
    NOVEL_DIMER = "novel_dimer"
    FORMULATION = "formulation"
    COMBINATION = "combination"
    BOTANICAL_EXTRACT = "botanical_extract"
    SYNTHETIC = "synthetic"
    
    @property
    def display_name(self) -> str:
        names = {
            ProductType.ISOLATED_CANNABINOID: "Isolated Cannabinoid",
            ProductType.NOVEL_DIMER: "Novel Dimeric Cannabinoid",
            ProductType.FORMULATION: "Cannabinoid Formulation",
            ProductType.COMBINATION: "Combination Product",
            ProductType.BOTANICAL_EXTRACT: "Botanical Extract",
            ProductType.SYNTHETIC: "Synthetic Cannabinoid"
        }
        return names.get(self, self.value)


class NoveltyLevel(Enum):
    """Novelty level for regulatory assessment."""
    LOW = "low"
    MODERATE = "moderate"
    HIGH = "high"
    
    @property
    def description(self) -> str:
        descs = {
            NoveltyLevel.LOW: "Well-characterized compound with existing safety data",
            NoveltyLevel.MODERATE: "Known scaffold with novel modifications",
            NoveltyLevel.HIGH: "Novel chemical entity with limited prior data"
        }
        return descs.get(self, "")


class RegulatoryPathway(Enum):
    """FDA regulatory pathways."""
    IND = "IND"
    NDA = "NDA"
    ANDA = "ANDA"
    NDA_505B2 = "505(b)(2)"
    BLA = "BLA"
    OTC_MONOGRAPH = "OTC_Monograph"
    DIETARY_SUPPLEMENT = "Dietary_Supplement"
    
    @property
    def full_name(self) -> str:
        names = {
            RegulatoryPathway.IND: "Investigational New Drug Application",
            RegulatoryPathway.NDA: "New Drug Application",
            RegulatoryPathway.ANDA: "Abbreviated New Drug Application",
            RegulatoryPathway.NDA_505B2: "505(b)(2) New Drug Application",
            RegulatoryPathway.BLA: "Biologics License Application",
            RegulatoryPathway.OTC_MONOGRAPH: "OTC Drug Monograph",
            RegulatoryPathway.DIETARY_SUPPLEMENT: "Dietary Supplement Notification"
        }
        return names.get(self, self.value)
    
    @property
    def typical_timeline_months(self) -> Tuple[int, int]:
        """Return (min, max) typical timeline in months."""
        timelines = {
            RegulatoryPathway.IND: (12, 24),
            RegulatoryPathway.NDA: (36, 60),
            RegulatoryPathway.ANDA: (24, 48),
            RegulatoryPathway.NDA_505B2: (24, 42),
            RegulatoryPathway.BLA: (48, 72),
            RegulatoryPathway.OTC_MONOGRAPH: (18, 36),
            RegulatoryPathway.DIETARY_SUPPLEMENT: (3, 6)
        }
        return timelines.get(self, (24, 48))


class ChecklistCategory(Enum):
    """Readiness checklist categories."""
    CMC = "CMC"
    NONCLINICAL = "Nonclinical"
    CLINICAL = "Clinical"
    DOCUMENTATION = "Documentation"
    REGULATORY = "Regulatory"
    MANUFACTURING = "Manufacturing"


class MilestonePhase(Enum):
    """Timeline milestone phases."""
    PRECLINICAL = "Preclinical"
    IND_ENABLING = "IND-Enabling"
    PHASE_1 = "Phase 1"
    PHASE_2 = "Phase 2"
    PHASE_3 = "Phase 3"
    NDA_PREP = "NDA Preparation"
    FDA_REVIEW = "FDA Review"
    POST_APPROVAL = "Post-Approval"


@dataclass
class ProductProfile:
    """Product profile for regulatory assessment."""
    product_name: str
    product_type: str  # Will be converted to ProductType
    intended_use: str
    route: str  # oral/inhaled/sublingual/topical
    target_markets: List[str] = field(default_factory=lambda: ["US"])
    therapeutic_area: Optional[str] = None
    active_ingredients: Optional[List[str]] = None
    
    def to_dict(self) -> Dict:
        return {
            "product_name": self.product_name,
            "product_type": self.product_type,
            "intended_use": self.intended_use,
            "route": self.route,
            "target_markets": self.target_markets,
            "therapeutic_area": self.therapeutic_area,
            "active_ingredients": self.active_ingredients
        }


@dataclass
class EvidenceInputs:
    """References to ChemPath/ToxPath analyses."""
    chempath_job_id: Optional[str] = None
    toxpath_assessment_id: Optional[str] = None
    existing_literature_refs: Optional[List[str]] = None
    predicate_drug: Optional[str] = None  # For 505(b)(2)
    
    def to_dict(self) -> Dict:
        return {
            "chempath_job_id": self.chempath_job_id,
            "toxpath_assessment_id": self.toxpath_assessment_id,
            "existing_literature_refs": self.existing_literature_refs,
            "predicate_drug": self.predicate_drug
        }


@dataclass
class StrategyConstraints:
    """Constraints for strategy generation."""
    budget_range: Optional[str] = None  # e.g., "500K-2M"
    launch_window: Optional[str] = None  # e.g., "2026 Q4"
    risk_tolerance: str = "moderate"  # low/moderate/high
    preferred_pathway: Optional[str] = None
    excluded_pathways: List[str] = field(default_factory=list)
    
    def to_dict(self) -> Dict:
        return {
            "budget_range": self.budget_range,
            "launch_window": self.launch_window,
            "risk_tolerance": self.risk_tolerance,
            "preferred_pathway": self.preferred_pathway,
            "excluded_pathways": self.excluded_pathways
        }


@dataclass
class PathwayRecommendation:
    """Pathway recommendation with rationale."""
    pathway: RegulatoryPathway
    confidence: str  # high/moderate/low
    rationale: List[str]
    key_requirements: List[str]
    estimated_cost_range: str
    estimated_timeline_months: Tuple[int, int]
    
    def to_dict(self) -> Dict:
        return {
            "pathway": self.pathway.value,
            "pathway_full_name": self.pathway.full_name,
            "confidence": self.confidence,
            "rationale": self.rationale,
            "key_requirements": self.key_requirements,
            "estimated_cost_range": self.estimated_cost_range,
            "estimated_timeline_months": {
                "min": self.estimated_timeline_months[0],
                "max": self.estimated_timeline_months[1]
            }
        }


@dataclass
class ChecklistItem:
    """Readiness checklist item."""
    category: ChecklistCategory
    item: str
    description: str
    status: str = "not_started"  # not_started/in_progress/complete/na
    priority: str = "required"  # required/recommended/optional
    dependencies: List[str] = field(default_factory=list)
    estimated_effort: Optional[str] = None
    
    def to_dict(self) -> Dict:
        return {
            "category": self.category.value,
            "item": self.item,
            "description": self.description,
            "status": self.status,
            "priority": self.priority,
            "dependencies": self.dependencies,
            "estimated_effort": self.estimated_effort
        }


@dataclass
class TimelineMilestone:
    """Timeline milestone."""
    phase: MilestonePhase
    milestone: str
    description: str
    duration_months: int
    start_month: int  # Relative to project start
    dependencies: List[str] = field(default_factory=list)
    deliverables: List[str] = field(default_factory=list)
    
    def to_dict(self) -> Dict:
        return {
            "phase": self.phase.value,
            "milestone": self.milestone,
            "description": self.description,
            "duration_months": self.duration_months,
            "start_month": self.start_month,
            "end_month": self.start_month + self.duration_months,
            "dependencies": self.dependencies,
            "deliverables": self.deliverables
        }


@dataclass
class RegPathRequest:
    """Request model for RegPath strategy generation."""
    product_profile: ProductProfile
    evidence_inputs: Optional[EvidenceInputs] = None
    constraints: Optional[StrategyConstraints] = None
    
    def to_dict(self) -> Dict:
        return {
            "product_profile": self.product_profile.to_dict(),
            "evidence_inputs": self.evidence_inputs.to_dict() if self.evidence_inputs else None,
            "constraints": self.constraints.to_dict() if self.constraints else None
        }


@dataclass
class RegPathResponse:
    """Response model for RegPath strategy."""
    regpath_strategy_id: str
    product_profile: ProductProfile
    novelty_assessment: str
    evidence_strength: str
    primary_pathway: PathwayRecommendation
    fallback_pathway: Optional[PathwayRecommendation]
    gating_questions: List[str]
    readiness_checklist: List[ChecklistItem]
    timeline: List[TimelineMilestone]
    next_actions: List[str]
    key_assumptions: List[str]
    status: str = "completed"
    error: Optional[str] = None
    generated_at: str = field(default_factory=lambda: datetime.utcnow().isoformat())
    
    def to_dict(self) -> Dict:
        return {
            "regpath_strategy_id": self.regpath_strategy_id,
            "product_profile": self.product_profile.to_dict(),
            "novelty_assessment": self.novelty_assessment,
            "evidence_strength": self.evidence_strength,
            "primary_pathway": self.primary_pathway.to_dict(),
            "fallback_pathway": self.fallback_pathway.to_dict() if self.fallback_pathway else None,
            "gating_questions": self.gating_questions,
            "readiness_checklist": [item.to_dict() for item in self.readiness_checklist],
            "timeline": [m.to_dict() for m in self.timeline],
            "next_actions": self.next_actions,
            "key_assumptions": self.key_assumptions,
            "status": self.status,
            "error": self.error,
            "generated_at": self.generated_at
        }


class RegPathStrategist:
    """Regulatory Pathway Optimization Engine.
    
    MVP Implementation Notes:
    - Uses decision matrix for pathway selection
    - Generates readiness checklists per pathway
    - Builds timeline from templates
    
    Trade Secret Elements (not exposed in API):
    - Pathway decision matrix weights
    - Scoring algorithms
    - Timeline estimation formulas
    - Checklist prioritization logic
    """
    
    # =========================================================================
    # TRADE SECRET: Pathway Decision Matrix
    # Product type × Novelty × Evidence → Primary/Fallback pathway
    # =========================================================================
    
    _PATHWAY_MATRIX = {
        # (product_type, novelty, evidence_strength) -> (primary, fallback)
        ("isolated_cannabinoid", "low", "strong"): (RegulatoryPathway.NDA_505B2, RegulatoryPathway.NDA),
        ("isolated_cannabinoid", "low", "moderate"): (RegulatoryPathway.NDA_505B2, RegulatoryPathway.IND),
        ("isolated_cannabinoid", "low", "weak"): (RegulatoryPathway.IND, RegulatoryPathway.NDA_505B2),
        ("isolated_cannabinoid", "moderate", "strong"): (RegulatoryPathway.NDA_505B2, RegulatoryPathway.NDA),
        ("isolated_cannabinoid", "moderate", "moderate"): (RegulatoryPathway.IND, RegulatoryPathway.NDA_505B2),
        ("isolated_cannabinoid", "moderate", "weak"): (RegulatoryPathway.IND, None),
        ("isolated_cannabinoid", "high", "strong"): (RegulatoryPathway.IND, RegulatoryPathway.NDA_505B2),
        ("isolated_cannabinoid", "high", "moderate"): (RegulatoryPathway.IND, None),
        ("isolated_cannabinoid", "high", "weak"): (RegulatoryPathway.IND, None),
        
        ("novel_dimer", "low", "strong"): (RegulatoryPathway.IND, RegulatoryPathway.NDA_505B2),
        ("novel_dimer", "low", "moderate"): (RegulatoryPathway.IND, RegulatoryPathway.NDA_505B2),
        ("novel_dimer", "low", "weak"): (RegulatoryPathway.IND, None),
        ("novel_dimer", "moderate", "strong"): (RegulatoryPathway.IND, RegulatoryPathway.NDA_505B2),
        ("novel_dimer", "moderate", "moderate"): (RegulatoryPathway.IND, None),
        ("novel_dimer", "moderate", "weak"): (RegulatoryPathway.IND, None),
        ("novel_dimer", "high", "strong"): (RegulatoryPathway.IND, None),
        ("novel_dimer", "high", "moderate"): (RegulatoryPathway.IND, None),
        ("novel_dimer", "high", "weak"): (RegulatoryPathway.IND, None),
        
        ("formulation", "low", "strong"): (RegulatoryPathway.ANDA, RegulatoryPathway.NDA_505B2),
        ("formulation", "low", "moderate"): (RegulatoryPathway.NDA_505B2, RegulatoryPathway.NDA),
        ("formulation", "low", "weak"): (RegulatoryPathway.NDA_505B2, RegulatoryPathway.IND),
        ("formulation", "moderate", "strong"): (RegulatoryPathway.NDA_505B2, RegulatoryPathway.NDA),
        ("formulation", "moderate", "moderate"): (RegulatoryPathway.NDA_505B2, RegulatoryPathway.IND),
        ("formulation", "moderate", "weak"): (RegulatoryPathway.IND, RegulatoryPathway.NDA_505B2),
        ("formulation", "high", "strong"): (RegulatoryPathway.IND, RegulatoryPathway.NDA_505B2),
        ("formulation", "high", "moderate"): (RegulatoryPathway.IND, None),
        ("formulation", "high", "weak"): (RegulatoryPathway.IND, None),
        
        ("combination", "low", "strong"): (RegulatoryPathway.NDA_505B2, RegulatoryPathway.NDA),
        ("combination", "low", "moderate"): (RegulatoryPathway.NDA_505B2, RegulatoryPathway.IND),
        ("combination", "low", "weak"): (RegulatoryPathway.IND, RegulatoryPathway.NDA_505B2),
        ("combination", "moderate", "strong"): (RegulatoryPathway.NDA_505B2, RegulatoryPathway.NDA),
        ("combination", "moderate", "moderate"): (RegulatoryPathway.IND, RegulatoryPathway.NDA_505B2),
        ("combination", "moderate", "weak"): (RegulatoryPathway.IND, None),
        ("combination", "high", "strong"): (RegulatoryPathway.IND, RegulatoryPathway.NDA_505B2),
        ("combination", "high", "moderate"): (RegulatoryPathway.IND, None),
        ("combination", "high", "weak"): (RegulatoryPathway.IND, None),
        
        ("botanical_extract", "low", "strong"): (RegulatoryPathway.NDA_505B2, RegulatoryPathway.DIETARY_SUPPLEMENT),
        ("botanical_extract", "low", "moderate"): (RegulatoryPathway.NDA_505B2, RegulatoryPathway.IND),
        ("botanical_extract", "low", "weak"): (RegulatoryPathway.IND, RegulatoryPathway.DIETARY_SUPPLEMENT),
        ("botanical_extract", "moderate", "strong"): (RegulatoryPathway.IND, RegulatoryPathway.NDA_505B2),
        ("botanical_extract", "moderate", "moderate"): (RegulatoryPathway.IND, None),
        ("botanical_extract", "moderate", "weak"): (RegulatoryPathway.IND, None),
        ("botanical_extract", "high", "strong"): (RegulatoryPathway.IND, None),
        ("botanical_extract", "high", "moderate"): (RegulatoryPathway.IND, None),
        ("botanical_extract", "high", "weak"): (RegulatoryPathway.IND, None),
        
        ("synthetic", "low", "strong"): (RegulatoryPathway.NDA_505B2, RegulatoryPathway.NDA),
        ("synthetic", "low", "moderate"): (RegulatoryPathway.NDA_505B2, RegulatoryPathway.IND),
        ("synthetic", "low", "weak"): (RegulatoryPathway.IND, RegulatoryPathway.NDA_505B2),
        ("synthetic", "moderate", "strong"): (RegulatoryPathway.IND, RegulatoryPathway.NDA_505B2),
        ("synthetic", "moderate", "moderate"): (RegulatoryPathway.IND, None),
        ("synthetic", "moderate", "weak"): (RegulatoryPathway.IND, None),
        ("synthetic", "high", "strong"): (RegulatoryPathway.IND, None),
        ("synthetic", "high", "moderate"): (RegulatoryPathway.IND, None),
        ("synthetic", "high", "weak"): (RegulatoryPathway.IND, None),
    }
    
    # TRADE SECRET: Cost estimation ranges by pathway
    _COST_RANGES = {
        RegulatoryPathway.IND: "$2M - $10M",
        RegulatoryPathway.NDA: "$50M - $200M+",
        RegulatoryPathway.ANDA: "$3M - $15M",
        RegulatoryPathway.NDA_505B2: "$15M - $75M",
        RegulatoryPathway.BLA: "$100M - $500M+",
        RegulatoryPathway.OTC_MONOGRAPH: "$5M - $25M",
        RegulatoryPathway.DIETARY_SUPPLEMENT: "$50K - $500K"
    }
    
    def __init__(self):
        """Initialize RegPath strategist."""
        pass
    
    def generate_strategy(self, request: RegPathRequest) -> RegPathResponse:
        """Generate full regulatory strategy.
        
        Pipeline steps:
        1. Classify product novelty
        2. Assess evidence strength
        3. Select pathways via decision matrix
        4. Build readiness checklist
        5. Generate timeline
        6. Compile next actions
        """
        strategy_id = str(uuid.uuid4())
        
        # Step 1: Assess novelty
        novelty = self._assess_novelty(request.product_profile)
        
        # Step 2: Assess evidence strength
        evidence_strength = self._assess_evidence(request.evidence_inputs)
        
        # Step 3: Select pathways
        primary, fallback = self._select_pathways(
            request.product_profile,
            novelty,
            evidence_strength,
            request.constraints
        )
        
        # Step 4: Generate gating questions
        gating_questions = self._generate_gating_questions(
            request.product_profile,
            primary,
            fallback
        )
        
        # Step 5: Build readiness checklist
        checklist = self._build_checklist(
            primary.pathway,
            request.product_profile,
            request.evidence_inputs
        )
        
        # Step 6: Generate timeline
        timeline = self._build_timeline(
            primary.pathway,
            request.product_profile
        )
        
        # Step 7: Generate next actions
        next_actions = self._generate_next_actions(
            primary,
            checklist,
            request.evidence_inputs
        )
        
        # Step 8: Document assumptions
        assumptions = self._document_assumptions(
            request.product_profile,
            request.constraints
        )
        
        return RegPathResponse(
            regpath_strategy_id=strategy_id,
            product_profile=request.product_profile,
            novelty_assessment=novelty,
            evidence_strength=evidence_strength,
            primary_pathway=primary,
            fallback_pathway=fallback,
            gating_questions=gating_questions,
            readiness_checklist=checklist,
            timeline=timeline,
            next_actions=next_actions,
            key_assumptions=assumptions
        )
    
    def _assess_novelty(self, profile: ProductProfile) -> str:
        """Assess novelty level of product.
        
        TRADE SECRET: Novelty scoring algorithm.
        """
        product_type = profile.product_type.lower()
        
        # Novel dimers are always high novelty
        if "novel" in product_type or "dimer" in product_type:
            return "high"
        
        # Synthetic cannabinoids typically high novelty
        if "synthetic" in product_type:
            return "moderate"
        
        # Formulations of known compounds are lower novelty
        if "formulation" in product_type:
            return "low"
        
        # Isolated cannabinoids depend on the specific compound
        if "isolated" in product_type:
            # Check if it's a well-known cannabinoid
            known_cannabinoids = ["thc", "cbd", "cbg", "cbn", "cbda", "thca"]
            ingredients = profile.active_ingredients or []
            for ing in ingredients:
                if any(known in ing.lower() for known in known_cannabinoids):
                    return "low"
            return "moderate"
        
        # Botanical extracts moderate by default
        if "botanical" in product_type or "extract" in product_type:
            return "moderate"
        
        # Default to moderate
        return "moderate"
    
    def _assess_evidence(self, evidence: Optional[EvidenceInputs]) -> str:
        """Assess evidence strength.
        
        TRADE SECRET: Evidence scoring weights.
        """
        if not evidence:
            return "weak"
        
        score = 0
        
        # ChemPath analysis adds confidence
        if evidence.chempath_job_id:
            score += 2
        
        # ToxPath assessment adds confidence
        if evidence.toxpath_assessment_id:
            score += 2
        
        # Literature references add confidence
        if evidence.existing_literature_refs:
            score += min(len(evidence.existing_literature_refs), 3)
        
        # Predicate drug is strong evidence for 505(b)(2)
        if evidence.predicate_drug:
            score += 3
        
        if score >= 6:
            return "strong"
        elif score >= 3:
            return "moderate"
        else:
            return "weak"
    
    def _select_pathways(
        self,
        profile: ProductProfile,
        novelty: str,
        evidence: str,
        constraints: Optional[StrategyConstraints]
    ) -> Tuple[PathwayRecommendation, Optional[PathwayRecommendation]]:
        """Select primary and fallback pathways.
        
        TRADE SECRET: Decision matrix and scoring logic.
        """
        product_type = profile.product_type.lower().replace(" ", "_")
        
        # Normalize product type
        type_map = {
            "isolated": "isolated_cannabinoid",
            "cannabinoid": "isolated_cannabinoid",
            "dimer": "novel_dimer",
            "dimeric": "novel_dimer",
            "formula": "formulation",
            "combo": "combination",
            "botanical": "botanical_extract",
            "extract": "botanical_extract",
            "synth": "synthetic"
        }
        
        for key, val in type_map.items():
            if key in product_type:
                product_type = val
                break
        
        # Lookup in matrix
        key = (product_type, novelty, evidence)
        if key not in self._PATHWAY_MATRIX:
            # Default to IND for unknown combinations
            key = ("isolated_cannabinoid", "moderate", "moderate")
        
        primary_pathway, fallback_pathway = self._PATHWAY_MATRIX[key]
        
        # Apply constraints
        if constraints:
            if constraints.excluded_pathways:
                excluded = [p.lower() for p in constraints.excluded_pathways]
                if primary_pathway.value.lower() in excluded:
                    primary_pathway, fallback_pathway = fallback_pathway, None
        
        # Build recommendation objects
        primary = self._build_pathway_recommendation(
            primary_pathway, 
            profile, 
            novelty, 
            evidence,
            is_primary=True
        )
        
        fallback = None
        if fallback_pathway:
            fallback = self._build_pathway_recommendation(
                fallback_pathway,
                profile,
                novelty,
                evidence,
                is_primary=False
            )
        
        return primary, fallback
    
    def _build_pathway_recommendation(
        self,
        pathway: RegulatoryPathway,
        profile: ProductProfile,
        novelty: str,
        evidence: str,
        is_primary: bool
    ) -> PathwayRecommendation:
        """Build detailed pathway recommendation."""
        rationale = self._generate_rationale(pathway, profile, novelty, evidence, is_primary)
        requirements = self._get_pathway_requirements(pathway, profile)
        
        # Adjust confidence based on factors
        if evidence == "strong" and novelty == "low":
            confidence = "high"
        elif evidence == "weak" or novelty == "high":
            confidence = "low"
        else:
            confidence = "moderate"
        
        return PathwayRecommendation(
            pathway=pathway,
            confidence=confidence,
            rationale=rationale,
            key_requirements=requirements,
            estimated_cost_range=self._COST_RANGES.get(pathway, "$10M - $50M"),
            estimated_timeline_months=pathway.typical_timeline_months
        )
    
    def _generate_rationale(
        self,
        pathway: RegulatoryPathway,
        profile: ProductProfile,
        novelty: str,
        evidence: str,
        is_primary: bool
    ) -> List[str]:
        """Generate rationale for pathway selection."""
        rationale = []
        
        if pathway == RegulatoryPathway.IND:
            rationale.append("IND pathway required for investigational studies in humans")
            if novelty == "high":
                rationale.append("Novel compound requires comprehensive safety evaluation")
            rationale.append("Enables early clinical data generation to support future applications")
        
        elif pathway == RegulatoryPathway.NDA:
            rationale.append("Full NDA provides most comprehensive approval pathway")
            rationale.append("Appropriate for novel drugs without suitable predicate")
            if evidence == "strong":
                rationale.append("Strong evidence base supports full NDA submission")
        
        elif pathway == RegulatoryPathway.NDA_505B2:
            rationale.append("505(b)(2) leverages existing literature and approved drug data")
            if novelty == "low":
                rationale.append("Well-characterized compound allows reference to prior approvals")
            rationale.append("Can reduce development timeline and costs vs full NDA")
        
        elif pathway == RegulatoryPathway.ANDA:
            rationale.append("ANDA appropriate for generic equivalent of approved product")
            rationale.append("Requires demonstration of bioequivalence")
            rationale.append("Most cost-effective pathway for non-novel formulations")
        
        elif pathway == RegulatoryPathway.DIETARY_SUPPLEMENT:
            rationale.append("Dietary supplement pathway for non-drug claims")
            rationale.append("Suitable for botanical extracts with structure/function claims")
            rationale.append("Note: Cannot make disease treatment claims")
        
        if not is_primary:
            rationale.insert(0, "Fallback option if primary pathway proves infeasible")
        
        return rationale
    
    def _get_pathway_requirements(
        self,
        pathway: RegulatoryPathway,
        profile: ProductProfile
    ) -> List[str]:
        """Get key requirements for pathway."""
        requirements = []
        
        if pathway == RegulatoryPathway.IND:
            requirements = [
                "CMC data (identity, purity, stability)",
                "Pharmacology studies",
                "Toxicology studies (acute, repeat-dose)",
                "Clinical protocol and investigator brochure",
                "FDA Form 1571"
            ]
        
        elif pathway == RegulatoryPathway.NDA:
            requirements = [
                "Complete CMC package (manufacturing, controls, stability)",
                "Full nonclinical package (pharmacology, toxicology, carcinogenicity)",
                "Clinical efficacy data (typically Phase 3)",
                "Clinical safety database",
                "Risk management plan"
            ]
        
        elif pathway == RegulatoryPathway.NDA_505B2:
            requirements = [
                "CMC data for new formulation/use",
                "Bridging studies to reference product",
                "Literature review and right of reference",
                "May require limited clinical studies",
                "Bioavailability/bioequivalence data"
            ]
        
        elif pathway == RegulatoryPathway.ANDA:
            requirements = [
                "Pharmaceutical equivalence demonstration",
                "Bioequivalence studies",
                "CMC data",
                "Patent certification (Paragraph IV if applicable)"
            ]
        
        elif pathway == RegulatoryPathway.DIETARY_SUPPLEMENT:
            requirements = [
                "New Dietary Ingredient notification (if applicable)",
                "Manufacturing compliance with cGMP",
                "Labeling compliance",
                "Structure/function claim substantiation"
            ]
        
        # Route-specific additions
        if profile.route == "inhaled":
            requirements.append("Pulmonary safety assessment")
            requirements.append("Device characterization (if applicable)")
        
        return requirements
    
    def _generate_gating_questions(
        self,
        profile: ProductProfile,
        primary: PathwayRecommendation,
        fallback: Optional[PathwayRecommendation]
    ) -> List[str]:
        """Generate critical gating questions."""
        questions = []
        
        # General questions
        questions.append("Is the active ingredient well-characterized with validated analytical methods?")
        questions.append("Is there adequate manufacturing capability for clinical supply?")
        
        # Pathway-specific
        if primary.pathway == RegulatoryPathway.NDA_505B2:
            questions.append("Is there an appropriate reference listed drug (RLD) for 505(b)(2)?")
            questions.append("Can right of reference be obtained for key studies?")
        
        if primary.pathway == RegulatoryPathway.IND:
            questions.append("Are IND-enabling toxicology studies completed or planned?")
            questions.append("Is there a Phase 1 clinical site identified?")
        
        # Product-specific
        if "dimer" in profile.product_type.lower():
            questions.append("Is the synthetic route scalable and reproducible?")
            questions.append("What is the target indication and unmet medical need?")
        
        if profile.route == "inhaled":
            questions.append("Is the delivery device approved or development needed?")
        
        # Evidence gaps
        questions.append("What are the key gaps in the current data package?")
        
        return questions[:6]  # Limit to most critical
    
    def _build_checklist(
        self,
        pathway: RegulatoryPathway,
        profile: ProductProfile,
        evidence: Optional[EvidenceInputs]
    ) -> List[ChecklistItem]:
        """Build readiness checklist for pathway."""
        checklist = []
        
        # CMC items (all pathways)
        checklist.extend([
            ChecklistItem(
                category=ChecklistCategory.CMC,
                item="Drug Substance Specification",
                description="Establish identity, purity, and potency specifications",
                priority="required",
                estimated_effort="2-4 weeks"
            ),
            ChecklistItem(
                category=ChecklistCategory.CMC,
                item="Analytical Method Validation",
                description="Validate methods for identity, purity, potency testing",
                priority="required",
                estimated_effort="4-8 weeks"
            ),
            ChecklistItem(
                category=ChecklistCategory.CMC,
                item="Stability Studies",
                description="Initiate ICH stability studies",
                priority="required",
                estimated_effort="Ongoing (6-24 months)"
            ),
            ChecklistItem(
                category=ChecklistCategory.CMC,
                item="Manufacturing Process Description",
                description="Document manufacturing process and controls",
                priority="required",
                estimated_effort="4-6 weeks"
            )
        ])
        
        # Nonclinical items
        if pathway in [RegulatoryPathway.IND, RegulatoryPathway.NDA, RegulatoryPathway.NDA_505B2]:
            checklist.extend([
                ChecklistItem(
                    category=ChecklistCategory.NONCLINICAL,
                    item="Pharmacology Studies",
                    description="In vitro and in vivo pharmacology characterization",
                    priority="required",
                    estimated_effort="3-6 months"
                ),
                ChecklistItem(
                    category=ChecklistCategory.NONCLINICAL,
                    item="GLP Toxicology Studies",
                    description="IND-enabling toxicology package",
                    priority="required",
                    dependencies=["Drug Substance Specification"],
                    estimated_effort="6-12 months"
                ),
                ChecklistItem(
                    category=ChecklistCategory.NONCLINICAL,
                    item="ADME Studies",
                    description="Absorption, distribution, metabolism, excretion characterization",
                    priority="required",
                    estimated_effort="3-6 months"
                )
            ])
        
        # Clinical items
        if pathway == RegulatoryPathway.IND:
            checklist.extend([
                ChecklistItem(
                    category=ChecklistCategory.CLINICAL,
                    item="Clinical Protocol",
                    description="Draft Phase 1 clinical protocol",
                    priority="required",
                    estimated_effort="4-8 weeks"
                ),
                ChecklistItem(
                    category=ChecklistCategory.CLINICAL,
                    item="Investigator's Brochure",
                    description="Compile IB with available data",
                    priority="required",
                    estimated_effort="4-6 weeks"
                )
            ])
        
        if pathway == RegulatoryPathway.NDA_505B2:
            checklist.extend([
                ChecklistItem(
                    category=ChecklistCategory.CLINICAL,
                    item="Literature Review",
                    description="Comprehensive review of reference product literature",
                    priority="required",
                    estimated_effort="4-8 weeks"
                ),
                ChecklistItem(
                    category=ChecklistCategory.CLINICAL,
                    item="Bridging Study Protocol",
                    description="Design BA/BE or bridging studies as needed",
                    priority="required",
                    estimated_effort="2-4 weeks"
                )
            ])
        
        # Documentation items
        checklist.extend([
            ChecklistItem(
                category=ChecklistCategory.DOCUMENTATION,
                item="Regulatory Strategy Document",
                description="Document pathway strategy and rationale",
                priority="required",
                estimated_effort="1-2 weeks"
            ),
            ChecklistItem(
                category=ChecklistCategory.DOCUMENTATION,
                item="Pre-IND/Pre-NDA Meeting Request",
                description="Prepare briefing document for FDA meeting",
                priority="recommended",
                estimated_effort="4-6 weeks"
            )
        ])
        
        # Update status based on evidence
        if evidence:
            if evidence.chempath_job_id:
                for item in checklist:
                    if "Specification" in item.item:
                        item.status = "in_progress"
            if evidence.toxpath_assessment_id:
                for item in checklist:
                    if "Toxicology" in item.item:
                        item.status = "in_progress"
        
        return checklist
    
    def _build_timeline(
        self,
        pathway: RegulatoryPathway,
        profile: ProductProfile
    ) -> List[TimelineMilestone]:
        """Build development timeline."""
        timeline = []
        
        if pathway == RegulatoryPathway.IND:
            timeline = [
                TimelineMilestone(
                    phase=MilestonePhase.PRECLINICAL,
                    milestone="CMC Package Development",
                    description="Complete drug substance/product CMC work",
                    duration_months=6,
                    start_month=0,
                    deliverables=["Drug substance specification", "Analytical methods"]
                ),
                TimelineMilestone(
                    phase=MilestonePhase.IND_ENABLING,
                    milestone="GLP Toxicology Studies",
                    description="Complete IND-enabling toxicology package",
                    duration_months=9,
                    start_month=3,
                    dependencies=["CMC Package Development"],
                    deliverables=["GLP tox reports", "ADME data"]
                ),
                TimelineMilestone(
                    phase=MilestonePhase.IND_ENABLING,
                    milestone="IND Submission",
                    description="Submit IND application to FDA",
                    duration_months=2,
                    start_month=10,
                    dependencies=["GLP Toxicology Studies"],
                    deliverables=["Form 1571", "IND modules"]
                ),
                TimelineMilestone(
                    phase=MilestonePhase.PHASE_1,
                    milestone="Phase 1 Clinical Trial",
                    description="First-in-human safety/PK study",
                    duration_months=9,
                    start_month=12,
                    dependencies=["IND Submission"],
                    deliverables=["Phase 1 clinical report"]
                )
            ]
        
        elif pathway == RegulatoryPathway.NDA_505B2:
            timeline = [
                TimelineMilestone(
                    phase=MilestonePhase.PRECLINICAL,
                    milestone="Literature/Data Package",
                    description="Compile reference product data and literature",
                    duration_months=4,
                    start_month=0,
                    deliverables=["Literature review", "Right of reference agreements"]
                ),
                TimelineMilestone(
                    phase=MilestonePhase.PRECLINICAL,
                    milestone="CMC Development",
                    description="Complete CMC package for new formulation",
                    duration_months=8,
                    start_month=0,
                    deliverables=["CMC package", "Stability data"]
                ),
                TimelineMilestone(
                    phase=MilestonePhase.IND_ENABLING,
                    milestone="Bridging Studies",
                    description="BA/BE or other bridging studies",
                    duration_months=6,
                    start_month=6,
                    dependencies=["CMC Development"],
                    deliverables=["BE study report"]
                ),
                TimelineMilestone(
                    phase=MilestonePhase.NDA_PREP,
                    milestone="505(b)(2) Submission",
                    description="Submit 505(b)(2) NDA",
                    duration_months=3,
                    start_month=14,
                    dependencies=["Bridging Studies"],
                    deliverables=["eCTD submission"]
                ),
                TimelineMilestone(
                    phase=MilestonePhase.FDA_REVIEW,
                    milestone="FDA Review",
                    description="FDA review period",
                    duration_months=12,
                    start_month=17,
                    dependencies=["505(b)(2) Submission"],
                    deliverables=["FDA approval/CRL"]
                )
            ]
        
        elif pathway == RegulatoryPathway.NDA:
            timeline = [
                TimelineMilestone(
                    phase=MilestonePhase.IND_ENABLING,
                    milestone="IND-Enabling Package",
                    description="Complete IND-enabling studies",
                    duration_months=12,
                    start_month=0,
                    deliverables=["Tox package", "CMC package"]
                ),
                TimelineMilestone(
                    phase=MilestonePhase.PHASE_1,
                    milestone="Phase 1",
                    description="Safety and PK studies",
                    duration_months=12,
                    start_month=12,
                    deliverables=["Phase 1 report"]
                ),
                TimelineMilestone(
                    phase=MilestonePhase.PHASE_2,
                    milestone="Phase 2",
                    description="Proof of concept efficacy",
                    duration_months=18,
                    start_month=24,
                    deliverables=["Phase 2 report"]
                ),
                TimelineMilestone(
                    phase=MilestonePhase.PHASE_3,
                    milestone="Phase 3",
                    description="Pivotal efficacy trials",
                    duration_months=24,
                    start_month=42,
                    deliverables=["Phase 3 reports"]
                ),
                TimelineMilestone(
                    phase=MilestonePhase.NDA_PREP,
                    milestone="NDA Submission",
                    description="Submit NDA",
                    duration_months=4,
                    start_month=66,
                    deliverables=["eCTD submission"]
                ),
                TimelineMilestone(
                    phase=MilestonePhase.FDA_REVIEW,
                    milestone="FDA Review",
                    description="FDA review and approval",
                    duration_months=12,
                    start_month=70,
                    deliverables=["FDA approval"]
                )
            ]
        
        else:
            # Default timeline for other pathways
            timeline = [
                TimelineMilestone(
                    phase=MilestonePhase.PRECLINICAL,
                    milestone="Development Package",
                    description="Complete required development work",
                    duration_months=12,
                    start_month=0,
                    deliverables=["Development package"]
                ),
                TimelineMilestone(
                    phase=MilestonePhase.NDA_PREP,
                    milestone="Submission Preparation",
                    description="Prepare regulatory submission",
                    duration_months=6,
                    start_month=12,
                    deliverables=["Submission package"]
                ),
                TimelineMilestone(
                    phase=MilestonePhase.FDA_REVIEW,
                    milestone="Agency Review",
                    description="Regulatory review period",
                    duration_months=12,
                    start_month=18,
                    deliverables=["Approval/authorization"]
                )
            ]
        
        return timeline
    
    def _generate_next_actions(
        self,
        primary: PathwayRecommendation,
        checklist: List[ChecklistItem],
        evidence: Optional[EvidenceInputs]
    ) -> List[str]:
        """Generate prioritized next actions."""
        actions = []
        
        # Check for data gaps
        if not evidence or not evidence.chempath_job_id:
            actions.append("Complete chemical characterization (ChemPath analysis)")
        
        if not evidence or not evidence.toxpath_assessment_id:
            actions.append("Conduct toxicity risk assessment (ToxPath analysis)")
        
        # Pathway-specific actions
        if primary.pathway == RegulatoryPathway.IND:
            actions.append("Initiate IND-enabling toxicology studies")
            actions.append("Prepare pre-IND meeting request for FDA")
        
        elif primary.pathway == RegulatoryPathway.NDA_505B2:
            actions.append("Identify reference listed drug (RLD)")
            actions.append("Assess need for bridging studies")
            actions.append("Request Type B pre-submission meeting with FDA")
        
        # Common actions
        actions.append("Engage regulatory affairs consultant")
        actions.append("Develop detailed project timeline and budget")
        
        # Prioritize incomplete checklist items
        incomplete = [item for item in checklist 
                     if item.status == "not_started" and item.priority == "required"]
        if incomplete:
            first_item = incomplete[0]
            actions.insert(0, f"Priority: {first_item.item}")
        
        return actions[:6]  # Top 6 actions
    
    def _document_assumptions(
        self,
        profile: ProductProfile,
        constraints: Optional[StrategyConstraints]
    ) -> List[str]:
        """Document key assumptions."""
        assumptions = [
            "Product will meet quality specifications for regulatory submission",
            "Manufacturing can be scaled to support clinical and commercial supply",
            "No significant unexpected safety signals will emerge",
            "Regulatory pathway guidance remains consistent with current FDA policy"
        ]
        
        if profile.route == "oral":
            assumptions.append("Oral bioavailability is adequate for intended use")
        elif profile.route == "inhaled":
            assumptions.append("Inhalation delivery system is compatible with product")
        
        if constraints and constraints.budget_range:
            assumptions.append(f"Budget of {constraints.budget_range} is available for development")
        
        assumptions.append("Timeline estimates assume no significant delays in study conduct or regulatory review")
        
        return assumptions
    
    def get_pathways(self) -> List[Dict[str, str]]:
        """Get available regulatory pathways."""
        return [
            {
                "pathway": p.value,
                "full_name": p.full_name,
                "typical_timeline_min": p.typical_timeline_months[0],
                "typical_timeline_max": p.typical_timeline_months[1]
            }
            for p in RegulatoryPathway
        ]
    
    def get_product_types(self) -> List[Dict[str, str]]:
        """Get supported product types."""
        return [
            {"type": pt.value, "display_name": pt.display_name}
            for pt in ProductType
        ]
    
    def get_checklist_categories(self) -> List[str]:
        """Get checklist categories."""
        return [c.value for c in ChecklistCategory]
