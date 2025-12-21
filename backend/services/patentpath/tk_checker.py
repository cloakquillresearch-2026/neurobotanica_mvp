"""
Traditional Knowledge Attribution Checker - PatentPath Lite Feature 5

This module checks if compounds require traditional knowledge attribution
under Nagoya Protocol and national ABS (Access and Benefit Sharing) laws.

MVP Implementation:
- Keyword-based detection for TK derivation indicators
- Sacred knowledge detection (absolute bar to patenting)
- Attribution requirements generation
- EquiPath integration placeholders for benefit-sharing

Full PatentPath Version (Future):
- ML classifier with 91.8% accuracy
- TKDL (Traditional Knowledge Digital Library) integration
- Real-time ABS law compliance checking
- Automated community consent verification
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional
from enum import Enum
from datetime import datetime
import logging

logger = logging.getLogger(__name__)


class TKIndicatorStrength(Enum):
    """Strength of traditional knowledge indicators."""
    NONE = "none"
    WEAK = "weak"
    MODERATE = "moderate"
    STRONG = "strong"
    DEFINITIVE = "definitive"


class SacredKnowledgeStatus(Enum):
    """Sacred knowledge detection status."""
    NOT_DETECTED = "not_detected"
    POSSIBLE = "possible"
    CONFIRMED = "confirmed"


@dataclass
class AttributionRequirement:
    """Single attribution requirement."""
    requirement: str
    category: str  # documentation, consent, benefit_sharing, legal
    priority: str  # required, recommended, optional
    deadline: str  # before_filing, during_prosecution, before_commercialization
    description: str


@dataclass
class TKCheckResult:
    """Result of traditional knowledge attribution check."""
    tk_derived: bool
    sacred_knowledge_detected: bool
    indicator_strength: TKIndicatorStrength
    attribution_required: bool
    can_proceed: bool
    blocking_reason: Optional[str] = None
    detected_indicators: List[str] = field(default_factory=list)
    requirements: List[AttributionRequirement] = field(default_factory=list)
    next_steps: List[str] = field(default_factory=list)
    equipath_integration: Optional[Dict] = None
    specification_requirements: Optional[Dict] = None
    legal_compliance: Optional[str] = None
    recommendation: str = ""
    checked_at: str = field(default_factory=lambda: datetime.utcnow().isoformat())
    
    def to_dict(self) -> Dict:
        """Convert result to dictionary."""
        return {
            "tk_derived": self.tk_derived,
            "sacred_knowledge_detected": self.sacred_knowledge_detected,
            "indicator_strength": self.indicator_strength.value,
            "attribution_required": self.attribution_required,
            "can_proceed": self.can_proceed,
            "blocking_reason": self.blocking_reason,
            "detected_indicators": self.detected_indicators,
            "requirements": [
                {
                    "requirement": r.requirement,
                    "category": r.category,
                    "priority": r.priority,
                    "deadline": r.deadline,
                    "description": r.description
                }
                for r in self.requirements
            ],
            "next_steps": self.next_steps,
            "equipath_integration": self.equipath_integration,
            "specification_requirements": self.specification_requirements,
            "legal_compliance": self.legal_compliance,
            "recommendation": self.recommendation,
            "checked_at": self.checked_at
        }


class TKAttributionChecker:
    """Check if compound requires traditional knowledge attribution.
    
    MVP version uses keyword detection.
    Full version (PatentPath Full) uses ML classifier (91.8% accuracy) + TKDL integration.
    
    Key concepts:
    - Traditional Knowledge (TK): Knowledge, innovations and practices of indigenous 
      and local communities relating to biodiversity and genetic resources
    - Nagoya Protocol: International agreement governing access and benefit-sharing
    - ABS Laws: National laws implementing access and benefit-sharing requirements
    - Sacred Knowledge: Traditional knowledge with spiritual/religious significance
      that cannot be patented under any circumstances
    """
    
    def __init__(self):
        # Keywords indicating TK derivation (weighted by strength)
        self.tk_keywords = {
            # Strong indicators (high confidence TK-derived)
            "strong": [
                "traditional knowledge",
                "indigenous knowledge",
                "ethnobotanical",
                "traditional medicine",
                "folk medicine",
                "traditional preparation",
                "traditional use",
                "ancestral knowledge",
                "indigenous peoples",
                "local community knowledge"
            ],
            # Moderate indicators (likely TK-derived)
            "moderate": [
                "ayurvedic",
                "traditional chinese medicine",
                "tcm",
                "kampo",
                "unani",
                "siddha",
                "aboriginal",
                "tribal medicine",
                "native healing",
                "ethnopharmacology"
            ],
            # Weak indicators (possible TK connection)
            "weak": [
                "traditional",
                "indigenous",
                "native",
                "ancient",
                "historical preparation",
                "historical use",
                "folk remedy"
            ]
        }
        
        # Sacred knowledge indicators (absolute bar to patenting)
        self.sacred_keywords = {
            # Definitive sacred indicators
            "definitive": [
                "sacred knowledge",
                "sacred medicine",
                "ceremonial use only",
                "sacred ritual",
                "sacred plant"
            ],
            # Strong sacred indicators
            "strong": [
                "sacred",
                "ceremonial",
                "spiritual healing",
                "religious ceremony",
                "ritual use",
                "holy medicine"
            ],
            # Possible sacred indicators (need context)
            "possible": [
                "spiritual",
                "religious",
                "ritual",
                "blessed",
                "consecrated"
            ]
        }
        
        # TK derivation source types
        self.tk_source_types = [
            "ethnobotanical_literature",
            "traditional_medicine_text",
            "indigenous_community",
            "folk_medicine_database",
            "ethnopharmacological_survey",
            "traditional_healer_interview"
        ]
        
        # Regions with strong ABS requirements
        self.high_abs_regions = [
            "india",
            "brazil",
            "south africa",
            "peru",
            "kenya",
            "philippines",
            "costa rica",
            "china"
        ]
    
    def check_tk_attribution(
        self,
        compound_description: str,
        derivation_source: Optional[str] = None,
        source_description: Optional[str] = None,
        source_region: Optional[str] = None,
        community_name: Optional[str] = None
    ) -> TKCheckResult:
        """Check if compound requires TK attribution.
        
        Args:
            compound_description: Description of compound and its development
            derivation_source: Source type (e.g., "ethnobotanical_literature")
            source_description: Detailed description of knowledge source
            source_region: Geographic region of knowledge source
            community_name: Name of indigenous/local community (if known)
            
        Returns:
            TKCheckResult with attribution requirements and next steps
        """
        logger.info(f"Checking TK attribution for compound")
        
        # Combine all text for analysis
        all_text = self._combine_text(
            compound_description,
            derivation_source,
            source_description,
            community_name
        )
        
        # Check for sacred knowledge first (blocking condition)
        sacred_status, sacred_indicators = self._detect_sacred_knowledge(all_text)
        
        if sacred_status == SacredKnowledgeStatus.CONFIRMED:
            return TKCheckResult(
                tk_derived=True,
                sacred_knowledge_detected=True,
                indicator_strength=TKIndicatorStrength.DEFINITIVE,
                attribution_required=False,
                can_proceed=False,
                blocking_reason="Sacred traditional knowledge detected - cannot be patented under any circumstances",
                detected_indicators=sacred_indicators,
                recommendation="Do not proceed with patent filing. Consider alternative protection strategies or community partnership without IP claims.",
                next_steps=[
                    "Consult with community leaders about appropriate use",
                    "Consider non-patent protection alternatives",
                    "Respect community wishes regarding sacred knowledge"
                ]
            )
        
        # Check for TK indicators
        tk_detected, tk_strength, tk_indicators = self._detect_tk_indicators(
            all_text,
            derivation_source
        )
        
        # Check source type
        is_tk_source = derivation_source in self.tk_source_types if derivation_source else False
        if is_tk_source and not tk_detected:
            tk_detected = True
            tk_strength = TKIndicatorStrength.MODERATE
            tk_indicators.append(f"Derivation source type: {derivation_source}")
        
        # Check high ABS region
        in_high_abs_region = False
        if source_region:
            in_high_abs_region = any(
                region in source_region.lower() 
                for region in self.high_abs_regions
            )
            if in_high_abs_region:
                tk_indicators.append(f"Source region with strong ABS laws: {source_region}")
        
        if tk_detected:
            # TK-derived but not sacred - attribution required
            requirements = self._generate_attribution_requirements(
                community_name,
                in_high_abs_region
            )
            
            return TKCheckResult(
                tk_derived=True,
                sacred_knowledge_detected=sacred_status == SacredKnowledgeStatus.POSSIBLE,
                indicator_strength=tk_strength,
                attribution_required=True,
                can_proceed=True,
                detected_indicators=tk_indicators,
                requirements=requirements,
                next_steps=self._generate_tk_next_steps(community_name, in_high_abs_region),
                equipath_integration=self._generate_equipath_config(community_name),
                specification_requirements=self._generate_specification_requirements(community_name),
                legal_compliance="Required under Nagoya Protocol and potential national ABS laws",
                recommendation="TK attribution required. Proceed with patent filing after establishing consent and benefit-sharing."
            )
        else:
            # Not TK-derived
            return TKCheckResult(
                tk_derived=False,
                sacred_knowledge_detected=False,
                indicator_strength=TKIndicatorStrength.NONE,
                attribution_required=False,
                can_proceed=True,
                recommendation="No TK attribution required. Proceed with standard patent filing."
            )
    
    def _combine_text(
        self,
        compound_description: str,
        derivation_source: Optional[str],
        source_description: Optional[str],
        community_name: Optional[str]
    ) -> str:
        """Combine all text sources for analysis."""
        parts = [
            compound_description or "",
            derivation_source or "",
            source_description or "",
            community_name or ""
        ]
        return " ".join(filter(None, parts)).lower()
    
    def _detect_sacred_knowledge(
        self,
        text: str
    ) -> tuple[SacredKnowledgeStatus, List[str]]:
        """Detect if knowledge is sacred (absolute bar to patenting)."""
        detected = []
        
        # Check definitive sacred indicators
        for keyword in self.sacred_keywords["definitive"]:
            if keyword in text:
                detected.append(f"Sacred indicator (definitive): {keyword}")
        
        if detected:
            return SacredKnowledgeStatus.CONFIRMED, detected
        
        # Check strong sacred indicators
        for keyword in self.sacred_keywords["strong"]:
            if keyword in text:
                detected.append(f"Sacred indicator (strong): {keyword}")
        
        if len(detected) >= 2:
            return SacredKnowledgeStatus.CONFIRMED, detected
        elif detected:
            return SacredKnowledgeStatus.POSSIBLE, detected
        
        # Check possible sacred indicators (need multiple)
        possible = []
        for keyword in self.sacred_keywords["possible"]:
            if keyword in text:
                possible.append(f"Sacred indicator (possible): {keyword}")
        
        if len(possible) >= 3:
            return SacredKnowledgeStatus.POSSIBLE, possible
        
        return SacredKnowledgeStatus.NOT_DETECTED, []
    
    def _detect_tk_indicators(
        self,
        text: str,
        derivation_source: Optional[str]
    ) -> tuple[bool, TKIndicatorStrength, List[str]]:
        """Detect if compound is derived from traditional knowledge."""
        detected = []
        
        # Check strong indicators
        for keyword in self.tk_keywords["strong"]:
            if keyword in text:
                detected.append(f"TK indicator (strong): {keyword}")
        
        if detected:
            return True, TKIndicatorStrength.STRONG, detected
        
        # Check moderate indicators
        for keyword in self.tk_keywords["moderate"]:
            if keyword in text:
                detected.append(f"TK indicator (moderate): {keyword}")
        
        if detected:
            return True, TKIndicatorStrength.MODERATE, detected
        
        # Check weak indicators (need multiple)
        weak = []
        for keyword in self.tk_keywords["weak"]:
            if keyword in text:
                weak.append(f"TK indicator (weak): {keyword}")
        
        if len(weak) >= 2:
            return True, TKIndicatorStrength.WEAK, weak
        
        return False, TKIndicatorStrength.NONE, []
    
    def _generate_attribution_requirements(
        self,
        community_name: Optional[str],
        in_high_abs_region: bool
    ) -> List[AttributionRequirement]:
        """Generate requirements for TK attribution."""
        requirements = [
            AttributionRequirement(
                requirement="Document traditional knowledge source",
                category="documentation",
                priority="required",
                deadline="before_filing",
                description="Record community name, geographic location, and nature of traditional knowledge"
            ),
            AttributionRequirement(
                requirement="Establish informed consent",
                category="consent",
                priority="required",
                deadline="before_filing",
                description="Obtain Prior Informed Consent (PIC) from knowledge holders via EthnoPath system"
            ),
            AttributionRequirement(
                requirement="Configure benefit-sharing agreement",
                category="benefit_sharing",
                priority="required",
                deadline="before_filing",
                description="Set up EquiPath benefit-sharing (recommended: 70% community, 25% STEM education, 5% infrastructure)"
            ),
            AttributionRequirement(
                requirement="Include attribution statement in specification",
                category="documentation",
                priority="required",
                deadline="before_filing",
                description="Draft attribution language acknowledging TK contribution in patent specification"
            ),
            AttributionRequirement(
                requirement="Cite traditional knowledge sources",
                category="documentation",
                priority="required",
                deadline="before_filing",
                description="Properly cite ethnobotanical literature and community sources in patent application"
            ),
            AttributionRequirement(
                requirement="Establish community partnership agreement",
                category="legal",
                priority="required",
                deadline="before_filing",
                description="Execute formal partnership agreement with community representatives"
            )
        ]
        
        if in_high_abs_region:
            requirements.append(
                AttributionRequirement(
                    requirement="National ABS permit",
                    category="legal",
                    priority="required",
                    deadline="before_filing",
                    description="Obtain national Access and Benefit Sharing permit from relevant authority"
                )
            )
        
        if community_name:
            requirements.append(
                AttributionRequirement(
                    requirement=f"Verify consent from {community_name}",
                    category="consent",
                    priority="required",
                    deadline="before_filing",
                    description=f"Confirm Prior Informed Consent specifically from {community_name} representatives"
                )
            )
        
        return requirements
    
    def _generate_tk_next_steps(
        self,
        community_name: Optional[str],
        in_high_abs_region: bool
    ) -> List[str]:
        """Generate actionable next steps for TK attribution."""
        steps = [
            "Contact community leaders to establish consent",
            "Register community benefit-sharing wallet in EquiPath",
            "Draft attribution language for patent specification",
            "Consult attorney experienced in TK and Nagoya Protocol compliance",
            "Delay patent filing until consent and benefit-sharing established"
        ]
        
        if in_high_abs_region:
            steps.insert(0, "Apply for national ABS permit in source country")
        
        if community_name:
            steps.insert(0, f"Schedule consultation with {community_name} leadership")
        
        return steps
    
    def _generate_equipath_config(
        self,
        community_name: Optional[str]
    ) -> Dict:
        """Generate EquiPath integration configuration."""
        return {
            "community_wallet_setup": "Required before filing patent",
            "community_name": community_name or "To be specified",
            "estimated_revenue_share": {
                "community": "70%",
                "stem_education": "25%",
                "infrastructure": "5%"
            },
            "attribution_blockchain_id": "pending",
            "consent_verification_status": "pending",
            "benefit_tracking_enabled": True
        }
    
    def _generate_specification_requirements(
        self,
        community_name: Optional[str]
    ) -> Dict:
        """Generate patent specification requirements for TK attribution."""
        community = community_name or "[COMMUNITY NAME]"
        
        return {
            "attribution_statement_template": (
                f"The compound disclosed herein was developed based on traditional knowledge "
                f"from {community}, who have historically used [PLANT/PREPARATION METHOD] "
                f"for [TRADITIONAL USE]. The inventors acknowledge the contribution of this "
                f"traditional knowledge and have established a benefit-sharing agreement whereby "
                f"[XX]% of commercialization revenues will be shared with the {community} "
                f"pursuant to the Nagoya Protocol."
            ),
            "citation_requirements": [
                "Include community name and geographic location",
                "Cite ethnobotanical literature documenting traditional use",
                "Reference benefit-sharing agreement",
                "Include EquiPath blockchain manifest ID"
            ],
            "background_section": "Dedicate Background section subsection to traditional knowledge source",
            "inventor_declaration": f"Inventors must declare TK contribution from {community}"
        }
    
    def get_tk_source_types(self) -> List[str]:
        """Get list of recognized TK source types."""
        return self.tk_source_types.copy()
    
    def get_high_abs_regions(self) -> List[str]:
        """Get list of regions with strong ABS requirements."""
        return self.high_abs_regions.copy()
    
    def validate_attribution_completeness(
        self,
        community_name: str,
        consent_obtained: bool,
        benefit_sharing_configured: bool,
        attribution_drafted: bool,
        abs_permit_obtained: bool = False,
        source_region: Optional[str] = None
    ) -> Dict:
        """Validate that all attribution requirements are met.
        
        Returns:
            Validation result with missing requirements
        """
        missing = []
        
        if not community_name:
            missing.append("Community name not specified")
        
        if not consent_obtained:
            missing.append("Prior Informed Consent not obtained")
        
        if not benefit_sharing_configured:
            missing.append("Benefit-sharing agreement not configured")
        
        if not attribution_drafted:
            missing.append("Attribution statement not drafted")
        
        # Check if ABS permit required
        needs_abs_permit = False
        if source_region:
            needs_abs_permit = any(
                region in source_region.lower()
                for region in self.high_abs_regions
            )
        
        if needs_abs_permit and not abs_permit_obtained:
            missing.append(f"National ABS permit required for {source_region}")
        
        return {
            "complete": len(missing) == 0,
            "missing_requirements": missing,
            "can_file_patent": len(missing) == 0,
            "recommendation": (
                "All TK attribution requirements met. Proceed with filing."
                if len(missing) == 0
                else f"Complete {len(missing)} requirements before filing patent."
            )
        }
