"""
FTO (Freedom-to-Operate) Checker - PatentPath Lite Week 7
Assess freedom to operate for cannabinoid compounds.

Features:
- Active patent filtering
- Blocking patent analysis
- Risk assessment
- Licensing cost estimation
- Mitigation strategies

DISCLAIMER: This is decision support only. Not legal advice.
Consult qualified patent counsel before filing or commercial decisions.
"""
from typing import Dict, List, Optional, Any
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
import logging

from .prior_art import PriorArtSearcher, PriorArtResult, SearchQuery

logger = logging.getLogger(__name__)

# Required disclaimer for all FTO outputs
FTO_DISCLAIMER = (
    "DISCLAIMER: This is decision support only. Not legal advice. "
    "Consult qualified patent counsel before filing or commercial decisions."
)


class RiskLevel(Enum):
    """FTO risk level classification."""
    LOW = "low"
    MODERATE = "moderate"
    HIGH = "high"
    CRITICAL = "critical"


@dataclass
class BlockingPatent:
    """Analysis of potentially blocking patent."""
    patent_number: str
    title: str
    relevance_score: float
    is_blocking: bool
    risk_level: RiskLevel
    expiration_year: int
    years_until_expiration: int
    blocking_reason: str
    mitigation_options: List[str]
    url: str
    
    def to_dict(self) -> Dict:
        return {
            "patent_number": self.patent_number,
            "title": self.title,
            "relevance_score": round(self.relevance_score, 3),
            "is_blocking": self.is_blocking,
            "risk_level": self.risk_level.value,
            "expiration_year": self.expiration_year,
            "years_until_expiration": self.years_until_expiration,
            "blocking_reason": self.blocking_reason,
            "mitigation_options": self.mitigation_options,
            "url": self.url
        }


@dataclass
class LicensingCostEstimate:
    """Estimated licensing costs for blocking patents."""
    total_estimated_cost: str
    per_patent_average: str
    breakdown: List[Dict]
    confidence: str
    note: str
    
    def to_dict(self) -> Dict:
        return {
            "total_estimated_cost": self.total_estimated_cost,
            "per_patent_average": self.per_patent_average,
            "breakdown": self.breakdown,
            "confidence": self.confidence,
            "note": self.note
        }


@dataclass
class FTOReport:
    """Complete FTO analysis report."""
    fto_risk_level: RiskLevel
    risk_score: float
    target_market: str
    intended_use: Optional[str]
    blocking_patents: List[BlockingPatent]
    num_blocking_patents: int
    num_patents_analyzed: int
    estimated_licensing_costs: LicensingCostEstimate
    recommendation: str
    mitigation_strategies: List[str]
    next_steps: List[str]
    analysis_date: str
    disclaimer: str = FTO_DISCLAIMER
    
    def to_dict(self) -> Dict:
        return {
            "fto_risk_level": self.fto_risk_level.value,
            "risk_score": round(self.risk_score, 3),
            "target_market": self.target_market,
            "intended_use": self.intended_use,
            "blocking_patents": [bp.to_dict() for bp in self.blocking_patents],
            "num_blocking_patents": self.num_blocking_patents,
            "num_patents_analyzed": self.num_patents_analyzed,
            "estimated_licensing_costs": self.estimated_licensing_costs.to_dict(),
            "recommendation": self.recommendation,
            "mitigation_strategies": self.mitigation_strategies,
            "next_steps": self.next_steps,
            "analysis_date": self.analysis_date,
            "disclaimer": self.disclaimer
        }


class FTOChecker:
    """Freedom-to-Operate analysis for cannabinoid compounds.
    
    MVP version checks:
    - Active patents in US only
    - Basic blocking patent identification
    - Estimated licensing costs
    
    Full version (PatentPath Full) includes:
    - Multi-jurisdiction analysis
    - Detailed infringement probability
    - Design-around strategies
    """
    
    # Industry average licensing costs
    LICENSING_COSTS = {
        RiskLevel.CRITICAL: 2000000,    # $2M for critical blocking patent
        RiskLevel.HIGH: 1000000,        # $1M for high-risk
        RiskLevel.MODERATE: 250000,     # $250K for moderate
        RiskLevel.LOW: 50000            # $50K for low-risk
    }
    
    def __init__(self, searcher: Optional[PriorArtSearcher] = None):
        """Initialize FTO checker.
        
        Args:
            searcher: PriorArtSearcher instance (creates new if None)
        """
        self.searcher = searcher or PriorArtSearcher()
    
    async def check_fto(
        self,
        compound_description: str,
        keywords: List[str],
        target_market: str = "US",
        intended_use: Optional[str] = None,
        smiles: Optional[str] = None
    ) -> FTOReport:
        """Perform FTO analysis for a compound.
        
        Args:
            compound_description: Technical description of compound
            keywords: Search keywords
            target_market: Geographic market (US only for MVP)
            intended_use: Therapeutic indication
            smiles: Optional SMILES string
            
        Returns:
            FTOReport with risk assessment and recommendations
        """
        # Search for potentially blocking patents
        search_results = await self.searcher.search(
            keywords=keywords,
            max_results=50
        )
        
        # Filter for active patents only
        active_patents = self._filter_active_patents(search_results)
        
        # Analyze each patent for blocking potential
        blocking_patents = []
        for patent in active_patents:
            blocking_analysis = self._analyze_blocking_potential(
                patent,
                compound_description,
                intended_use
            )
            
            if blocking_analysis.is_blocking:
                blocking_patents.append(blocking_analysis)
        
        # Calculate overall FTO risk
        risk_level, risk_score = self._calculate_fto_risk(blocking_patents)
        
        # Estimate licensing costs
        licensing_costs = self._estimate_licensing_costs(blocking_patents)
        
        # Generate recommendations
        recommendation = self._generate_fto_recommendation(risk_level, blocking_patents)
        mitigation = self._generate_mitigation_strategies(blocking_patents)
        next_steps = self._generate_fto_next_steps(risk_level)
        
        return FTOReport(
            fto_risk_level=risk_level,
            risk_score=risk_score,
            target_market=target_market,
            intended_use=intended_use,
            blocking_patents=blocking_patents,
            num_blocking_patents=len(blocking_patents),
            num_patents_analyzed=len(search_results),
            estimated_licensing_costs=licensing_costs,
            recommendation=recommendation,
            mitigation_strategies=mitigation,
            next_steps=next_steps,
            analysis_date=datetime.now().isoformat()
        )
    
    def _filter_active_patents(
        self,
        patents: List[PriorArtResult]
    ) -> List[PriorArtResult]:
        """Filter for active (not expired/abandoned) patents."""
        active = []
        current_year = datetime.now().year
        
        for patent in patents:
            # Extract year from publication date
            patent_year = self._extract_year(patent.publication_date)
            
            if patent_year:
                # Patents expire 20 years from filing
                # Use publication date as approximation (filing is typically 1-2 years before)
                estimated_expiration = patent_year + 18  # Conservative estimate
                
                if estimated_expiration > current_year:
                    # Still active
                    patent.filing_date = str(patent_year)
                    active.append(patent)
        
        return active
    
    def _extract_year(self, date_str: str) -> Optional[int]:
        """Extract year from date string."""
        if not date_str:
            return None
        try:
            return int(date_str[:4])
        except (ValueError, IndexError):
            return None
    
    def _analyze_blocking_potential(
        self,
        patent: PriorArtResult,
        compound_description: str,
        intended_use: Optional[str]
    ) -> BlockingPatent:
        """Analyze if a patent could block commercialization."""
        relevance = patent.relevance_score
        current_year = datetime.now().year
        
        # Estimate expiration
        patent_year = self._extract_year(patent.publication_date) or current_year
        estimated_expiration = patent_year + 20
        years_until_expiration = max(0, estimated_expiration - current_year)
        
        # Determine if blocking and risk level
        if relevance > 0.8:
            is_blocking = True
            risk_level = RiskLevel.CRITICAL if years_until_expiration > 10 else RiskLevel.HIGH
        elif relevance > 0.6:
            is_blocking = True
            risk_level = RiskLevel.HIGH if years_until_expiration > 10 else RiskLevel.MODERATE
        elif relevance > 0.4:
            is_blocking = relevance > 0.5
            risk_level = RiskLevel.MODERATE
        else:
            is_blocking = False
            risk_level = RiskLevel.LOW
        
        return BlockingPatent(
            patent_number=patent.patent_number,
            title=patent.title,
            relevance_score=relevance,
            is_blocking=is_blocking,
            risk_level=risk_level,
            expiration_year=estimated_expiration,
            years_until_expiration=years_until_expiration,
            blocking_reason=self._generate_blocking_reason(relevance, intended_use),
            mitigation_options=self._generate_mitigation_options(risk_level),
            url=patent.url
        )
    
    def _generate_blocking_reason(
        self,
        relevance: float,
        intended_use: Optional[str]
    ) -> str:
        """Generate explanation for why patent may be blocking."""
        use_text = intended_use or "therapeutic use"
        
        if relevance > 0.8:
            return f"Highly similar compound structure and overlapping {use_text}"
        elif relevance > 0.6:
            return f"Similar compound class with potential overlap in {use_text}"
        elif relevance > 0.4:
            return f"Related technology area, some overlap in {use_text}"
        else:
            return "Low blocking risk - significant structural or use differences"
    
    def _generate_mitigation_options(self, risk_level: RiskLevel) -> List[str]:
        """Generate mitigation strategies for blocking patents."""
        if risk_level == RiskLevel.CRITICAL:
            return [
                "Immediate patent attorney consultation required",
                "Negotiate licensing agreement",
                "Challenge patent validity (IPR/PGR)",
                "Major structural redesign required",
                "Consider alternative therapeutic targets"
            ]
        elif risk_level == RiskLevel.HIGH:
            return [
                "Engage patent attorney for detailed claim analysis",
                "Explore licensing negotiations",
                "Challenge patent validity (prior art search)",
                "Design-around with structural modifications",
                "Wait for patent expiration if timeline permits"
            ]
        elif risk_level == RiskLevel.MODERATE:
            return [
                "Conduct detailed claim analysis with attorney",
                "Explore design-around opportunities",
                "Consider licensing if commercialization imminent",
                "Monitor patent prosecution history"
            ]
        else:
            return [
                "Monitor patent status",
                "Proceed with caution",
                "Standard patent due diligence"
            ]
    
    def _calculate_fto_risk(
        self,
        blocking_patents: List[BlockingPatent]
    ) -> tuple[RiskLevel, float]:
        """Calculate overall FTO risk level and score."""
        if not blocking_patents:
            return RiskLevel.LOW, 0.1
        
        # Count patents by risk level
        critical_count = sum(1 for p in blocking_patents if p.risk_level == RiskLevel.CRITICAL)
        high_count = sum(1 for p in blocking_patents if p.risk_level == RiskLevel.HIGH)
        moderate_count = sum(1 for p in blocking_patents if p.risk_level == RiskLevel.MODERATE)
        
        # Risk score (0-1)
        risk_score = min(1.0, (
            critical_count * 0.4 +
            high_count * 0.25 +
            moderate_count * 0.1
        ))
        
        # Risk level classification
        if critical_count > 0 or risk_score > 0.7:
            level = RiskLevel.CRITICAL
        elif high_count > 0 or risk_score > 0.5:
            level = RiskLevel.HIGH
        elif risk_score > 0.2:
            level = RiskLevel.MODERATE
        else:
            level = RiskLevel.LOW
        
        return level, risk_score
    
    def _estimate_licensing_costs(
        self,
        blocking_patents: List[BlockingPatent]
    ) -> LicensingCostEstimate:
        """Estimate licensing costs for blocking patents."""
        if not blocking_patents:
            return LicensingCostEstimate(
                total_estimated_cost="$0",
                per_patent_average="$0",
                breakdown=[],
                confidence="High",
                note="No blocking patents identified"
            )
        
        breakdown = []
        total_cost = 0
        
        for patent in blocking_patents:
            cost = self.LICENSING_COSTS.get(patent.risk_level, 0)
            total_cost += cost
            
            breakdown.append({
                "patent_number": patent.patent_number,
                "risk_level": patent.risk_level.value,
                "estimated_cost": f"${cost:,}",
                "expiration_year": patent.expiration_year
            })
        
        avg_cost = total_cost // len(blocking_patents)
        
        return LicensingCostEstimate(
            total_estimated_cost=f"${total_cost:,}",
            per_patent_average=f"${avg_cost:,}",
            breakdown=breakdown,
            confidence="Low to Moderate",
            note="Estimates based on industry averages. Actual costs vary significantly based on negotiation, patent strength, and market factors."
        )
    
    def _generate_fto_recommendation(
        self,
        risk_level: RiskLevel,
        blocking_patents: List[BlockingPatent]
    ) -> str:
        """Generate FTO recommendation."""
        if risk_level == RiskLevel.CRITICAL:
            return (
                "CRITICAL RISK: Commercialization highly likely blocked without licensing or major redesign. "
                "Engage patent attorney immediately before any further development investment."
            )
        elif risk_level == RiskLevel.HIGH:
            return (
                "HIGH RISK: Significant blocking patents identified. "
                "Attorney review required before proceeding. Consider licensing or design-around strategies."
            )
        elif risk_level == RiskLevel.MODERATE:
            return (
                "MODERATE RISK: Some potentially blocking patents identified. "
                "Attorney consultation recommended before significant commercial investment."
            )
        else:
            return (
                "LOW RISK: No major blocking patents identified. "
                "Proceed with standard patent due diligence and periodic monitoring."
            )
    
    def _generate_mitigation_strategies(
        self,
        blocking_patents: List[BlockingPatent]
    ) -> List[str]:
        """Generate overall mitigation strategies."""
        strategies = []
        
        if not blocking_patents:
            return ["Continue development with standard IP monitoring"]
        
        # Collect unique strategies
        strategy_set = set()
        for patent in blocking_patents:
            for option in patent.mitigation_options[:2]:  # Top 2 from each
                strategy_set.add(option)
        
        strategies = list(strategy_set)[:6]  # Limit to 6
        
        # Add general strategies
        if any(p.years_until_expiration < 5 for p in blocking_patents):
            strategies.append("Consider delayed market entry to avoid licensing costs on soon-to-expire patents")
        
        return strategies
    
    def _generate_fto_next_steps(self, risk_level: RiskLevel) -> List[str]:
        """Generate actionable next steps."""
        if risk_level == RiskLevel.CRITICAL:
            return [
                "STOP: Do not proceed with further development investment",
                "Engage patent attorney for formal FTO opinion",
                "Conduct detailed claim-by-claim infringement analysis",
                "Evaluate licensing feasibility with patent holders",
                "Assess design-around possibilities with R&D team"
            ]
        elif risk_level == RiskLevel.HIGH:
            return [
                "Engage patent attorney for claim mapping",
                "Obtain formal freedom-to-operate opinion letter",
                "Initiate licensing discussions with key patent holders",
                "Investigate design-around opportunities",
                "File provisional patent to establish priority date"
            ]
        elif risk_level == RiskLevel.MODERATE:
            return [
                "Consult patent attorney for preliminary assessment",
                "File provisional patent application to secure priority",
                "Monitor patent prosecution history for changes",
                "Document design choices to support non-infringement position"
            ]
        else:
            return [
                "Proceed with compound development",
                "File provisional patent application",
                "Set up periodic patent landscape monitoring",
                "Document development timeline for prior art defense"
            ]


async def quick_fto_check(
    keywords: List[str],
    intended_use: Optional[str] = None
) -> Dict:
    """Quick FTO assessment (convenience function).
    
    Args:
        keywords: Search keywords for compound
        intended_use: Therapeutic indication
        
    Returns:
        Simplified FTO assessment
    """
    checker = FTOChecker()
    report = await checker.check_fto(
        compound_description=" ".join(keywords),
        keywords=keywords,
        intended_use=intended_use
    )
    
    return {
        "risk_level": report.fto_risk_level.value,
        "risk_score": round(report.risk_score, 2),
        "blocking_patents": report.num_blocking_patents,
        "patents_analyzed": report.num_patents_analyzed,
        "estimated_licensing_cost": report.estimated_licensing_costs.total_estimated_cost,
        "recommendation": report.recommendation,
        "top_next_step": report.next_steps[0] if report.next_steps else None,
        "disclaimer": FTO_DISCLAIMER
    }
