"""
Patent Cost Estimator - PatentPath Lite Feature 6

This module estimates patent filing, prosecution, and maintenance costs
for cannabinoid compound IP protection strategies.

MVP Implementation:
- USPTO fee schedule (2025 rates)
- Attorney cost estimates (industry averages)
- Patent vs. Trade Secret comparison
- 20-year cost projections

Full PatentPath Version (Future):
- Multi-jurisdiction cost analysis (PCT, EPO, CNIPA, etc.)
- Portfolio optimization algorithms
- Real-time fee schedule updates
- Attorney network rate negotiation
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional
from enum import Enum
from datetime import datetime
import logging

logger = logging.getLogger(__name__)


class ProtectionStrategy(Enum):
    """IP protection strategy options."""
    PATENT = "patent"
    TRADE_SECRET = "trade_secret"
    HYBRID = "hybrid"  # Patent core innovations, trade secret for processes


class ComplexityLevel(Enum):
    """Patent application complexity levels."""
    SIMPLE = "simple"
    MODERATE = "moderate"
    COMPLEX = "complex"


class EntitySize(Enum):
    """USPTO entity size for fee calculation."""
    LARGE = "large"
    SMALL = "small"  # 50% discount
    MICRO = "micro"  # 75% discount


@dataclass
class USPTOFeeSchedule:
    """2025 USPTO fee schedule."""
    # Filing fees (large entity)
    filing_fee: int = 1280
    search_fee: int = 2640
    examination_fee: int = 3040
    issue_fee: int = 1000
    
    # Claim fees
    per_claim_over_20: int = 100
    per_independent_claim_over_3: int = 460
    
    # Maintenance fees
    maintenance_3_5_years: int = 2000
    maintenance_7_5_years: int = 3760
    maintenance_11_5_years: int = 7060
    
    # Other fees
    petition_fee: int = 200
    rce_fee: int = 1360  # Request for Continued Examination
    appeal_fee: int = 840
    
    def get_maintenance_total(self) -> int:
        """Get total maintenance fees over patent life."""
        return (
            self.maintenance_3_5_years +
            self.maintenance_7_5_years +
            self.maintenance_11_5_years
        )
    
    def apply_entity_discount(self, entity_size: EntitySize) -> 'USPTOFeeSchedule':
        """Apply entity size discount to fees."""
        if entity_size == EntitySize.LARGE:
            return self
        
        multiplier = 0.5 if entity_size == EntitySize.SMALL else 0.25
        
        return USPTOFeeSchedule(
            filing_fee=int(self.filing_fee * multiplier),
            search_fee=int(self.search_fee * multiplier),
            examination_fee=int(self.examination_fee * multiplier),
            issue_fee=int(self.issue_fee * multiplier),
            per_claim_over_20=int(self.per_claim_over_20 * multiplier),
            per_independent_claim_over_3=int(self.per_independent_claim_over_3 * multiplier),
            maintenance_3_5_years=int(self.maintenance_3_5_years * multiplier),
            maintenance_7_5_years=int(self.maintenance_7_5_years * multiplier),
            maintenance_11_5_years=int(self.maintenance_11_5_years * multiplier),
            petition_fee=int(self.petition_fee * multiplier),
            rce_fee=int(self.rce_fee * multiplier),
            appeal_fee=int(self.appeal_fee * multiplier)
        )


@dataclass
class AttorneyCostEstimates:
    """Attorney cost estimates by complexity."""
    # Drafting costs
    simple_drafting: int = 5000
    moderate_drafting: int = 10000
    complex_drafting: int = 15000
    
    # Prosecution costs
    simple_prosecution: int = 3000
    moderate_prosecution: int = 7500
    complex_prosecution: int = 12000
    
    # Per-action costs
    office_action_response: int = 2500
    rce_response: int = 3500
    appeal_brief: int = 5000
    oral_argument: int = 3000
    
    def get_drafting_cost(self, complexity: ComplexityLevel) -> int:
        """Get drafting cost by complexity."""
        costs = {
            ComplexityLevel.SIMPLE: self.simple_drafting,
            ComplexityLevel.MODERATE: self.moderate_drafting,
            ComplexityLevel.COMPLEX: self.complex_drafting
        }
        return costs[complexity]
    
    def get_prosecution_cost(self, complexity: ComplexityLevel) -> int:
        """Get prosecution cost by complexity."""
        costs = {
            ComplexityLevel.SIMPLE: self.simple_prosecution,
            ComplexityLevel.MODERATE: self.moderate_prosecution,
            ComplexityLevel.COMPLEX: self.complex_prosecution
        }
        return costs[complexity]


@dataclass
class CostBreakdown:
    """Detailed cost breakdown."""
    filing_costs: int
    prosecution_costs: int
    maintenance_costs: int
    total_cost: int
    
    def to_dict(self) -> Dict:
        return {
            "filing_costs": f"${self.filing_costs:,}",
            "prosecution_costs": f"${self.prosecution_costs:,}",
            "maintenance_fees_20yr": f"${self.maintenance_costs:,}",
            "total_investment": f"${self.total_cost:,}"
        }


@dataclass
class CostScheduleEvent:
    """Single event in cost schedule."""
    year: float
    event: str
    cost: int
    cumulative: int
    
    def to_dict(self) -> Dict:
        return {
            "year": self.year,
            "event": self.event,
            "cost": f"${self.cost:,}",
            "cumulative": f"${self.cumulative:,}"
        }


@dataclass
class CostEstimate:
    """Complete cost estimate result."""
    strategy: ProtectionStrategy
    num_compounds: int
    jurisdictions: List[str]
    entity_size: EntitySize
    complexity: ComplexityLevel
    
    total_cost_estimate: CostBreakdown
    per_patent_breakdown: CostBreakdown
    
    timeline: Dict[str, str]
    cost_schedule: List[CostScheduleEvent]
    
    advantages: List[str] = field(default_factory=list)
    disadvantages: List[str] = field(default_factory=list)
    recommendation: str = ""
    
    comparison_savings: Optional[int] = None  # vs. alternative strategy
    estimated_at: str = field(default_factory=lambda: datetime.utcnow().isoformat())
    
    def to_dict(self) -> Dict:
        return {
            "strategy": self.strategy.value,
            "num_compounds": self.num_compounds,
            "jurisdictions": self.jurisdictions,
            "entity_size": self.entity_size.value,
            "complexity": self.complexity.value,
            "total_cost_estimate": self.total_cost_estimate.to_dict(),
            "per_patent_breakdown": self.per_patent_breakdown.to_dict(),
            "timeline": self.timeline,
            "cost_schedule": [e.to_dict() for e in self.cost_schedule],
            "advantages": self.advantages,
            "disadvantages": self.disadvantages,
            "recommendation": self.recommendation,
            "comparison_savings": f"${self.comparison_savings:,}" if self.comparison_savings else None,
            "estimated_at": self.estimated_at
        }


class CostEstimator:
    """Estimate patent filing and maintenance costs.
    
    MVP version uses USPTO fee schedule and industry averages.
    Full version includes multi-jurisdiction analysis and portfolio optimization.
    """
    
    def __init__(self):
        self.uspto_fees = USPTOFeeSchedule()
        self.attorney_costs = AttorneyCostEstimates()
        
        # Trade secret costs
        self.trade_secret_costs = {
            "documentation_per_compound": 2000,
            "security_annual": 5000,
            "nda_management_annual": 1500,
            "training_annual": 2000
        }
        
        # International filing costs (estimates)
        self.international_costs = {
            "PCT": {
                "filing": 4000,
                "search": 2500,
                "national_phase_per_country": 3000
            },
            "EPO": {
                "filing": 3500,
                "examination": 2800,
                "grant": 1500,
                "maintenance_annual_avg": 1200
            },
            "CNIPA": {
                "filing": 800,
                "examination": 600,
                "maintenance_annual_avg": 400
            }
        }
    
    def estimate_costs(
        self,
        num_compounds: int,
        jurisdictions: List[str] = None,
        claim_count_avg: int = 15,
        independent_claims_avg: int = 3,
        complexity: str = "moderate",
        strategy: str = "patent",
        entity_size: str = "large",
        expected_office_actions: int = 2
    ) -> CostEstimate:
        """Estimate total IP protection costs.
        
        Args:
            num_compounds: Number of compounds to protect
            jurisdictions: List of countries (US only for MVP)
            claim_count_avg: Average claims per patent
            independent_claims_avg: Average independent claims per patent
            complexity: "simple", "moderate", or "complex"
            strategy: "patent", "trade_secret", or "hybrid"
            entity_size: "large", "small", or "micro"
            expected_office_actions: Expected number of office actions
            
        Returns:
            Detailed cost estimate with schedule and recommendations
        """
        jurisdictions = jurisdictions or ["US"]
        complexity_level = ComplexityLevel(complexity)
        entity = EntitySize(entity_size)
        strat = ProtectionStrategy(strategy)
        
        logger.info(f"Estimating costs for {num_compounds} compounds, strategy: {strategy}")
        
        if strat == ProtectionStrategy.PATENT:
            return self._estimate_patent_costs(
                num_compounds,
                jurisdictions,
                claim_count_avg,
                independent_claims_avg,
                complexity_level,
                entity,
                expected_office_actions
            )
        elif strat == ProtectionStrategy.TRADE_SECRET:
            return self._estimate_trade_secret_costs(
                num_compounds,
                complexity_level,
                entity
            )
        else:  # HYBRID
            return self._estimate_hybrid_costs(
                num_compounds,
                jurisdictions,
                claim_count_avg,
                independent_claims_avg,
                complexity_level,
                entity,
                expected_office_actions
            )
    
    def _estimate_patent_costs(
        self,
        num_compounds: int,
        jurisdictions: List[str],
        claim_count_avg: int,
        independent_claims_avg: int,
        complexity: ComplexityLevel,
        entity: EntitySize,
        expected_office_actions: int
    ) -> CostEstimate:
        """Estimate patent costs."""
        # Apply entity discount
        fees = self.uspto_fees.apply_entity_discount(entity)
        
        # Per-patent USPTO fees
        filing_cost_per = self._calculate_filing_cost(
            fees, claim_count_avg, independent_claims_avg
        )
        
        # Per-patent attorney costs
        prosecution_cost_per = self._calculate_prosecution_cost(
            complexity, expected_office_actions
        )
        
        # Maintenance costs
        maintenance_cost_per = fees.get_maintenance_total()
        
        # Total per patent
        total_per_patent = filing_cost_per + prosecution_cost_per + maintenance_cost_per
        
        # Total costs
        filing_total = filing_cost_per * num_compounds
        prosecution_total = prosecution_cost_per * num_compounds
        maintenance_total = maintenance_cost_per * num_compounds
        total_20_year = total_per_patent * num_compounds
        
        # Generate cost schedule
        schedule = self._generate_cost_schedule(
            filing_cost_per, prosecution_cost_per, fees, num_compounds
        )
        
        return CostEstimate(
            strategy=ProtectionStrategy.PATENT,
            num_compounds=num_compounds,
            jurisdictions=jurisdictions,
            entity_size=entity,
            complexity=complexity,
            total_cost_estimate=CostBreakdown(
                filing_costs=filing_total,
                prosecution_costs=prosecution_total,
                maintenance_costs=maintenance_total,
                total_cost=total_20_year
            ),
            per_patent_breakdown=CostBreakdown(
                filing_costs=filing_cost_per,
                prosecution_costs=prosecution_cost_per,
                maintenance_costs=maintenance_cost_per,
                total_cost=total_per_patent
            ),
            timeline={
                "filing_to_examination": "12-18 months",
                "examination_to_grant": "18-36 months",
                "total_to_grant": "2.5-4.5 years",
                "protection_duration": "20 years from filing"
            },
            cost_schedule=schedule,
            advantages=[
                "Legal exclusivity for 20 years",
                "Can license to others for revenue",
                "Blocking power against competitors",
                "Valuable business asset",
                "Public recognition of innovation"
            ],
            disadvantages=[
                "Public disclosure required",
                "Fixed 20-year term",
                "Ongoing maintenance fees",
                "Prosecution uncertainty",
                "Can be designed around"
            ],
            recommendation=self._generate_patent_recommendation(
                num_compounds, total_20_year
            )
        )
    
    def _estimate_trade_secret_costs(
        self,
        num_compounds: int,
        complexity: ComplexityLevel,
        entity: EntitySize
    ) -> CostEstimate:
        """Estimate trade secret protection costs."""
        # Documentation costs
        doc_per_compound = self.trade_secret_costs["documentation_per_compound"]
        if complexity == ComplexityLevel.COMPLEX:
            doc_per_compound *= 1.5
        elif complexity == ComplexityLevel.SIMPLE:
            doc_per_compound *= 0.75
        
        documentation_total = int(doc_per_compound * num_compounds)
        
        # Annual security costs (20 years)
        annual_security = (
            self.trade_secret_costs["security_annual"] +
            self.trade_secret_costs["nda_management_annual"] +
            self.trade_secret_costs["training_annual"]
        )
        security_20_year = annual_security * 20
        
        total_20_year = documentation_total + security_20_year
        
        # Generate schedule
        schedule = [
            CostScheduleEvent(0, "Initial documentation", documentation_total, documentation_total),
            CostScheduleEvent(1, "Year 1 security/NDA", annual_security, documentation_total + annual_security),
        ]
        
        cumulative = documentation_total + annual_security
        for year in [5, 10, 15, 20]:
            years_cost = annual_security * (year - (1 if year == 5 else (year - 5)))
            cumulative += years_cost
            schedule.append(
                CostScheduleEvent(year, f"Security through year {year}", years_cost, cumulative)
            )
        
        return CostEstimate(
            strategy=ProtectionStrategy.TRADE_SECRET,
            num_compounds=num_compounds,
            jurisdictions=["Global"],
            entity_size=entity,
            complexity=complexity,
            total_cost_estimate=CostBreakdown(
                filing_costs=documentation_total,
                prosecution_costs=0,
                maintenance_costs=security_20_year,
                total_cost=total_20_year
            ),
            per_patent_breakdown=CostBreakdown(
                filing_costs=int(doc_per_compound),
                prosecution_costs=0,
                maintenance_costs=int(security_20_year / num_compounds),
                total_cost=int(total_20_year / num_compounds)
            ),
            timeline={
                "protection_start": "Immediately upon documentation",
                "protection_duration": "Indefinite (while secret maintained)",
                "loss_of_protection": "If independently discovered or reverse-engineered"
            },
            cost_schedule=schedule,
            advantages=[
                "No expiration (vs. 20 years for patents)",
                "No public disclosure required",
                "Lower upfront costs",
                "No maintenance fees to USPTO",
                "Global protection without filing"
            ],
            disadvantages=[
                "No legal exclusivity if reverse-engineered",
                "Cannot prevent independent discovery",
                "Difficult to enforce",
                "No blocking power against competitors",
                "Value lost if employee leaves",
                "Cannot license without disclosure risk"
            ],
            recommendation=(
                "Trade secrets best for: manufacturing processes, formulations "
                "difficult to reverse-engineer, and proprietary methods. "
                "Patents best for: novel compound structures, therapeutic uses, "
                "and inventions that will be publicly disclosed anyway."
            )
        )
    
    def _estimate_hybrid_costs(
        self,
        num_compounds: int,
        jurisdictions: List[str],
        claim_count_avg: int,
        independent_claims_avg: int,
        complexity: ComplexityLevel,
        entity: EntitySize,
        expected_office_actions: int
    ) -> CostEstimate:
        """Estimate hybrid strategy costs (patent core, trade secret for processes)."""
        # Patent half the compounds (core innovations)
        patent_compounds = max(1, num_compounds // 2)
        trade_secret_compounds = num_compounds - patent_compounds
        
        # Get patent estimate
        patent_estimate = self._estimate_patent_costs(
            patent_compounds, jurisdictions, claim_count_avg,
            independent_claims_avg, complexity, entity, expected_office_actions
        )
        
        # Get trade secret estimate
        ts_estimate = self._estimate_trade_secret_costs(
            trade_secret_compounds, complexity, entity
        )
        
        # Combine
        total_filing = patent_estimate.total_cost_estimate.filing_costs + ts_estimate.total_cost_estimate.filing_costs
        total_prosecution = patent_estimate.total_cost_estimate.prosecution_costs
        total_maintenance = patent_estimate.total_cost_estimate.maintenance_costs + ts_estimate.total_cost_estimate.maintenance_costs
        total_cost = total_filing + total_prosecution + total_maintenance
        
        # Calculate savings vs. all-patent
        all_patent = self._estimate_patent_costs(
            num_compounds, jurisdictions, claim_count_avg,
            independent_claims_avg, complexity, entity, expected_office_actions
        )
        savings = all_patent.total_cost_estimate.total_cost - total_cost
        
        return CostEstimate(
            strategy=ProtectionStrategy.HYBRID,
            num_compounds=num_compounds,
            jurisdictions=jurisdictions,
            entity_size=entity,
            complexity=complexity,
            total_cost_estimate=CostBreakdown(
                filing_costs=total_filing,
                prosecution_costs=total_prosecution,
                maintenance_costs=total_maintenance,
                total_cost=total_cost
            ),
            per_patent_breakdown=CostBreakdown(
                filing_costs=total_filing // num_compounds,
                prosecution_costs=total_prosecution // num_compounds,
                maintenance_costs=total_maintenance // num_compounds,
                total_cost=total_cost // num_compounds
            ),
            timeline={
                "patent_portion": f"{patent_compounds} compounds patented",
                "trade_secret_portion": f"{trade_secret_compounds} compounds as trade secrets",
                "patent_timeline": "2.5-4.5 years to grant",
                "trade_secret_timeline": "Immediate protection"
            },
            cost_schedule=patent_estimate.cost_schedule,  # Simplified
            advantages=[
                "Cost savings vs. all-patent approach",
                "Legal exclusivity for core innovations",
                "Perpetual protection for processes",
                "Flexible portfolio management",
                "Optimized IP spend"
            ],
            disadvantages=[
                "More complex IP management",
                "Trade secret portions still vulnerable",
                "Requires careful decision on what to patent"
            ],
            recommendation=(
                f"Hybrid approach saves ${savings:,} vs. all-patent strategy. "
                f"Patent the {patent_compounds} most novel compounds, keep "
                f"{trade_secret_compounds} manufacturing/formulation trade secrets."
            ),
            comparison_savings=savings
        )
    
    def _calculate_filing_cost(
        self,
        fees: USPTOFeeSchedule,
        claim_count: int,
        independent_claims: int
    ) -> int:
        """Calculate filing cost including claim fees."""
        base = fees.filing_fee + fees.search_fee + fees.examination_fee + fees.issue_fee
        
        # Extra claim fees
        if claim_count > 20:
            base += (claim_count - 20) * fees.per_claim_over_20
        
        if independent_claims > 3:
            base += (independent_claims - 3) * fees.per_independent_claim_over_3
        
        return base
    
    def _calculate_prosecution_cost(
        self,
        complexity: ComplexityLevel,
        expected_office_actions: int
    ) -> int:
        """Calculate attorney prosecution costs."""
        drafting = self.attorney_costs.get_drafting_cost(complexity)
        prosecution = self.attorney_costs.get_prosecution_cost(complexity)
        
        # Office action responses
        oa_cost = expected_office_actions * self.attorney_costs.office_action_response
        
        return drafting + prosecution + oa_cost
    
    def _generate_cost_schedule(
        self,
        filing_per: int,
        prosecution_per: int,
        fees: USPTOFeeSchedule,
        num_compounds: int
    ) -> List[CostScheduleEvent]:
        """Generate payment schedule over time."""
        schedule = []
        cumulative = 0
        
        # Year 0: Filing
        filing_total = filing_per * num_compounds
        cumulative += filing_total
        schedule.append(CostScheduleEvent(0, "Filing", filing_total, cumulative))
        
        # Year 1-2: Prosecution
        prosecution_total = prosecution_per * num_compounds
        cumulative += prosecution_total
        schedule.append(CostScheduleEvent(1.5, "Prosecution", prosecution_total, cumulative))
        
        # Maintenance fees
        maint_3_5 = fees.maintenance_3_5_years * num_compounds
        cumulative += maint_3_5
        schedule.append(CostScheduleEvent(3.5, "1st maintenance fee", maint_3_5, cumulative))
        
        maint_7_5 = fees.maintenance_7_5_years * num_compounds
        cumulative += maint_7_5
        schedule.append(CostScheduleEvent(7.5, "2nd maintenance fee", maint_7_5, cumulative))
        
        maint_11_5 = fees.maintenance_11_5_years * num_compounds
        cumulative += maint_11_5
        schedule.append(CostScheduleEvent(11.5, "3rd maintenance fee", maint_11_5, cumulative))
        
        return schedule
    
    def _generate_patent_recommendation(
        self,
        num_compounds: int,
        total_cost: int
    ) -> str:
        """Generate cost-based recommendation."""
        cost_per_compound = total_cost / num_compounds
        
        if cost_per_compound > 30000:
            return (
                f"HIGH COST: ${cost_per_compound:,.0f}/compound. Consider prioritizing "
                f"top {num_compounds // 2} compounds and keeping others as trade secrets."
            )
        elif cost_per_compound > 20000:
            return (
                f"MODERATE COST: ${cost_per_compound:,.0f}/compound. Consider filing "
                f"provisionals first to defer costs while validating commercial potential."
            )
        else:
            return (
                f"REASONABLE COST: ${cost_per_compound:,.0f}/compound. Proceed with "
                f"patent filing for maximum IP protection."
            )
    
    def compare_strategies(
        self,
        num_compounds: int,
        complexity: str = "moderate",
        entity_size: str = "large"
    ) -> Dict:
        """Compare patent vs. trade secret costs.
        
        Returns:
            Side-by-side comparison of both strategies
        """
        patent = self.estimate_costs(
            num_compounds=num_compounds,
            complexity=complexity,
            strategy="patent",
            entity_size=entity_size
        )
        
        trade_secret = self.estimate_costs(
            num_compounds=num_compounds,
            complexity=complexity,
            strategy="trade_secret",
            entity_size=entity_size
        )
        
        hybrid = self.estimate_costs(
            num_compounds=num_compounds,
            complexity=complexity,
            strategy="hybrid",
            entity_size=entity_size
        )
        
        patent_cost = patent.total_cost_estimate.total_cost
        ts_cost = trade_secret.total_cost_estimate.total_cost
        hybrid_cost = hybrid.total_cost_estimate.total_cost
        
        return {
            "patent": patent.to_dict(),
            "trade_secret": trade_secret.to_dict(),
            "hybrid": hybrid.to_dict(),
            "comparison": {
                "patent_total": f"${patent_cost:,}",
                "trade_secret_total": f"${ts_cost:,}",
                "hybrid_total": f"${hybrid_cost:,}",
                "patent_vs_ts_difference": f"${patent_cost - ts_cost:,}",
                "hybrid_savings_vs_patent": f"${patent_cost - hybrid_cost:,}",
                "recommended_strategy": self._recommend_strategy(
                    patent_cost, ts_cost, hybrid_cost, num_compounds
                )
            }
        }
    
    def _recommend_strategy(
        self,
        patent_cost: int,
        ts_cost: int,
        hybrid_cost: int,
        num_compounds: int
    ) -> str:
        """Recommend optimal strategy based on costs."""
        if num_compounds == 1:
            return "PATENT: Single compound benefits most from full patent protection"
        elif num_compounds <= 3:
            return "PATENT: Small portfolio warrants full patent coverage"
        elif patent_cost / ts_cost > 3:
            return "HYBRID: Significant cost savings with balanced protection"
        else:
            return "PATENT: Full portfolio protection justified by commercial potential"
    
    def get_entity_discount_info(self) -> Dict:
        """Get information about entity size discounts."""
        return {
            "large_entity": {
                "discount": "0%",
                "qualification": "More than 500 employees or licensed to large entity"
            },
            "small_entity": {
                "discount": "50%",
                "qualification": "Fewer than 500 employees, universities, nonprofits"
            },
            "micro_entity": {
                "discount": "75%",
                "qualification": "Fewer than 5 patents, income below $229,907 (2025)"
            }
        }
