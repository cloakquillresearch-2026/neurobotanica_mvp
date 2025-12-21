"""
Comparative Efficacy Analyzer
Patent Claim 1(j) - FDA Schedule III Compliance Support

Analyzes clinical evidence across the 320-study database to generate
comparative efficacy reports for FDA regulatory submissions.

Features:
- Cross-condition efficacy comparison
- Effect size meta-analysis
- Evidence grade weighting
- FDA-approved drug benchmarking
- Statistical confidence intervals
"""
import json
from typing import List, Dict, Any, Optional
from dataclasses import dataclass
from sqlalchemy.orm import Session
from sqlalchemy import func, and_
import statistics


@dataclass
class EffectSizeStats:
    """Statistical summary of effect sizes."""
    mean: float
    median: float
    std_dev: float
    min_value: float
    max_value: float
    n_studies: int
    confidence_95_lower: float
    confidence_95_upper: float


@dataclass
class EfficacyComparison:
    """Comparison result between cannabinoid and FDA-approved drug."""
    cannabinoid: str
    fda_comparator: str
    condition: str
    cannabinoid_effect: float
    comparator_effect: float
    relative_efficacy: float  # ratio
    statistical_significance: str
    evidence_level: str
    recommendation: str


class ComparativeEfficacyAnalyzer:
    """
    Analyzes comparative efficacy across cannabinoids and conditions.
    
    Supports FDA Schedule III documentation by:
    - Computing weighted effect sizes from RCT data
    - Comparing against FDA-approved drug benchmarks
    - Generating regulatory-ready efficacy summaries
    """
    
    # FDA-approved cannabinoid benchmarks (from literature)
    FDA_BENCHMARKS = {
        "Epidiolex": {
            "compound": "CBD",
            "indications": ["epilepsy", "seizures", "Dravet syndrome", "Lennox-Gastaut"],
            "effect_sizes": {
                "epilepsy": 0.85,
                "seizures": 0.82
            }
        },
        "Marinol": {
            "compound": "Dronabinol (synthetic THC)",
            "indications": ["nausea", "chemotherapy-induced nausea", "appetite", "AIDS wasting"],
            "effect_sizes": {
                "nausea": 0.75,
                "appetite": 0.68,
                "chemotherapy": 0.72
            }
        },
        "Cesamet": {
            "compound": "Nabilone (synthetic THC analog)",
            "indications": ["nausea", "chemotherapy-induced nausea"],
            "effect_sizes": {
                "nausea": 0.71,
                "chemotherapy": 0.69
            }
        },
        "Sativex": {
            "compound": "THC:CBD (1:1)",
            "indications": ["multiple sclerosis", "spasticity", "neuropathic pain"],
            "effect_sizes": {
                "multiple_sclerosis": 0.65,
                "spasticity": 0.62,
                "neuropathic_pain": 0.58
            }
        }
    }
    
    # Evidence grade weights for meta-analysis
    EVIDENCE_WEIGHTS = {
        "Level 1": 1.0,   # High-quality RCT
        "Level 2": 0.8,   # Lower-quality RCT
        "Level 3": 0.5,   # Observational
        "Level 4": 0.3,   # Case series
        "Level 5": 0.1    # Expert opinion
    }
    
    def __init__(self, db: Session):
        self.db = db
    
    def compute_effect_size_stats(
        self, 
        effect_sizes: List[float],
        weights: Optional[List[float]] = None
    ) -> EffectSizeStats:
        """Compute statistical summary of effect sizes."""
        if not effect_sizes:
            return EffectSizeStats(
                mean=0, median=0, std_dev=0, min_value=0, max_value=0,
                n_studies=0, confidence_95_lower=0, confidence_95_upper=0
            )
        
        n = len(effect_sizes)
        
        # Weighted mean if weights provided
        if weights and len(weights) == n:
            weighted_sum = sum(e * w for e, w in zip(effect_sizes, weights))
            weight_sum = sum(weights)
            mean = weighted_sum / weight_sum if weight_sum > 0 else 0
        else:
            mean = statistics.mean(effect_sizes)
        
        median = statistics.median(effect_sizes)
        std_dev = statistics.stdev(effect_sizes) if n > 1 else 0
        
        # 95% confidence interval (assuming normal distribution)
        se = std_dev / (n ** 0.5) if n > 0 else 0
        ci_lower = mean - 1.96 * se
        ci_upper = mean + 1.96 * se
        
        return EffectSizeStats(
            mean=round(mean, 3),
            median=round(median, 3),
            std_dev=round(std_dev, 3),
            min_value=round(min(effect_sizes), 3),
            max_value=round(max(effect_sizes), 3),
            n_studies=n,
            confidence_95_lower=round(ci_lower, 3),
            confidence_95_upper=round(ci_upper, 3)
        )
    
    def analyze_condition(
        self, 
        condition: str,
        cannabinoid: Optional[str] = None
    ) -> Dict[str, Any]:
        """
        Analyze efficacy evidence for a specific condition.
        
        Returns comprehensive efficacy analysis including:
        - Effect size statistics
        - Study quality breakdown
        - Comparisons with FDA-approved drugs
        """
        from backend.models.study import ClinicalStudy
        
        # Build query
        query = self.db.query(ClinicalStudy).filter(
            ClinicalStudy.condition.ilike(f"%{condition}%")
        )
        
        if cannabinoid:
            query = query.filter(
                ClinicalStudy.cannabinoid.ilike(f"%{cannabinoid}%")
            )
        
        studies = query.all()
        
        if not studies:
            return {
                "condition": condition,
                "cannabinoid": cannabinoid,
                "error": "No studies found",
                "studies_analyzed": 0
            }
        
        # Extract effect sizes with evidence weights
        effect_sizes = []
        weights = []
        evidence_breakdown = {"Level 1": 0, "Level 2": 0, "Level 3": 0, "Level 4": 0, "Level 5": 0}
        pivotal_studies = []
        
        for study in studies:
            if study.effect_size_numeric:
                effect_sizes.append(study.effect_size_numeric)
                weight = self.EVIDENCE_WEIGHTS.get(study.evidence_grade, 0.5)
                weights.append(weight)
            
            if study.evidence_grade in evidence_breakdown:
                evidence_breakdown[study.evidence_grade] += 1
            
            if study.is_pivotal_trial:
                pivotal_studies.append({
                    "study_id": study.study_id,
                    "title": study.study_title[:100] if study.study_title else "N/A",
                    "effect_size": study.effect_size_numeric,
                    "evidence_grade": study.evidence_grade
                })
        
        # Compute statistics
        stats = self.compute_effect_size_stats(effect_sizes, weights)
        
        # Find FDA comparator
        comparator = self._find_fda_comparator(condition)
        comparison = None
        
        if comparator and stats.mean > 0:
            comparison = self._compare_to_fda(
                condition=condition,
                cannabinoid_effect=stats.mean,
                fda_drug=comparator["drug"],
                fda_effect=comparator["effect_size"]
            )
        
        return {
            "condition": condition,
            "cannabinoid": cannabinoid or "All cannabinoids",
            "studies_analyzed": len(studies),
            "studies_with_effect_size": len(effect_sizes),
            "effect_size_statistics": {
                "weighted_mean": stats.mean,
                "median": stats.median,
                "std_deviation": stats.std_dev,
                "range": [stats.min_value, stats.max_value],
                "95_ci": [stats.confidence_95_lower, stats.confidence_95_upper]
            },
            "evidence_quality": evidence_breakdown,
            "pivotal_trials": pivotal_studies[:5],  # Top 5
            "fda_comparison": comparison,
            "regulatory_assessment": self._generate_regulatory_assessment(
                stats, evidence_breakdown, comparison
            )
        }
    
    def cross_condition_analysis(
        self, 
        cannabinoid: str
    ) -> Dict[str, Any]:
        """
        Analyze a cannabinoid's efficacy across all conditions.
        
        Useful for identifying therapeutic potential and
        supporting expanded indication applications.
        """
        from backend.models.study import ClinicalStudy
        
        # Get all conditions for this cannabinoid
        conditions = self.db.query(
            ClinicalStudy.condition,
            func.count(ClinicalStudy.id).label('count'),
            func.avg(ClinicalStudy.effect_size_numeric).label('avg_effect')
        ).filter(
            ClinicalStudy.cannabinoid.ilike(f"%{cannabinoid}%")
        ).group_by(
            ClinicalStudy.condition
        ).all()
        
        results = []
        for cond, count, avg_effect in conditions:
            analysis = self.analyze_condition(cond, cannabinoid)
            results.append({
                "condition": cond,
                "study_count": count,
                "average_effect": round(avg_effect, 3) if avg_effect else None,
                "evidence_strength": self._classify_evidence_strength(analysis)
            })
        
        # Sort by evidence strength
        results.sort(key=lambda x: x.get("average_effect") or 0, reverse=True)
        
        return {
            "cannabinoid": cannabinoid,
            "total_conditions_studied": len(results),
            "conditions": results,
            "strongest_indications": results[:3] if results else [],
            "fda_pathway_recommendation": self._recommend_fda_pathway(results)
        }
    
    def head_to_head_comparison(
        self,
        cannabinoid_1: str,
        cannabinoid_2: str,
        condition: str
    ) -> Dict[str, Any]:
        """
        Direct comparison between two cannabinoids for a condition.
        
        Useful for formulation optimization and treatment selection.
        """
        analysis_1 = self.analyze_condition(condition, cannabinoid_1)
        analysis_2 = self.analyze_condition(condition, cannabinoid_2)
        
        effect_1 = analysis_1.get("effect_size_statistics", {}).get("weighted_mean", 0)
        effect_2 = analysis_2.get("effect_size_statistics", {}).get("weighted_mean", 0)
        
        if effect_1 == 0 and effect_2 == 0:
            winner = "Insufficient data"
            difference = 0
        elif effect_1 > effect_2:
            winner = cannabinoid_1
            difference = effect_1 - effect_2
        elif effect_2 > effect_1:
            winner = cannabinoid_2
            difference = effect_2 - effect_1
        else:
            winner = "Equivalent"
            difference = 0
        
        return {
            "condition": condition,
            "comparison": {
                cannabinoid_1: {
                    "effect_size": effect_1,
                    "study_count": analysis_1.get("studies_analyzed", 0),
                    "evidence_quality": analysis_1.get("evidence_quality", {})
                },
                cannabinoid_2: {
                    "effect_size": effect_2,
                    "study_count": analysis_2.get("studies_analyzed", 0),
                    "evidence_quality": analysis_2.get("evidence_quality", {})
                }
            },
            "superior_cannabinoid": winner,
            "effect_difference": round(difference, 3),
            "clinical_significance": "Clinically significant" if difference > 0.2 else "Not clinically significant",
            "recommendation": self._generate_head_to_head_recommendation(
                cannabinoid_1, cannabinoid_2, effect_1, effect_2, condition
            )
        }
    
    def _find_fda_comparator(self, condition: str) -> Optional[Dict]:
        """Find the most relevant FDA-approved comparator for a condition."""
        condition_lower = condition.lower()
        
        for drug, info in self.FDA_BENCHMARKS.items():
            for indication in info["indications"]:
                if indication.lower() in condition_lower or condition_lower in indication.lower():
                    # Find best matching effect size
                    for eff_condition, effect in info["effect_sizes"].items():
                        if eff_condition.lower() in condition_lower or condition_lower in eff_condition.lower():
                            return {
                                "drug": drug,
                                "compound": info["compound"],
                                "indication": indication,
                                "effect_size": effect
                            }
        return None
    
    def _compare_to_fda(
        self,
        condition: str,
        cannabinoid_effect: float,
        fda_drug: str,
        fda_effect: float
    ) -> Dict[str, Any]:
        """Compare cannabinoid efficacy to FDA-approved drug."""
        if fda_effect == 0:
            relative = 0
            assessment = "No comparison available"
        else:
            relative = cannabinoid_effect / fda_effect
            if relative >= 1.0:
                assessment = "Non-inferior to FDA-approved treatment"
            elif relative >= 0.8:
                assessment = "Potentially non-inferior (within 80% threshold)"
            else:
                assessment = "May require additional evidence for non-inferiority"
        
        return {
            "fda_comparator": fda_drug,
            "fda_effect_size": fda_effect,
            "cannabinoid_effect_size": cannabinoid_effect,
            "relative_efficacy": round(relative, 3),
            "non_inferiority_assessment": assessment,
            "regulatory_implication": self._get_regulatory_implication(relative)
        }
    
    def _get_regulatory_implication(self, relative_efficacy: float) -> str:
        """Determine regulatory pathway implication based on relative efficacy."""
        if relative_efficacy >= 1.1:
            return "Strong case for efficacy claims; consider superiority trial design"
        elif relative_efficacy >= 1.0:
            return "Supports non-inferiority; standard 505(b)(2) pathway appropriate"
        elif relative_efficacy >= 0.8:
            return "Borderline non-inferiority; may need larger confirmatory studies"
        else:
            return "Additional efficacy evidence needed; consider dose optimization"
    
    def _generate_regulatory_assessment(
        self,
        stats: EffectSizeStats,
        evidence_breakdown: Dict,
        comparison: Optional[Dict]
    ) -> Dict[str, Any]:
        """Generate regulatory assessment summary."""
        # Calculate evidence quality score
        total_studies = sum(evidence_breakdown.values())
        if total_studies == 0:
            quality_score = 0
        else:
            weighted_score = (
                evidence_breakdown["Level 1"] * 1.0 +
                evidence_breakdown["Level 2"] * 0.8 +
                evidence_breakdown["Level 3"] * 0.5 +
                evidence_breakdown["Level 4"] * 0.3 +
                evidence_breakdown["Level 5"] * 0.1
            ) / total_studies
            quality_score = round(weighted_score * 100, 1)
        
        # Determine regulatory readiness
        if stats.n_studies >= 5 and evidence_breakdown["Level 1"] >= 2 and stats.mean >= 0.5:
            readiness = "High - Sufficient for NDA submission"
        elif stats.n_studies >= 3 and (evidence_breakdown["Level 1"] + evidence_breakdown["Level 2"]) >= 1:
            readiness = "Moderate - May support 505(b)(2) pathway"
        else:
            readiness = "Low - Additional controlled trials recommended"
        
        return {
            "evidence_quality_score": quality_score,
            "regulatory_readiness": readiness,
            "rct_count": evidence_breakdown["Level 1"] + evidence_breakdown["Level 2"],
            "effect_size_adequate": stats.mean >= 0.3,
            "sample_size_adequate": stats.n_studies >= 3,
            "recommended_next_steps": self._get_next_steps(stats, evidence_breakdown)
        }
    
    def _get_next_steps(
        self,
        stats: EffectSizeStats,
        evidence_breakdown: Dict
    ) -> List[str]:
        """Recommend next steps for regulatory pathway."""
        steps = []
        
        if evidence_breakdown["Level 1"] < 2:
            steps.append("Conduct additional Phase III RCTs")
        
        if stats.n_studies < 5:
            steps.append("Expand evidence base with controlled studies")
        
        if stats.mean < 0.5:
            steps.append("Optimize dosing regimen for improved efficacy")
        
        if stats.std_dev > 0.5:
            steps.append("Investigate sources of efficacy variability")
        
        if not steps:
            steps.append("Evidence base supports regulatory submission")
            steps.append("Prepare Chemistry Manufacturing and Controls (CMC) documentation")
            steps.append("Schedule pre-NDA meeting with FDA")
        
        return steps
    
    def _classify_evidence_strength(self, analysis: Dict) -> str:
        """Classify overall evidence strength."""
        quality = analysis.get("evidence_quality", {})
        stats = analysis.get("effect_size_statistics", {})
        
        rct_count = quality.get("Level 1", 0) + quality.get("Level 2", 0)
        mean_effect = stats.get("weighted_mean", 0)
        
        if rct_count >= 3 and mean_effect >= 0.6:
            return "Strong"
        elif rct_count >= 1 and mean_effect >= 0.4:
            return "Moderate"
        elif mean_effect >= 0.2:
            return "Weak"
        else:
            return "Insufficient"
    
    def _recommend_fda_pathway(self, results: List[Dict]) -> str:
        """Recommend FDA approval pathway based on evidence."""
        if not results:
            return "Insufficient evidence for regulatory pathway recommendation"
        
        strong_indications = [r for r in results if r.get("evidence_strength") == "Strong"]
        moderate_indications = [r for r in results if r.get("evidence_strength") == "Moderate"]
        
        if len(strong_indications) >= 2:
            return "505(b)(2) NDA with multiple indications - Strong evidence supports expedited pathway"
        elif len(strong_indications) >= 1:
            return "505(b)(2) NDA for primary indication - Consider breakthrough therapy designation"
        elif len(moderate_indications) >= 2:
            return "Phase III trials recommended - Evidence supports IND application"
        else:
            return "Additional preclinical/Phase II studies needed before IND submission"
    
    def _generate_head_to_head_recommendation(
        self,
        cannabinoid_1: str,
        cannabinoid_2: str,
        effect_1: float,
        effect_2: float,
        condition: str
    ) -> str:
        """Generate recommendation for head-to-head comparison."""
        if effect_1 == 0 and effect_2 == 0:
            return f"Insufficient efficacy data for {condition}. Additional controlled trials needed."
        
        winner = cannabinoid_1 if effect_1 > effect_2 else cannabinoid_2
        difference = abs(effect_1 - effect_2)
        
        if difference < 0.1:
            return f"Both {cannabinoid_1} and {cannabinoid_2} show similar efficacy for {condition}. Consider combination therapy or patient-specific factors for treatment selection."
        elif difference < 0.3:
            return f"{winner} shows modest advantage for {condition}. Clinical significance may depend on individual patient response and tolerability."
        else:
            return f"{winner} demonstrates substantial efficacy advantage for {condition}. Recommend as primary treatment option pending safety profile review."
