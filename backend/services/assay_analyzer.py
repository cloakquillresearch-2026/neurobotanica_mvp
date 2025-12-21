"""
Assay Analyzer Service
Detects heterogeneity in receptor binding assay data

Supports NeuroBotanica patent claims:
- Quality-weighted evidence aggregation
- Cross-source validation
- Assay standardization and normalization

Reference: NeuroBotanica MVP Development Plan - Week 3 Task 3.3
"""
import math
from typing import Optional, List, Dict, Any, Tuple
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
import logging

logger = logging.getLogger(__name__)


class HeterogeneitySeverity(str, Enum):
    """Severity levels for heterogeneity issues."""
    CRITICAL = "critical"
    HIGH = "high"
    MEDIUM = "medium"
    LOW = "low"
    NONE = "none"


@dataclass
class HeterogeneitIssue:
    """Single heterogeneity issue detected."""
    issue_type: str
    severity: HeterogeneitySeverity
    description: str
    affected_measurements: List[int] = field(default_factory=list)
    remediation: Optional[str] = None


@dataclass
class HeterogeneitReport:
    """Full heterogeneity analysis report."""
    compound_name: str
    receptor: str
    total_measurements: int
    issues: List[HeterogeneitIssue]
    heterogeneity_score: float  # 0.0 (homogeneous) to 1.0 (highly heterogeneous)
    can_aggregate: bool
    aggregation_warnings: List[str]
    recommended_actions: List[str]
    analysis_timestamp: datetime = field(default_factory=datetime.utcnow)
    
    @property
    def overall_severity(self) -> HeterogeneitySeverity:
        """Get overall severity from all issues."""
        if not self.issues:
            return HeterogeneitySeverity.NONE
        
        severity_order = [
            HeterogeneitySeverity.CRITICAL,
            HeterogeneitySeverity.HIGH,
            HeterogeneitySeverity.MEDIUM,
            HeterogeneitySeverity.LOW
        ]
        
        for severity in severity_order:
            if any(i.severity == severity for i in self.issues):
                return severity
        
        return HeterogeneitySeverity.NONE


class AssayAnalyzer:
    """Analyzer for detecting heterogeneity in receptor binding data.
    
    Implements:
    - Unit consistency checking
    - Assay type variance analysis
    - Statistical outlier detection
    - Source agreement scoring
    - Quality-weighted aggregation
    """
    
    # Standard unit conversions to nM
    UNIT_CONVERSIONS = {
        "pM": 0.001,
        "nM": 1.0,
        "µM": 1000.0,
        "uM": 1000.0,
        "mM": 1000000.0,
        "M": 1000000000.0
    }
    
    # Assay type groupings for compatibility
    COMPATIBLE_ASSAY_GROUPS = {
        "binding": ["radioligand_binding", "SPR", "ITC"],
        "functional_camp": ["functional_cAMP", "GTPγS"],
        "functional_signaling": ["beta_arrestin", "functional_calcium"],
        "biophysical": ["BRET", "FRET"]
    }
    
    # Thresholds for heterogeneity detection
    CV_THRESHOLD_HIGH = 1.0  # 100% coefficient of variation
    CV_THRESHOLD_MEDIUM = 0.5  # 50%
    OUTLIER_THRESHOLD = 2.5  # Standard deviations from mean
    
    def __init__(self, strict_mode: bool = False):
        """Initialize analyzer.
        
        Args:
            strict_mode: If True, use stricter thresholds for heterogeneity
        """
        self.strict_mode = strict_mode
        if strict_mode:
            self.CV_THRESHOLD_HIGH = 0.75
            self.CV_THRESHOLD_MEDIUM = 0.35
    
    def analyze(
        self,
        measurements: List[Dict[str, Any]],
        compound_name: str,
        receptor: str
    ) -> HeterogeneitReport:
        """Perform full heterogeneity analysis on measurements.
        
        Args:
            measurements: List of measurement dictionaries with structure:
                {
                    "id": int,
                    "value": float,
                    "unit": str,
                    "assay_type": str,
                    "source_type": str,
                    "confidence": float,
                    ...
                }
            compound_name: Name of compound being analyzed
            receptor: Target receptor
            
        Returns:
            HeterogeneitReport with all detected issues and recommendations
        """
        if not measurements:
            return HeterogeneitReport(
                compound_name=compound_name,
                receptor=receptor,
                total_measurements=0,
                issues=[],
                heterogeneity_score=0.0,
                can_aggregate=False,
                aggregation_warnings=["No measurements to analyze"],
                recommended_actions=["Collect receptor binding data from literature or databases"]
            )
        
        issues = []
        warnings = []
        
        # 1. Unit consistency analysis
        unit_issues = self._analyze_units(measurements)
        issues.extend(unit_issues)
        
        # 2. Assay type variance analysis
        assay_issues = self._analyze_assay_types(measurements)
        issues.extend(assay_issues)
        
        # 3. Statistical variance analysis
        variance_issues = self._analyze_variance(measurements)
        issues.extend(variance_issues)
        
        # 4. Outlier detection
        outlier_issues = self._detect_outliers(measurements)
        issues.extend(outlier_issues)
        
        # 5. Source agreement analysis
        source_issues = self._analyze_source_agreement(measurements)
        issues.extend(source_issues)
        
        # 6. Organism/cell type consistency
        organism_issues = self._analyze_organism_consistency(measurements)
        issues.extend(organism_issues)
        
        # Calculate heterogeneity score
        het_score = self._calculate_heterogeneity_score(issues, len(measurements))
        
        # Determine if aggregation is safe
        can_aggregate = self._can_safely_aggregate(issues, het_score)
        
        # Generate warnings and recommendations
        if not can_aggregate:
            warnings.append("High heterogeneity detected - aggregation not recommended")
        
        recommendations = self._generate_recommendations(issues, measurements)
        
        return HeterogeneitReport(
            compound_name=compound_name,
            receptor=receptor,
            total_measurements=len(measurements),
            issues=issues,
            heterogeneity_score=het_score,
            can_aggregate=can_aggregate,
            aggregation_warnings=warnings,
            recommended_actions=recommendations
        )
    
    def _analyze_units(self, measurements: List[Dict[str, Any]]) -> List[HeterogeneitIssue]:
        """Check for unit inconsistencies."""
        issues = []
        
        units = {}
        for i, m in enumerate(measurements):
            unit = m.get("unit", "nM")
            if unit not in units:
                units[unit] = []
            units[unit].append(m.get("id", i))
        
        if len(units) > 1:
            # Check if units are convertible
            convertible = all(u in self.UNIT_CONVERSIONS for u in units.keys())
            
            if not convertible:
                issues.append(HeterogeneitIssue(
                    issue_type="unit_mismatch",
                    severity=HeterogeneitySeverity.CRITICAL,
                    description=f"Incompatible units detected: {list(units.keys())}",
                    affected_measurements=[id for ids in units.values() for id in ids],
                    remediation="Manual unit conversion or exclusion required"
                ))
            else:
                issues.append(HeterogeneitIssue(
                    issue_type="unit_variance",
                    severity=HeterogeneitySeverity.LOW,
                    description=f"Multiple convertible units: {list(units.keys())}",
                    affected_measurements=[id for ids in units.values() for id in ids],
                    remediation="Automatic unit normalization available"
                ))
        
        return issues
    
    def _analyze_assay_types(self, measurements: List[Dict[str, Any]]) -> List[HeterogeneitIssue]:
        """Analyze assay type consistency."""
        issues = []
        
        assay_types = {}
        for i, m in enumerate(measurements):
            atype = m.get("assay_type", "unknown")
            if atype not in assay_types:
                assay_types[atype] = []
            assay_types[atype].append(m.get("id", i))
        
        if len(assay_types) == 1:
            return issues
        
        # Check if assay types are compatible
        detected_groups = set()
        for atype in assay_types.keys():
            for group_name, group_types in self.COMPATIBLE_ASSAY_GROUPS.items():
                if atype in group_types:
                    detected_groups.add(group_name)
                    break
            else:
                detected_groups.add(f"other_{atype}")
        
        if len(detected_groups) > 1:
            issues.append(HeterogeneitIssue(
                issue_type="assay_incompatibility",
                severity=HeterogeneitySeverity.HIGH,
                description=f"Incompatible assay types: {list(assay_types.keys())} span groups: {detected_groups}",
                affected_measurements=[id for ids in assay_types.values() for id in ids],
                remediation="Stratify analysis by assay type or use assay-specific weighting"
            ))
        elif len(assay_types) > 2:
            issues.append(HeterogeneitIssue(
                issue_type="assay_variance",
                severity=HeterogeneitySeverity.MEDIUM,
                description=f"Multiple assay types: {list(assay_types.keys())}",
                affected_measurements=[id for ids in assay_types.values() for id in ids],
                remediation="Consider assay-weighted aggregation"
            ))
        
        return issues
    
    def _analyze_variance(self, measurements: List[Dict[str, Any]]) -> List[HeterogeneitIssue]:
        """Analyze statistical variance in measurements."""
        issues = []
        
        # Normalize values to nM
        normalized = []
        for m in measurements:
            value = m.get("value")
            unit = m.get("unit", "nM")
            if value and unit in self.UNIT_CONVERSIONS:
                normalized.append(value * self.UNIT_CONVERSIONS[unit])
        
        if len(normalized) < 2:
            return issues
        
        # Calculate statistics
        mean_val = sum(normalized) / len(normalized)
        variance = sum((v - mean_val) ** 2 for v in normalized) / len(normalized)
        std_dev = math.sqrt(variance) if variance > 0 else 0
        cv = std_dev / mean_val if mean_val > 0 else 0
        
        # Log-transformed CV (more appropriate for binding affinities)
        log_values = [math.log10(v) for v in normalized if v > 0]
        if len(log_values) > 1:
            log_mean = sum(log_values) / len(log_values)
            log_variance = sum((v - log_mean) ** 2 for v in log_values) / len(log_values)
            log_std = math.sqrt(log_variance) if log_variance > 0 else 0
            
            # Log-transformed values typically have lower CV
            if log_std > 1.0:  # >1 log unit spread
                issues.append(HeterogeneitIssue(
                    issue_type="high_variance",
                    severity=HeterogeneitySeverity.HIGH,
                    description=f"High log-scale variance: ±{log_std:.2f} log units (values span {10**log_std:.1f}-fold range)",
                    affected_measurements=[m.get("id", i) for i, m in enumerate(measurements)],
                    remediation="Review individual measurements for outliers or stratify by assay conditions"
                ))
        
        if cv > self.CV_THRESHOLD_HIGH:
            issues.append(HeterogeneitIssue(
                issue_type="extreme_cv",
                severity=HeterogeneitySeverity.HIGH,
                description=f"Extreme coefficient of variation: {cv:.2%}",
                affected_measurements=[m.get("id", i) for i, m in enumerate(measurements)],
                remediation="Use median or geometric mean instead of arithmetic mean"
            ))
        elif cv > self.CV_THRESHOLD_MEDIUM:
            issues.append(HeterogeneitIssue(
                issue_type="high_cv",
                severity=HeterogeneitySeverity.MEDIUM,
                description=f"High coefficient of variation: {cv:.2%}",
                affected_measurements=[m.get("id", i) for i, m in enumerate(measurements)],
                remediation="Consider confidence-weighted aggregation"
            ))
        
        return issues
    
    def _detect_outliers(self, measurements: List[Dict[str, Any]]) -> List[HeterogeneitIssue]:
        """Detect statistical outliers using modified Z-score."""
        issues = []
        
        # Normalize values
        values_with_ids = []
        for i, m in enumerate(measurements):
            value = m.get("value")
            unit = m.get("unit", "nM")
            if value and unit in self.UNIT_CONVERSIONS:
                norm_val = value * self.UNIT_CONVERSIONS[unit]
                values_with_ids.append((norm_val, m.get("id", i)))
        
        if len(values_with_ids) < 3:
            return issues
        
        values = [v[0] for v in values_with_ids]
        
        # Use log-transformed values for outlier detection
        log_values = [math.log10(v) for v in values if v > 0]
        if len(log_values) < 3:
            return issues
        
        # Calculate median and MAD (Median Absolute Deviation)
        sorted_logs = sorted(log_values)
        median_log = sorted_logs[len(sorted_logs) // 2]
        
        abs_deviations = sorted([abs(v - median_log) for v in log_values])
        mad = abs_deviations[len(abs_deviations) // 2]
        
        if mad == 0:
            return issues
        
        # Modified Z-scores
        outliers = []
        for val, mid in values_with_ids:
            if val > 0:
                log_val = math.log10(val)
                modified_z = 0.6745 * (log_val - median_log) / mad
                if abs(modified_z) > self.OUTLIER_THRESHOLD:
                    outliers.append(mid)
        
        if outliers:
            issues.append(HeterogeneitIssue(
                issue_type="outliers",
                severity=HeterogeneitySeverity.MEDIUM,
                description=f"Detected {len(outliers)} statistical outlier(s)",
                affected_measurements=outliers,
                remediation="Review outliers for data entry errors or exclude from aggregation"
            ))
        
        return issues
    
    def _analyze_source_agreement(self, measurements: List[Dict[str, Any]]) -> List[HeterogeneitIssue]:
        """Analyze agreement between different data sources."""
        issues = []
        
        # Group by source type
        by_source = {}
        for i, m in enumerate(measurements):
            source = m.get("source_type", "unknown")
            if source not in by_source:
                by_source[source] = []
            
            value = m.get("value")
            unit = m.get("unit", "nM")
            if value and unit in self.UNIT_CONVERSIONS:
                by_source[source].append(value * self.UNIT_CONVERSIONS[unit])
        
        if len(by_source) < 2:
            return issues
        
        # Calculate mean per source
        source_means = {}
        for source, values in by_source.items():
            if values:
                source_means[source] = sum(values) / len(values)
        
        if len(source_means) < 2:
            return issues
        
        # Check for disagreement between sources
        mean_values = list(source_means.values())
        overall_mean = sum(mean_values) / len(mean_values)
        
        disagreements = []
        for source, mean in source_means.items():
            fold_diff = max(mean, overall_mean) / min(mean, overall_mean) if min(mean, overall_mean) > 0 else float('inf')
            if fold_diff > 3:  # >3-fold difference
                disagreements.append(f"{source}: {mean:.2f} nM")
        
        if len(disagreements) > 0:
            issues.append(HeterogeneitIssue(
                issue_type="source_disagreement",
                severity=HeterogeneitySeverity.MEDIUM,
                description=f"Source disagreement detected: {', '.join(disagreements)}",
                affected_measurements=[m.get("id", i) for i, m in enumerate(measurements)],
                remediation="Prioritize peer-reviewed sources or use source-weighted averaging"
            ))
        
        return issues
    
    def _analyze_organism_consistency(self, measurements: List[Dict[str, Any]]) -> List[HeterogeneitIssue]:
        """Check for organism/species consistency."""
        issues = []
        
        organisms = {}
        for i, m in enumerate(measurements):
            organism = m.get("organism", m.get("target_organism", "unknown"))
            if organism not in organisms:
                organisms[organism] = []
            organisms[organism].append(m.get("id", i))
        
        if len(organisms) > 1:
            # Human vs non-human is most important
            human_present = any(
                "homo" in o.lower() or "human" in o.lower()
                for o in organisms.keys()
            )
            non_human = [o for o in organisms.keys() 
                        if "homo" not in o.lower() and "human" not in o.lower()]
            
            if human_present and non_human:
                issues.append(HeterogeneitIssue(
                    issue_type="species_variance",
                    severity=HeterogeneitySeverity.MEDIUM,
                    description=f"Mixed species data: human + {non_human}",
                    affected_measurements=[id for ids in organisms.values() for id in ids],
                    remediation="Prioritize human data or stratify analysis by species"
                ))
        
        return issues
    
    def _calculate_heterogeneity_score(
        self,
        issues: List[HeterogeneitIssue],
        measurement_count: int
    ) -> float:
        """Calculate overall heterogeneity score from issues."""
        if not issues:
            return 0.0
        
        # Weight by severity
        severity_weights = {
            HeterogeneitySeverity.CRITICAL: 1.0,
            HeterogeneitySeverity.HIGH: 0.7,
            HeterogeneitySeverity.MEDIUM: 0.4,
            HeterogeneitySeverity.LOW: 0.1,
            HeterogeneitySeverity.NONE: 0.0
        }
        
        total_weight = sum(severity_weights[i.severity] for i in issues)
        
        # Normalize by number of possible issue types (6 types)
        max_possible = 6 * severity_weights[HeterogeneitySeverity.CRITICAL]
        
        score = min(total_weight / max_possible, 1.0)
        
        # Adjust for measurement count (more measurements = more tolerance)
        if measurement_count > 10:
            score *= 0.9
        elif measurement_count > 20:
            score *= 0.8
        
        return round(score, 3)
    
    def _can_safely_aggregate(
        self,
        issues: List[HeterogeneitIssue],
        het_score: float
    ) -> bool:
        """Determine if aggregation is safe given detected issues."""
        # Never aggregate with critical issues
        if any(i.severity == HeterogeneitySeverity.CRITICAL for i in issues):
            return False
        
        # Check heterogeneity score threshold
        if het_score > 0.7:
            return False
        
        # High issues should be reviewed
        high_issues = [i for i in issues if i.severity == HeterogeneitySeverity.HIGH]
        if len(high_issues) > 2:
            return False
        
        return True
    
    def _generate_recommendations(
        self,
        issues: List[HeterogeneitIssue],
        measurements: List[Dict[str, Any]]
    ) -> List[str]:
        """Generate actionable recommendations based on issues."""
        recommendations = []
        
        if not issues:
            recommendations.append("Data appears homogeneous - standard aggregation recommended")
            return recommendations
        
        # Unit issues
        unit_issues = [i for i in issues if "unit" in i.issue_type]
        if unit_issues:
            recommendations.append("Normalize all values to standard units (nM) before aggregation")
        
        # Variance issues
        variance_issues = [i for i in issues if "variance" in i.issue_type or "cv" in i.issue_type]
        if variance_issues:
            recommendations.append("Use geometric mean instead of arithmetic mean for aggregation")
            recommendations.append("Apply confidence weighting to reduce impact of uncertain values")
        
        # Outlier issues
        if any(i.issue_type == "outliers" for i in issues):
            recommendations.append("Review and consider excluding statistical outliers")
        
        # Assay issues
        assay_issues = [i for i in issues if "assay" in i.issue_type]
        if assay_issues:
            recommendations.append("Stratify analysis by assay type for accurate comparison")
        
        # Source issues
        if any(i.issue_type == "source_disagreement" for i in issues):
            recommendations.append("Prioritize peer-reviewed sources in aggregation")
            recommendations.append("Use source quality weighting in confidence calculation")
        
        # Species issues
        if any(i.issue_type == "species_variance" for i in issues):
            recommendations.append("Prioritize human data for clinical relevance")
        
        return recommendations
    
    def aggregate_with_weights(
        self,
        measurements: List[Dict[str, Any]],
        weighting: str = "confidence"
    ) -> Dict[str, Any]:
        """Aggregate measurements with quality weighting.
        
        Args:
            measurements: List of measurement dictionaries
            weighting: Weighting method - "confidence", "source", "equal"
            
        Returns:
            Aggregated result with statistics
        """
        if not measurements:
            return {"error": "No measurements to aggregate"}
        
        # Normalize to nM
        normalized = []
        for m in measurements:
            value = m.get("value")
            unit = m.get("unit", "nM")
            if value and unit in self.UNIT_CONVERSIONS:
                norm_val = value * self.UNIT_CONVERSIONS[unit]
                confidence = m.get("confidence", 0.5)
                
                # Adjust confidence by source quality
                source = m.get("source_type", "unknown")
                source_weights = {
                    "peer_reviewed": 1.0,
                    "database_curated": 0.9,
                    "preprint": 0.7,
                    "vendor_data": 0.5,
                    "predicted": 0.3
                }
                
                if weighting == "confidence":
                    weight = confidence
                elif weighting == "source":
                    weight = source_weights.get(source, 0.5)
                else:
                    weight = 1.0
                
                normalized.append({
                    "value": norm_val,
                    "weight": weight,
                    "log_value": math.log10(norm_val) if norm_val > 0 else None
                })
        
        if not normalized:
            return {"error": "No valid measurements after normalization"}
        
        values = [n["value"] for n in normalized]
        weights = [n["weight"] for n in normalized]
        log_values = [n["log_value"] for n in normalized if n["log_value"] is not None]
        
        # Weighted arithmetic mean
        total_weight = sum(weights)
        weighted_mean = sum(v * w for v, w in zip(values, weights)) / total_weight
        
        # Geometric mean (weighted)
        if log_values:
            weighted_log_mean = sum(
                lv * w for lv, w in zip(log_values, weights) 
                if lv is not None
            ) / total_weight
            geometric_mean = 10 ** weighted_log_mean
        else:
            geometric_mean = None
        
        # Statistics
        sorted_values = sorted(values)
        median = sorted_values[len(sorted_values) // 2]
        
        return {
            "weighted_mean": round(weighted_mean, 3),
            "geometric_mean": round(geometric_mean, 3) if geometric_mean else None,
            "median": round(median, 3),
            "min": round(min(values), 3),
            "max": round(max(values), 3),
            "n_measurements": len(values),
            "total_weight": round(total_weight, 3),
            "unit": "nM",
            "weighting_method": weighting
        }
