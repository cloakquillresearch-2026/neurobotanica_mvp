"""
ToxPath Assessor - Toxicity Risk Assessment Engine

MVP Implementation for toxicity risk evaluation.
Trade secret: Risk scoring weights, alert thresholds, and decision trees are proprietary.
"""

from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
from typing import Optional, List, Dict, Any, Tuple
import uuid
import re


class RiskTier(Enum):
    """Risk tier classification."""
    LOW = 1
    MODERATE = 2
    HIGH = 3
    VERY_HIGH = 4
    
    @property
    def label(self) -> str:
        labels = {
            RiskTier.LOW: "Low Risk",
            RiskTier.MODERATE: "Moderate Risk", 
            RiskTier.HIGH: "High Risk",
            RiskTier.VERY_HIGH: "Very High Risk"
        }
        return labels[self]
    
    @property
    def testing_depth(self) -> str:
        depths = {
            RiskTier.LOW: "Minimal",
            RiskTier.MODERATE: "Standard panel",
            RiskTier.HIGH: "Extended + consult",
            RiskTier.VERY_HIGH: "Full GLP recommended"
        }
        return depths[self]


class AlertSeverity(Enum):
    """Structural alert severity levels."""
    INFO = "info"
    CAUTION = "caution"
    WARNING = "warning"
    CRITICAL = "critical"


class AdministrationRoute(Enum):
    """Routes of administration."""
    ORAL = "oral"
    INHALED = "inhaled"
    SUBLINGUAL = "sublingual"
    TOPICAL = "topical"
    INTRAVENOUS = "intravenous"
    UNKNOWN = "unknown"


@dataclass
class ExposureProfile:
    """Exposure profile for risk assessment."""
    dose_mg: Optional[float] = None
    frequency: Optional[str] = None  # once, daily, multiple_daily
    duration: Optional[str] = None  # acute, subacute, chronic
    
    def to_dict(self) -> Dict:
        return {
            "dose_mg": self.dose_mg,
            "frequency": self.frequency,
            "duration": self.duration
        }


@dataclass
class ImpurityEntry:
    """Known impurity entry."""
    name: str
    ppm: float
    cas_number: Optional[str] = None
    
    def to_dict(self) -> Dict:
        return {
            "name": self.name,
            "ppm": self.ppm,
            "cas_number": self.cas_number
        }


@dataclass
class StructuralAlert:
    """Structural toxicology alert."""
    code: str
    name: str
    severity: AlertSeverity
    description: str
    mechanism: Optional[str] = None
    affected_organs: List[str] = field(default_factory=list)
    
    def to_dict(self) -> Dict:
        return {
            "code": self.code,
            "name": self.name,
            "severity": self.severity.value,
            "description": self.description,
            "mechanism": self.mechanism,
            "affected_organs": self.affected_organs
        }


@dataclass
class TestingStep:
    """Recommended testing step."""
    order: int
    test_name: str
    test_type: str  # in_silico, in_vitro, in_vivo, clinical
    rationale: str
    estimated_cost_range: str
    timeline: str
    required: bool = True
    gmp_glp_required: bool = False
    
    def to_dict(self) -> Dict:
        return {
            "order": self.order,
            "test_name": self.test_name,
            "test_type": self.test_type,
            "rationale": self.rationale,
            "estimated_cost_range": self.estimated_cost_range,
            "timeline": self.timeline,
            "required": self.required,
            "gmp_glp_required": self.gmp_glp_required
        }


@dataclass
class RiskSummary:
    """Risk assessment summary."""
    overall_tier: RiskTier
    tier_rationale: str
    top_risks: List[str]
    key_unknowns: List[str]
    key_assumptions: List[str]
    route_specific_concerns: List[str] = field(default_factory=list)
    
    def to_dict(self) -> Dict:
        return {
            "overall_tier": self.overall_tier.value,
            "tier_label": self.overall_tier.label,
            "testing_depth": self.overall_tier.testing_depth,
            "tier_rationale": self.tier_rationale,
            "top_risks": self.top_risks,
            "key_unknowns": self.key_unknowns,
            "key_assumptions": self.key_assumptions,
            "route_specific_concerns": self.route_specific_concerns
        }


@dataclass
class ToxPathRequest:
    """Request model for ToxPath assessment."""
    compound_ref: str  # canonical_smiles or compound_id
    compound_name: Optional[str] = None
    route: Optional[str] = None
    exposure: Optional[ExposureProfile] = None
    known_impurities: Optional[List[ImpurityEntry]] = None
    chempath_job_id: Optional[str] = None  # Link to ChemPath analysis
    
    def to_dict(self) -> Dict:
        return {
            "compound_ref": self.compound_ref,
            "compound_name": self.compound_name,
            "route": self.route,
            "exposure": self.exposure.to_dict() if self.exposure else None,
            "known_impurities": [i.to_dict() for i in self.known_impurities] if self.known_impurities else None,
            "chempath_job_id": self.chempath_job_id
        }


@dataclass
class ToxPathResponse:
    """Response model for ToxPath assessment."""
    toxpath_assessment_id: str
    compound_ref: str
    compound_name: Optional[str]
    route: str
    risk_summary: RiskSummary
    alerts: List[StructuralAlert]
    testing_plan: List[TestingStep]
    consultation_required: bool
    consultation_flags: List[str]
    properties_snapshot: Dict[str, Any]
    status: str = "completed"
    error: Optional[str] = None
    assessed_at: str = field(default_factory=lambda: datetime.utcnow().isoformat())
    
    def to_dict(self) -> Dict:
        return {
            "toxpath_assessment_id": self.toxpath_assessment_id,
            "compound_ref": self.compound_ref,
            "compound_name": self.compound_name,
            "route": self.route,
            "risk_summary": self.risk_summary.to_dict(),
            "alerts": [a.to_dict() for a in self.alerts],
            "testing_plan": [t.to_dict() for t in self.testing_plan],
            "consultation_required": self.consultation_required,
            "consultation_flags": self.consultation_flags,
            "properties_snapshot": self.properties_snapshot,
            "status": self.status,
            "error": self.error,
            "assessed_at": self.assessed_at
        }


class ToxPathAssessor:
    """Toxicity Risk Assessment Engine.
    
    MVP Implementation Notes:
    - Uses rule-based structural alert detection
    - Integrates with ChemPath for descriptor data
    - Route-aware risk framing
    - Generates testing recommendations
    
    Trade Secret Elements (not exposed in API):
    - Alert detection SMARTS patterns
    - Risk scoring weights
    - Tier classification thresholds
    - Testing plan decision trees
    """
    
    # =========================================================================
    # TRADE SECRET: Structural Alert Patterns
    # These SMARTS patterns and thresholds are proprietary
    # =========================================================================
    
    # Alert codes (exposed) - patterns (secret)
    ALERT_CODES = {
        "ALERT_REACTIVE_METABOLITE": "Potential reactive metabolite formation",
        "ALERT_HEPATOTOXICITY": "Structural features associated with hepatotoxicity",
        "ALERT_CARDIOTOXICITY": "Structural features associated with cardiotoxicity",
        "ALERT_MUTAGENICITY": "Potential mutagenic structural alert",
        "ALERT_GENOTOXICITY": "Potential genotoxic structural feature",
        "ALERT_MICHAEL_ACCEPTOR": "Michael acceptor moiety detected",
        "ALERT_QUINONE": "Quinone/quinone-like structure",
        "ALERT_NITRO": "Nitro group detected",
        "ALERT_EPOXIDE": "Epoxide or epoxide precursor",
        "ALERT_ALDEHYDE": "Aldehyde group detected",
        "ALERT_ACYL_HALIDE": "Acyl halide detected",
        "ALERT_THIOUREA": "Thiourea moiety detected",
        "ALERT_HYDRAZINE": "Hydrazine derivative detected",
        "ALERT_HIGH_LOGP": "High lipophilicity may affect distribution",
        "ALERT_PAINS": "Pan-assay interference compound features",
        "ALERT_CYP_INHIBITION": "Potential CYP450 inhibition",
        "ALERT_PGLYCOPROTEIN": "P-glycoprotein substrate/inhibitor features",
        "ALERT_BBB_PENETRATION": "Blood-brain barrier penetration likely"
    }
    
    # TRADE SECRET: SMARTS patterns for alerts (simplified for MVP)
    _ALERT_SMARTS = {
        "ALERT_MICHAEL_ACCEPTOR": "[C,c]=[C,c][C,c,N,n,O,o]=[O,o,N,n,S,s]",
        "ALERT_QUINONE": "O=C1C=CC(=O)C=C1",
        "ALERT_NITRO": "[N+](=O)[O-]",
        "ALERT_EPOXIDE": "C1OC1",
        "ALERT_ALDEHYDE": "[CH]=O",
        "ALERT_ACYL_HALIDE": "[CX3](=[OX1])[F,Cl,Br,I]",
        "ALERT_THIOUREA": "[NX3][CX3](=[SX1])[NX3]",
        "ALERT_HYDRAZINE": "[NX3][NX3]"
    }
    
    # Route-specific risk modifiers (TRADE SECRET)
    _ROUTE_RISK_MODIFIERS = {
        "oral": {"hepatotoxicity": 1.5, "gi_toxicity": 1.3, "bioavailability": 0.8},
        "inhaled": {"pulmonary": 2.0, "systemic": 1.2, "onset_speed": 2.0},
        "sublingual": {"mucosal": 1.3, "systemic": 1.0, "hepatic_bypass": 0.7},
        "topical": {"dermal": 1.5, "systemic": 0.3, "local_effects": 2.0},
        "intravenous": {"systemic": 2.0, "cardiotoxicity": 1.5, "onset_speed": 3.0}
    }
    
    # Risk tier thresholds (TRADE SECRET)
    _TIER_THRESHOLDS = {
        "critical_alerts": 1,  # Any critical → Very High
        "warning_alerts": 2,   # 2+ warnings → High
        "caution_alerts": 3,   # 3+ cautions → Moderate
        "novelty_penalty": 0.5,  # Unknown structure multiplier
        "characterization_bonus": -0.3  # Good ChemPath data bonus
    }
    
    def __init__(self):
        """Initialize ToxPath assessor."""
        self._rdkit_available = self._check_rdkit()
    
    def _check_rdkit(self) -> bool:
        """Check if RDKit is available."""
        try:
            from rdkit import Chem
            return True
        except ImportError:
            return False
    
    def assess(self, request: ToxPathRequest) -> ToxPathResponse:
        """Run full ToxPath assessment pipeline.
        
        Pipeline steps:
        1. Parse compound structure
        2. Get/compute molecular properties
        3. Run structural alert detection
        4. Calculate risk tier
        5. Generate route-aware risk summary
        6. Build testing plan
        """
        assessment_id = str(uuid.uuid4())
        
        # Step 1: Parse compound
        smiles = request.compound_ref
        properties = self._get_compound_properties(smiles)
        
        # Step 2: Detect structural alerts
        alerts = self._detect_alerts(smiles, properties)
        
        # Step 3: Calculate risk tier
        route = self._parse_route(request.route)
        risk_summary = self._calculate_risk(
            alerts, 
            properties, 
            route,
            request.exposure,
            request.known_impurities
        )
        
        # Step 4: Generate testing plan
        testing_plan = self._generate_testing_plan(
            risk_summary.overall_tier,
            alerts,
            route,
            properties
        )
        
        # Step 5: Determine consultation needs
        consultation_required, consultation_flags = self._check_consultation_needs(
            risk_summary,
            alerts
        )
        
        return ToxPathResponse(
            toxpath_assessment_id=assessment_id,
            compound_ref=smiles,
            compound_name=request.compound_name,
            route=route.value,
            risk_summary=risk_summary,
            alerts=alerts,
            testing_plan=testing_plan,
            consultation_required=consultation_required,
            consultation_flags=consultation_flags,
            properties_snapshot=properties
        )
    
    def _parse_route(self, route: Optional[str]) -> AdministrationRoute:
        """Parse administration route."""
        if not route:
            return AdministrationRoute.UNKNOWN
        
        route_lower = route.lower()
        for r in AdministrationRoute:
            if r.value == route_lower:
                return r
        
        return AdministrationRoute.UNKNOWN
    
    def _get_compound_properties(self, smiles: str) -> Dict[str, Any]:
        """Get compound properties for assessment."""
        properties = {
            "smiles": smiles,
            "is_valid": False,
            "molecular_weight": None,
            "logp": None,
            "tpsa": None,
            "hbd": None,
            "hba": None,
            "rotatable_bonds": None,
            "num_rings": None,
            "num_aromatic_rings": None
        }
        
        if not self._rdkit_available:
            return properties
        
        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors, rdMolDescriptors
            
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return properties
            
            properties["is_valid"] = True
            properties["canonical_smiles"] = Chem.MolToSmiles(mol, canonical=True)
            properties["molecular_weight"] = round(Descriptors.MolWt(mol), 2)
            properties["logp"] = round(Descriptors.MolLogP(mol), 2)
            properties["tpsa"] = round(Descriptors.TPSA(mol), 2)
            properties["hbd"] = Descriptors.NumHDonors(mol)
            properties["hba"] = Descriptors.NumHAcceptors(mol)
            properties["rotatable_bonds"] = rdMolDescriptors.CalcNumRotatableBonds(mol)
            properties["num_rings"] = rdMolDescriptors.CalcNumRings(mol)
            properties["num_aromatic_rings"] = rdMolDescriptors.CalcNumAromaticRings(mol)
            
        except Exception:
            pass
        
        return properties
    
    def _detect_alerts(
        self, 
        smiles: str, 
        properties: Dict[str, Any]
    ) -> List[StructuralAlert]:
        """Detect structural toxicology alerts.
        
        TRADE SECRET: Specific SMARTS patterns and scoring logic.
        """
        alerts: List[StructuralAlert] = []
        
        if not self._rdkit_available or not properties.get("is_valid"):
            return alerts
        
        from rdkit import Chem
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return alerts
        
        # Run SMARTS-based alerts
        for alert_code, smarts in self._ALERT_SMARTS.items():
            try:
                pattern = Chem.MolFromSmarts(smarts)
                if pattern and mol.HasSubstructMatch(pattern):
                    alerts.append(self._create_alert(alert_code))
            except:
                pass
        
        # Property-based alerts
        logp = properties.get("logp", 0)
        if logp and logp > 5:
            alerts.append(StructuralAlert(
                code="ALERT_HIGH_LOGP",
                name="High Lipophilicity",
                severity=AlertSeverity.CAUTION,
                description=f"LogP of {logp} may affect distribution and accumulation",
                mechanism="Increased membrane partitioning and tissue accumulation",
                affected_organs=["liver", "adipose tissue"]
            ))
        
        # Check for BBB penetration potential (cannabinoid-relevant)
        mw = properties.get("molecular_weight", 0)
        tpsa = properties.get("tpsa", 0)
        if mw and mw < 450 and tpsa and tpsa < 90:
            alerts.append(StructuralAlert(
                code="ALERT_BBB_PENETRATION",
                name="CNS Penetration Likely",
                severity=AlertSeverity.INFO,
                description="Compound properties suggest blood-brain barrier penetration",
                mechanism="Low MW and TPSA favor CNS distribution",
                affected_organs=["central nervous system"]
            ))
        
        return alerts
    
    def _create_alert(self, alert_code: str) -> StructuralAlert:
        """Create structural alert from code."""
        alert_info = {
            "ALERT_MICHAEL_ACCEPTOR": {
                "name": "Michael Acceptor",
                "severity": AlertSeverity.WARNING,
                "mechanism": "Covalent binding to proteins via nucleophilic attack",
                "affected_organs": ["liver", "kidney"]
            },
            "ALERT_QUINONE": {
                "name": "Quinone Structure",
                "severity": AlertSeverity.WARNING,
                "mechanism": "Redox cycling and oxidative stress",
                "affected_organs": ["liver", "bone marrow"]
            },
            "ALERT_NITRO": {
                "name": "Nitro Group",
                "severity": AlertSeverity.WARNING,
                "mechanism": "Metabolic reduction to reactive nitroso intermediates",
                "affected_organs": ["liver"]
            },
            "ALERT_EPOXIDE": {
                "name": "Epoxide Moiety",
                "severity": AlertSeverity.CRITICAL,
                "mechanism": "Electrophilic attack on DNA and proteins",
                "affected_organs": ["liver", "systemic"]
            },
            "ALERT_ALDEHYDE": {
                "name": "Aldehyde Group",
                "severity": AlertSeverity.CAUTION,
                "mechanism": "Protein cross-linking and Schiff base formation",
                "affected_organs": ["respiratory tract", "eyes"]
            },
            "ALERT_ACYL_HALIDE": {
                "name": "Acyl Halide",
                "severity": AlertSeverity.CRITICAL,
                "mechanism": "Highly reactive acylating agent",
                "affected_organs": ["skin", "respiratory tract", "eyes"]
            },
            "ALERT_THIOUREA": {
                "name": "Thiourea Moiety",
                "severity": AlertSeverity.WARNING,
                "mechanism": "Thyroid disruption and hepatotoxicity",
                "affected_organs": ["thyroid", "liver"]
            },
            "ALERT_HYDRAZINE": {
                "name": "Hydrazine Derivative",
                "severity": AlertSeverity.WARNING,
                "mechanism": "Oxidative stress and hepatotoxicity",
                "affected_organs": ["liver", "lung"]
            }
        }
        
        info = alert_info.get(alert_code, {
            "name": alert_code,
            "severity": AlertSeverity.CAUTION,
            "mechanism": "Unknown mechanism",
            "affected_organs": []
        })
        
        return StructuralAlert(
            code=alert_code,
            name=info["name"],
            severity=info["severity"],
            description=self.ALERT_CODES.get(alert_code, "Structural alert detected"),
            mechanism=info["mechanism"],
            affected_organs=info["affected_organs"]
        )
    
    def _calculate_risk(
        self,
        alerts: List[StructuralAlert],
        properties: Dict[str, Any],
        route: AdministrationRoute,
        exposure: Optional[ExposureProfile],
        impurities: Optional[List[ImpurityEntry]]
    ) -> RiskSummary:
        """Calculate overall risk tier.
        
        TRADE SECRET: Scoring weights and tier thresholds.
        """
        # Count alerts by severity
        critical_count = len([a for a in alerts if a.severity == AlertSeverity.CRITICAL])
        warning_count = len([a for a in alerts if a.severity == AlertSeverity.WARNING])
        caution_count = len([a for a in alerts if a.severity == AlertSeverity.CAUTION])
        
        # Base tier determination
        if critical_count >= self._TIER_THRESHOLDS["critical_alerts"]:
            tier = RiskTier.VERY_HIGH
            rationale = f"Critical structural alert(s) detected ({critical_count})"
        elif warning_count >= self._TIER_THRESHOLDS["warning_alerts"]:
            tier = RiskTier.HIGH
            rationale = f"Multiple warning-level alerts ({warning_count})"
        elif caution_count >= self._TIER_THRESHOLDS["caution_alerts"]:
            tier = RiskTier.MODERATE
            rationale = f"Several caution-level alerts ({caution_count})"
        elif len(alerts) > 0:
            tier = RiskTier.MODERATE
            rationale = "Structural alerts present but low severity"
        else:
            tier = RiskTier.LOW
            rationale = "No significant structural alerts detected"
        
        # Adjust for characterization quality
        if not properties.get("is_valid"):
            if tier.value < RiskTier.HIGH.value:
                tier = RiskTier(tier.value + 1)
                rationale += "; structure validation issues"
        
        # Route-specific concerns
        route_concerns = self._get_route_concerns(route, alerts, properties)
        
        # Generate top risks
        top_risks = []
        if critical_count > 0:
            top_risks.append("Critical structural alerts requiring GLP toxicology studies")
        if warning_count > 0:
            top_risks.append("Potential for reactive metabolite formation")
        for alert in alerts[:3]:  # Top 3 alerts
            if alert.mechanism:
                top_risks.append(f"{alert.name}: {alert.mechanism}")
        
        # Handle exposure-based risks
        if exposure and exposure.duration == "chronic":
            top_risks.append("Chronic exposure may increase cumulative toxicity risk")
        
        # Handle impurity risks
        if impurities and len(impurities) > 0:
            high_ppm_impurities = [i for i in impurities if i.ppm > 100]
            if high_ppm_impurities:
                top_risks.append(f"Known impurities at elevated levels ({len(high_ppm_impurities)})")
        
        # Key unknowns
        key_unknowns = []
        if not properties.get("is_valid"):
            key_unknowns.append("Structure could not be fully validated")
        if route == AdministrationRoute.UNKNOWN:
            key_unknowns.append("Route of administration not specified")
        if not exposure:
            key_unknowns.append("Exposure profile not provided")
        key_unknowns.append("Long-term safety data not available without studies")
        key_unknowns.append("Individual patient variability not assessable in silico")
        
        # Key assumptions
        key_assumptions = [
            "Assessment based on structural features only",
            "Does not account for formulation effects",
            "Assumes typical metabolic pathways",
            "Human relevance extrapolated from structural alerts"
        ]
        
        return RiskSummary(
            overall_tier=tier,
            tier_rationale=rationale,
            top_risks=top_risks[:5],  # Limit to 5
            key_unknowns=key_unknowns[:5],
            key_assumptions=key_assumptions,
            route_specific_concerns=route_concerns
        )
    
    def _get_route_concerns(
        self,
        route: AdministrationRoute,
        alerts: List[StructuralAlert],
        properties: Dict[str, Any]
    ) -> List[str]:
        """Get route-specific toxicity concerns."""
        concerns = []
        
        if route == AdministrationRoute.ORAL:
            concerns.append("First-pass hepatic metabolism may generate reactive metabolites")
            concerns.append("GI irritation potential should be evaluated")
            if properties.get("logp", 0) < 0:
                concerns.append("Low lipophilicity may limit oral absorption")
        
        elif route == AdministrationRoute.INHALED:
            concerns.append("Pulmonary irritation and deposition must be assessed")
            concerns.append("Rapid systemic exposure via pulmonary route")
            concerns.append("Particle size and formulation critical for safety")
        
        elif route == AdministrationRoute.SUBLINGUAL:
            concerns.append("Mucosal irritation potential")
            concerns.append("Partial hepatic bypass affects metabolite profile")
        
        elif route == AdministrationRoute.TOPICAL:
            concerns.append("Dermal sensitization potential should be assessed")
            concerns.append("Penetration depth affects systemic exposure")
        
        elif route == AdministrationRoute.INTRAVENOUS:
            concerns.append("Direct systemic exposure - no first-pass protection")
            concerns.append("Acute cardiovascular effects priority")
            concerns.append("Sterility and pyrogen testing required")
        
        return concerns
    
    def _generate_testing_plan(
        self,
        tier: RiskTier,
        alerts: List[StructuralAlert],
        route: AdministrationRoute,
        properties: Dict[str, Any]
    ) -> List[TestingStep]:
        """Generate recommended testing sequence.
        
        TRADE SECRET: Testing decision tree and prioritization logic.
        """
        plan: List[TestingStep] = []
        order = 1
        
        # Tier 1-4 all need basic characterization
        plan.append(TestingStep(
            order=order,
            test_name="Structural Confirmation",
            test_type="analytical",
            rationale="Confirm identity and purity before toxicology studies",
            estimated_cost_range="$500-2,000",
            timeline="1-2 weeks",
            required=True,
            gmp_glp_required=False
        ))
        order += 1
        
        # In silico predictions (all tiers)
        plan.append(TestingStep(
            order=order,
            test_name="In Silico ADMET Prediction",
            test_type="in_silico",
            rationale="Computational prediction of absorption, distribution, metabolism, excretion, toxicity",
            estimated_cost_range="$500-1,500",
            timeline="1-3 days",
            required=True,
            gmp_glp_required=False
        ))
        order += 1
        
        if tier.value >= RiskTier.MODERATE.value:
            # In vitro cytotoxicity
            plan.append(TestingStep(
                order=order,
                test_name="In Vitro Cytotoxicity Panel",
                test_type="in_vitro",
                rationale="Assess general cytotoxicity in relevant cell lines",
                estimated_cost_range="$3,000-8,000",
                timeline="2-4 weeks",
                required=True,
                gmp_glp_required=False
            ))
            order += 1
            
            # Ames test for mutagenicity
            has_mutagen_alert = any(a.code in ["ALERT_MUTAGENICITY", "ALERT_GENOTOXICITY", "ALERT_NITRO"] 
                                   for a in alerts)
            plan.append(TestingStep(
                order=order,
                test_name="Ames Mutagenicity Assay",
                test_type="in_vitro",
                rationale="Standard mutagenicity screen" + 
                         (" - prioritized due to structural alerts" if has_mutagen_alert else ""),
                estimated_cost_range="$5,000-15,000",
                timeline="4-6 weeks",
                required=has_mutagen_alert or tier.value >= RiskTier.HIGH.value,
                gmp_glp_required=tier.value >= RiskTier.HIGH.value
            ))
            order += 1
        
        if tier.value >= RiskTier.HIGH.value:
            # Metabolite ID
            plan.append(TestingStep(
                order=order,
                test_name="Metabolite Identification",
                test_type="in_vitro",
                rationale="Identify major metabolites for reactive intermediate assessment",
                estimated_cost_range="$15,000-30,000",
                timeline="4-8 weeks",
                required=True,
                gmp_glp_required=False
            ))
            order += 1
            
            # Hepatotoxicity panel
            has_hepato_alert = any("liver" in a.affected_organs or "hepato" in a.code.lower() 
                                  for a in alerts)
            plan.append(TestingStep(
                order=order,
                test_name="Hepatotoxicity Panel (HepaRG/PHH)",
                test_type="in_vitro",
                rationale="Primary hepatocyte panel for liver toxicity" +
                         (" - prioritized due to hepatotoxicity alerts" if has_hepato_alert else ""),
                estimated_cost_range="$10,000-25,000",
                timeline="3-6 weeks",
                required=has_hepato_alert,
                gmp_glp_required=False
            ))
            order += 1
            
            # Route-specific testing
            if route == AdministrationRoute.INHALED:
                plan.append(TestingStep(
                    order=order,
                    test_name="Pulmonary Cytotoxicity Panel",
                    test_type="in_vitro",
                    rationale="Lung cell line panel for inhalation route",
                    estimated_cost_range="$8,000-20,000",
                    timeline="3-5 weeks",
                    required=True,
                    gmp_glp_required=False
                ))
                order += 1
        
        if tier.value >= RiskTier.VERY_HIGH.value:
            # GLP acute toxicity
            plan.append(TestingStep(
                order=order,
                test_name="GLP Acute Toxicity Study (Rodent)",
                test_type="in_vivo",
                rationale="Full GLP acute toxicity for IND-enabling package",
                estimated_cost_range="$50,000-100,000",
                timeline="8-12 weeks",
                required=True,
                gmp_glp_required=True
            ))
            order += 1
            
            # Genetic toxicology battery
            plan.append(TestingStep(
                order=order,
                test_name="GLP Genetic Toxicology Battery",
                test_type="in_vivo",
                rationale="Full genetic tox package (Ames, chromosomal aberration, micronucleus)",
                estimated_cost_range="$80,000-150,000",
                timeline="12-16 weeks",
                required=True,
                gmp_glp_required=True
            ))
            order += 1
        
        # Consultation recommendation for high/very high
        if tier.value >= RiskTier.HIGH.value:
            plan.append(TestingStep(
                order=order,
                test_name="Toxicology Consultation",
                test_type="consultation",
                rationale="Expert toxicologist review recommended before proceeding",
                estimated_cost_range="$2,000-10,000",
                timeline="1-2 weeks",
                required=True,
                gmp_glp_required=False
            ))
        
        return plan
    
    def _check_consultation_needs(
        self,
        risk_summary: RiskSummary,
        alerts: List[StructuralAlert]
    ) -> Tuple[bool, List[str]]:
        """Determine if expert consultation is required."""
        flags = []
        
        if risk_summary.overall_tier.value >= RiskTier.HIGH.value:
            flags.append("High risk tier requires toxicologist review")
        
        critical_alerts = [a for a in alerts if a.severity == AlertSeverity.CRITICAL]
        if critical_alerts:
            flags.append(f"Critical structural alert(s): {', '.join(a.name for a in critical_alerts)}")
        
        if len(alerts) >= 5:
            flags.append("Multiple structural alerts warrant expert interpretation")
        
        return len(flags) > 0, flags
    
    def get_alert_codes(self) -> Dict[str, str]:
        """Get available alert codes and descriptions."""
        return self.ALERT_CODES.copy()
    
    def get_risk_tiers(self) -> List[Dict[str, str]]:
        """Get risk tier information."""
        return [
            {
                "tier": tier.value,
                "label": tier.label,
                "testing_depth": tier.testing_depth
            }
            for tier in RiskTier
        ]
    
    def get_routes(self) -> List[str]:
        """Get supported administration routes."""
        return [r.value for r in AdministrationRoute if r != AdministrationRoute.UNKNOWN]
