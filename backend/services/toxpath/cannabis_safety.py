"""
ToxPath Cannabis Module - Cannabis Safety Assessor

Schedule III abuse liability assessment, drug-drug interaction prediction,
and safety profile generation for cannabis formulations.

Trade Secret: Abuse liability algorithms, interaction matrices,
safety scoring algorithms derived from clinical evidence.

Regulatory Compliance: DEA Schedule III requirements, FDA safety standards.
"""

from enum import Enum
from typing import Dict, List, Optional, Tuple, Set
from pydantic import BaseModel, Field
from datetime import datetime
import uuid


class AbuseClass(str, Enum):
    """DEA abuse potential classification."""
    SCHEDULE_I = "Schedule I"       # High abuse, no medical use
    SCHEDULE_II = "Schedule II"     # High abuse, accepted medical use
    SCHEDULE_III = "Schedule III"   # Moderate abuse, accepted medical use
    SCHEDULE_IV = "Schedule IV"     # Low abuse
    SCHEDULE_V = "Schedule V"       # Lowest abuse
    UNSCHEDULED = "Unscheduled"     # No DEA scheduling


class RiskLevel(str, Enum):
    """Risk severity levels."""
    MINIMAL = "minimal"
    LOW = "low"
    MODERATE = "moderate"
    HIGH = "high"
    SEVERE = "severe"
    CONTRAINDICATED = "contraindicated"


class InteractionType(str, Enum):
    """Drug interaction types."""
    PHARMACOKINETIC = "pharmacokinetic"  # Affects drug metabolism
    PHARMACODYNAMIC = "pharmacodynamic"  # Affects drug action
    ADDITIVE = "additive"                # Combined effects
    SYNERGISTIC = "synergistic"          # Enhanced effects
    ANTAGONISTIC = "antagonistic"        # Reduced effects
    UNKNOWN = "unknown"


class InteractionMechanism(str, Enum):
    """Specific interaction mechanisms."""
    CYP3A4_INHIBITION = "CYP3A4 inhibition"
    CYP2C9_INHIBITION = "CYP2C9 inhibition"
    CYP2C19_INHIBITION = "CYP2C19 inhibition"
    CYP2D6_INHIBITION = "CYP2D6 inhibition"
    CYP1A2_INDUCTION = "CYP1A2 induction"
    P_GLYCOPROTEIN = "P-glycoprotein interaction"
    CNS_DEPRESSION = "CNS depression"
    BLEEDING_RISK = "bleeding risk"
    HYPOTENSION = "hypotension"
    SEDATION = "additive sedation"
    CARDIAC = "cardiac effects"


class VulnerablePopulation(str, Enum):
    """Populations requiring special consideration."""
    PEDIATRIC = "pediatric"
    GERIATRIC = "geriatric"
    PREGNANT = "pregnant"
    LACTATING = "lactating"
    HEPATIC_IMPAIRED = "hepatic_impaired"
    RENAL_IMPAIRED = "renal_impaired"
    CARDIAC_DISEASE = "cardiac_disease"
    PSYCHIATRIC = "psychiatric"
    SUBSTANCE_USE_DISORDER = "substance_use_disorder"


class AdverseEvent(BaseModel):
    """Adverse event information."""
    event_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    name: str
    severity: RiskLevel
    frequency: str  # "common", "uncommon", "rare", "very_rare"
    onset: str  # "immediate", "delayed", "chronic"
    reversible: bool = True
    management: str
    monitoring_required: bool = False


class DrugInteraction(BaseModel):
    """Drug-drug interaction details."""
    interaction_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    drug_name: str
    drug_class: str
    interaction_type: InteractionType
    mechanism: InteractionMechanism
    risk_level: RiskLevel
    clinical_significance: str
    recommendation: str
    evidence_quality: str  # "strong", "moderate", "limited"
    monitoring_parameters: List[str] = Field(default_factory=list)


class AbuseLiabilityProfile(BaseModel):
    """Abuse liability assessment."""
    profile_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    
    # Overall classification
    proposed_schedule: AbuseClass
    schedule_justification: str
    
    # Component scores (0-100 scale)
    dependence_potential: float = 0.0     # Physical/psychological dependence
    reinforcement_potential: float = 0.0  # Likelihood of repeated use
    euphoria_potential: float = 0.0       # Subjective "high"
    withdrawal_severity: float = 0.0      # Withdrawal syndrome severity
    tolerance_development: float = 0.0    # Rate of tolerance
    
    # Overall abuse liability (0-100)
    overall_abuse_liability: float = 0.0
    
    # Risk factors
    risk_factors: List[str] = Field(default_factory=list)
    protective_factors: List[str] = Field(default_factory=list)
    
    # DEA requirements
    dea_quota_required: bool = False
    controlled_substance_act_requirements: List[str] = Field(default_factory=list)


class SafetyProfile(BaseModel):
    """Complete safety profile for a formulation."""
    profile_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    formulation_name: str
    
    # Abuse liability
    abuse_liability: AbuseLiabilityProfile
    
    # Drug interactions
    drug_interactions: List[DrugInteraction] = Field(default_factory=list)
    high_risk_interactions: int = 0
    
    # Adverse events
    adverse_events: List[AdverseEvent] = Field(default_factory=list)
    
    # Population-specific risks
    population_risks: Dict[str, RiskLevel] = Field(default_factory=dict)
    
    # Contraindications
    absolute_contraindications: List[str] = Field(default_factory=list)
    relative_contraindications: List[str] = Field(default_factory=list)
    
    # Monitoring requirements
    monitoring_parameters: List[str] = Field(default_factory=list)
    monitoring_frequency: str = ""
    
    # Overall safety score (0-100, higher = safer)
    overall_safety_score: float = 0.0
    
    # Metadata
    analysis_timestamp: datetime = Field(default_factory=datetime.utcnow)
    confidence_score: float = 0.0


class SafetyRecommendation(BaseModel):
    """Safety-based recommendation."""
    recommendation_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    category: str  # "dosing", "monitoring", "contraindication", "interaction"
    priority: str  # "critical", "important", "advisory"
    recommendation: str
    rationale: str
    supporting_evidence: List[str] = Field(default_factory=list)


class CannabisSafetyAssessor:
    """
    Cannabis safety assessment for Schedule III compliance.
    
    Trade Secret: Abuse liability algorithms, interaction matrices,
    safety scoring derived from clinical evidence and traditional knowledge.
    """
    
    # ==========================================================================
    # TRADE SECRET: Cannabinoid Abuse Liability Coefficients
    # Values represent abuse potential (0-1 scale)
    # ==========================================================================
    _CANNABINOID_ABUSE_COEFFICIENTS: Dict[str, Dict[str, float]] = {
        "THC": {
            "dependence": 0.30,
            "reinforcement": 0.55,
            "euphoria": 0.70,
            "withdrawal": 0.25,
            "tolerance": 0.45,
        },
        "delta-9-THC": {
            "dependence": 0.30,
            "reinforcement": 0.55,
            "euphoria": 0.70,
            "withdrawal": 0.25,
            "tolerance": 0.45,
        },
        "CBD": {
            "dependence": 0.02,
            "reinforcement": 0.05,
            "euphoria": 0.03,
            "withdrawal": 0.01,
            "tolerance": 0.10,
        },
        "CBG": {
            "dependence": 0.03,
            "reinforcement": 0.05,
            "euphoria": 0.02,
            "withdrawal": 0.01,
            "tolerance": 0.08,
        },
        "CBN": {
            "dependence": 0.15,
            "reinforcement": 0.20,
            "euphoria": 0.15,
            "withdrawal": 0.10,
            "tolerance": 0.25,
        },
        "delta-8-THC": {
            "dependence": 0.25,
            "reinforcement": 0.45,
            "euphoria": 0.55,
            "withdrawal": 0.20,
            "tolerance": 0.40,
        },
        "THCV": {
            "dependence": 0.15,
            "reinforcement": 0.25,
            "euphoria": 0.35,
            "withdrawal": 0.10,
            "tolerance": 0.30,
        },
    }
    
    # ==========================================================================
    # TRADE SECRET: Drug Interaction Database
    # ==========================================================================
    _DRUG_INTERACTIONS: Dict[str, Dict] = {
        "warfarin": {
            "class": "Anticoagulant",
            "interaction_type": InteractionType.PHARMACOKINETIC,
            "mechanism": InteractionMechanism.CYP2C9_INHIBITION,
            "risk_level": RiskLevel.HIGH,
            "significance": "CBD inhibits CYP2C9, increasing warfarin levels and bleeding risk",
            "recommendation": "Monitor INR closely, consider 25-50% warfarin dose reduction",
            "evidence": "strong",
            "monitoring": ["INR", "bleeding signs"],
            "affected_by": ["CBD", "CBG"],
        },
        "clobazam": {
            "class": "Benzodiazepine anticonvulsant",
            "interaction_type": InteractionType.PHARMACOKINETIC,
            "mechanism": InteractionMechanism.CYP2C19_INHIBITION,
            "risk_level": RiskLevel.HIGH,
            "significance": "CBD inhibits CYP2C19, increasing N-desmethylclobazam levels 3-5x",
            "recommendation": "Reduce clobazam dose by 50%, monitor for sedation",
            "evidence": "strong",
            "monitoring": ["sedation", "clobazam levels"],
            "affected_by": ["CBD"],
        },
        "opioids": {
            "class": "Opioid analgesic",
            "interaction_type": InteractionType.PHARMACODYNAMIC,
            "mechanism": InteractionMechanism.CNS_DEPRESSION,
            "risk_level": RiskLevel.MODERATE,
            "significance": "Additive CNS depression, but may allow opioid dose reduction",
            "recommendation": "Monitor for excessive sedation, consider opioid dose adjustment",
            "evidence": "moderate",
            "monitoring": ["respiratory rate", "sedation level"],
            "affected_by": ["THC", "CBN"],
        },
        "benzodiazepines": {
            "class": "Anxiolytic/sedative",
            "interaction_type": InteractionType.PHARMACODYNAMIC,
            "mechanism": InteractionMechanism.SEDATION,
            "risk_level": RiskLevel.MODERATE,
            "significance": "Additive sedation and CNS depression",
            "recommendation": "Avoid concurrent high doses, monitor sedation",
            "evidence": "moderate",
            "monitoring": ["sedation", "respiratory status"],
            "affected_by": ["THC", "CBN", "CBD"],
        },
        "ssri_antidepressants": {
            "class": "Antidepressant",
            "interaction_type": InteractionType.PHARMACOKINETIC,
            "mechanism": InteractionMechanism.CYP2D6_INHIBITION,
            "risk_level": RiskLevel.LOW,
            "significance": "CBD may inhibit CYP2D6 metabolism of some SSRIs",
            "recommendation": "Monitor for SSRI side effects",
            "evidence": "limited",
            "monitoring": ["serotonin syndrome signs", "SSRI side effects"],
            "affected_by": ["CBD"],
        },
        "immunosuppressants": {
            "class": "Immunosuppressant",
            "interaction_type": InteractionType.PHARMACOKINETIC,
            "mechanism": InteractionMechanism.CYP3A4_INHIBITION,
            "risk_level": RiskLevel.HIGH,
            "significance": "CBD inhibits CYP3A4, may increase tacrolimus/cyclosporine levels",
            "recommendation": "Monitor drug levels closely, consider dose adjustment",
            "evidence": "moderate",
            "monitoring": ["drug levels", "toxicity signs"],
            "affected_by": ["CBD"],
        },
        "antihypertensives": {
            "class": "Cardiovascular",
            "interaction_type": InteractionType.PHARMACODYNAMIC,
            "mechanism": InteractionMechanism.HYPOTENSION,
            "risk_level": RiskLevel.LOW,
            "significance": "Additive hypotensive effects",
            "recommendation": "Monitor blood pressure, especially at initiation",
            "evidence": "moderate",
            "monitoring": ["blood pressure", "orthostatic symptoms"],
            "affected_by": ["THC", "CBD"],
        },
        "antiplatelet_agents": {
            "class": "Antiplatelet",
            "interaction_type": InteractionType.PHARMACODYNAMIC,
            "mechanism": InteractionMechanism.BLEEDING_RISK,
            "risk_level": RiskLevel.MODERATE,
            "significance": "Increased bleeding risk with CBD",
            "recommendation": "Monitor for bleeding, use caution with procedures",
            "evidence": "limited",
            "monitoring": ["bleeding signs", "platelet function"],
            "affected_by": ["CBD"],
        },
    }
    
    # ==========================================================================
    # TRADE SECRET: Adverse Event Database
    # ==========================================================================
    _ADVERSE_EVENTS: Dict[str, Dict] = {
        "somnolence": {
            "severity": RiskLevel.LOW,
            "frequency": "common",
            "onset": "immediate",
            "reversible": True,
            "management": "Dose reduction, avoid evening doses",
            "monitoring": False,
            "associated_compounds": ["THC", "CBN", "CBD"],
        },
        "dizziness": {
            "severity": RiskLevel.LOW,
            "frequency": "common",
            "onset": "immediate",
            "reversible": True,
            "management": "Slow titration, avoid sudden position changes",
            "monitoring": False,
            "associated_compounds": ["THC", "CBD"],
        },
        "elevated_liver_enzymes": {
            "severity": RiskLevel.MODERATE,
            "frequency": "uncommon",
            "onset": "delayed",
            "reversible": True,
            "management": "Monitor LFTs, dose reduction or discontinuation",
            "monitoring": True,
            "associated_compounds": ["CBD"],
        },
        "tachycardia": {
            "severity": RiskLevel.MODERATE,
            "frequency": "common",
            "onset": "immediate",
            "reversible": True,
            "management": "Dose reduction, avoid in cardiac patients",
            "monitoring": True,
            "associated_compounds": ["THC"],
        },
        "anxiety_paranoia": {
            "severity": RiskLevel.MODERATE,
            "frequency": "uncommon",
            "onset": "immediate",
            "reversible": True,
            "management": "Reduce THC dose, increase CBD ratio",
            "monitoring": False,
            "associated_compounds": ["THC"],
        },
        "dry_mouth": {
            "severity": RiskLevel.MINIMAL,
            "frequency": "common",
            "onset": "immediate",
            "reversible": True,
            "management": "Adequate hydration, artificial saliva",
            "monitoring": False,
            "associated_compounds": ["THC", "CBD"],
        },
        "appetite_changes": {
            "severity": RiskLevel.MINIMAL,
            "frequency": "common",
            "onset": "immediate",
            "reversible": True,
            "management": "Monitor weight, dietary counseling",
            "monitoring": False,
            "associated_compounds": ["THC", "THCV"],
        },
        "cannabinoid_hyperemesis": {
            "severity": RiskLevel.HIGH,
            "frequency": "rare",
            "onset": "chronic",
            "reversible": True,
            "management": "Discontinuation, hot water bathing for acute relief",
            "monitoring": True,
            "associated_compounds": ["THC"],
        },
    }
    
    # ==========================================================================
    # TRADE SECRET: Population Risk Modifiers
    # ==========================================================================
    _POPULATION_RISK_MODIFIERS: Dict[str, Dict] = {
        VulnerablePopulation.PEDIATRIC.value: {
            "risk_modifier": 1.5,
            "concerns": ["neurodevelopment", "dosing uncertainty"],
            "thc_max_recommendation": 0.3,  # % THC
            "special_monitoring": ["developmental milestones", "cognition"],
        },
        VulnerablePopulation.GERIATRIC.value: {
            "risk_modifier": 1.3,
            "concerns": ["falls", "cognitive impairment", "polypharmacy"],
            "starting_dose_modifier": 0.5,
            "special_monitoring": ["fall risk", "cognition", "drug interactions"],
        },
        VulnerablePopulation.PREGNANT.value: {
            "risk_modifier": 2.0,
            "concerns": ["fetal neurodevelopment", "low birth weight"],
            "thc_max_recommendation": 0.0,  # Avoid
            "special_monitoring": ["fetal growth", "developmental outcomes"],
        },
        VulnerablePopulation.LACTATING.value: {
            "risk_modifier": 1.8,
            "concerns": ["infant exposure", "neurodevelopment"],
            "thc_max_recommendation": 0.0,  # Avoid
            "special_monitoring": ["infant development"],
        },
        VulnerablePopulation.HEPATIC_IMPAIRED.value: {
            "risk_modifier": 1.5,
            "concerns": ["reduced metabolism", "increased toxicity"],
            "starting_dose_modifier": 0.5,
            "special_monitoring": ["LFTs", "drug levels"],
        },
        VulnerablePopulation.CARDIAC_DISEASE.value: {
            "risk_modifier": 1.4,
            "concerns": ["tachycardia", "orthostatic hypotension"],
            "thc_max_recommendation": 5.0,
            "special_monitoring": ["heart rate", "blood pressure", "ECG"],
        },
        VulnerablePopulation.PSYCHIATRIC.value: {
            "risk_modifier": 1.6,
            "concerns": ["psychosis risk", "anxiety", "depression"],
            "thc_max_recommendation": 5.0,
            "special_monitoring": ["psychiatric symptoms", "psychosis screening"],
        },
        VulnerablePopulation.SUBSTANCE_USE_DISORDER.value: {
            "risk_modifier": 2.0,
            "concerns": ["addiction", "abuse potential", "relapse"],
            "special_monitoring": ["substance use patterns", "craving"],
        },
    }
    
    def __init__(self):
        self._assessment_cache: Dict[str, SafetyProfile] = {}
    
    def assess_abuse_liability(
        self,
        cannabinoids: Dict[str, float],
        formulation_type: str = "oral",
    ) -> AbuseLiabilityProfile:
        """
        Assess abuse liability for DEA Schedule III compliance.
        
        Trade Secret: Abuse liability scoring algorithm derived from
        clinical data and historical scheduling decisions.
        
        Args:
            cannabinoids: Dict of cannabinoid concentrations (%)
            formulation_type: Delivery method affecting abuse potential
            
        Returns:
            AbuseLiabilityProfile with scheduling recommendation
        """
        profile = AbuseLiabilityProfile(
            proposed_schedule=AbuseClass.UNSCHEDULED,
            schedule_justification="",
        )
        
        # Calculate component scores
        total_weight = 0.0
        
        for cannabinoid, concentration in cannabinoids.items():
            if cannabinoid in self._CANNABINOID_ABUSE_COEFFICIENTS:
                coefficients = self._CANNABINOID_ABUSE_COEFFICIENTS[cannabinoid]
                weight = min(concentration / 30.0, 1.0)  # Normalize to 30%
                total_weight += weight
                
                profile.dependence_potential += coefficients["dependence"] * weight * 100
                profile.reinforcement_potential += coefficients["reinforcement"] * weight * 100
                profile.euphoria_potential += coefficients["euphoria"] * weight * 100
                profile.withdrawal_severity += coefficients["withdrawal"] * weight * 100
                profile.tolerance_development += coefficients["tolerance"] * weight * 100
        
        # Normalize scores
        if total_weight > 0:
            profile.dependence_potential /= total_weight
            profile.reinforcement_potential /= total_weight
            profile.euphoria_potential /= total_weight
            profile.withdrawal_severity /= total_weight
            profile.tolerance_development /= total_weight
        
        # Apply formulation modifier
        formulation_modifiers = {
            "oral": 1.0,
            "sublingual": 1.1,
            "inhalation": 1.4,  # Faster onset = higher abuse potential
            "topical": 0.5,
            "transdermal": 0.7,
            "suppository": 0.8,
        }
        modifier = formulation_modifiers.get(formulation_type.lower(), 1.0)
        
        # Calculate overall abuse liability
        profile.overall_abuse_liability = (
            profile.dependence_potential * 0.20 +
            profile.reinforcement_potential * 0.30 +
            profile.euphoria_potential * 0.25 +
            profile.withdrawal_severity * 0.10 +
            profile.tolerance_development * 0.15
        ) * modifier
        
        # Determine proposed schedule
        profile = self._determine_schedule(profile, cannabinoids)
        
        # Identify risk and protective factors
        profile.risk_factors = self._identify_abuse_risk_factors(cannabinoids, formulation_type)
        profile.protective_factors = self._identify_protective_factors(cannabinoids)
        
        # DEA requirements
        if profile.proposed_schedule in [AbuseClass.SCHEDULE_II, AbuseClass.SCHEDULE_III]:
            profile.dea_quota_required = True
            profile.controlled_substance_act_requirements = [
                "DEA registration required",
                "Prescription limitations apply",
                "Record-keeping requirements",
                "Security requirements",
            ]
            if profile.proposed_schedule == AbuseClass.SCHEDULE_III:
                profile.controlled_substance_act_requirements.extend([
                    "Up to 6 refills within 6 months permitted",
                    "Telephone prescriptions permitted",
                ])
        
        return profile
    
    def _determine_schedule(
        self,
        profile: AbuseLiabilityProfile,
        cannabinoids: Dict[str, float],
    ) -> AbuseLiabilityProfile:
        """Determine appropriate DEA schedule based on abuse liability."""
        
        thc_content = cannabinoids.get("THC", cannabinoids.get("delta-9-THC", 0))
        cbd_content = cannabinoids.get("CBD", 0)
        
        # CBD-dominant products (CBD:THC > 2:1)
        if cbd_content > 0 and thc_content > 0:
            cbd_thc_ratio = cbd_content / thc_content
            if cbd_thc_ratio >= 2.0 and thc_content <= 5.0:
                profile.proposed_schedule = AbuseClass.UNSCHEDULED
                profile.schedule_justification = (
                    f"CBD-dominant formulation (CBD:THC ratio {cbd_thc_ratio:.1f}:1) "
                    f"with low THC content ({thc_content:.1f}%) demonstrates minimal abuse potential"
                )
                return profile
        
        # Very low THC content
        if thc_content < 0.3:
            profile.proposed_schedule = AbuseClass.UNSCHEDULED
            profile.schedule_justification = (
                f"THC content ({thc_content:.2f}%) below 0.3% threshold, "
                "qualifies as hemp under 2018 Farm Bill"
            )
            return profile
        
        # Schedule determination based on overall abuse liability
        if profile.overall_abuse_liability < 20:
            profile.proposed_schedule = AbuseClass.SCHEDULE_V
            profile.schedule_justification = (
                f"Very low abuse liability ({profile.overall_abuse_liability:.1f}%) "
                "with demonstrated medical utility"
            )
        elif profile.overall_abuse_liability < 35:
            profile.proposed_schedule = AbuseClass.SCHEDULE_IV
            profile.schedule_justification = (
                f"Low abuse liability ({profile.overall_abuse_liability:.1f}%) "
                "relative to Schedule III substances"
            )
        elif profile.overall_abuse_liability < 55:
            profile.proposed_schedule = AbuseClass.SCHEDULE_III
            profile.schedule_justification = (
                f"Moderate abuse liability ({profile.overall_abuse_liability:.1f}%) "
                "consistent with Schedule III classification. "
                "Comparable to approved dronabinol (Marinol®)"
            )
        elif profile.overall_abuse_liability < 70:
            profile.proposed_schedule = AbuseClass.SCHEDULE_II
            profile.schedule_justification = (
                f"Significant abuse liability ({profile.overall_abuse_liability:.1f}%) "
                "requires Schedule II classification"
            )
        else:
            profile.proposed_schedule = AbuseClass.SCHEDULE_I
            profile.schedule_justification = (
                f"High abuse liability ({profile.overall_abuse_liability:.1f}%) "
                "with characteristics of Schedule I substances"
            )
        
        return profile
    
    def _identify_abuse_risk_factors(
        self,
        cannabinoids: Dict[str, float],
        formulation_type: str,
    ) -> List[str]:
        """Identify factors that increase abuse potential."""
        risk_factors = []
        
        thc = cannabinoids.get("THC", cannabinoids.get("delta-9-THC", 0))
        
        if thc > 20:
            risk_factors.append(f"High THC concentration ({thc:.1f}%)")
        
        if formulation_type.lower() == "inhalation":
            risk_factors.append("Inhalation route provides rapid onset (higher reinforcement)")
        
        if "delta-8-THC" in cannabinoids and cannabinoids["delta-8-THC"] > 5:
            risk_factors.append("Significant delta-8-THC content")
        
        if "THCP" in cannabinoids:
            risk_factors.append("THCP present (33x THC potency at CB1)")
        
        cbd = cannabinoids.get("CBD", 0)
        if cbd < thc * 0.1 and thc > 5:
            risk_factors.append("Minimal CBD to modulate THC effects")
        
        return risk_factors
    
    def _identify_protective_factors(
        self,
        cannabinoids: Dict[str, float],
    ) -> List[str]:
        """Identify factors that reduce abuse potential."""
        protective_factors = []
        
        thc = cannabinoids.get("THC", cannabinoids.get("delta-9-THC", 0))
        cbd = cannabinoids.get("CBD", 0)
        
        if cbd > thc:
            protective_factors.append(f"CBD-dominant ratio ({cbd/thc:.1f}:1)")
        
        if cbd > 5:
            protective_factors.append("CBD modulates THC psychoactivity")
        
        if "CBG" in cannabinoids and cannabinoids["CBG"] > 1:
            protective_factors.append("CBG provides non-intoxicating therapeutic effects")
        
        if thc < 10:
            protective_factors.append("Moderate THC concentration")
        
        if "THCV" in cannabinoids and cannabinoids["THCV"] > 1:
            protective_factors.append("THCV acts as CB1 antagonist at low doses")
        
        return protective_factors
    
    def assess_drug_interactions(
        self,
        cannabinoids: Dict[str, float],
        concomitant_medications: List[str],
    ) -> List[DrugInteraction]:
        """
        Assess drug-drug interactions with concomitant medications.
        
        Trade Secret: Interaction prediction algorithm and severity scoring.
        
        Args:
            cannabinoids: Dict of cannabinoid concentrations
            concomitant_medications: List of medication names/classes
            
        Returns:
            List of DrugInteraction objects
        """
        interactions = []
        
        for medication in concomitant_medications:
            medication_lower = medication.lower()
            
            for drug_key, interaction_data in self._DRUG_INTERACTIONS.items():
                if drug_key in medication_lower or medication_lower in drug_key:
                    # Check if relevant cannabinoids are present
                    relevant_cannabinoids = [
                        c for c in interaction_data["affected_by"]
                        if c in cannabinoids or c.replace("-", "_") in cannabinoids
                    ]
                    
                    if relevant_cannabinoids:
                        interaction = DrugInteraction(
                            drug_name=medication,
                            drug_class=interaction_data["class"],
                            interaction_type=interaction_data["interaction_type"],
                            mechanism=interaction_data["mechanism"],
                            risk_level=interaction_data["risk_level"],
                            clinical_significance=interaction_data["significance"],
                            recommendation=interaction_data["recommendation"],
                            evidence_quality=interaction_data["evidence"],
                            monitoring_parameters=interaction_data["monitoring"],
                        )
                        interactions.append(interaction)
        
        return interactions
    
    def assess_adverse_events(
        self,
        cannabinoids: Dict[str, float],
    ) -> List[AdverseEvent]:
        """
        Predict likely adverse events based on formulation.
        
        Trade Secret: Adverse event prediction algorithm.
        """
        events = []
        
        for event_name, event_data in self._ADVERSE_EVENTS.items():
            # Check if relevant cannabinoids are present
            relevant = any(
                c in cannabinoids or c.replace("-", "_") in cannabinoids
                for c in event_data["associated_compounds"]
            )
            
            if relevant:
                event = AdverseEvent(
                    name=event_name,
                    severity=event_data["severity"],
                    frequency=event_data["frequency"],
                    onset=event_data["onset"],
                    reversible=event_data["reversible"],
                    management=event_data["management"],
                    monitoring_required=event_data["monitoring"],
                )
                events.append(event)
        
        return events
    
    def assess_population_risks(
        self,
        cannabinoids: Dict[str, float],
        populations: List[str],
    ) -> Dict[str, RiskLevel]:
        """
        Assess risks for specific populations.
        
        Trade Secret: Population risk algorithm.
        """
        risks = {}
        
        thc = cannabinoids.get("THC", cannabinoids.get("delta-9-THC", 0))
        
        for population in populations:
            if population in self._POPULATION_RISK_MODIFIERS:
                modifier_data = self._POPULATION_RISK_MODIFIERS[population]
                
                # Base risk on THC content and modifier
                base_risk = thc / 30.0  # Normalize to 30%
                modified_risk = base_risk * modifier_data["risk_modifier"]
                
                # Check THC max recommendation
                if "thc_max_recommendation" in modifier_data:
                    max_thc = modifier_data["thc_max_recommendation"]
                    if thc > max_thc:
                        modified_risk += 0.3
                
                # Classify risk level
                if modified_risk >= 0.8:
                    risks[population] = RiskLevel.CONTRAINDICATED
                elif modified_risk >= 0.6:
                    risks[population] = RiskLevel.SEVERE
                elif modified_risk >= 0.4:
                    risks[population] = RiskLevel.HIGH
                elif modified_risk >= 0.25:
                    risks[population] = RiskLevel.MODERATE
                elif modified_risk >= 0.1:
                    risks[population] = RiskLevel.LOW
                else:
                    risks[population] = RiskLevel.MINIMAL
        
        return risks
    
    def generate_safety_profile(
        self,
        cannabinoids: Dict[str, float],
        formulation_name: str,
        formulation_type: str = "oral",
        concomitant_medications: Optional[List[str]] = None,
        populations: Optional[List[str]] = None,
    ) -> SafetyProfile:
        """
        Generate comprehensive safety profile for a formulation.
        
        Args:
            cannabinoids: Dict of cannabinoid concentrations
            formulation_name: Name of the formulation
            formulation_type: Delivery method
            concomitant_medications: List of medications to check interactions
            populations: List of populations to assess
            
        Returns:
            Complete SafetyProfile
        """
        profile = SafetyProfile(
            formulation_name=formulation_name,
            abuse_liability=self.assess_abuse_liability(cannabinoids, formulation_type),
        )
        
        # Drug interactions
        if concomitant_medications:
            profile.drug_interactions = self.assess_drug_interactions(
                cannabinoids, concomitant_medications
            )
            profile.high_risk_interactions = sum(
                1 for i in profile.drug_interactions
                if i.risk_level in [RiskLevel.HIGH, RiskLevel.SEVERE]
            )
        
        # Adverse events
        profile.adverse_events = self.assess_adverse_events(cannabinoids)
        
        # Population risks
        if populations:
            profile.population_risks = self.assess_population_risks(
                cannabinoids, populations
            )
        
        # Contraindications
        profile.absolute_contraindications = self._get_contraindications(
            cannabinoids, "absolute"
        )
        profile.relative_contraindications = self._get_contraindications(
            cannabinoids, "relative"
        )
        
        # Monitoring requirements
        profile.monitoring_parameters = self._get_monitoring_requirements(
            cannabinoids, profile.drug_interactions
        )
        profile.monitoring_frequency = self._determine_monitoring_frequency(profile)
        
        # Calculate overall safety score
        profile.overall_safety_score = self._calculate_safety_score(profile)
        
        # Confidence
        profile.confidence_score = 0.75  # Base confidence for safety assessment
        
        # Cache
        self._assessment_cache[profile.profile_id] = profile
        
        return profile
    
    def _get_contraindications(
        self,
        cannabinoids: Dict[str, float],
        contraindication_type: str,
    ) -> List[str]:
        """Get contraindications based on formulation."""
        absolute = [
            "History of cannabis-induced psychosis",
            "Hypersensitivity to cannabis or cannabinoids",
        ]
        
        relative = [
            "Active psychiatric disorder",
            "Cardiovascular disease",
            "Hepatic impairment",
            "History of substance use disorder",
            "Concurrent use of CNS depressants",
        ]
        
        thc = cannabinoids.get("THC", cannabinoids.get("delta-9-THC", 0))
        
        if thc > 10:
            relative.append("Occupations requiring alertness")
        
        if contraindication_type == "absolute":
            return absolute
        return relative
    
    def _get_monitoring_requirements(
        self,
        cannabinoids: Dict[str, float],
        interactions: List[DrugInteraction],
    ) -> List[str]:
        """Determine monitoring requirements."""
        monitoring = set()
        
        cbd = cannabinoids.get("CBD", 0)
        if cbd > 5:
            monitoring.add("Liver function tests (baseline and periodic)")
        
        thc = cannabinoids.get("THC", cannabinoids.get("delta-9-THC", 0))
        if thc > 5:
            monitoring.add("Psychiatric symptom screening")
            monitoring.add("Heart rate monitoring")
        
        for interaction in interactions:
            monitoring.update(interaction.monitoring_parameters)
        
        return list(monitoring)
    
    def _determine_monitoring_frequency(self, profile: SafetyProfile) -> str:
        """Determine appropriate monitoring frequency."""
        if profile.high_risk_interactions > 0:
            return "Weekly for first month, then monthly"
        if len(profile.monitoring_parameters) > 3:
            return "Every 2 weeks for first month, then monthly"
        return "Monthly for first 3 months, then quarterly"
    
    def _calculate_safety_score(self, profile: SafetyProfile) -> float:
        """
        Calculate overall safety score (0-100, higher = safer).
        Trade Secret: Safety scoring algorithm.
        """
        score = 100.0
        
        # Deduct for abuse liability
        score -= profile.abuse_liability.overall_abuse_liability * 0.3
        
        # Deduct for interactions
        for interaction in profile.drug_interactions:
            if interaction.risk_level == RiskLevel.SEVERE:
                score -= 15
            elif interaction.risk_level == RiskLevel.HIGH:
                score -= 10
            elif interaction.risk_level == RiskLevel.MODERATE:
                score -= 5
            else:
                score -= 2
        
        # Deduct for severe adverse events
        for event in profile.adverse_events:
            if event.severity == RiskLevel.SEVERE:
                score -= 10
            elif event.severity == RiskLevel.HIGH:
                score -= 5
            elif event.severity == RiskLevel.MODERATE:
                score -= 3
        
        # Deduct for contraindications
        score -= len(profile.absolute_contraindications) * 5
        score -= len(profile.relative_contraindications) * 2
        
        return max(0, min(100, score))
    
    def get_safety_recommendations(
        self,
        profile: SafetyProfile,
    ) -> List[SafetyRecommendation]:
        """Generate safety recommendations based on profile."""
        recommendations = []
        
        # Abuse liability recommendations
        if profile.abuse_liability.overall_abuse_liability > 40:
            recommendations.append(SafetyRecommendation(
                category="abuse_prevention",
                priority="important",
                recommendation="Implement prescriber monitoring program",
                rationale=f"Abuse liability score of {profile.abuse_liability.overall_abuse_liability:.1f}% warrants enhanced monitoring",
            ))
        
        # Interaction recommendations
        for interaction in profile.drug_interactions:
            if interaction.risk_level in [RiskLevel.HIGH, RiskLevel.SEVERE]:
                recommendations.append(SafetyRecommendation(
                    category="interaction",
                    priority="critical",
                    recommendation=interaction.recommendation,
                    rationale=interaction.clinical_significance,
                ))
        
        # Population-specific recommendations
        for population, risk in profile.population_risks.items():
            if risk == RiskLevel.CONTRAINDICATED:
                recommendations.append(SafetyRecommendation(
                    category="contraindication",
                    priority="critical",
                    recommendation=f"Contraindicated in {population} population",
                    rationale="Risk outweighs potential benefit",
                ))
        
        return recommendations
    
    def get_assessment(self, profile_id: str) -> Optional[SafetyProfile]:
        """Retrieve cached assessment by ID."""
        return self._assessment_cache.get(profile_id)
