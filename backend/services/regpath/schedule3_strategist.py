"""
RegPath Cannabis Module - Schedule III Strategist

DEA registration planning, FDA botanical drug pathway design,
and regulatory compliance strategy for cannabis products.

Trade Secret: Regulatory pathway optimization algorithms,
DEA quota calculation, FDA submission timeline generators.

Regulatory Compliance: DEA 21 CFR 1301, FDA 21 CFR 312/314.
"""

from enum import Enum
from typing import Dict, List, Optional, Tuple, Set
from pydantic import BaseModel, Field
from datetime import datetime, date, timedelta
import uuid


class RegulatoryPathway(str, Enum):
    """FDA regulatory pathways for cannabis products."""
    NDA_505B1 = "NDA 505(b)(1)"           # Full NDA
    NDA_505B2 = "NDA 505(b)(2)"           # Abbreviated NDA
    ANDA = "ANDA"                          # Generic drug
    BLA = "BLA"                            # Biologic
    BOTANICAL_DRUG_DEVELOPMENT = "Botanical IND/NDA"  # Section 201(ff)
    OTC_MONOGRAPH = "OTC Monograph"        # Over-the-counter
    DIETARY_SUPPLEMENT = "DSHEA"           # Dietary supplement
    COSMETIC = "Cosmetic (FD&C Act)"       # Cosmetic


class DEASchedule(str, Enum):
    """DEA scheduling classifications."""
    SCHEDULE_I = "Schedule I"
    SCHEDULE_II = "Schedule II"
    SCHEDULE_III = "Schedule III"
    SCHEDULE_IV = "Schedule IV"
    SCHEDULE_V = "Schedule V"
    UNSCHEDULED = "Unscheduled"


class DEALicenseType(str, Enum):
    """DEA registration categories."""
    MANUFACTURER = "Manufacturer"
    DISTRIBUTOR = "Distributor"
    RESEARCHER = "Researcher"
    PHARMACY = "Pharmacy"
    HOSPITAL = "Hospital"
    PRACTITIONER = "Practitioner"
    IMPORTER = "Importer"
    EXPORTER = "Exporter"


class SubmissionType(str, Enum):
    """FDA submission types."""
    PRE_IND = "Pre-IND Meeting Request"
    IND = "IND Application"
    IND_AMENDMENT = "IND Amendment"
    IND_SAFETY_REPORT = "IND Safety Report"
    END_OF_PHASE_1 = "End-of-Phase 1 Meeting"
    END_OF_PHASE_2 = "End-of-Phase 2 Meeting"
    PRE_NDA = "Pre-NDA Meeting"
    NDA = "NDA Application"
    NDA_AMENDMENT = "NDA Amendment"
    BPD_GUIDANCE = "Botanical Product Development Guidance"


class MilestoneStatus(str, Enum):
    """Regulatory milestone status."""
    NOT_STARTED = "not_started"
    IN_PROGRESS = "in_progress"
    COMPLETED = "completed"
    DELAYED = "delayed"
    BLOCKED = "blocked"


class DEARequirement(BaseModel):
    """DEA registration requirement."""
    requirement_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    category: str
    requirement: str
    cfr_reference: str
    timeline: str
    complexity: str  # "simple", "moderate", "complex"
    estimated_cost: Tuple[int, int]  # Min-max range
    documentation_needed: List[str] = Field(default_factory=list)


class FDASubmissionMilestone(BaseModel):
    """FDA submission milestone."""
    milestone_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    name: str
    submission_type: SubmissionType
    target_date: Optional[date] = None
    dependencies: List[str] = Field(default_factory=list)
    status: MilestoneStatus = MilestoneStatus.NOT_STARTED
    estimated_duration_days: int = 30
    required_documents: List[str] = Field(default_factory=list)
    estimated_cost: Tuple[int, int] = (0, 0)


class ClinicalTrialDesign(BaseModel):
    """Clinical trial design specification."""
    trial_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    phase: str  # "Phase 1", "Phase 2", "Phase 3"
    design_type: str  # "randomized", "double-blind", "placebo-controlled"
    primary_endpoint: str
    secondary_endpoints: List[str] = Field(default_factory=list)
    sample_size: int
    duration_weeks: int
    estimated_cost: Tuple[int, int]
    regulatory_considerations: List[str] = Field(default_factory=list)


class BotanicalDrugRequirement(BaseModel):
    """Botanical drug-specific requirement."""
    requirement_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    category: str  # "CMC", "Nonclinical", "Clinical", "Quality"
    requirement: str
    guidance_reference: str
    rationale: str
    special_considerations: List[str] = Field(default_factory=list)


class DEARegistrationPlan(BaseModel):
    """Complete DEA registration plan."""
    plan_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    proposed_schedule: DEASchedule
    license_types_needed: List[DEALicenseType] = Field(default_factory=list)
    
    # Requirements
    requirements: List[DEARequirement] = Field(default_factory=list)
    
    # Timeline
    total_timeline_months: int = 0
    critical_path_items: List[str] = Field(default_factory=list)
    
    # Costs
    estimated_total_cost: Tuple[int, int] = (0, 0)
    annual_compliance_cost: Tuple[int, int] = (0, 0)
    
    # Quota
    requires_quota: bool = False
    quota_application_deadline: Optional[date] = None
    estimated_annual_quota_kg: float = 0.0
    
    # Security
    security_requirements: List[str] = Field(default_factory=list)
    
    # Metadata
    generated_date: datetime = Field(default_factory=datetime.utcnow)


class BotanicalDrugPathway(BaseModel):
    """FDA botanical drug development pathway."""
    pathway_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    pathway_type: RegulatoryPathway = RegulatoryPathway.BOTANICAL_DRUG_DEVELOPMENT
    
    # Product characterization
    product_description: str = ""
    botanical_source: str = ""
    active_constituents: List[str] = Field(default_factory=list)
    standardization_markers: List[str] = Field(default_factory=list)
    
    # Milestones
    milestones: List[FDASubmissionMilestone] = Field(default_factory=list)
    current_phase: str = "Pre-clinical"
    
    # Clinical trials
    clinical_trials: List[ClinicalTrialDesign] = Field(default_factory=list)
    
    # Botanical-specific
    botanical_requirements: List[BotanicalDrugRequirement] = Field(default_factory=list)
    
    # Timeline & Cost
    total_development_years: float = 0.0
    estimated_total_cost_millions: Tuple[float, float] = (0.0, 0.0)
    probability_of_success: float = 0.0
    
    # Risk factors
    risk_factors: List[str] = Field(default_factory=list)
    mitigation_strategies: List[str] = Field(default_factory=list)
    
    # Metadata
    generated_date: datetime = Field(default_factory=datetime.utcnow)


class Schedule3Strategist:
    """
    Strategic planning for Schedule III cannabis products.
    
    Trade Secret: Regulatory pathway optimization, timeline prediction,
    cost estimation algorithms for DEA and FDA compliance.
    """
    
    # ==========================================================================
    # TRADE SECRET: DEA Requirements Database
    # ==========================================================================
    _DEA_REQUIREMENTS: Dict[str, Dict] = {
        "registration": {
            "category": "Registration",
            "requirement": "DEA Form 225 for manufacturers/distributors",
            "cfr_reference": "21 CFR 1301.13",
            "timeline": "3-6 months",
            "complexity": "moderate",
            "cost": (5000, 15000),
            "documents": [
                "DEA Form 225",
                "State license documentation",
                "Proof of business registration",
                "Background check authorization",
                "Physical security assessment",
            ],
        },
        "security": {
            "category": "Physical Security",
            "requirement": "Schedule III security requirements",
            "cfr_reference": "21 CFR 1301.71-76",
            "timeline": "1-3 months",
            "complexity": "complex",
            "cost": (50000, 200000),
            "documents": [
                "Security plan",
                "Vault specifications",
                "Alarm system documentation",
                "Access control procedures",
                "Employee background checks",
            ],
        },
        "record_keeping": {
            "category": "Record Keeping",
            "requirement": "Controlled substance inventory and records",
            "cfr_reference": "21 CFR 1304.03-04",
            "timeline": "Ongoing",
            "complexity": "moderate",
            "cost": (10000, 30000),
            "documents": [
                "Inventory management SOP",
                "Electronic record keeping system",
                "Biennial inventory procedures",
                "DEA Form 222/CSOS documentation",
            ],
        },
        "reporting": {
            "category": "Reporting",
            "requirement": "ARCOS reporting system registration",
            "cfr_reference": "21 CFR 1304.33",
            "timeline": "1-2 months",
            "complexity": "simple",
            "cost": (2000, 5000),
            "documents": [
                "ARCOS registration",
                "Reporting procedures",
                "Quarterly report templates",
            ],
        },
        "quota": {
            "category": "Quota Management",
            "requirement": "Annual production quota application",
            "cfr_reference": "21 CFR 1303.11-13",
            "timeline": "Annual (May 1 deadline)",
            "complexity": "moderate",
            "cost": (3000, 10000),
            "documents": [
                "DEA Form 189",
                "Production forecasts",
                "Sales projections",
                "Inventory reports",
            ],
        },
        "manufacturing": {
            "category": "Manufacturing",
            "requirement": "cGMP compliance for controlled substances",
            "cfr_reference": "21 CFR 211 + 1301.74",
            "timeline": "6-12 months",
            "complexity": "complex",
            "cost": (500000, 2000000),
            "documents": [
                "Master batch records",
                "Process validation protocols",
                "Quality assurance program",
                "Personnel training records",
                "Equipment qualification",
            ],
        },
    }
    
    # ==========================================================================
    # TRADE SECRET: FDA Botanical Drug Requirements
    # ==========================================================================
    _BOTANICAL_REQUIREMENTS: Dict[str, Dict] = {
        "cmc_identity": {
            "category": "CMC",
            "requirement": "Botanical drug substance identity testing",
            "guidance": "Botanical Drug Development Guidance for Industry (2016)",
            "rationale": "Ensure consistent botanical identity throughout development",
            "considerations": [
                "Macroscopic and microscopic examination",
                "Chromatographic fingerprinting (HPLC, GC)",
                "Marker compound quantification",
                "DNA barcoding for species verification",
            ],
        },
        "cmc_quality": {
            "category": "CMC",
            "requirement": "Quality control specifications",
            "guidance": "ICH Q6B, Botanical Guidance",
            "rationale": "Establish acceptance criteria for botanical consistency",
            "considerations": [
                "Cannabinoid content specifications",
                "Terpene profile specifications",
                "Residual solvent limits",
                "Pesticide residue limits",
                "Heavy metal limits",
                "Microbial limits",
            ],
        },
        "cmc_manufacturing": {
            "category": "CMC",
            "requirement": "Manufacturing process controls",
            "guidance": "21 CFR 211, Botanical Guidance",
            "rationale": "Ensure batch-to-batch consistency",
            "considerations": [
                "Raw material qualification",
                "In-process controls",
                "Process validation",
                "Hold time studies",
                "Container closure qualification",
            ],
        },
        "nonclinical_safety": {
            "category": "Nonclinical",
            "requirement": "Nonclinical safety assessment",
            "guidance": "ICH M3(R2), Botanical Guidance",
            "rationale": "Support human exposure in clinical trials",
            "considerations": [
                "Acute toxicity studies",
                "Repeat-dose toxicity (28-day minimum)",
                "Genotoxicity battery",
                "Reproductive toxicity (if applicable)",
                "Carcinogenicity (if long-term use)",
            ],
        },
        "clinical_phase1": {
            "category": "Clinical",
            "requirement": "Phase 1 clinical trial design",
            "guidance": "21 CFR 312, Botanical Guidance",
            "rationale": "Establish safety and PK in humans",
            "considerations": [
                "Dose-escalation design",
                "PK sampling strategy",
                "Safety monitoring plan",
                "Dose selection rationale",
            ],
        },
        "clinical_phase2": {
            "category": "Clinical",
            "requirement": "Phase 2 proof-of-concept",
            "guidance": "21 CFR 312, Botanical Guidance",
            "rationale": "Establish preliminary efficacy",
            "considerations": [
                "Indication-specific endpoints",
                "Dose-ranging design",
                "Patient population selection",
                "Efficacy signals to support Phase 3",
            ],
        },
        "clinical_phase3": {
            "category": "Clinical",
            "requirement": "Phase 3 pivotal trials",
            "guidance": "21 CFR 312, ICH E9",
            "rationale": "Confirm efficacy for NDA submission",
            "considerations": [
                "Adequate and well-controlled studies",
                "Pre-specified endpoints",
                "Statistical analysis plan",
                "Multi-center design",
            ],
        },
    }
    
    # ==========================================================================
    # TRADE SECRET: Timeline and Cost Matrices
    # ==========================================================================
    _TIMELINE_ESTIMATES: Dict[str, Dict] = {
        "Pre-IND": {"duration_months": 3, "cost_range": (50000, 150000)},
        "IND_preparation": {"duration_months": 6, "cost_range": (500000, 1500000)},
        "Phase_1": {"duration_months": 12, "cost_range": (2000000, 5000000)},
        "Phase_2": {"duration_months": 24, "cost_range": (10000000, 30000000)},
        "Phase_3": {"duration_months": 36, "cost_range": (50000000, 200000000)},
        "NDA_preparation": {"duration_months": 12, "cost_range": (5000000, 15000000)},
        "FDA_review": {"duration_months": 12, "cost_range": (2500000, 5000000)},
    }
    
    # Success probability by phase (historical botanical drugs)
    _SUCCESS_PROBABILITIES: Dict[str, float] = {
        "Pre-IND_to_IND": 0.85,
        "IND_to_Phase_1": 0.90,
        "Phase_1_to_Phase_2": 0.65,
        "Phase_2_to_Phase_3": 0.35,
        "Phase_3_to_NDA": 0.65,
        "NDA_to_Approval": 0.85,
    }
    
    def __init__(self):
        self._plans_cache: Dict[str, DEARegistrationPlan] = {}
        self._pathways_cache: Dict[str, BotanicalDrugPathway] = {}
    
    def generate_dea_registration_plan(
        self,
        proposed_schedule: DEASchedule = DEASchedule.SCHEDULE_III,
        license_types: Optional[List[DEALicenseType]] = None,
        estimated_annual_production_kg: float = 100.0,
    ) -> DEARegistrationPlan:
        """
        Generate comprehensive DEA registration plan for Schedule III cannabis.
        
        Trade Secret: Regulatory timeline optimization, cost estimation algorithms.
        
        Args:
            proposed_schedule: Target DEA schedule
            license_types: Types of licenses needed
            estimated_annual_production_kg: Projected annual production
            
        Returns:
            Complete DEARegistrationPlan
        """
        if license_types is None:
            license_types = [DEALicenseType.MANUFACTURER, DEALicenseType.DISTRIBUTOR]
        
        plan = DEARegistrationPlan(
            proposed_schedule=proposed_schedule,
            license_types_needed=license_types,
        )
        
        # Generate requirements
        plan.requirements = self._generate_dea_requirements(proposed_schedule, license_types)
        
        # Calculate timeline
        plan.total_timeline_months = self._calculate_dea_timeline(plan.requirements)
        
        # Identify critical path
        plan.critical_path_items = self._identify_critical_path(plan.requirements)
        
        # Calculate costs
        plan.estimated_total_cost = self._calculate_total_cost(plan.requirements)
        plan.annual_compliance_cost = self._calculate_annual_compliance_cost(
            proposed_schedule, license_types
        )
        
        # Quota information
        if proposed_schedule in [DEASchedule.SCHEDULE_II, DEASchedule.SCHEDULE_III]:
            plan.requires_quota = True
            plan.quota_application_deadline = self._calculate_quota_deadline()
            plan.estimated_annual_quota_kg = estimated_annual_production_kg
        
        # Security requirements
        plan.security_requirements = self._generate_security_requirements(
            proposed_schedule
        )
        
        # Cache
        self._plans_cache[plan.plan_id] = plan
        
        return plan
    
    def _generate_dea_requirements(
        self,
        schedule: DEASchedule,
        license_types: List[DEALicenseType],
    ) -> List[DEARequirement]:
        """Generate applicable DEA requirements."""
        requirements = []
        
        for req_key, req_data in self._DEA_REQUIREMENTS.items():
            requirement = DEARequirement(
                category=req_data["category"],
                requirement=req_data["requirement"],
                cfr_reference=req_data["cfr_reference"],
                timeline=req_data["timeline"],
                complexity=req_data["complexity"],
                estimated_cost=req_data["cost"],
                documentation_needed=req_data["documents"],
            )
            requirements.append(requirement)
        
        # Add manufacturer-specific requirements
        if DEALicenseType.MANUFACTURER in license_types:
            requirements.append(DEARequirement(
                category="Manufacturing",
                requirement="Schedule III manufacturing facility approval",
                cfr_reference="21 CFR 1301.74",
                timeline="6-12 months",
                complexity="complex",
                estimated_cost=(200000, 500000),
                documentation_needed=[
                    "Facility floor plans",
                    "Equipment specifications",
                    "Environmental controls documentation",
                    "Quality management system",
                ],
            ))
        
        # Add researcher-specific requirements
        if DEALicenseType.RESEARCHER in license_types:
            requirements.append(DEARequirement(
                category="Research",
                requirement="Research protocol approval",
                cfr_reference="21 CFR 1301.32",
                timeline="2-4 months",
                complexity="moderate",
                estimated_cost=(5000, 15000),
                documentation_needed=[
                    "Research protocol",
                    "Institutional review",
                    "Security plan for research quantities",
                ],
            ))
        
        return requirements
    
    def _calculate_dea_timeline(self, requirements: List[DEARequirement]) -> int:
        """Calculate total DEA registration timeline in months."""
        # Parse timelines and find critical path
        max_months = 0
        
        for req in requirements:
            # Parse timeline string to months
            timeline = req.timeline.lower()
            if "month" in timeline:
                # Extract number range
                parts = timeline.split("-")
                if len(parts) >= 2:
                    max_months = max(max_months, int(parts[1].split()[0]))
                else:
                    max_months = max(max_months, int(parts[0].split()[0]))
        
        # Add buffer for sequential dependencies
        return max_months + 3  # 3-month buffer
    
    def _identify_critical_path(
        self,
        requirements: List[DEARequirement],
    ) -> List[str]:
        """Identify critical path items."""
        critical = []
        
        for req in requirements:
            if req.complexity == "complex":
                critical.append(f"{req.category}: {req.requirement}")
        
        return critical
    
    def _calculate_total_cost(
        self,
        requirements: List[DEARequirement],
    ) -> Tuple[int, int]:
        """Calculate total estimated cost range."""
        min_total = sum(req.estimated_cost[0] for req in requirements)
        max_total = sum(req.estimated_cost[1] for req in requirements)
        return (min_total, max_total)
    
    def _calculate_annual_compliance_cost(
        self,
        schedule: DEASchedule,
        license_types: List[DEALicenseType],
    ) -> Tuple[int, int]:
        """Calculate annual compliance cost estimate."""
        base_cost = 50000  # Base compliance staff/systems
        
        if schedule == DEASchedule.SCHEDULE_III:
            base_cost += 30000
        elif schedule == DEASchedule.SCHEDULE_II:
            base_cost += 75000
        
        if DEALicenseType.MANUFACTURER in license_types:
            base_cost += 100000  # Quality systems, audits
        
        return (int(base_cost * 0.8), int(base_cost * 1.3))
    
    def _calculate_quota_deadline(self) -> date:
        """Calculate next quota application deadline (May 1)."""
        today = date.today()
        year = today.year if today.month < 5 else today.year + 1
        return date(year, 5, 1)
    
    def _generate_security_requirements(
        self,
        schedule: DEASchedule,
    ) -> List[str]:
        """Generate security requirements based on schedule."""
        base_requirements = [
            "Controlled access to storage areas",
            "Alarm system with 24/7 monitoring",
            "Video surveillance with 90-day retention",
            "Employee background checks",
            "Visitor log and escort procedures",
            "Key control procedures",
        ]
        
        if schedule in [DEASchedule.SCHEDULE_II, DEASchedule.SCHEDULE_III]:
            base_requirements.extend([
                "Vault or safe meeting DEA specifications",
                "Access limited to authorized personnel",
                "Physical barrier requirements",
            ])
        
        return base_requirements
    
    def design_botanical_drug_pathway(
        self,
        product_name: str,
        botanical_source: str = "Cannabis sativa L.",
        active_constituents: Optional[List[str]] = None,
        target_indication: str = "chronic pain",
        start_date: Optional[date] = None,
    ) -> BotanicalDrugPathway:
        """
        Design FDA botanical drug development pathway.
        
        Trade Secret: Pathway optimization, timeline prediction,
        cost estimation for botanical drug development.
        
        Args:
            product_name: Name of the botanical drug product
            botanical_source: Scientific name of botanical source
            active_constituents: Key active compounds
            target_indication: Primary therapeutic indication
            start_date: Development start date
            
        Returns:
            Complete BotanicalDrugPathway
        """
        if active_constituents is None:
            active_constituents = ["THC", "CBD", "terpenes"]
        if start_date is None:
            start_date = date.today()
        
        pathway = BotanicalDrugPathway(
            product_description=f"{product_name} for {target_indication}",
            botanical_source=botanical_source,
            active_constituents=active_constituents,
            standardization_markers=["THC", "CBD", "total cannabinoids", "myrcene"],
        )
        
        # Generate milestones
        pathway.milestones = self._generate_fda_milestones(start_date)
        
        # Design clinical trials
        pathway.clinical_trials = self._design_clinical_program(target_indication)
        
        # Generate botanical-specific requirements
        pathway.botanical_requirements = self._generate_botanical_requirements()
        
        # Calculate timeline and costs
        pathway.total_development_years = self._calculate_total_development_time()
        pathway.estimated_total_cost_millions = self._calculate_total_development_cost()
        pathway.probability_of_success = self._calculate_cumulative_success_probability()
        
        # Identify risks
        pathway.risk_factors = self._identify_development_risks(target_indication)
        pathway.mitigation_strategies = self._generate_mitigation_strategies(
            pathway.risk_factors
        )
        
        # Cache
        self._pathways_cache[pathway.pathway_id] = pathway
        
        return pathway
    
    def _generate_fda_milestones(
        self,
        start_date: date,
    ) -> List[FDASubmissionMilestone]:
        """Generate FDA submission milestones."""
        milestones = []
        current_date = start_date
        
        # Pre-IND
        milestones.append(FDASubmissionMilestone(
            name="Pre-IND Meeting",
            submission_type=SubmissionType.PRE_IND,
            target_date=current_date + timedelta(days=90),
            estimated_duration_days=90,
            required_documents=[
                "Briefing document",
                "Preliminary CMC information",
                "Nonclinical summary",
                "Proposed clinical plan",
            ],
            estimated_cost=(50000, 150000),
        ))
        current_date += timedelta(days=90)
        
        # IND Submission
        milestones.append(FDASubmissionMilestone(
            name="IND Submission",
            submission_type=SubmissionType.IND,
            target_date=current_date + timedelta(days=180),
            dependencies=["Pre-IND Meeting"],
            estimated_duration_days=180,
            required_documents=[
                "Form FDA 1571",
                "CMC section (Module 3)",
                "Nonclinical section (Module 4)",
                "Clinical protocol (Module 5)",
                "Investigator's Brochure",
            ],
            estimated_cost=(500000, 1500000),
        ))
        current_date += timedelta(days=180 + 30)  # +30 for FDA review
        
        # End of Phase 1
        milestones.append(FDASubmissionMilestone(
            name="End-of-Phase 1 Meeting",
            submission_type=SubmissionType.END_OF_PHASE_1,
            target_date=current_date + timedelta(days=365),
            dependencies=["IND Submission"],
            estimated_duration_days=60,
            required_documents=[
                "Phase 1 data summary",
                "Updated IB",
                "Phase 2 protocol proposal",
            ],
            estimated_cost=(30000, 100000),
        ))
        current_date += timedelta(days=365 + 60)
        
        # End of Phase 2
        milestones.append(FDASubmissionMilestone(
            name="End-of-Phase 2 Meeting",
            submission_type=SubmissionType.END_OF_PHASE_2,
            target_date=current_date + timedelta(days=730),
            dependencies=["End-of-Phase 1 Meeting"],
            estimated_duration_days=60,
            required_documents=[
                "Phase 2 data summary",
                "Phase 3 design proposal",
                "Special Protocol Assessment request",
            ],
            estimated_cost=(50000, 150000),
        ))
        current_date += timedelta(days=730 + 60)
        
        # Pre-NDA Meeting
        milestones.append(FDASubmissionMilestone(
            name="Pre-NDA Meeting",
            submission_type=SubmissionType.PRE_NDA,
            target_date=current_date + timedelta(days=1095),
            dependencies=["End-of-Phase 2 Meeting"],
            estimated_duration_days=60,
            required_documents=[
                "Phase 3 data summary",
                "NDA content proposal",
                "CMC update",
            ],
            estimated_cost=(50000, 150000),
        ))
        current_date += timedelta(days=1095 + 60)
        
        # NDA Submission
        milestones.append(FDASubmissionMilestone(
            name="NDA Submission",
            submission_type=SubmissionType.NDA,
            target_date=current_date + timedelta(days=365),
            dependencies=["Pre-NDA Meeting"],
            estimated_duration_days=365,
            required_documents=[
                "Complete CTD (Modules 1-5)",
                "All clinical study reports",
                "Integrated safety summary",
                "Integrated efficacy summary",
                "Risk Evaluation and Mitigation Strategy (REMS)",
            ],
            estimated_cost=(5000000, 15000000),
        ))
        
        return milestones
    
    def _design_clinical_program(
        self,
        indication: str,
    ) -> List[ClinicalTrialDesign]:
        """Design clinical trial program."""
        trials = []
        
        # Phase 1
        trials.append(ClinicalTrialDesign(
            phase="Phase 1",
            design_type="open-label, dose-escalation",
            primary_endpoint="Safety and tolerability",
            secondary_endpoints=[
                "Pharmacokinetics",
                "Preliminary efficacy signals",
                "Dose-response relationship",
            ],
            sample_size=60,
            duration_weeks=8,
            estimated_cost=(2000000, 5000000),
            regulatory_considerations=[
                "Sentinel dosing required",
                "Schedule III handling procedures",
                "Abuse liability assessment",
            ],
        ))
        
        # Phase 2
        trials.append(ClinicalTrialDesign(
            phase="Phase 2",
            design_type="randomized, double-blind, placebo-controlled",
            primary_endpoint=self._get_primary_endpoint(indication),
            secondary_endpoints=[
                "Patient-reported outcomes",
                "Quality of life measures",
                "Dose optimization",
            ],
            sample_size=200,
            duration_weeks=12,
            estimated_cost=(10000000, 30000000),
            regulatory_considerations=[
                "Enrichment design considerations",
                "Adaptive design elements",
                "Schedule III compliance at sites",
            ],
        ))
        
        # Phase 3 (two pivotal trials typically required)
        trials.append(ClinicalTrialDesign(
            phase="Phase 3a",
            design_type="randomized, double-blind, placebo-controlled, multi-center",
            primary_endpoint=self._get_primary_endpoint(indication),
            secondary_endpoints=[
                "Durability of response",
                "Safety in broader population",
                "Health economic outcomes",
            ],
            sample_size=400,
            duration_weeks=24,
            estimated_cost=(25000000, 75000000),
            regulatory_considerations=[
                "Pre-specified statistical analysis plan",
                "Independent Data Monitoring Committee",
                "Special Protocol Assessment",
            ],
        ))
        
        trials.append(ClinicalTrialDesign(
            phase="Phase 3b",
            design_type="randomized, double-blind, active-controlled, multi-center",
            primary_endpoint=self._get_primary_endpoint(indication),
            secondary_endpoints=[
                "Non-inferiority to standard of care",
                "Long-term safety",
                "Abuse potential assessment",
            ],
            sample_size=500,
            duration_weeks=52,
            estimated_cost=(30000000, 100000000),
            regulatory_considerations=[
                "Active comparator selection",
                "Non-inferiority margin justification",
                "Post-marketing commitments",
            ],
        ))
        
        return trials
    
    def _get_primary_endpoint(self, indication: str) -> str:
        """Get indication-specific primary endpoint."""
        endpoints = {
            "chronic pain": "Change in average pain intensity (11-point NRS)",
            "epilepsy": "Percent change in seizure frequency",
            "anxiety": "Change in HAM-A score from baseline",
            "nausea": "Complete response rate (no emesis, no rescue)",
            "spasticity": "Change in Modified Ashworth Scale",
            "ptsd": "Change in CAPS-5 total severity score",
            "insomnia": "Change in subjective sleep quality (ISI)",
        }
        return endpoints.get(indication.lower(), "Change in disease-specific primary measure")
    
    def _generate_botanical_requirements(self) -> List[BotanicalDrugRequirement]:
        """Generate botanical drug-specific requirements."""
        requirements = []
        
        for req_key, req_data in self._BOTANICAL_REQUIREMENTS.items():
            requirement = BotanicalDrugRequirement(
                category=req_data["category"],
                requirement=req_data["requirement"],
                guidance_reference=req_data["guidance"],
                rationale=req_data["rationale"],
                special_considerations=req_data["considerations"],
            )
            requirements.append(requirement)
        
        return requirements
    
    def _calculate_total_development_time(self) -> float:
        """Calculate total development time in years."""
        total_months = sum(
            data["duration_months"] 
            for data in self._TIMELINE_ESTIMATES.values()
        )
        return total_months / 12.0
    
    def _calculate_total_development_cost(self) -> Tuple[float, float]:
        """Calculate total development cost in millions."""
        min_total = sum(
            data["cost_range"][0] 
            for data in self._TIMELINE_ESTIMATES.values()
        ) / 1_000_000
        max_total = sum(
            data["cost_range"][1] 
            for data in self._TIMELINE_ESTIMATES.values()
        ) / 1_000_000
        return (min_total, max_total)
    
    def _calculate_cumulative_success_probability(self) -> float:
        """Calculate cumulative probability of success."""
        probability = 1.0
        for phase_prob in self._SUCCESS_PROBABILITIES.values():
            probability *= phase_prob
        return probability
    
    def _identify_development_risks(self, indication: str) -> List[str]:
        """Identify development risks specific to cannabis botanicals."""
        risks = [
            "Batch-to-batch variability in cannabinoid/terpene content",
            "Regulatory uncertainty around Schedule III reclassification timing",
            "Competition from state-legal cannabis markets",
            "Difficulty recruiting cannabis-naïve patients for trials",
            "Placebo response rates in pain/anxiety indications",
            "DEA quota limitations affecting clinical supply",
            "Standardization challenges for complex botanical mixture",
            "Potential for abuse liability findings affecting scheduling",
            "Long-term safety data requirements",
            "REMS requirement likelihood",
        ]
        
        if indication.lower() == "chronic pain":
            risks.append("Opioid comparator trial requirements")
        
        return risks
    
    def _generate_mitigation_strategies(
        self,
        risks: List[str],
    ) -> List[str]:
        """Generate mitigation strategies for identified risks."""
        return [
            "Implement robust agricultural and manufacturing controls for consistency",
            "Engage early with DEA and FDA through Pre-IND meetings",
            "Design clinical program to demonstrate clear differentiation from recreational use",
            "Include experienced cannabis researchers as investigators",
            "Use validated, FDA-accepted endpoints with historical context",
            "Secure adequate DEA quota early with contingency planning",
            "Develop comprehensive fingerprinting and standardization methods",
            "Conduct thorough human abuse potential studies early",
            "Plan for 12-month Phase 3 safety extension studies",
            "Design REMS-ready risk management program from Phase 2",
        ]
    
    def get_registration_plan(self, plan_id: str) -> Optional[DEARegistrationPlan]:
        """Retrieve cached DEA registration plan."""
        return self._plans_cache.get(plan_id)
    
    def get_botanical_pathway(self, pathway_id: str) -> Optional[BotanicalDrugPathway]:
        """Retrieve cached botanical drug pathway."""
        return self._pathways_cache.get(pathway_id)
    
    def estimate_time_to_market(
        self,
        current_phase: str = "Pre-clinical",
    ) -> Dict[str, any]:
        """
        Estimate time to market from current phase.
        
        Trade Secret: Timeline optimization algorithm.
        """
        phase_order = [
            "Pre-clinical",
            "Pre-IND",
            "IND_preparation", 
            "Phase_1",
            "Phase_2",
            "Phase_3",
            "NDA_preparation",
            "FDA_review",
        ]
        
        try:
            start_idx = phase_order.index(current_phase)
        except ValueError:
            start_idx = 0
        
        remaining_phases = phase_order[start_idx:]
        
        total_months = sum(
            self._TIMELINE_ESTIMATES.get(phase, {"duration_months": 0})["duration_months"]
            for phase in remaining_phases
            if phase in self._TIMELINE_ESTIMATES
        )
        
        min_cost = sum(
            self._TIMELINE_ESTIMATES.get(phase, {"cost_range": (0, 0)})["cost_range"][0]
            for phase in remaining_phases
            if phase in self._TIMELINE_ESTIMATES
        )
        max_cost = sum(
            self._TIMELINE_ESTIMATES.get(phase, {"cost_range": (0, 0)})["cost_range"][1]
            for phase in remaining_phases
            if phase in self._TIMELINE_ESTIMATES
        )
        
        # Calculate success probability from current phase
        success_prob = 1.0
        for i, phase in enumerate(remaining_phases[:-1]):
            next_phase = remaining_phases[i + 1]
            transition = f"{phase}_to_{next_phase}"
            if transition in self._SUCCESS_PROBABILITIES:
                success_prob *= self._SUCCESS_PROBABILITIES[transition]
        
        return {
            "current_phase": current_phase,
            "estimated_months_to_approval": total_months,
            "estimated_years_to_approval": total_months / 12.0,
            "remaining_phases": remaining_phases,
            "estimated_cost_range_millions": (min_cost / 1_000_000, max_cost / 1_000_000),
            "success_probability": success_prob,
            "target_approval_date": date.today() + timedelta(days=total_months * 30),
        }
