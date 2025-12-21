"""
Schedule III Documentation Templates
Patent Claim 1(j) - FDA Cannabis Rescheduling Support

Generates comprehensive documentation templates for:
- Schedule III transition applications
- DEA registration requirements
- State-federal compliance bridging
- Quality assurance frameworks

Reference: NeuroBotanica Provisional Patent [0033], [0141], Claim 1(j)
"""
from typing import Dict, Any, List, Optional
from datetime import datetime
from dataclasses import dataclass, asdict
import json


@dataclass
class RegulatoryContact:
    """Regulatory agency contact information."""
    agency: str
    division: str
    address: str
    phone: str
    submission_portal: str


@dataclass
class ScheduleIIIRequirement:
    """Individual Schedule III requirement."""
    requirement_id: str
    category: str
    description: str
    regulatory_citation: str
    evidence_needed: List[str]
    deadline_type: str  # "pre-submission", "ongoing", "annual"
    compliance_status: str


class ScheduleIIIDocumentationGenerator:
    """
    Generates Schedule III compliance documentation templates.
    
    Supports the cannabis rescheduling transition by providing:
    - Pre-formatted regulatory templates
    - Evidence organization frameworks
    - Compliance checklists
    - State-federal bridging documentation
    """
    
    # DEA Schedule III requirements for cannabis
    SCHEDULE_III_REQUIREMENTS = [
        {
            "requirement_id": "S3-001",
            "category": "Abuse Potential",
            "description": "Demonstrate moderate to low potential for physical and psychological dependence",
            "regulatory_citation": "21 USC § 812(b)(3)",
            "evidence_needed": [
                "Clinical dependence studies",
                "Withdrawal syndrome characterization",
                "Comparative abuse liability data",
                "Post-marketing surveillance data"
            ],
            "deadline_type": "pre-submission"
        },
        {
            "requirement_id": "S3-002",
            "category": "Accepted Medical Use",
            "description": "Establish currently accepted medical use in treatment in the United States",
            "regulatory_citation": "21 USC § 812(b)(3)(A)",
            "evidence_needed": [
                "FDA-approved product labeling",
                "Peer-reviewed efficacy studies",
                "Clinical practice guidelines",
                "Expert consensus statements"
            ],
            "deadline_type": "pre-submission"
        },
        {
            "requirement_id": "S3-003",
            "category": "Safety Profile",
            "description": "Document safety profile supporting medical use under supervision",
            "regulatory_citation": "21 USC § 812(b)(3)(B)",
            "evidence_needed": [
                "Adverse event database analysis",
                "Drug-drug interaction studies",
                "Special population safety data",
                "Overdose and toxicity reports"
            ],
            "deadline_type": "pre-submission"
        },
        {
            "requirement_id": "S3-004",
            "category": "Manufacturing Controls",
            "description": "Implement Schedule III manufacturing and distribution controls",
            "regulatory_citation": "21 CFR Part 1301-1321",
            "evidence_needed": [
                "DEA registration application",
                "Facility security assessment",
                "Inventory control procedures",
                "Distribution chain documentation"
            ],
            "deadline_type": "ongoing"
        },
        {
            "requirement_id": "S3-005",
            "category": "Prescription Requirements",
            "description": "Establish prescription and dispensing protocols",
            "regulatory_citation": "21 CFR § 1306.04",
            "evidence_needed": [
                "Prescriber registration procedures",
                "Prescription monitoring integration",
                "Refill limitation protocols",
                "Patient counseling requirements"
            ],
            "deadline_type": "ongoing"
        },
        {
            "requirement_id": "S3-006",
            "category": "Record Keeping",
            "description": "Maintain comprehensive Schedule III record keeping",
            "regulatory_citation": "21 CFR § 1304",
            "evidence_needed": [
                "Inventory records (initial and biennial)",
                "Dispensing records",
                "Distribution records",
                "Theft/loss reports"
            ],
            "deadline_type": "ongoing"
        },
        {
            "requirement_id": "S3-007",
            "category": "Reporting Requirements",
            "description": "Submit required reports to DEA and FDA",
            "regulatory_citation": "21 CFR § 1301.74-76",
            "evidence_needed": [
                "Annual production quotas",
                "ARCOS transaction reports",
                "Suspicious order monitoring",
                "Adverse event reports"
            ],
            "deadline_type": "annual"
        }
    ]
    
    # FDA regulatory contacts
    REGULATORY_CONTACTS = {
        "FDA_CDER": RegulatoryContact(
            agency="FDA Center for Drug Evaluation and Research",
            division="Division of Drug Information",
            address="10001 New Hampshire Ave, Silver Spring, MD 20993",
            phone="1-855-543-3784",
            submission_portal="https://www.fda.gov/drugs/development-approval-process-drugs"
        ),
        "DEA_Diversion": RegulatoryContact(
            agency="DEA Diversion Control Division",
            division="Drug & Chemical Evaluation Section",
            address="8701 Morrissette Drive, Springfield, VA 22152",
            phone="1-202-307-7165",
            submission_portal="https://www.deadiversion.usdoj.gov/"
        ),
        "FDA_Cannabis": RegulatoryContact(
            agency="FDA Cannabis Product Committee",
            division="Office of Regulatory Affairs",
            address="10903 New Hampshire Ave, Silver Spring, MD 20993",
            phone="1-888-INFO-FDA",
            submission_portal="https://www.fda.gov/news-events/public-health-focus/fda-and-cannabis"
        )
    }
    
    def __init__(self, company_name: str = "NeuroBotanica Client"):
        self.company_name = company_name
        self.generation_date = datetime.now()
    
    def generate_full_documentation_package(
        self,
        compound_name: str,
        indication: str,
        clinical_evidence: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        Generate complete Schedule III documentation package.
        
        Returns all templates and checklists needed for
        Schedule III transition compliance.
        """
        return {
            "metadata": {
                "document_type": "Schedule III Transition Documentation Package",
                "compound": compound_name,
                "primary_indication": indication,
                "generated_date": self.generation_date.isoformat(),
                "prepared_for": self.company_name,
                "patent_reference": "NeuroBotanica Claim 1(j)"
            },
            "executive_summary": self._generate_executive_summary(
                compound_name, indication, clinical_evidence
            ),
            "requirements_checklist": self._generate_requirements_checklist(),
            "evidence_matrix": self._generate_evidence_matrix(clinical_evidence),
            "dea_registration_template": self._generate_dea_registration(),
            "state_compliance_bridge": self._generate_state_bridge(),
            "quality_assurance_framework": self._generate_qa_framework(),
            "timeline_milestones": self._generate_timeline(),
            "regulatory_contacts": self._get_contacts_dict(),
            "appendices": {
                "appendix_a": "Clinical Evidence Summary Tables",
                "appendix_b": "Manufacturing SOPs Template",
                "appendix_c": "Inventory Control Procedures",
                "appendix_d": "Prescriber Education Materials"
            }
        }
    
    def _generate_executive_summary(
        self,
        compound: str,
        indication: str,
        evidence: Dict[str, Any]
    ) -> Dict[str, Any]:
        """Generate executive summary for regulatory submission."""
        study_count = evidence.get("studies_analyzed", 0)
        effect_size = evidence.get("effect_size_statistics", {}).get("weighted_mean", 0)
        
        return {
            "title": f"Schedule III Transition Application: {compound} for {indication}",
            "summary": f"""
This documentation package supports the Schedule III classification of {compound} 
for the treatment of {indication}. The evidence base includes {study_count} clinical 
studies demonstrating efficacy (weighted mean effect size: {effect_size:.2f}).

Key findings supporting Schedule III classification:
1. Established medical use with FDA-recognized therapeutic applications
2. Moderate abuse potential relative to Schedule I/II substances
3. Acceptable safety profile with appropriate medical supervision
4. Comparable efficacy to existing Schedule III medications

This package includes all documentation required under 21 USC § 812(b)(3) and 
21 CFR Parts 1301-1321 for Schedule III controlled substance compliance.
            """.strip(),
            "recommendation": "Proceed with Schedule III application",
            "confidence_level": "High" if study_count >= 5 and effect_size >= 0.5 else "Moderate"
        }
    
    def _generate_requirements_checklist(self) -> List[Dict[str, Any]]:
        """Generate compliance requirements checklist."""
        checklist = []
        for req in self.SCHEDULE_III_REQUIREMENTS:
            checklist.append({
                **req,
                "compliance_status": "Pending Review",
                "assigned_to": "",
                "target_completion": "",
                "notes": ""
            })
        return checklist
    
    def _generate_evidence_matrix(
        self,
        clinical_evidence: Dict[str, Any]
    ) -> Dict[str, Any]:
        """Generate evidence-to-requirement mapping matrix."""
        matrix = {
            "abuse_potential": {
                "requirement": "S3-001",
                "evidence_available": [],
                "gaps_identified": [],
                "action_items": []
            },
            "medical_use": {
                "requirement": "S3-002",
                "evidence_available": [],
                "gaps_identified": [],
                "action_items": []
            },
            "safety_profile": {
                "requirement": "S3-003",
                "evidence_available": [],
                "gaps_identified": [],
                "action_items": []
            }
        }
        
        # Map clinical evidence to requirements
        studies = clinical_evidence.get("studies_analyzed", 0)
        evidence_quality = clinical_evidence.get("evidence_quality", {})
        
        # Abuse potential evidence
        if studies > 0:
            matrix["abuse_potential"]["evidence_available"].append(
                f"{studies} clinical studies with safety/tolerability data"
            )
        
        # Medical use evidence
        rcts = evidence_quality.get("Level 1", 0) + evidence_quality.get("Level 2", 0)
        if rcts > 0:
            matrix["medical_use"]["evidence_available"].append(
                f"{rcts} randomized controlled trials"
            )
        
        pivotal = clinical_evidence.get("pivotal_trials", [])
        if pivotal:
            matrix["medical_use"]["evidence_available"].append(
                f"{len(pivotal)} pivotal trials identified"
            )
        
        # Safety profile evidence
        matrix["safety_profile"]["evidence_available"].append(
            "Adverse event data from clinical studies"
        )
        
        # Identify gaps
        if rcts < 2:
            matrix["medical_use"]["gaps_identified"].append(
                "Additional Phase III RCTs may strengthen application"
            )
            matrix["medical_use"]["action_items"].append(
                "Consider conducting confirmatory Phase III trial"
            )
        
        return matrix
    
    def _generate_dea_registration(self) -> Dict[str, Any]:
        """Generate DEA registration application template."""
        return {
            "form_type": "DEA Form 224 (New Registration)",
            "registration_category": "Schedule III Manufacturer/Distributor",
            "required_sections": [
                {
                    "section": "1. Applicant Information",
                    "fields": [
                        "Legal Business Name",
                        "Trade/DBA Name",
                        "Physical Address (no P.O. Box)",
                        "Mailing Address",
                        "DUNS Number",
                        "Tax ID (EIN)"
                    ]
                },
                {
                    "section": "2. Business Activity",
                    "fields": [
                        "Primary Business Activity Code",
                        "Schedule III Substances to be Handled",
                        "Estimated Annual Volume",
                        "State License Numbers"
                    ]
                },
                {
                    "section": "3. Facility Security",
                    "fields": [
                        "Physical Security Description",
                        "Alarm System Details",
                        "Safe/Vault Specifications",
                        "Employee Background Check Procedures"
                    ]
                },
                {
                    "section": "4. Responsible Personnel",
                    "fields": [
                        "Responsible Person Name and Title",
                        "DEA Registration (if existing)",
                        "Professional License Information",
                        "Background Check Consent"
                    ]
                }
            ],
            "submission_instructions": {
                "online": "https://www.deadiversion.usdoj.gov/online_forms.html",
                "fee": "$3,047 (Manufacturer) / $888 (Distributor)",
                "processing_time": "4-6 weeks",
                "renewal": "Every 3 years"
            }
        }
    
    def _generate_state_bridge(self) -> Dict[str, Any]:
        """Generate state-federal compliance bridging documentation."""
        return {
            "purpose": "Bridge state medical cannabis programs with federal Schedule III framework",
            "key_considerations": [
                {
                    "area": "Licensing Harmonization",
                    "federal_requirement": "DEA Schedule III registration",
                    "state_requirement": "State cannabis license",
                    "bridge_strategy": "Dual licensure with cross-referenced compliance programs"
                },
                {
                    "area": "Product Standards",
                    "federal_requirement": "FDA cGMP compliance",
                    "state_requirement": "State-specific testing and labeling",
                    "bridge_strategy": "Implement unified quality system meeting both requirements"
                },
                {
                    "area": "Distribution",
                    "federal_requirement": "DEA-registered distributors only",
                    "state_requirement": "State-licensed dispensaries",
                    "bridge_strategy": "Establish DEA-registered distribution hub model"
                },
                {
                    "area": "Prescribing",
                    "federal_requirement": "DEA-registered prescriber with valid prescription",
                    "state_requirement": "State medical cannabis card/recommendation",
                    "bridge_strategy": "Transition to prescription model with state integration"
                },
                {
                    "area": "Record Keeping",
                    "federal_requirement": "DEA ARCOS reporting",
                    "state_requirement": "State seed-to-sale tracking",
                    "bridge_strategy": "Unified tracking system with dual reporting capability"
                }
            ],
            "implementation_phases": [
                {
                    "phase": 1,
                    "timeline": "Months 1-3",
                    "activities": [
                        "Gap analysis of current state compliance vs federal requirements",
                        "DEA registration application submission",
                        "Quality system upgrade planning"
                    ]
                },
                {
                    "phase": 2,
                    "timeline": "Months 4-6",
                    "activities": [
                        "Facility security upgrades",
                        "Staff training on federal requirements",
                        "Inventory system integration"
                    ]
                },
                {
                    "phase": 3,
                    "timeline": "Months 7-12",
                    "activities": [
                        "Dual compliance audit",
                        "Prescriber network development",
                        "Full Schedule III operation launch"
                    ]
                }
            ]
        }
    
    def _generate_qa_framework(self) -> Dict[str, Any]:
        """Generate quality assurance framework for Schedule III compliance."""
        return {
            "framework_name": "NeuroBotanica Schedule III Quality System",
            "standard_reference": "21 CFR Parts 210, 211 (cGMP); 21 CFR Part 1300 (DEA)",
            "quality_elements": [
                {
                    "element": "Quality Management",
                    "requirements": [
                        "Quality policy and objectives",
                        "Management responsibility",
                        "Quality manual",
                        "Document control"
                    ],
                    "sop_templates": ["QM-001: Quality Manual", "QM-002: Document Control"]
                },
                {
                    "element": "Personnel",
                    "requirements": [
                        "Qualification requirements",
                        "Training programs",
                        "Background checks",
                        "Competency assessments"
                    ],
                    "sop_templates": ["HR-001: Personnel Qualification", "HR-002: Training Program"]
                },
                {
                    "element": "Facilities and Equipment",
                    "requirements": [
                        "Facility design and security",
                        "Equipment qualification",
                        "Calibration program",
                        "Maintenance procedures"
                    ],
                    "sop_templates": ["FE-001: Facility Security", "FE-002: Equipment Qualification"]
                },
                {
                    "element": "Materials Management",
                    "requirements": [
                        "Supplier qualification",
                        "Incoming material testing",
                        "Inventory control",
                        "Storage conditions"
                    ],
                    "sop_templates": ["MM-001: Supplier Qualification", "MM-002: Inventory Control"]
                },
                {
                    "element": "Production",
                    "requirements": [
                        "Batch records",
                        "In-process controls",
                        "Contamination prevention",
                        "Yield reconciliation"
                    ],
                    "sop_templates": ["PR-001: Batch Production", "PR-002: In-Process Testing"]
                },
                {
                    "element": "Laboratory Controls",
                    "requirements": [
                        "Testing procedures",
                        "Specifications",
                        "Stability program",
                        "Reference standards"
                    ],
                    "sop_templates": ["LC-001: Laboratory Testing", "LC-002: Stability Program"]
                },
                {
                    "element": "Packaging and Labeling",
                    "requirements": [
                        "Label control",
                        "Child-resistant packaging",
                        "Serialization",
                        "Label accuracy verification"
                    ],
                    "sop_templates": ["PL-001: Label Control", "PL-002: Packaging Operations"]
                },
                {
                    "element": "Controlled Substance Controls",
                    "requirements": [
                        "DEA registration maintenance",
                        "Physical security",
                        "Inventory accountability",
                        "Theft/loss reporting"
                    ],
                    "sop_templates": ["CS-001: Controlled Substance Handling", "CS-002: Security Procedures"]
                }
            ],
            "audit_schedule": {
                "internal_audits": "Quarterly",
                "management_review": "Annual",
                "dea_inspection_prep": "As required",
                "fda_inspection_prep": "As required"
            }
        }
    
    def _generate_timeline(self) -> Dict[str, Any]:
        """Generate Schedule III transition timeline."""
        return {
            "project_name": "Schedule III Compliance Transition",
            "total_duration": "12-18 months",
            "milestones": [
                {
                    "milestone": "M1: Project Initiation",
                    "target_month": 1,
                    "deliverables": [
                        "Project charter approval",
                        "Resource allocation",
                        "Gap analysis completion"
                    ],
                    "status": "Not Started"
                },
                {
                    "milestone": "M2: DEA Registration",
                    "target_month": 3,
                    "deliverables": [
                        "DEA Form 224 submission",
                        "Background check completion",
                        "Security system installation"
                    ],
                    "status": "Not Started"
                },
                {
                    "milestone": "M3: Quality System Implementation",
                    "target_month": 6,
                    "deliverables": [
                        "All SOPs approved",
                        "Staff training complete",
                        "Equipment qualified"
                    ],
                    "status": "Not Started"
                },
                {
                    "milestone": "M4: Pilot Operations",
                    "target_month": 9,
                    "deliverables": [
                        "First Schedule III batch produced",
                        "Distribution network validated",
                        "Prescriber onboarding initiated"
                    ],
                    "status": "Not Started"
                },
                {
                    "milestone": "M5: Full Commercial Launch",
                    "target_month": 12,
                    "deliverables": [
                        "Full-scale operations",
                        "Compliance audit passed",
                        "Post-market surveillance active"
                    ],
                    "status": "Not Started"
                }
            ],
            "critical_path_items": [
                "DEA registration approval",
                "cGMP facility certification",
                "State-federal license reconciliation"
            ]
        }
    
    def _get_contacts_dict(self) -> Dict[str, Dict]:
        """Convert regulatory contacts to dictionary format."""
        return {
            name: asdict(contact) 
            for name, contact in self.REGULATORY_CONTACTS.items()
        }
    
    def generate_compliance_report(
        self,
        requirements_status: List[Dict[str, Any]]
    ) -> Dict[str, Any]:
        """Generate compliance status report from checklist data."""
        total = len(requirements_status)
        compliant = sum(1 for r in requirements_status if r.get("compliance_status") == "Compliant")
        in_progress = sum(1 for r in requirements_status if r.get("compliance_status") == "In Progress")
        non_compliant = sum(1 for r in requirements_status if r.get("compliance_status") == "Non-Compliant")
        pending = total - compliant - in_progress - non_compliant
        
        compliance_percentage = (compliant / total * 100) if total > 0 else 0
        
        return {
            "report_date": datetime.now().isoformat(),
            "overall_status": "Ready" if compliance_percentage >= 90 else "In Progress",
            "compliance_summary": {
                "total_requirements": total,
                "compliant": compliant,
                "in_progress": in_progress,
                "non_compliant": non_compliant,
                "pending_review": pending,
                "compliance_percentage": round(compliance_percentage, 1)
            },
            "risk_areas": [
                r for r in requirements_status 
                if r.get("compliance_status") in ["Non-Compliant", "Pending Review"]
            ],
            "next_actions": self._determine_next_actions(requirements_status),
            "estimated_completion": self._estimate_completion(compliance_percentage)
        }
    
    def _determine_next_actions(
        self,
        requirements_status: List[Dict[str, Any]]
    ) -> List[str]:
        """Determine next actions based on compliance status."""
        actions = []
        
        non_compliant = [r for r in requirements_status if r.get("compliance_status") == "Non-Compliant"]
        if non_compliant:
            for req in non_compliant[:3]:  # Top 3 priority
                actions.append(f"Address {req['requirement_id']}: {req['category']}")
        
        pending = [r for r in requirements_status if r.get("compliance_status") == "Pending Review"]
        if pending:
            actions.append(f"Review {len(pending)} pending requirements")
        
        if not actions:
            actions.append("Prepare for final compliance audit")
            actions.append("Schedule DEA inspection")
        
        return actions
    
    def _estimate_completion(self, current_percentage: float) -> str:
        """Estimate completion timeline based on current progress."""
        if current_percentage >= 90:
            return "1-2 months"
        elif current_percentage >= 70:
            return "3-4 months"
        elif current_percentage >= 50:
            return "5-6 months"
        else:
            return "6+ months"
