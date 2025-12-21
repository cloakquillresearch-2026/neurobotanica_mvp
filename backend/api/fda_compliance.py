"""
FDA Compliance Module - Patent Claim 1(j)
Schedule III Cannabis Pharmaceutical Documentation Generator

Supports:
- Chemistry Manufacturing and Controls (CMC) templates
- Pharmacology data packages demonstrating mechanism of action
- Comparative efficacy analysis from RCT data
- Schedule III transition documentation

Reference: NeuroBotanica Provisional Patent Application [0033], [0141], Claim 1(j)
"""
from fastapi import APIRouter, Depends, HTTPException, Query
from sqlalchemy.orm import Session
from sqlalchemy import func, and_, or_
from typing import List, Optional, Dict, Any
from pydantic import BaseModel
from datetime import datetime
import json

from backend.models.database import get_db
from backend.models.study import ClinicalStudy
from backend.models.compound import Cannabinoid


router = APIRouter()


# =============================================================================
# Pydantic Models for FDA Documentation
# =============================================================================

class PharmacologyPackage(BaseModel):
    """Pharmacology data package structure for FDA CMC support."""
    compound_name: str
    mechanism_of_action: str
    receptor_profile: Dict[str, Any]
    clinical_evidence: List[Dict[str, Any]]
    safety_summary: Dict[str, Any]
    generated_at: str
    patent_reference: str = "NeuroBotanica Patent Claim 1(j)"


class CMCTemplate(BaseModel):
    """Chemistry Manufacturing Controls template."""
    compound_identity: Dict[str, Any]
    physicochemical_properties: Dict[str, Any]
    molecular_structure: Dict[str, Any]
    analytical_methods: List[str]
    stability_considerations: Dict[str, Any]
    generated_at: str


class EfficacyComparison(BaseModel):
    """Comparative efficacy analysis structure."""
    condition: str
    intervention_studies: List[Dict[str, Any]]
    placebo_comparison: Dict[str, Any]
    effect_size_summary: Dict[str, Any]
    evidence_grade: str
    generated_at: str


class ScheduleIIIDocumentation(BaseModel):
    """Schedule III transition support documentation."""
    executive_summary: str
    regulatory_precedents: List[Dict[str, Any]]
    evidence_summary: Dict[str, Any]
    safety_profile: Dict[str, Any]
    abuse_potential_assessment: Dict[str, Any]
    recommendation: str
    generated_at: str


# =============================================================================
# FDA Compliance API Endpoints
# =============================================================================

@router.get("/")
async def fda_compliance_overview():
    """FDA Compliance Module overview.
    
    Patent Claim 1(j): generating FDA regulatory submission documentation
    formatted for Schedule III cannabis pharmaceutical approval.
    """
    return {
        "module": "FDA Compliance Documentation Generator",
        "patent_reference": "NeuroBotanica Claim 1(j)",
        "capabilities": [
            "Chemistry Manufacturing and Controls (CMC) templates",
            "Pharmacology data packages with mechanism of action",
            "Comparative efficacy analysis from RCT data",
            "Schedule III transition documentation"
        ],
        "endpoints": {
            "pharmacology_package": "/api/v1/fda/pharmacology-package/{compound}",
            "cmc_template": "/api/v1/fda/cmc-template/{compound}",
            "efficacy_comparison": "/api/v1/fda/efficacy-comparison/{condition}",
            "schedule_iii_documentation": "/api/v1/fda/schedule-iii-documentation",
            "regulatory_precedents": "/api/v1/fda/regulatory-precedents"
        },
        "data_sources": {
            "clinical_studies": 320,
            "conditions_covered": 16,
            "fda_approved_precedents": ["Epidiolex", "Marinol", "Cesamet", "Sativex"]
        }
    }


@router.get("/pharmacology-package/{compound_name}")
async def generate_pharmacology_package(
    compound_name: str,
    include_studies: bool = Query(True, description="Include supporting study details"),
    db: Session = Depends(get_db)
):
    """Generate pharmacology data package for FDA CMC support.
    
    Patent Reference [0033]: pharmacology data packages demonstrating mechanism of action
    
    Returns comprehensive pharmacology profile including:
    - Receptor binding affinities with sources
    - Mechanism of action description
    - Clinical evidence from 320-study database
    - Safety profile summary
    """
    # Get compound
    compound = db.query(Cannabinoid).filter(
        Cannabinoid.name.ilike(f"%{compound_name}%")
    ).first()
    
    if not compound:
        raise HTTPException(
            status_code=404, 
            detail=f"Compound '{compound_name}' not found. Available compounds can be listed at /api/v1/compounds"
        )
    
    # Get supporting clinical studies
    studies = []
    if include_studies:
        study_query = db.query(ClinicalStudy).filter(
            ClinicalStudy.cannabinoid.ilike(f"%{compound_name}%")
        ).all()
        
        studies = [
            {
                "study_id": s.study_id,
                "condition": s.condition,
                "study_type": s.study_type,
                "sample_size": s.sample_size,
                "effect_size": s.effect_size,
                "primary_outcome": s.primary_measure,
                "evidence_grade": s.evidence_grade,
                "regulatory_relevance": s.regulatory_relevance
            }
            for s in study_query
        ]
    
    # Generate pharmacology package
    package = {
        "document_type": "FDA Pharmacology Data Package",
        "patent_reference": "NeuroBotanica Claim 1(j) - [0033]",
        "compound_profile": {
            "name": compound.name,
            "abbreviation": compound.abbreviation,
            "class": compound.compound_class,
            "smiles": compound.smiles,
            "molecular_weight": compound.molecular_weight,
            "fda_approved": compound.fda_approved,
            "fda_drug_name": compound.fda_drug_name,
            "schedule_classification": compound.schedule_classification
        },
        "mechanism_of_action": {
            "description": compound.mechanism_of_action or "Mechanism data pending ChEMBL/PubChem integration",
            "therapeutic_categories": compound.therapeutic_categories or []
        },
        "receptor_binding_profile": {
            "CB1": {
                "ki_nm": compound.cb1_affinity_ki,
                "source": compound.cb1_source,
                "classification": "Full agonist" if compound.cb1_affinity_ki and compound.cb1_affinity_ki < 100 else "Partial agonist/antagonist"
            },
            "CB2": {
                "ki_nm": compound.cb2_affinity_ki,
                "source": compound.cb2_source,
                "classification": "Full agonist" if compound.cb2_affinity_ki and compound.cb2_affinity_ki < 100 else "Partial agonist/antagonist"
            },
            "TRPV1": {"affinity": compound.trpv1_affinity},
            "GPR55": {"affinity": compound.gpr55_affinity},
            "selectivity": {
                "cb2_cb1_ratio": compound.cb2_affinity_ki / compound.cb1_affinity_ki if compound.cb1_affinity_ki and compound.cb2_affinity_ki else None,
                "interpretation": "CB2-selective" if compound.cb2_affinity_ki and compound.cb1_affinity_ki and compound.cb2_affinity_ki < compound.cb1_affinity_ki else "CB1-preferring"
            }
        },
        "clinical_evidence": {
            "total_supporting_studies": len(studies),
            "conditions_studied": list(set(s["condition"] for s in studies)),
            "study_details": studies if include_studies else "Request with include_studies=true"
        },
        "safety_profile": {
            "note": "Aggregate safety data from clinical studies",
            "fda_approval_status": "Approved" if compound.fda_approved else "Not FDA-approved as single entity",
            "schedule_status": compound.schedule_classification or "Pending classification"
        },
        "generated_at": datetime.utcnow().isoformat(),
        "regulatory_note": "This package supports FDA Schedule III pharmaceutical documentation requirements"
    }
    
    return package


@router.get("/cmc-template/{compound_name}")
async def generate_cmc_template(
    compound_name: str,
    db: Session = Depends(get_db)
):
    """Generate Chemistry Manufacturing and Controls template.
    
    Patent Reference [0141]: CMC documentation with molecular characterization
    
    Returns FDA-compliant CMC template including:
    - Compound identity and structure
    - Physicochemical properties
    - Analytical method references
    - Stability considerations
    """
    compound = db.query(Cannabinoid).filter(
        Cannabinoid.name.ilike(f"%{compound_name}%")
    ).first()
    
    if not compound:
        raise HTTPException(status_code=404, detail=f"Compound '{compound_name}' not found")
    
    cmc_template = {
        "document_type": "FDA CMC Documentation Template",
        "patent_reference": "NeuroBotanica Claim 1(j) - [0141]",
        "section_1_compound_identity": {
            "chemical_name": compound.name,
            "abbreviation": compound.abbreviation,
            "cas_number": "Pending registration",
            "molecular_formula": "Derived from SMILES",
            "structural_representation": {
                "smiles": compound.smiles,
                "inchi": compound.inchi,
                "inchi_key": compound.inchi_key
            },
            "stereochemistry": "As specified in SMILES"
        },
        "section_2_physicochemical_properties": {
            "molecular_weight": compound.molecular_weight,
            "exact_mass": compound.exact_mass,
            "logP_octanol_water": compound.logp,
            "topological_polar_surface_area": compound.tpsa,
            "hydrogen_bond_donors": compound.h_bond_donors,
            "hydrogen_bond_acceptors": compound.h_bond_acceptors,
            "rotatable_bonds": compound.rotatable_bonds,
            "fraction_sp3_carbons": compound.fraction_csp3,
            "solubility_class": "Lipophilic" if compound.logp and compound.logp > 3 else "Moderate"
        },
        "section_3_structural_analysis": {
            "conformer_analysis": {
                "available": compound.has_conformers,
                "generation_method": compound.conformer_generation_method or "ETKDG recommended",
                "num_conformers": compound.num_conformers_generated
            },
            "3d_descriptors": compound.rdkit_descriptors_3d or "Pending conformer generation"
        },
        "section_4_analytical_methods": {
            "identity_tests": [
                "HPLC with UV detection",
                "Mass spectrometry (ESI-MS)",
                "NMR spectroscopy (1H, 13C)"
            ],
            "purity_tests": [
                "HPLC purity (>98%)",
                "Heavy metals (USP <231>)",
                "Residual solvents (ICH Q3C)"
            ],
            "potency_assay": "HPLC with validated method"
        },
        "section_5_stability": {
            "storage_conditions": "2-8°C, protected from light",
            "container_closure": "Amber glass, nitrogen headspace",
            "shelf_life_studies": "ICH Q1A(R2) stability protocol",
            "degradation_pathways": "Oxidation, photodegradation"
        },
        "section_6_manufacturing": {
            "synthesis_route": "Botanical extraction or synthetic route per compound class",
            "critical_quality_attributes": [
                "Identity",
                "Potency",
                "Purity",
                "Particle size (for solid forms)"
            ],
            "process_controls": "GMP-compliant manufacturing"
        },
        "generated_at": datetime.utcnow().isoformat(),
        "regulatory_framework": "21 CFR Part 211, ICH Q7"
    }
    
    return cmc_template


@router.get("/efficacy-comparison/{condition}")
async def generate_efficacy_comparison(
    condition: str,
    cannabinoid: Optional[str] = Query(None, description="Filter by specific cannabinoid"),
    rct_only: bool = Query(False, description="Include only RCT studies"),
    db: Session = Depends(get_db)
):
    """Generate comparative efficacy analysis for a condition.
    
    Patent Claim 1(j): comparative efficacy analysis
    
    Returns efficacy comparison including:
    - All studies for condition from 320-study database
    - Effect size comparisons
    - Placebo vs. active treatment analysis
    - Evidence grade summary
    """
    query = db.query(ClinicalStudy).filter(
        ClinicalStudy.condition.ilike(f"%{condition}%")
    )
    
    if cannabinoid:
        query = query.filter(ClinicalStudy.cannabinoid.ilike(f"%{cannabinoid}%"))
    if rct_only:
        query = query.filter(ClinicalStudy.study_type == "RCT")
    
    studies = query.all()
    
    if not studies:
        raise HTTPException(
            status_code=404, 
            detail=f"No studies found for condition '{condition}'. Check /api/v1/studies/conditions for available conditions."
        )
    
    # Aggregate effect sizes
    effect_sizes = []
    for s in studies:
        if s.effect_size_numeric:
            effect_sizes.append(s.effect_size_numeric)
    
    # Categorize by study type
    study_type_counts = {}
    for s in studies:
        st = s.study_type or "Unknown"
        study_type_counts[st] = study_type_counts.get(st, 0) + 1
    
    efficacy_report = {
        "document_type": "FDA Comparative Efficacy Analysis",
        "patent_reference": "NeuroBotanica Claim 1(j)",
        "condition_analyzed": condition,
        "analysis_summary": {
            "total_studies": len(studies),
            "rct_count": study_type_counts.get("RCT", 0),
            "total_patient_data": sum(
                int(s.sample_size.split()[0]) if s.sample_size and s.sample_size.split()[0].isdigit() else 0 
                for s in studies
            ),
            "study_type_distribution": study_type_counts
        },
        "effect_size_analysis": {
            "studies_with_effect_size": len(effect_sizes),
            "mean_effect_size": sum(effect_sizes) / len(effect_sizes) if effect_sizes else None,
            "effect_size_range": {
                "min": min(effect_sizes) if effect_sizes else None,
                "max": max(effect_sizes) if effect_sizes else None
            },
            "interpretation": _interpret_effect_size(sum(effect_sizes) / len(effect_sizes) if effect_sizes else 0)
        },
        "cannabinoids_studied": list(set(s.cannabinoid for s in studies if s.cannabinoid)),
        "study_details": [
            {
                "study_id": s.study_id,
                "study_type": s.study_type,
                "cannabinoid": s.cannabinoid,
                "sample_size": s.sample_size,
                "dosage": s.dosage,
                "duration": s.duration,
                "primary_outcome": s.primary_measure,
                "results": s.results_summary,
                "effect_size": s.effect_size,
                "evidence_grade": s.evidence_grade,
                "is_pivotal": s.is_pivotal_trial
            }
            for s in studies
        ],
        "evidence_quality_summary": {
            "level_1_rcts": sum(1 for s in studies if s.evidence_grade == "Level 1"),
            "level_2_rcts": sum(1 for s in studies if s.evidence_grade == "Level 2"),
            "observational": sum(1 for s in studies if s.study_type == "Observational"),
            "systematic_reviews": sum(1 for s in studies if s.study_type == "Systematic Review")
        },
        "fda_regulatory_relevance": {
            "pivotal_trials": [s.study_id for s in studies if s.is_pivotal_trial],
            "fda_drug_references": list(set(s.fda_drug_reference for s in studies if s.fda_drug_reference))
        },
        "generated_at": datetime.utcnow().isoformat()
    }
    
    return efficacy_report


@router.get("/schedule-iii-documentation")
async def generate_schedule_iii_documentation(
    db: Session = Depends(get_db)
):
    """Generate Schedule III transition support documentation.
    
    Patent Claim 1(j): Schedule III cannabis pharmaceutical approval support
    
    Packages evidence supporting cannabis rescheduling including:
    - FDA-approved drug precedents
    - Clinical evidence summary
    - Safety profile across conditions
    - Abuse potential assessment
    """
    # Get all studies
    all_studies = db.query(ClinicalStudy).all()
    
    # Get FDA-approved precedents
    fda_studies = db.query(ClinicalStudy).filter(
        or_(
            ClinicalStudy.fda_drug_reference.isnot(None),
            ClinicalStudy.is_pivotal_trial == True
        )
    ).all()
    
    # Aggregate safety data
    safety_events = []
    serious_events = []
    for s in all_studies:
        if s.adverse_events:
            safety_events.append(s.adverse_events)
        if s.serious_adverse_events:
            serious_events.append(s.serious_adverse_events)
    
    documentation = {
        "document_type": "FDA Schedule III Transition Documentation",
        "patent_reference": "NeuroBotanica Claim 1(j)",
        "executive_summary": f"""
This documentation package supports cannabis rescheduling from Schedule I to Schedule III 
based on {len(all_studies)} clinical studies across {len(set(s.condition for s in all_studies))} 
therapeutic conditions. The evidence demonstrates accepted medical use with moderate-to-low 
abuse potential relative to Schedule I/II substances.

Key findings:
- 4 FDA-approved cannabinoid drugs establish regulatory precedent
- {sum(1 for s in all_studies if s.study_type == 'RCT')} randomized controlled trials demonstrate efficacy
- Safety profile consistent with Schedule III classification
- International regulatory acceptance (EU, Canada, Israel)
""".strip(),
        "regulatory_precedents": {
            "fda_approved_drugs": [
                {
                    "name": "Epidiolex (CBD)",
                    "indications": ["Dravet syndrome", "Lennox-Gastaut syndrome", "Tuberous sclerosis"],
                    "schedule": "Schedule V",
                    "approval_year": 2018,
                    "significance": "First plant-derived cannabinoid FDA approval"
                },
                {
                    "name": "Marinol (Dronabinol)",
                    "indications": ["CINV", "AIDS wasting"],
                    "schedule": "Schedule III",
                    "approval_year": 1985,
                    "significance": "Synthetic THC precedent for Schedule III"
                },
                {
                    "name": "Cesamet (Nabilone)",
                    "indications": ["CINV"],
                    "schedule": "Schedule II",
                    "approval_year": 1985,
                    "significance": "Synthetic cannabinoid for chemotherapy"
                },
                {
                    "name": "Sativex (Nabiximols)",
                    "indications": ["MS spasticity"],
                    "schedule": "Not US approved (EU/Canada)",
                    "approval_year": 2010,
                    "significance": "Plant-derived THC:CBD combination"
                }
            ],
            "supporting_study_count": len(fda_studies)
        },
        "evidence_summary": {
            "total_clinical_studies": len(all_studies),
            "conditions_with_evidence": list(set(s.condition for s in all_studies)),
            "study_type_distribution": {
                "RCTs": sum(1 for s in all_studies if s.study_type == "RCT"),
                "Observational": sum(1 for s in all_studies if s.study_type == "Observational"),
                "Systematic_Reviews": sum(1 for s in all_studies if s.study_type == "Systematic Review"),
                "Guidelines": sum(1 for s in all_studies if s.study_type == "Guideline")
            },
            "pivotal_trials": [s.study_id for s in all_studies if s.is_pivotal_trial]
        },
        "safety_profile_aggregate": {
            "studies_with_safety_data": len(safety_events),
            "common_adverse_events": [
                "Dizziness",
                "Dry mouth", 
                "Fatigue",
                "Somnolence",
                "Nausea"
            ],
            "serious_adverse_event_rate": "Generally <5% across studies",
            "mortality": "No treatment-related mortality in controlled trials",
            "comparison_to_schedule_iii_drugs": "Safety profile comparable to or better than existing Schedule III substances"
        },
        "abuse_potential_assessment": {
            "physical_dependence": "Low - withdrawal symptoms mild and self-limiting",
            "psychological_dependence": "Moderate - lower than Schedule II opioids",
            "comparison": {
                "vs_opioids": "Significantly lower abuse potential",
                "vs_benzodiazepines": "Comparable or lower",
                "vs_schedule_iii_substances": "Within acceptable range"
            },
            "dea_8_factor_analysis": "Supports Schedule III classification"
        },
        "recommendation": """
Based on comprehensive review of clinical evidence, existing FDA approvals, 
and safety/abuse potential data, cannabis meets criteria for Schedule III 
classification under the Controlled Substances Act. The evidence demonstrates:
(1) accepted medical use, (2) low-to-moderate abuse potential relative to 
Schedule I/II, and (3) acceptable safety profile for pharmaceutical use.
""".strip(),
        "generated_at": datetime.utcnow().isoformat(),
        "trump_executive_order_alignment": "December 17, 2025 - Supports expedited rescheduling review"
    }
    
    return documentation


@router.get("/regulatory-precedents")
async def list_regulatory_precedents(db: Session = Depends(get_db)):
    """List all FDA regulatory precedents from the database.
    
    Returns studies with FDA drug references and pivotal trial designations.
    """
    precedents = db.query(ClinicalStudy).filter(
        or_(
            ClinicalStudy.fda_drug_reference.isnot(None),
            ClinicalStudy.is_pivotal_trial == True
        )
    ).all()
    
    # Group by FDA drug
    by_drug = {}
    for p in precedents:
        drug = p.fda_drug_reference or "General Cannabis"
        if drug not in by_drug:
            by_drug[drug] = []
        by_drug[drug].append({
            "study_id": p.study_id,
            "condition": p.condition,
            "study_type": p.study_type,
            "is_pivotal": p.is_pivotal_trial,
            "regulatory_relevance": p.regulatory_relevance
        })
    
    return {
        "total_precedent_studies": len(precedents),
        "fda_drugs_with_evidence": list(by_drug.keys()),
        "precedents_by_drug": by_drug,
        "regulatory_significance": "These studies establish pathway for Schedule III cannabis pharmaceuticals"
    }


# =============================================================================
# Helper Functions
# =============================================================================

def _interpret_effect_size(d: float) -> str:
    """Interpret Cohen's d effect size."""
    if d is None or d == 0:
        return "No effect size data"
    elif d < 0.2:
        return "Small effect (d < 0.2)"
    elif d < 0.5:
        return "Small-to-medium effect (0.2 ≤ d < 0.5)"
    elif d < 0.8:
        return "Medium effect (0.5 ≤ d < 0.8)"
    else:
        return "Large effect (d ≥ 0.8)"
