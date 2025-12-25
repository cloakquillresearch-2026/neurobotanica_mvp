"""
Adjuvant Optimization Service - NeuroBotanica

Implements patent claims 14-17: Adjuvant Enhancement for therapeutic optimization.

Features:
- Receptor priming score calculation
- Signal-to-noise optimization
- Temporal phasing (optimal timing)
- Drug interaction checking
- Evidence-tiered recommendations

Reference: NeuroBotanica Patent - Section 3: Adjuvant Enhancement
"""

from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
import logging

logger = logging.getLogger(__name__)


class TherapeuticTarget(Enum):
    """Supported therapeutic targets for adjuvant optimization."""
    CHRONIC_PAIN = "chronic_pain"
    ANXIETY = "anxiety"
    INSOMNIA = "insomnia"
    INFLAMMATION = "inflammation"
    NAUSEA = "nausea"
    EPILEPSY = "epilepsy"
    PTSD = "ptsd"
    NEUROPROTECTION = "neuroprotection"
    APPETITE = "appetite"
    MIGRAINE = "migraine"
    METABOLIC_SUPPORT = "metabolic_support"
    WEIGHT_MANAGEMENT = "weight_management"


@dataclass
class PatientProfile:
    """Patient profile for adjuvant optimization."""
    age: int
    weight_kg: float
    conditions: List[str] = field(default_factory=list)
    current_medications: List[str] = field(default_factory=list)
    supplements: List[str] = field(default_factory=list)
    allergies: List[str] = field(default_factory=list)
    metabolizer_status: str = "normal"  # slow, normal, fast
    experience_level: str = "regular"  # naive, occasional, regular, experienced


@dataclass
class AdjuvantRecommendation:
    """Recommended adjuvant compound with dosing details."""
    adjuvant_name: str
    dosage_mg: float
    timing_offset_minutes: int
    expected_enhancement_percent: float
    confidence_interval: Tuple[float, float]
    mechanism: str
    evidence_tier: int
    citations: List[str] = field(default_factory=list)
    contraindications: List[str] = field(default_factory=list)
    synergy_note: Optional[str] = None


@dataclass
class OptimizationResult:
    """Complete adjuvant optimization result."""
    primary_compound: str
    therapeutic_target: str
    recommendations: List[AdjuvantRecommendation]
    total_expected_enhancement: float
    protocol_summary: str
    warnings: List[str] = field(default_factory=list)
    created_at: datetime = field(default_factory=datetime.utcnow)


class AdjuvantOptimizer:
    """
    Adjuvant optimization engine for cannabinoid therapeutic enhancement.
    
    Implements the magnesium glycinate + cannabinoid synergy model
    and extends to other adjuvant-cannabinoid combinations.
    """
    
    # Evidence-based adjuvant database
    # Evidence tiers: 1=RCT/Meta, 2=Clinical, 3=In vivo, 4=In vitro, 5=Theoretical
    ADJUVANT_DATABASE: Dict[str, Dict[str, Any]] = {
        "magnesium_glycinate": {
            "name": "Magnesium Glycinate",
            "category": "Receptor Priming",
            "mechanism": "NMDA antagonism + GlyR potentiation",
            "target_receptors": ["NMDA", "GlyR"],
            "recommended_dose_mg": 400,
            "dose_range": (200, 600),
            "timing_offset_minutes": -45,
            "timing_flexibility": 15,
            "evidence_tier": 3,
            "citations": ["PMID:25765142", "PMID:28394950"],
            "therapeutic_targets": ["insomnia", "chronic_pain", "anxiety", "migraine"],
            "enhancement_by_target": {
                "insomnia": 50,
                "chronic_pain": 35,
                "anxiety": 30,
                "migraine": 40,
            },
            "contraindications": ["kidney_disease", "myasthenia_gravis"],
            "drug_interactions": ["bisphosphonates", "antibiotics_quinolone"],
            "is_otc": True,
        },
        "l_theanine": {
            "name": "L-Theanine",
            "category": "Synergistic",
            "mechanism": "GABA enhancement + glutamate modulation",
            "target_receptors": ["GABA-A", "NMDA"],
            "recommended_dose_mg": 200,
            "dose_range": (100, 400),
            "timing_offset_minutes": -30,
            "timing_flexibility": 15,
            "evidence_tier": 2,
            "citations": ["PMID:16930802", "PMID:21040626"],
            "therapeutic_targets": ["anxiety", "insomnia", "neuroprotection"],
            "enhancement_by_target": {
                "anxiety": 40,
                "insomnia": 35,
                "neuroprotection": 25,
            },
            "contraindications": [],
            "drug_interactions": ["blood_pressure_medications"],
            "is_otc": True,
        },
        "omega_3_fatty_acids": {
            "name": "Omega-3 Fatty Acids (EPA/DHA)",
            "category": "Bioavailability",
            "mechanism": "Lipid carrier enhancement + anti-inflammatory",
            "target_receptors": ["PPAR-gamma", "CB1"],
            "recommended_dose_mg": 2000,
            "dose_range": (1000, 4000),
            "timing_offset_minutes": 0,
            "timing_flexibility": 60,
            "evidence_tier": 2,
            "citations": ["PMID:26068424", "PMID:29387462"],
            "therapeutic_targets": ["inflammation", "neuroprotection", "chronic_pain"],
            "enhancement_by_target": {
                "inflammation": 45,
                "neuroprotection": 40,
                "chronic_pain": 30,
            },
            "contraindications": ["fish_allergy"],
            "drug_interactions": ["blood_thinners"],
            "is_otc": True,
        },
        "curcumin": {
            "name": "Curcumin (with Piperine)",
            "category": "Synergistic",
            "mechanism": "COX-2 inhibition + NF-κB modulation",
            "target_receptors": ["COX-2", "NF-κB"],
            "recommended_dose_mg": 500,
            "dose_range": (250, 1000),
            "timing_offset_minutes": 0,
            "timing_flexibility": 30,
            "evidence_tier": 2,
            "citations": ["PMID:17569207", "PMID:27213821"],
            "therapeutic_targets": ["inflammation", "chronic_pain", "neuroprotection"],
            "enhancement_by_target": {
                "inflammation": 55,
                "chronic_pain": 40,
                "neuroprotection": 35,
            },
            "contraindications": ["gallbladder_disease"],
            "drug_interactions": ["blood_thinners", "diabetes_medications"],
            "is_otc": True,
            "bioavailability_note": "Take with black pepper (piperine) for 2000% absorption increase",
        },
        "vitamin_d3": {
            "name": "Vitamin D3",
            "category": "Receptor Priming",
            "mechanism": "CB1/CB2 receptor expression enhancement",
            "target_receptors": ["CB1", "CB2", "VDR"],
            "recommended_dose_mg": 50,  # 2000 IU
            "dose_range": (25, 125),  # 1000-5000 IU
            "timing_offset_minutes": 0,
            "timing_flexibility": 120,
            "evidence_tier": 3,
            "citations": ["PMID:28648215", "PMID:27030531"],
            "therapeutic_targets": ["chronic_pain", "inflammation", "neuroprotection"],
            "enhancement_by_target": {
                "chronic_pain": 25,
                "inflammation": 30,
                "neuroprotection": 35,
            },
            "contraindications": ["hypercalcemia", "kidney_stones"],
            "drug_interactions": ["thiazide_diuretics", "steroids"],
            "is_otc": True,
        },
        "glycine": {
            "name": "Glycine",
            "category": "Receptor Priming",
            "mechanism": "Glycine receptor agonism + NMDA co-agonism",
            "target_receptors": ["GlyR", "NMDA"],
            "recommended_dose_mg": 3000,
            "dose_range": (2000, 5000),
            "timing_offset_minutes": -30,
            "timing_flexibility": 15,
            "evidence_tier": 3,
            "citations": ["PMID:25533534", "PMID:22293292"],
            "therapeutic_targets": ["insomnia", "anxiety", "neuroprotection"],
            "enhancement_by_target": {
                "insomnia": 45,
                "anxiety": 25,
                "neuroprotection": 30,
            },
            "contraindications": [],
            "drug_interactions": ["clozapine"],
            "is_otc": True,
        },
        "black_pepper_extract": {
            "name": "Black Pepper Extract (Piperine)",
            "category": "Bioavailability",
            "mechanism": "CYP450 inhibition → increased cannabinoid bioavailability",
            "target_receptors": ["CYP3A4", "CYP2C9"],
            "recommended_dose_mg": 20,
            "dose_range": (5, 30),
            "timing_offset_minutes": 0,
            "timing_flexibility": 15,
            "evidence_tier": 3,
            "citations": ["PMID:9619120", "PMID:20578802"],
            "therapeutic_targets": ["all"],  # Universal bioavailability enhancer
            "enhancement_by_target": {
                "chronic_pain": 30,
                "anxiety": 25,
                "inflammation": 30,
                "insomnia": 25,
            },
            "contraindications": [],
            "drug_interactions": ["many_medications"],  # CYP450 inhibitor
            "is_otc": True,
            "caution": "May increase blood levels of many medications",
        },
        "palmitoylethanolamide": {
            "name": "Palmitoylethanolamide (PEA)",
            "category": "Synergistic",
            "mechanism": "PPAR-alpha agonism + FAAH inhibition (entourage effect)",
            "target_receptors": ["PPAR-alpha", "CB1", "CB2"],
            "recommended_dose_mg": 600,
            "dose_range": (300, 1200),
            "timing_offset_minutes": 0,
            "timing_flexibility": 30,
            "evidence_tier": 2,
            "citations": ["PMID:23850343", "PMID:26068424"],
            "therapeutic_targets": ["chronic_pain", "inflammation", "neuroprotection"],
            "enhancement_by_target": {
                "chronic_pain": 50,
                "inflammation": 45,
                "neuroprotection": 40,
            },
            "contraindications": [],
            "drug_interactions": [],
            "is_otc": True,
        },
        "melatonin": {
            "name": "Melatonin",
            "category": "Synergistic",
            "mechanism": "MT1/MT2 receptor activation + antioxidant",
            "target_receptors": ["MT1", "MT2"],
            "recommended_dose_mg": 3,
            "dose_range": (0.5, 10),
            "timing_offset_minutes": -60,
            "timing_flexibility": 30,
            "evidence_tier": 1,
            "citations": ["PMID:28648359", "PMID:23814343"],
            "therapeutic_targets": ["insomnia"],
            "enhancement_by_target": {
                "insomnia": 60,
            },
            "contraindications": ["autoimmune_disease", "depression"],
            "drug_interactions": ["immunosuppressants", "blood_thinners", "diabetes_medications"],
            "is_otc": True,
        },
        "nac": {
            "name": "N-Acetyl Cysteine (NAC)",
            "category": "Neuroprotection",
            "mechanism": "Glutathione precursor + glutamate modulation + reduces oxidative stress",
            "target_receptors": ["NMDA", "mGluR5"],
            "recommended_dose_mg": 600,
            "dose_range": (600, 1800),
            "timing_offset_minutes": -60,
            "timing_flexibility": 30,
            "evidence_tier": 2,
            "citations": ["PMID:28011714", "PMID:27137659", "PMID:23369637"],
            "therapeutic_targets": ["ptsd", "depression", "neuroprotection", "anxiety"],
            "enhancement_by_target": {
                "ptsd": 45,
                "depression": 40,
                "neuroprotection": 50,
                "anxiety": 35,
            },
            "contraindications": ["asthma"],
            "drug_interactions": ["nitroglycerin", "activated_charcoal"],
            "is_otc": True,
            "special_note": "May reduce cannabis tolerance and dependence",
        },
        "taurine": {
            "name": "Taurine",
            "category": "Synergistic",
            "mechanism": "GABA-A receptor modulation + glycine receptor agonism + neuroprotection",
            "target_receptors": ["GABA-A", "GlyR"],
            "recommended_dose_mg": 1000,
            "dose_range": (500, 3000),
            "timing_offset_minutes": -30,
            "timing_flexibility": 30,
            "evidence_tier": 3,
            "citations": ["PMID:24126426", "PMID:23170060"],
            "therapeutic_targets": ["anxiety", "epilepsy", "neuroprotection"],
            "enhancement_by_target": {
                "anxiety": 35,
                "epilepsy": 40,
                "neuroprotection": 30,
            },
            "contraindications": [],
            "drug_interactions": ["lithium"],
            "is_otc": True,
        },
        "ashwagandha": {
            "name": "Ashwagandha (KSM-66)",
            "category": "Adaptogen",
            "mechanism": "Cortisol reduction + GABA mimetic + anxiolytic + withanolides",
            "target_receptors": ["GABA-A", "5-HT1A", "cortisol"],
            "recommended_dose_mg": 300,
            "dose_range": (300, 600),
            "timing_offset_minutes": 0,
            "timing_flexibility": 60,
            "evidence_tier": 2,
            "citations": ["PMID:23439798", "PMID:32021735", "PMID:28471731"],
            "therapeutic_targets": ["anxiety", "chronic_pain", "insomnia", "depression"],
            "enhancement_by_target": {
                "anxiety": 50,
                "chronic_pain": 30,
                "insomnia": 40,
                "depression": 35,
            },
            "contraindications": ["thyroid_disorders", "pregnancy"],
            "drug_interactions": ["thyroid_medications", "immunosuppressants", "sedatives"],
            "is_otc": True,
        },
        "coq10": {
            "name": "Coenzyme Q10 (Ubiquinol)",
            "category": "Mitochondrial Support",
            "mechanism": "Mitochondrial electron transport + antioxidant + energy production",
            "target_receptors": ["mitochondrial_complex_I_III"],
            "recommended_dose_mg": 200,
            "dose_range": (100, 400),
            "timing_offset_minutes": 0,
            "timing_flexibility": 120,
            "evidence_tier": 2,
            "citations": ["PMID:29225069", "PMID:26792398"],
            "therapeutic_targets": ["migraine", "neuroprotection", "chronic_pain"],
            "enhancement_by_target": {
                "migraine": 50,
                "neuroprotection": 40,
                "chronic_pain": 25,
            },
            "contraindications": [],
            "drug_interactions": ["blood_thinners", "chemotherapy"],
            "is_otc": True,
            "bioavailability_note": "Ubiquinol form is 2-8x more bioavailable than ubiquinone",
        },
        "alpha_lipoic_acid": {
            "name": "Alpha-Lipoic Acid (ALA)",
            "category": "Neuroprotection",
            "mechanism": "Antioxidant + regenerates other antioxidants + nerve growth factor",
            "target_receptors": ["NRF2", "insulin_receptor"],
            "recommended_dose_mg": 600,
            "dose_range": (300, 1200),
            "timing_offset_minutes": -30,
            "timing_flexibility": 30,
            "evidence_tier": 2,
            "citations": ["PMID:22891897", "PMID:16464129"],
            "therapeutic_targets": ["chronic_pain", "neuroprotection", "inflammation"],
            "enhancement_by_target": {
                "chronic_pain": 45,
                "neuroprotection": 50,
                "inflammation": 35,
            },
            "contraindications": ["thiamine_deficiency"],
            "drug_interactions": ["diabetes_medications", "thyroid_medications"],
            "is_otc": True,
            "special_note": "R-ALA form is more bioavailable than racemic ALA",
        },
        "phosphatidylserine": {
            "name": "Phosphatidylserine (PS)",
            "category": "Cognitive Support",
            "mechanism": "Cortisol modulation + cell membrane support + neurotransmitter release",
            "target_receptors": ["cortisol", "dopamine", "acetylcholine"],
            "recommended_dose_mg": 100,
            "dose_range": (100, 300),
            "timing_offset_minutes": 0,
            "timing_flexibility": 60,
            "evidence_tier": 2,
            "citations": ["PMID:24449470", "PMID:21103034"],
            "therapeutic_targets": ["ptsd", "neuroprotection", "anxiety", "depression"],
            "enhancement_by_target": {
                "ptsd": 40,
                "neuroprotection": 45,
                "anxiety": 35,
                "depression": 30,
            },
            "contraindications": ["blood_clotting_disorders"],
            "drug_interactions": ["blood_thinners", "anticholinergics"],
            "is_otc": True,
            "special_note": "Blunts cortisol response to stress",
        },
        "berberine": {
            "name": "Berberine",
            "category": "Metabolic Support",
            "mechanism": "AMPK activation + improves insulin sensitivity + GLP-1 pathway modulation",
            "target_receptors": ["AMPK", "insulin_receptor", "GLP-1"],
            "recommended_dose_mg": 500,
            "dose_range": (500, 1500),
            "timing_offset_minutes": -30,
            "timing_flexibility": 30,
            "evidence_tier": 1,
            "citations": ["PMID:18442638", "PMID:23118793", "PMID:32268516"],
            "therapeutic_targets": ["metabolic_support", "weight_management", "inflammation"],
            "enhancement_by_target": {
                "metabolic_support": 55,
                "weight_management": 50,
                "inflammation": 35,
            },
            "contraindications": ["pregnancy", "breastfeeding"],
            "drug_interactions": ["metformin", "diabetes_medications", "cyclosporine"],
            "is_otc": True,
            "special_note": "Often called 'nature's metformin' - comparable efficacy to prescription diabetes drugs",
        },
        "inositol": {
            "name": "Myo-Inositol",
            "category": "Metabolic Support",
            "mechanism": "Insulin signaling enhancement + PCOS support + mood regulation",
            "target_receptors": ["insulin_receptor", "5-HT2A"],
            "recommended_dose_mg": 2000,
            "dose_range": (2000, 4000),
            "timing_offset_minutes": 0,
            "timing_flexibility": 60,
            "evidence_tier": 1,
            "citations": ["PMID:22774396", "PMID:29498933", "PMID:27568069"],
            "therapeutic_targets": ["metabolic_support", "anxiety", "weight_management"],
            "enhancement_by_target": {
                "metabolic_support": 50,
                "anxiety": 40,
                "weight_management": 45,
            },
            "contraindications": [],
            "drug_interactions": [],
            "is_otc": True,
            "special_note": "Especially effective for PCOS and metabolic syndrome",
        },
        "chromium_picolinate": {
            "name": "Chromium Picolinate",
            "category": "Metabolic Support",
            "mechanism": "Enhances insulin receptor sensitivity + glucose transport + reduces cravings",
            "target_receptors": ["insulin_receptor", "GLUT4"],
            "recommended_dose_mg": 200,
            "dose_range": (200, 1000),
            "timing_offset_minutes": 0,
            "timing_flexibility": 60,
            "evidence_tier": 2,
            "citations": ["PMID:17109600", "PMID:25439135"],
            "therapeutic_targets": ["metabolic_support", "weight_management"],
            "enhancement_by_target": {
                "metabolic_support": 35,
                "weight_management": 40,
            },
            "contraindications": ["kidney_disease"],
            "drug_interactions": ["insulin", "diabetes_medications", "thyroid_medications"],
            "is_otc": True,
            "special_note": "May reduce carbohydrate cravings",
        },
        "egcg": {
            "name": "Green Tea Extract (EGCG)",
            "category": "Metabolic Support",
            "mechanism": "Thermogenesis + fat oxidation + catechin-mediated metabolism boost",
            "target_receptors": ["AMPK", "catechol-O-methyltransferase"],
            "recommended_dose_mg": 400,
            "dose_range": (250, 500),
            "timing_offset_minutes": -30,
            "timing_flexibility": 30,
            "evidence_tier": 1,
            "citations": ["PMID:19597519", "PMID:20156466", "PMID:21115335"],
            "therapeutic_targets": ["metabolic_support", "weight_management", "neuroprotection"],
            "enhancement_by_target": {
                "metabolic_support": 45,
                "weight_management": 50,
                "neuroprotection": 30,
            },
            "contraindications": ["liver_disease"],
            "drug_interactions": ["blood_thinners", "stimulants", "beta_blockers"],
            "is_otc": True,
            "special_note": "Contains caffeine - take earlier in day; decaf versions available",
        },
    }
    
    # Cannabinoid-specific adjuvant preferences
    CANNABINOID_ADJUVANT_SYNERGIES: Dict[str, List[str]] = {
        "CBD": ["magnesium_glycinate", "l_theanine", "omega_3_fatty_acids", "curcumin", "ashwagandha"],
        "THC": ["black_pepper_extract", "vitamin_d3", "omega_3_fatty_acids", "nac"],
        "CBN": ["magnesium_glycinate", "glycine", "melatonin"],
        "CBG": ["curcumin", "palmitoylethanolamide", "omega_3_fatty_acids", "coq10"],
        "CBC": ["omega_3_fatty_acids", "curcumin", "alpha_lipoic_acid"],
        "THCA": ["curcumin", "palmitoylethanolamide"],
        "CBDA": ["curcumin", "l_theanine", "taurine"],
        "THCV": ["berberine", "egcg", "chromium_picolinate", "inositol"],
    }
    
    def __init__(self):
        self.logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")
    
    def optimize(
        self,
        primary_compound: str,
        therapeutic_target: str,
        patient_profile: Optional[PatientProfile] = None,
        max_adjuvants: int = 3,
    ) -> OptimizationResult:
        """
        Optimize adjuvant selection for a cannabinoid compound.
        
        Args:
            primary_compound: Primary cannabinoid (e.g., "CBD", "THC", dimer name)
            therapeutic_target: Target condition (e.g., "insomnia", "chronic_pain")
            patient_profile: Optional patient-specific data for personalization
            max_adjuvants: Maximum number of adjuvants to recommend
            
        Returns:
            OptimizationResult with ranked adjuvant recommendations
        """
        self.logger.info(f"Optimizing adjuvants for {primary_compound} targeting {therapeutic_target}")
        
        # Normalize inputs
        target = therapeutic_target.lower().replace(" ", "_")
        compound = primary_compound.upper()
        
        # Get suitable adjuvants for this target
        candidates = self._get_candidate_adjuvants(target, compound)
        
        # Score and rank candidates
        scored_candidates = []
        for adj_key, adj_data in candidates.items():
            score = self._score_adjuvant(adj_key, adj_data, target, compound, patient_profile)
            scored_candidates.append((adj_key, adj_data, score))
        
        # Sort by score and take top N
        scored_candidates.sort(key=lambda x: x[2], reverse=True)
        top_candidates = scored_candidates[:max_adjuvants]
        
        # Build recommendations
        recommendations = []
        warnings = []
        
        for adj_key, adj_data, score in top_candidates:
            # Check contraindications
            if patient_profile:
                contraindication_warnings = self._check_contraindications(
                    adj_data, patient_profile
                )
                if contraindication_warnings:
                    warnings.extend(contraindication_warnings)
                    continue  # Skip this adjuvant
            
            # Calculate personalized dose
            dose = self._calculate_dose(adj_data, patient_profile)
            
            # Get enhancement for this target
            enhancement = adj_data["enhancement_by_target"].get(target, 20)
            
            # Confidence interval based on evidence tier
            ci_range = self._get_confidence_interval(enhancement, adj_data["evidence_tier"])
            
            rec = AdjuvantRecommendation(
                adjuvant_name=adj_data["name"],
                dosage_mg=dose,
                timing_offset_minutes=adj_data["timing_offset_minutes"],
                expected_enhancement_percent=enhancement,
                confidence_interval=ci_range,
                mechanism=adj_data["mechanism"],
                evidence_tier=adj_data["evidence_tier"],
                citations=adj_data.get("citations", []),
                contraindications=adj_data.get("contraindications", []),
                synergy_note=adj_data.get("bioavailability_note"),
            )
            recommendations.append(rec)
        
        # Calculate total expected enhancement (diminishing returns)
        total_enhancement = self._calculate_total_enhancement(recommendations)
        
        # Generate protocol summary
        protocol = self._generate_protocol_summary(primary_compound, recommendations)
        
        return OptimizationResult(
            primary_compound=primary_compound,
            therapeutic_target=therapeutic_target,
            recommendations=recommendations,
            total_expected_enhancement=total_enhancement,
            protocol_summary=protocol,
            warnings=warnings,
        )
    
    def _get_candidate_adjuvants(
        self, target: str, compound: str
    ) -> Dict[str, Dict[str, Any]]:
        """Get adjuvants suitable for the target and compound."""
        candidates = {}
        
        # Get compound-specific synergies
        synergy_list = self.CANNABINOID_ADJUVANT_SYNERGIES.get(
            compound, list(self.ADJUVANT_DATABASE.keys())
        )
        
        for adj_key, adj_data in self.ADJUVANT_DATABASE.items():
            targets = adj_data["therapeutic_targets"]
            
            # Check if adjuvant is suitable for target
            if target in targets or "all" in targets:
                # Prioritize compound-specific synergies
                if adj_key in synergy_list:
                    candidates[adj_key] = adj_data
                elif len(candidates) < 10:  # Include some general options
                    candidates[adj_key] = adj_data
        
        return candidates
    
    def _score_adjuvant(
        self,
        adj_key: str,
        adj_data: Dict[str, Any],
        target: str,
        compound: str,
        patient_profile: Optional[PatientProfile],
    ) -> float:
        """Score an adjuvant based on relevance and evidence."""
        score = 0.0
        
        # Base enhancement for target (0-60 points)
        enhancement = adj_data["enhancement_by_target"].get(target, 0)
        score += enhancement
        
        # Evidence tier bonus (higher tier = more points, 0-25 points)
        evidence_score = (6 - adj_data["evidence_tier"]) * 5
        score += evidence_score
        
        # Compound synergy bonus (0-15 points)
        synergy_list = self.CANNABINOID_ADJUVANT_SYNERGIES.get(compound, [])
        if adj_key in synergy_list:
            score += 15
        
        # Patient profile adjustments
        if patient_profile:
            # Already taking this supplement (avoid duplication)
            if adj_data["name"].lower() in [s.lower() for s in patient_profile.supplements]:
                score += 10  # Bonus: patient already has it
            
            # Check for contraindications
            for condition in patient_profile.conditions:
                if condition.lower() in [c.lower() for c in adj_data.get("contraindications", [])]:
                    score -= 50  # Heavy penalty
        
        return score
    
    def _check_contraindications(
        self, adj_data: Dict[str, Any], patient_profile: PatientProfile
    ) -> List[str]:
        """Check for contraindications and drug interactions."""
        warnings = []
        
        # Check condition contraindications
        for condition in patient_profile.conditions:
            if condition.lower() in [c.lower() for c in adj_data.get("contraindications", [])]:
                warnings.append(
                    f"{adj_data['name']} contraindicated with {condition}"
                )
        
        # Check drug interactions
        for med in patient_profile.current_medications:
            med_lower = med.lower()
            for interaction in adj_data.get("drug_interactions", []):
                if interaction.lower() in med_lower or med_lower in interaction.lower():
                    warnings.append(
                        f"{adj_data['name']} may interact with {med}"
                    )
        
        return warnings
    
    def _calculate_dose(
        self, adj_data: Dict[str, Any], patient_profile: Optional[PatientProfile]
    ) -> float:
        """Calculate personalized dose based on patient profile."""
        base_dose = adj_data["recommended_dose_mg"]
        
        if not patient_profile:
            return base_dose
        
        dose = base_dose
        
        # Adjust for weight (use 70kg as reference)
        weight_factor = patient_profile.weight_kg / 70
        dose *= min(1.3, max(0.7, weight_factor))
        
        # Adjust for age (reduce for elderly)
        if patient_profile.age >= 65:
            dose *= 0.75
        
        # Adjust for metabolizer status
        if patient_profile.metabolizer_status == "slow":
            dose *= 0.75
        elif patient_profile.metabolizer_status == "fast":
            dose *= 1.25
        
        # Ensure within range
        min_dose, max_dose = adj_data["dose_range"]
        dose = max(min_dose, min(max_dose, dose))
        
        return round(dose, 0)
    
    def _get_confidence_interval(
        self, enhancement: float, evidence_tier: int
    ) -> Tuple[float, float]:
        """Calculate confidence interval based on evidence tier."""
        # Higher tier = tighter CI
        ci_width = {
            1: 0.1,   # ±10%
            2: 0.15,  # ±15%
            3: 0.2,   # ±20%
            4: 0.25,  # ±25%
            5: 0.35,  # ±35%
        }.get(evidence_tier, 0.25)
        
        lower = max(0, enhancement * (1 - ci_width))
        upper = enhancement * (1 + ci_width)
        
        return (round(lower, 1), round(upper, 1))
    
    def _calculate_total_enhancement(
        self, recommendations: List[AdjuvantRecommendation]
    ) -> float:
        """Calculate total enhancement with diminishing returns."""
        if not recommendations:
            return 0.0
        
        # First adjuvant: full effect
        # Second: 70% of stated effect
        # Third: 50% of stated effect
        # etc.
        
        multipliers = [1.0, 0.7, 0.5, 0.3, 0.2]
        total = 0.0
        
        for i, rec in enumerate(recommendations):
            mult = multipliers[i] if i < len(multipliers) else 0.1
            total += rec.expected_enhancement_percent * mult
        
        return round(total, 1)
    
    def _generate_protocol_summary(
        self,
        primary_compound: str,
        recommendations: List[AdjuvantRecommendation],
    ) -> str:
        """Generate human-readable protocol summary."""
        if not recommendations:
            return f"Take {primary_compound} as directed. No adjuvants recommended."
        
        # Sort by timing (most negative first = earliest)
        sorted_recs = sorted(recommendations, key=lambda r: r.timing_offset_minutes)
        
        steps = []
        for rec in sorted_recs:
            timing = rec.timing_offset_minutes
            if timing < 0:
                timing_str = f"{abs(timing)} minutes before"
            elif timing > 0:
                timing_str = f"{timing} minutes after"
            else:
                timing_str = "together with"
            
            steps.append(f"• {rec.adjuvant_name} {rec.dosage_mg}mg — {timing_str} {primary_compound}")
        
        # Add primary compound timing
        steps.append(f"• {primary_compound} — standard dose")
        
        return "\n".join([
            "PROTOCOL:",
            *steps,
            f"\nExpected enhancement: {self._calculate_total_enhancement(recommendations)}%"
        ])


# Singleton instance
_optimizer_instance: Optional[AdjuvantOptimizer] = None


def get_adjuvant_optimizer() -> AdjuvantOptimizer:
    """Get or create the AdjuvantOptimizer singleton."""
    global _optimizer_instance
    if _optimizer_instance is None:
        _optimizer_instance = AdjuvantOptimizer()
    return _optimizer_instance
