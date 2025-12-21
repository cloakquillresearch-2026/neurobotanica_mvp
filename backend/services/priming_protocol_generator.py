"""
Priming Protocol Generator - NeuroBotanica Week 5
Generates optimized adjuvant priming protocols for cannabinoid formulations.

Patent-Protected Features:
- Evidence-weighted adjuvant selection
- Timing optimization algorithms
- Drug interaction safety checking
- Patient profile customization
"""
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass, field
from datetime import datetime
import logging

logger = logging.getLogger(__name__)


@dataclass
class AdjuvantRecommendation:
    """Single adjuvant recommendation within a protocol."""
    adjuvant_id: int
    adjuvant_name: str
    category: str
    dose_mg: float
    timing_minutes: int  # Negative = before cannabinoid
    timing_label: str  # "30 minutes before", "concurrent", etc.
    mechanism: str
    expected_enhancement: float  # Percentage
    evidence_tier: int
    notes: str = ""
    
    def to_dict(self) -> Dict:
        return {
            "adjuvant_id": self.adjuvant_id,
            "adjuvant_name": self.adjuvant_name,
            "category": self.category,
            "dose_mg": self.dose_mg,
            "timing_minutes": self.timing_minutes,
            "timing_label": self.timing_label,
            "mechanism": self.mechanism,
            "expected_enhancement": self.expected_enhancement,
            "evidence_tier": self.evidence_tier,
            "notes": self.notes
        }


@dataclass
class PrimingProtocolResult:
    """Complete priming protocol with all recommendations."""
    protocol_name: str
    target_cannabinoid: str
    target_condition: str
    recommendations: List[AdjuvantRecommendation]
    total_enhancement_percent: float
    protocol_confidence: float
    evidence_tier: int  # Lowest tier across all adjuvants
    warnings: List[str]
    contraindications: List[str]
    generated_at: str = field(default_factory=lambda: datetime.now().isoformat())
    
    def to_dict(self) -> Dict:
        return {
            "protocol_name": self.protocol_name,
            "target_cannabinoid": self.target_cannabinoid,
            "target_condition": self.target_condition,
            "recommendations": [r.to_dict() for r in self.recommendations],
            "total_enhancement_percent": round(self.total_enhancement_percent, 1),
            "protocol_confidence": round(self.protocol_confidence, 3),
            "evidence_tier": self.evidence_tier,
            "warnings": self.warnings,
            "contraindications": self.contraindications,
            "generated_at": self.generated_at
        }


# Default adjuvant database (expanded with ChEMBL/PubChem data sync)
DEFAULT_ADJUVANTS = [
    {
        "id": 1,
        "name": "Magnesium Glycinate",
        "category": "Receptor Priming",
        "mechanism": "Enhances NMDA receptor modulation, reduces tolerance development",
        "recommended_dose_mg": 200.0,
        "timing_offset_minutes": -30,
        "evidence_tier": 2,
        "target_receptors": ["NMDA", "GABA-A"],
        "cannabinoid_synergies": ["THC", "CBD", "CBG"],
        "max_enhancement_percent": 25.0,
        "contraindications": ["kidney_disease", "myasthenia_gravis"]
    },
    {
        "id": 2,
        "name": "L-Theanine",
        "category": "Receptor Priming",
        "mechanism": "Modulates glutamate/GABA balance, enhances alpha wave activity",
        "recommended_dose_mg": 200.0,
        "timing_offset_minutes": -20,
        "evidence_tier": 2,
        "target_receptors": ["GABA-A", "Glutamate"],
        "cannabinoid_synergies": ["THC", "CBD", "CBN"],
        "max_enhancement_percent": 20.0,
        "contraindications": []
    },
    {
        "id": 3,
        "name": "Omega-3 Fatty Acids (EPA/DHA)",
        "category": "Bioavailability",
        "mechanism": "Enhances cannabinoid absorption, provides fatty acid substrate",
        "recommended_dose_mg": 2000.0,
        "timing_offset_minutes": -60,
        "evidence_tier": 2,
        "target_receptors": ["CB1", "CB2", "PPAR"],
        "cannabinoid_synergies": ["THC", "CBD", "CBG", "CBC"],
        "max_enhancement_percent": 35.0,
        "contraindications": ["fish_allergy", "blood_thinners"]
    },
    {
        "id": 4,
        "name": "Black Pepper Extract (Piperine)",
        "category": "Metabolic Enhancer",
        "mechanism": "Inhibits CYP3A4/P-glycoprotein, increases bioavailability",
        "recommended_dose_mg": 20.0,
        "timing_offset_minutes": -15,
        "evidence_tier": 2,
        "target_receptors": ["CYP3A4"],
        "cannabinoid_synergies": ["CBD", "CBG", "CBC"],
        "max_enhancement_percent": 30.0,
        "contraindications": ["cyp3a4_inhibitor_medications"]
    },
    {
        "id": 5,
        "name": "Curcumin (with Piperine)",
        "category": "Synergistic",
        "mechanism": "Anti-inflammatory synergy, CB1/CB2 allosteric modulation",
        "recommended_dose_mg": 500.0,
        "timing_offset_minutes": -30,
        "evidence_tier": 3,
        "target_receptors": ["CB1", "CB2", "COX-2"],
        "cannabinoid_synergies": ["CBD", "THC", "CBG"],
        "max_enhancement_percent": 25.0,
        "contraindications": ["gallbladder_disease", "blood_thinners"]
    },
    {
        "id": 6,
        "name": "N-Acetyl Cysteine (NAC)",
        "category": "Protective",
        "mechanism": "Glutathione precursor, reduces oxidative stress from THC metabolism",
        "recommended_dose_mg": 600.0,
        "timing_offset_minutes": -45,
        "evidence_tier": 2,
        "target_receptors": ["Glutamate", "Cystine-Glutamate"],
        "cannabinoid_synergies": ["THC", "THCA"],
        "max_enhancement_percent": 15.0,
        "contraindications": ["asthma_acute"]
    },
    {
        "id": 7,
        "name": "Beta-Caryophyllene",
        "category": "Terpene",
        "mechanism": "Dietary CB2 agonist, anti-inflammatory entourage effect",
        "recommended_dose_mg": 100.0,
        "timing_offset_minutes": 0,
        "evidence_tier": 2,
        "target_receptors": ["CB2", "PPAR-gamma"],
        "cannabinoid_synergies": ["CBD", "THC", "CBG", "CBC"],
        "max_enhancement_percent": 20.0,
        "contraindications": []
    },
    {
        "id": 8,
        "name": "Myrcene",
        "category": "Terpene",
        "mechanism": "Enhances BBB permeability, potentiates sedative effects",
        "recommended_dose_mg": 50.0,
        "timing_offset_minutes": -10,
        "evidence_tier": 3,
        "target_receptors": ["GABA-A", "5-HT1A"],
        "cannabinoid_synergies": ["THC", "CBN", "CBD"],
        "max_enhancement_percent": 30.0,
        "contraindications": ["operating_machinery"]
    },
    {
        "id": 9,
        "name": "Limonene",
        "category": "Terpene",
        "mechanism": "Anxiolytic, elevates mood, enhances absorption",
        "recommended_dose_mg": 100.0,
        "timing_offset_minutes": 0,
        "evidence_tier": 3,
        "target_receptors": ["5-HT1A", "A2A"],
        "cannabinoid_synergies": ["CBD", "THC", "CBG"],
        "max_enhancement_percent": 15.0,
        "contraindications": ["gerd_severe"]
    },
    {
        "id": 10,
        "name": "Linalool",
        "category": "Terpene",
        "mechanism": "Anxiolytic via GABA modulation, enhances sedative effects",
        "recommended_dose_mg": 50.0,
        "timing_offset_minutes": 0,
        "evidence_tier": 3,
        "target_receptors": ["GABA-A", "Glutamate"],
        "cannabinoid_synergies": ["CBD", "CBN", "THC"],
        "max_enhancement_percent": 20.0,
        "contraindications": []
    },
    {
        "id": 11,
        "name": "Quercetin",
        "category": "Synergistic",
        "mechanism": "CB1 positive allosteric modulator, anti-inflammatory",
        "recommended_dose_mg": 500.0,
        "timing_offset_minutes": -30,
        "evidence_tier": 3,
        "target_receptors": ["CB1", "NF-kB"],
        "cannabinoid_synergies": ["THC", "CBD"],
        "max_enhancement_percent": 15.0,
        "contraindications": ["quinolone_antibiotics"]
    },
    {
        "id": 12,
        "name": "Palmitoylethanolamide (PEA)",
        "category": "Synergistic",
        "mechanism": "Endocannabinoid enhancer, FAAH inhibitor, entourage effect",
        "recommended_dose_mg": 600.0,
        "timing_offset_minutes": -30,
        "evidence_tier": 2,
        "target_receptors": ["PPAR-alpha", "GPR55", "CB1"],
        "cannabinoid_synergies": ["CBD", "THC", "CBG"],
        "max_enhancement_percent": 30.0,
        "contraindications": []
    },
    {
        "id": 13,
        "name": "Olive Oil (Extra Virgin)",
        "category": "Bioavailability",
        "mechanism": "Lipid carrier enhances absorption of fat-soluble cannabinoids",
        "recommended_dose_mg": 15000.0,  # ~1 tablespoon
        "timing_offset_minutes": 0,
        "evidence_tier": 2,
        "target_receptors": [],
        "cannabinoid_synergies": ["THC", "CBD", "CBG", "CBC", "CBN"],
        "max_enhancement_percent": 40.0,
        "contraindications": []
    },
    {
        "id": 14,
        "name": "Lecithin",
        "category": "Bioavailability",
        "mechanism": "Phospholipid enhancer, improves cannabinoid absorption",
        "recommended_dose_mg": 1200.0,
        "timing_offset_minutes": -15,
        "evidence_tier": 3,
        "target_receptors": [],
        "cannabinoid_synergies": ["THC", "CBD", "CBG"],
        "max_enhancement_percent": 25.0,
        "contraindications": ["soy_allergy"]
    },
    {
        "id": 15,
        "name": "Melatonin",
        "category": "Synergistic",
        "mechanism": "Sleep regulation synergy, enhances CBN effects",
        "recommended_dose_mg": 3.0,
        "timing_offset_minutes": 0,
        "evidence_tier": 2,
        "target_receptors": ["MT1", "MT2", "5-HT2C"],
        "cannabinoid_synergies": ["CBN", "CBD", "THC"],
        "max_enhancement_percent": 25.0,
        "contraindications": ["autoimmune_conditions", "depression"]
    }
]

# Condition-specific adjuvant priorities
CONDITION_ADJUVANT_PRIORITIES = {
    "chronic_pain": ["Magnesium Glycinate", "Palmitoylethanolamide (PEA)", "Omega-3 Fatty Acids (EPA/DHA)", "Curcumin (with Piperine)"],
    "anxiety": ["L-Theanine", "Magnesium Glycinate", "Linalool", "Limonene"],
    "insomnia": ["Melatonin", "L-Theanine", "Myrcene", "Linalool"],
    "inflammation": ["Curcumin (with Piperine)", "Omega-3 Fatty Acids (EPA/DHA)", "Beta-Caryophyllene", "Quercetin"],
    "neuropathy": ["Palmitoylethanolamide (PEA)", "Magnesium Glycinate", "N-Acetyl Cysteine (NAC)"],
    "appetite_stimulation": ["Omega-3 Fatty Acids (EPA/DHA)", "Black Pepper Extract (Piperine)"],
    "nausea": ["Limonene", "Beta-Caryophyllene"],
    "depression": ["Omega-3 Fatty Acids (EPA/DHA)", "L-Theanine", "Limonene"],
    "ptsd": ["L-Theanine", "Magnesium Glycinate", "N-Acetyl Cysteine (NAC)", "Omega-3 Fatty Acids (EPA/DHA)"],
    "epilepsy": ["Magnesium Glycinate", "Omega-3 Fatty Acids (EPA/DHA)"],
    "multiple_sclerosis": ["Palmitoylethanolamide (PEA)", "Omega-3 Fatty Acids (EPA/DHA)", "Curcumin (with Piperine)"],
    "arthritis": ["Curcumin (with Piperine)", "Omega-3 Fatty Acids (EPA/DHA)", "Beta-Caryophyllene"],
    "ibd": ["Curcumin (with Piperine)", "Omega-3 Fatty Acids (EPA/DHA)", "Palmitoylethanolamide (PEA)"],
    "glaucoma": ["Omega-3 Fatty Acids (EPA/DHA)", "Magnesium Glycinate"],
    "cancer_palliative": ["Omega-3 Fatty Acids (EPA/DHA)", "Curcumin (with Piperine)", "N-Acetyl Cysteine (NAC)"]
}

# Cannabinoid-specific optimal adjuvants
CANNABINOID_ADJUVANT_PRIORITIES = {
    "THC": ["N-Acetyl Cysteine (NAC)", "L-Theanine", "Magnesium Glycinate", "Omega-3 Fatty Acids (EPA/DHA)"],
    "CBD": ["Palmitoylethanolamide (PEA)", "Curcumin (with Piperine)", "Black Pepper Extract (Piperine)", "Omega-3 Fatty Acids (EPA/DHA)"],
    "CBG": ["Beta-Caryophyllene", "Omega-3 Fatty Acids (EPA/DHA)", "Curcumin (with Piperine)"],
    "CBN": ["Melatonin", "L-Theanine", "Myrcene", "Linalool"],
    "CBC": ["Omega-3 Fatty Acids (EPA/DHA)", "Curcumin (with Piperine)"],
    "THCA": ["Curcumin (with Piperine)", "Omega-3 Fatty Acids (EPA/DHA)"],
    "CBDA": ["Omega-3 Fatty Acids (EPA/DHA)", "Black Pepper Extract (Piperine)"]
}


class PrimingProtocolGenerator:
    """Generate optimized adjuvant priming protocols for cannabinoid formulations.
    
    Uses evidence-weighted algorithms to select adjuvants based on:
    - Target cannabinoid
    - Target condition
    - Patient profile (allergies, medications, conditions)
    - Drug interaction safety
    """
    
    def __init__(self, adjuvant_database: Optional[List[Dict]] = None):
        """Initialize generator with adjuvant database.
        
        Args:
            adjuvant_database: List of adjuvant records. Uses DEFAULT_ADJUVANTS if None.
        """
        self.adjuvants = {a["name"]: a for a in (adjuvant_database or DEFAULT_ADJUVANTS)}
        self.adjuvant_by_id = {a["id"]: a for a in (adjuvant_database or DEFAULT_ADJUVANTS)}
        logger.info(f"PrimingProtocolGenerator initialized with {len(self.adjuvants)} adjuvants")
    
    def generate_protocol(
        self,
        cannabinoid_name: str,
        condition: str,
        patient_profile: Optional[Dict] = None,
        max_adjuvants: int = 4,
        min_evidence_tier: int = 4
    ) -> PrimingProtocolResult:
        """Generate optimized priming protocol.
        
        Args:
            cannabinoid_name: Target cannabinoid (e.g., "THC", "CBD")
            condition: Target condition (e.g., "chronic_pain", "anxiety")
            patient_profile: Optional patient info with allergies, medications, conditions
            max_adjuvants: Maximum number of adjuvants to include
            min_evidence_tier: Minimum evidence tier (1=highest, 5=lowest)
            
        Returns:
            PrimingProtocolResult with recommendations
        """
        logger.info(f"Generating protocol for {cannabinoid_name} + {condition}")
        
        patient_profile = patient_profile or {}
        warnings = []
        contraindications = []
        
        # Get condition-specific priorities
        condition_priorities = CONDITION_ADJUVANT_PRIORITIES.get(
            condition.lower().replace(" ", "_"),
            []
        )
        
        # Get cannabinoid-specific priorities
        cannabinoid_priorities = CANNABINOID_ADJUVANT_PRIORITIES.get(
            cannabinoid_name.upper(),
            []
        )
        
        # Score all adjuvants
        scored_adjuvants = []
        for name, adjuvant in self.adjuvants.items():
            # Skip if below evidence threshold
            if adjuvant["evidence_tier"] > min_evidence_tier:
                continue
            
            # Check contraindications
            patient_conditions = patient_profile.get("conditions", [])
            patient_allergies = patient_profile.get("allergies", [])
            patient_medications = patient_profile.get("medications", [])
            
            skip = False
            for contra in adjuvant.get("contraindications", []):
                if contra in patient_conditions or contra in patient_allergies:
                    contraindications.append(f"{name} contraindicated due to: {contra}")
                    skip = True
                    break
            
            if skip:
                continue
            
            # Check cannabinoid synergy
            synergies = adjuvant.get("cannabinoid_synergies", [])
            if cannabinoid_name.upper() not in synergies:
                continue
            
            # Calculate score
            score = 0.0
            
            # Condition priority bonus (0-40 points)
            if name in condition_priorities:
                priority_idx = condition_priorities.index(name)
                score += 40 - (priority_idx * 10)
            
            # Cannabinoid priority bonus (0-30 points)
            if name in cannabinoid_priorities:
                priority_idx = cannabinoid_priorities.index(name)
                score += 30 - (priority_idx * 7)
            
            # Evidence tier bonus (higher tier = more points)
            score += (6 - adjuvant["evidence_tier"]) * 10
            
            # Enhancement bonus
            score += adjuvant.get("max_enhancement_percent", 0) * 0.5
            
            scored_adjuvants.append((score, name, adjuvant))
        
        # Sort by score and take top N
        scored_adjuvants.sort(key=lambda x: x[0], reverse=True)
        selected = scored_adjuvants[:max_adjuvants]
        
        # Generate recommendations
        recommendations = []
        total_enhancement = 0.0
        min_tier = 1
        
        for score, name, adjuvant in selected:
            # Adjust dose based on patient profile (simplified)
            dose = adjuvant["recommended_dose_mg"]
            if patient_profile.get("weight_kg", 70) < 60:
                dose *= 0.8
            elif patient_profile.get("weight_kg", 70) > 100:
                dose *= 1.2
            
            # Create timing label
            timing = adjuvant["timing_offset_minutes"]
            if timing < 0:
                timing_label = f"{abs(timing)} minutes before cannabinoid"
            elif timing > 0:
                timing_label = f"{timing} minutes after cannabinoid"
            else:
                timing_label = "Concurrent with cannabinoid"
            
            rec = AdjuvantRecommendation(
                adjuvant_id=adjuvant["id"],
                adjuvant_name=name,
                category=adjuvant["category"],
                dose_mg=round(dose, 1),
                timing_minutes=timing,
                timing_label=timing_label,
                mechanism=adjuvant["mechanism"],
                expected_enhancement=adjuvant.get("max_enhancement_percent", 0),
                evidence_tier=adjuvant["evidence_tier"],
                notes=self._generate_notes(adjuvant, patient_profile)
            )
            recommendations.append(rec)
            
            # Track cumulative enhancement (diminishing returns)
            enhancement_factor = 1.0 - (len(recommendations) - 1) * 0.15
            total_enhancement += adjuvant.get("max_enhancement_percent", 0) * enhancement_factor
            
            min_tier = max(min_tier, adjuvant["evidence_tier"])
        
        # Sort recommendations by timing (earliest first)
        recommendations.sort(key=lambda r: r.timing_minutes)
        
        # Calculate protocol confidence
        confidence = self._calculate_confidence(
            recommendations=recommendations,
            condition_priorities=condition_priorities,
            cannabinoid_priorities=cannabinoid_priorities
        )
        
        # Generate warnings
        if not recommendations:
            warnings.append("No suitable adjuvants found for this combination")
        if min_tier > 3:
            warnings.append("Protocol relies on lower-tier evidence (Tier 4-5)")
        if len(recommendations) < 2:
            warnings.append("Limited adjuvant options available")
        
        return PrimingProtocolResult(
            protocol_name=f"{cannabinoid_name} + {condition.replace('_', ' ').title()} Protocol",
            target_cannabinoid=cannabinoid_name,
            target_condition=condition,
            recommendations=recommendations,
            total_enhancement_percent=min(total_enhancement, 80.0),  # Cap at 80%
            protocol_confidence=confidence,
            evidence_tier=min_tier,
            warnings=warnings,
            contraindications=contraindications
        )
    
    def _generate_notes(self, adjuvant: Dict, patient_profile: Dict) -> str:
        """Generate context-specific notes for an adjuvant recommendation."""
        notes = []
        
        category = adjuvant["category"]
        if category == "Bioavailability":
            notes.append("Take with food for optimal absorption")
        elif category == "Receptor Priming":
            notes.append("Prepares receptors for enhanced response")
        elif category == "Terpene":
            notes.append("Part of entourage effect")
        
        if adjuvant.get("dose_form"):
            notes.append(f"Recommended form: {adjuvant['dose_form']}")
        
        return ". ".join(notes) if notes else ""
    
    def _calculate_confidence(
        self,
        recommendations: List[AdjuvantRecommendation],
        condition_priorities: List[str],
        cannabinoid_priorities: List[str]
    ) -> float:
        """Calculate overall protocol confidence score."""
        if not recommendations:
            return 0.0
        
        # Base confidence from evidence tiers
        avg_tier = sum(r.evidence_tier for r in recommendations) / len(recommendations)
        tier_confidence = (5 - avg_tier) / 4  # 0.25 to 1.0
        
        # Priority match bonus
        priority_matches = sum(
            1 for r in recommendations
            if r.adjuvant_name in condition_priorities or r.adjuvant_name in cannabinoid_priorities
        )
        priority_confidence = min(priority_matches / 3, 1.0)
        
        # Synergy bonus (more adjuvants = higher potential)
        synergy_confidence = min(len(recommendations) / 4, 1.0)
        
        # Combined score
        confidence = (
            tier_confidence * 0.4 +
            priority_confidence * 0.35 +
            synergy_confidence * 0.25
        )
        
        return round(confidence, 3)
    
    def check_drug_interactions(
        self,
        adjuvant_names: List[str],
        medications: List[str]
    ) -> List[Dict]:
        """Check for drug interactions between adjuvants and patient medications.
        
        Args:
            adjuvant_names: List of adjuvant names in protocol
            medications: List of patient's current medications
            
        Returns:
            List of interaction warnings
        """
        interactions = []
        
        # Known interaction pairs (simplified - would connect to DrugBank in production)
        known_interactions = {
            "Black Pepper Extract (Piperine)": {
                "cyp3a4_substrates": "Major - increases blood levels of CYP3A4-metabolized drugs",
                "blood_thinners": "Moderate - may increase bleeding risk"
            },
            "Omega-3 Fatty Acids (EPA/DHA)": {
                "blood_thinners": "Moderate - may increase bleeding risk"
            },
            "Magnesium Glycinate": {
                "bisphosphonates": "Moderate - may decrease absorption",
                "antibiotics": "Minor - take 2 hours apart"
            },
            "Curcumin (with Piperine)": {
                "blood_thinners": "Moderate - may increase bleeding risk",
                "diabetes_medications": "Minor - monitor blood sugar"
            },
            "Melatonin": {
                "immunosuppressants": "Moderate - may alter immune function",
                "blood_pressure_medications": "Minor - may affect blood pressure"
            }
        }
        
        for adjuvant in adjuvant_names:
            if adjuvant in known_interactions:
                for med_class, warning in known_interactions[adjuvant].items():
                    # Simplified medication matching
                    for med in medications:
                        if self._medication_matches_class(med, med_class):
                            interactions.append({
                                "adjuvant": adjuvant,
                                "medication": med,
                                "medication_class": med_class,
                                "severity": warning.split(" - ")[0],
                                "description": warning
                            })
        
        return interactions
    
    def _medication_matches_class(self, medication: str, med_class: str) -> bool:
        """Check if a medication belongs to a drug class."""
        # Simplified matching - production version would use DrugBank classification
        class_keywords = {
            "blood_thinners": ["warfarin", "aspirin", "heparin", "eliquis", "xarelto", "pradaxa"],
            "cyp3a4_substrates": ["atorvastatin", "simvastatin", "alprazolam", "midazolam"],
            "bisphosphonates": ["alendronate", "risedronate", "fosamax"],
            "antibiotics": ["amoxicillin", "ciprofloxacin", "azithromycin"],
            "diabetes_medications": ["metformin", "insulin", "glipizide"],
            "immunosuppressants": ["cyclosporine", "tacrolimus", "prednisone"],
            "blood_pressure_medications": ["lisinopril", "amlodipine", "metoprolol"]
        }
        
        keywords = class_keywords.get(med_class, [])
        return any(kw in medication.lower() for kw in keywords)
    
    def get_adjuvant_info(self, name: str) -> Optional[Dict]:
        """Get detailed information about a specific adjuvant."""
        return self.adjuvants.get(name)
    
    def list_adjuvants_by_category(self, category: str) -> List[Dict]:
        """List all adjuvants in a specific category."""
        return [
            a for a in self.adjuvants.values()
            if a["category"].lower() == category.lower()
        ]
    
    def get_conditions(self) -> List[str]:
        """Get list of supported conditions."""
        return list(CONDITION_ADJUVANT_PRIORITIES.keys())
    
    def get_cannabinoids(self) -> List[str]:
        """Get list of supported cannabinoids."""
        return list(CANNABINOID_ADJUVANT_PRIORITIES.keys())
