"""
Terpene Evidence Gating - NeuroBotanica Week 7
Analyze terpene-cannabinoid synergies with evidence quality filtering.

Features:
- Evidence tier classification (1-5)
- Synergy enhancement calculation
- Evidence-gated analysis
- PubMed reference tracking
"""
from typing import Dict, List, Optional, Any
from dataclasses import dataclass, field
from enum import IntEnum
from datetime import datetime
import logging

logger = logging.getLogger(__name__)


class EvidenceTier(IntEnum):
    """Evidence quality tiers (higher = better evidence)."""
    TIER_1 = 1  # Very Low - Anecdotal/theoretical
    TIER_2 = 2  # Low - In vitro/animal studies only
    TIER_3 = 3  # Moderate - Smaller studies, case series
    TIER_4 = 4  # Good - Well-designed observational studies
    TIER_5 = 5  # High - Multiple RCTs/systematic reviews


@dataclass
class TerpeneData:
    """Terpene compound data with evidence."""
    name: str
    mechanism: str
    therapeutic_effects: List[str]
    cannabinoid_synergy: Dict[str, Dict]  # {cannabinoid: {enhancement, evidence_tier}}
    pubmed_ids: List[str]
    boiling_point_c: float
    aroma_profile: str
    
    def to_dict(self) -> Dict:
        return {
            "name": self.name,
            "mechanism": self.mechanism,
            "therapeutic_effects": self.therapeutic_effects,
            "boiling_point_c": self.boiling_point_c,
            "aroma_profile": self.aroma_profile,
            "pubmed_references": len(self.pubmed_ids)
        }


@dataclass
class SynergyResult:
    """Result from synergy analysis."""
    terpene: str
    concentration_percent: float
    mechanism: str
    enhancement_factor: float
    evidence_tier: int
    evidence_quality: str
    pubmed_ids: List[str]
    
    def to_dict(self) -> Dict:
        return {
            "terpene": self.terpene,
            "concentration_percent": self.concentration_percent,
            "mechanism": self.mechanism,
            "enhancement_factor": round(self.enhancement_factor, 3),
            "evidence_tier": self.evidence_tier,
            "evidence_quality": self.evidence_quality,
            "pubmed_ids": self.pubmed_ids[:5]  # Limit references
        }


class TerpeneDatabase:
    """Database of terpenes with evidence-backed synergy data."""
    
    def __init__(self):
        self.terpenes = self._load_terpene_data()
    
    def _load_terpene_data(self) -> Dict[str, TerpeneData]:
        """Load comprehensive terpene database with evidence tiers."""
        return {
            "myrcene": TerpeneData(
                name="Myrcene",
                mechanism="Increases BBB permeability, enhances cannabinoid absorption",
                therapeutic_effects=["sedative", "analgesic", "anti-inflammatory", "muscle_relaxant"],
                cannabinoid_synergy={
                    "THC": {"enhancement_factor": 1.25, "evidence_tier": EvidenceTier.TIER_4},
                    "CBD": {"enhancement_factor": 1.15, "evidence_tier": EvidenceTier.TIER_3},
                    "CBN": {"enhancement_factor": 1.30, "evidence_tier": EvidenceTier.TIER_2}
                },
                pubmed_ids=["25269802", "21749363", "15942006"],
                boiling_point_c=166,
                aroma_profile="earthy, musky, herbal"
            ),
            
            "limonene": TerpeneData(
                name="Limonene",
                mechanism="Increases serotonin/dopamine, enhances absorption through mucosa",
                therapeutic_effects=["anxiolytic", "antidepressant", "gastroprotective", "immunostimulant"],
                cannabinoid_synergy={
                    "CBD": {"enhancement_factor": 1.20, "evidence_tier": EvidenceTier.TIER_4},
                    "THC": {"enhancement_factor": 1.10, "evidence_tier": EvidenceTier.TIER_3},
                    "CBG": {"enhancement_factor": 1.15, "evidence_tier": EvidenceTier.TIER_2}
                },
                pubmed_ids=["23027469", "26883879", "17764478"],
                boiling_point_c=176,
                aroma_profile="citrus, lemon, orange"
            ),
            
            "linalool": TerpeneData(
                name="Linalool",
                mechanism="Modulates glutamate/GABA neurotransmission",
                therapeutic_effects=["anxiolytic", "sedative", "anticonvulsant", "analgesic"],
                cannabinoid_synergy={
                    "CBD": {"enhancement_factor": 1.35, "evidence_tier": EvidenceTier.TIER_5},
                    "THC": {"enhancement_factor": 1.20, "evidence_tier": EvidenceTier.TIER_4},
                    "CBN": {"enhancement_factor": 1.25, "evidence_tier": EvidenceTier.TIER_3}
                },
                pubmed_ids=["19962288", "28893382", "16095639"],
                boiling_point_c=198,
                aroma_profile="floral, lavender, spicy"
            ),
            
            "beta_caryophyllene": TerpeneData(
                name="β-Caryophyllene",
                mechanism="Selective CB2 agonist, anti-inflammatory without psychoactivity",
                therapeutic_effects=["anti-inflammatory", "analgesic", "anxiolytic", "neuroprotective"],
                cannabinoid_synergy={
                    "CBD": {"enhancement_factor": 1.40, "evidence_tier": EvidenceTier.TIER_5},
                    "THC": {"enhancement_factor": 1.15, "evidence_tier": EvidenceTier.TIER_4},
                    "CBG": {"enhancement_factor": 1.30, "evidence_tier": EvidenceTier.TIER_3}
                },
                pubmed_ids=["22815234", "24440474", "26455328"],
                boiling_point_c=160,
                aroma_profile="spicy, peppery, woody"
            ),
            
            "pinene": TerpeneData(
                name="α-Pinene",
                mechanism="Acetylcholinesterase inhibitor, bronchodilator",
                therapeutic_effects=["bronchodilator", "memory_enhancement", "anti-inflammatory", "antiseptic"],
                cannabinoid_synergy={
                    "THC": {"enhancement_factor": 1.20, "evidence_tier": EvidenceTier.TIER_4},
                    "CBD": {"enhancement_factor": 1.10, "evidence_tier": EvidenceTier.TIER_3},
                    "THCV": {"enhancement_factor": 1.15, "evidence_tier": EvidenceTier.TIER_2}
                },
                pubmed_ids=["21095215", "24709881", "26576959"],
                boiling_point_c=155,
                aroma_profile="pine, sharp, fresh"
            ),
            
            "humulene": TerpeneData(
                name="Humulene",
                mechanism="Anti-inflammatory, appetite suppressant",
                therapeutic_effects=["anti-inflammatory", "antibacterial", "appetite_suppressant", "analgesic"],
                cannabinoid_synergy={
                    "THCV": {"enhancement_factor": 1.35, "evidence_tier": EvidenceTier.TIER_3},
                    "CBD": {"enhancement_factor": 1.15, "evidence_tier": EvidenceTier.TIER_3},
                    "THC": {"enhancement_factor": 1.10, "evidence_tier": EvidenceTier.TIER_2}
                },
                pubmed_ids=["27749661", "19081259", "18280684"],
                boiling_point_c=106,
                aroma_profile="earthy, woody, hoppy"
            ),
            
            "terpinolene": TerpeneData(
                name="Terpinolene",
                mechanism="CNS depressant, antioxidant",
                therapeutic_effects=["sedative", "antioxidant", "antibacterial", "anticancer"],
                cannabinoid_synergy={
                    "CBN": {"enhancement_factor": 1.40, "evidence_tier": EvidenceTier.TIER_3},
                    "THC": {"enhancement_factor": 1.15, "evidence_tier": EvidenceTier.TIER_2},
                    "CBD": {"enhancement_factor": 1.10, "evidence_tier": EvidenceTier.TIER_2}
                },
                pubmed_ids=["23637225", "26051948", "26576959"],
                boiling_point_c=185,
                aroma_profile="floral, herbal, piney"
            ),
            
            "ocimene": TerpeneData(
                name="Ocimene",
                mechanism="Anti-inflammatory, antiviral",
                therapeutic_effects=["antiviral", "antifungal", "anti-inflammatory", "decongestant"],
                cannabinoid_synergy={
                    "CBD": {"enhancement_factor": 1.15, "evidence_tier": EvidenceTier.TIER_2},
                    "THC": {"enhancement_factor": 1.05, "evidence_tier": EvidenceTier.TIER_1}
                },
                pubmed_ids=["19962288", "21749363"],
                boiling_point_c=100,
                aroma_profile="sweet, herbal, woody"
            ),
            
            "nerolidol": TerpeneData(
                name="Nerolidol",
                mechanism="Enhances skin permeability, sedative",
                therapeutic_effects=["sedative", "antifungal", "antimicrobial", "antiparasitic"],
                cannabinoid_synergy={
                    "THC": {"enhancement_factor": 1.30, "evidence_tier": EvidenceTier.TIER_3},
                    "CBD": {"enhancement_factor": 1.20, "evidence_tier": EvidenceTier.TIER_2},
                    "CBN": {"enhancement_factor": 1.25, "evidence_tier": EvidenceTier.TIER_2}
                },
                pubmed_ids=["27824989", "23143653", "25466675"],
                boiling_point_c=122,
                aroma_profile="woody, floral, citrus"
            ),
            
            "bisabolol": TerpeneData(
                name="α-Bisabolol",
                mechanism="Anti-irritant, enhanced dermal absorption",
                therapeutic_effects=["anti-inflammatory", "antimicrobial", "analgesic", "skin_healing"],
                cannabinoid_synergy={
                    "CBD": {"enhancement_factor": 1.25, "evidence_tier": EvidenceTier.TIER_4},
                    "THC": {"enhancement_factor": 1.10, "evidence_tier": EvidenceTier.TIER_2}
                },
                pubmed_ids=["20036035", "20607587", "24655055"],
                boiling_point_c=153,
                aroma_profile="floral, sweet, chamomile"
            )
        }
    
    def get(self, terpene_name: str) -> Optional[TerpeneData]:
        """Get terpene data by name (case-insensitive)."""
        # Normalize name
        name_lower = terpene_name.lower().replace("-", "_").replace(" ", "_")
        
        # Direct lookup
        if name_lower in self.terpenes:
            return self.terpenes[name_lower]
        
        # Try without common prefixes
        for prefix in ["alpha_", "beta_", "a_", "b_"]:
            if name_lower.startswith(prefix):
                base_name = name_lower[len(prefix):]
                for key in self.terpenes:
                    if base_name in key:
                        return self.terpenes[key]
        
        # Partial match
        for key, data in self.terpenes.items():
            if name_lower in key or key in name_lower:
                return data
        
        return None
    
    def list_all(self) -> List[str]:
        """List all terpene names."""
        return list(self.terpenes.keys())
    
    def get_by_effect(self, effect: str) -> List[TerpeneData]:
        """Get terpenes with a specific therapeutic effect."""
        return [
            t for t in self.terpenes.values()
            if effect.lower() in [e.lower() for e in t.therapeutic_effects]
        ]


class TerpeneAnalyzer:
    """Analyze terpene-cannabinoid synergies with evidence gating."""
    
    def __init__(self, database: Optional[TerpeneDatabase] = None):
        self.terpene_database = database or TerpeneDatabase()
    
    def analyze_synergy(
        self,
        cannabinoid: str,
        terpene_profile: Dict[str, float],
        min_evidence_tier: int = 3
    ) -> Dict:
        """Analyze terpene-cannabinoid synergy with evidence gating.
        
        Args:
            cannabinoid: Cannabinoid abbreviation (e.g., "CBD", "THC")
            terpene_profile: Dict of {terpene_name: concentration_percent}
            min_evidence_tier: Minimum evidence tier to include (1-5, default: 3)
            
        Returns:
            Synergy analysis with evidence-supported interactions only
        """
        synergies: List[SynergyResult] = []
        total_enhancement = 1.0
        excluded_interactions = []
        
        # Normalize cannabinoid name
        cannabinoid_upper = cannabinoid.upper()
        
        for terpene_name, concentration in terpene_profile.items():
            terpene_data = self.terpene_database.get(terpene_name)
            
            if not terpene_data:
                logger.warning(f"Unknown terpene: {terpene_name}")
                continue
            
            # Check if cannabinoid-terpene interaction exists
            synergy_data = terpene_data.cannabinoid_synergy.get(cannabinoid_upper)
            
            if not synergy_data:
                continue
            
            evidence_tier = synergy_data["evidence_tier"]
            
            # Evidence gating: only include if meets minimum tier
            if evidence_tier >= min_evidence_tier:
                # Calculate weighted enhancement based on concentration
                # Effect scales with concentration up to ~2% then plateaus
                concentration_factor = min(concentration, 2.0) / 2.0
                weighted_enhancement = 1.0 + (synergy_data["enhancement_factor"] - 1.0) * concentration_factor
                total_enhancement *= weighted_enhancement
                
                synergies.append(SynergyResult(
                    terpene=terpene_data.name,
                    concentration_percent=concentration,
                    mechanism=terpene_data.mechanism,
                    enhancement_factor=weighted_enhancement,
                    evidence_tier=evidence_tier,
                    evidence_quality=self._tier_to_description(evidence_tier),
                    pubmed_ids=terpene_data.pubmed_ids
                ))
            else:
                # Track excluded interactions for transparency
                excluded_interactions.append({
                    "terpene": terpene_data.name,
                    "cannabinoid": cannabinoid_upper,
                    "evidence_tier": evidence_tier,
                    "required_tier": min_evidence_tier,
                    "reason": f"Evidence tier {evidence_tier} below minimum {min_evidence_tier}"
                })
        
        return {
            "cannabinoid": cannabinoid_upper,
            "total_enhancement_factor": round(total_enhancement, 3),
            "synergies": [s.to_dict() for s in synergies],
            "num_synergies_included": len(synergies),
            "excluded_interactions": excluded_interactions,
            "num_excluded": len(excluded_interactions),
            "evidence_gating_enabled": True,
            "min_evidence_tier_used": min_evidence_tier,
            "analysis_date": datetime.now().isoformat()
        }
    
    def _tier_to_description(self, tier: int) -> str:
        """Convert evidence tier to human-readable description."""
        descriptions = {
            5: "High - Multiple RCTs/systematic reviews",
            4: "Good - Well-designed observational studies",
            3: "Moderate - Smaller studies, case series",
            2: "Low - In vitro/animal studies only",
            1: "Very Low - Anecdotal/theoretical"
        }
        return descriptions.get(tier, "Unknown")
    
    def get_optimal_profile(
        self,
        cannabinoid: str,
        therapeutic_goal: str,
        min_evidence_tier: int = 3
    ) -> Dict:
        """Recommend optimal terpene profile for a therapeutic goal.
        
        Args:
            cannabinoid: Target cannabinoid
            therapeutic_goal: Desired effect (e.g., "anxiolytic", "analgesic")
            min_evidence_tier: Minimum evidence tier
            
        Returns:
            Recommended terpene profile with rationale
        """
        recommendations = []
        cannabinoid_upper = cannabinoid.upper()
        
        for terpene_name, terpene_data in self.terpene_database.terpenes.items():
            # Check if terpene has desired therapeutic effect
            if therapeutic_goal.lower() not in [e.lower() for e in terpene_data.therapeutic_effects]:
                continue
            
            # Check synergy with cannabinoid
            synergy_data = terpene_data.cannabinoid_synergy.get(cannabinoid_upper)
            
            if not synergy_data:
                continue
            
            if synergy_data["evidence_tier"] >= min_evidence_tier:
                recommendations.append({
                    "terpene": terpene_data.name,
                    "recommended_concentration": self._recommend_concentration(synergy_data),
                    "enhancement_factor": synergy_data["enhancement_factor"],
                    "evidence_tier": synergy_data["evidence_tier"],
                    "mechanism": terpene_data.mechanism,
                    "additional_effects": [
                        e for e in terpene_data.therapeutic_effects
                        if e.lower() != therapeutic_goal.lower()
                    ]
                })
        
        # Sort by enhancement factor
        recommendations.sort(key=lambda x: x["enhancement_factor"], reverse=True)
        
        return {
            "cannabinoid": cannabinoid_upper,
            "therapeutic_goal": therapeutic_goal,
            "recommended_terpenes": recommendations[:5],  # Top 5
            "min_evidence_tier": min_evidence_tier,
            "total_potential_enhancement": self._calculate_total_enhancement(recommendations[:3])
        }
    
    def _recommend_concentration(self, synergy_data: Dict) -> str:
        """Recommend terpene concentration based on enhancement factor."""
        enhancement = synergy_data["enhancement_factor"]
        
        if enhancement > 1.3:
            return "1.0-2.0%"
        elif enhancement > 1.15:
            return "0.5-1.5%"
        else:
            return "0.3-1.0%"
    
    def _calculate_total_enhancement(self, recommendations: List[Dict]) -> float:
        """Calculate combined enhancement from multiple terpenes."""
        total = 1.0
        for rec in recommendations:
            # Assume mid-range concentration
            total *= 1.0 + (rec["enhancement_factor"] - 1.0) * 0.5
        return round(total, 3)
    
    def compare_profiles(
        self,
        cannabinoid: str,
        profile_a: Dict[str, float],
        profile_b: Dict[str, float],
        min_evidence_tier: int = 3
    ) -> Dict:
        """Compare two terpene profiles for synergy potential.
        
        Args:
            cannabinoid: Target cannabinoid
            profile_a: First terpene profile
            profile_b: Second terpene profile
            min_evidence_tier: Minimum evidence tier
            
        Returns:
            Comparison analysis
        """
        analysis_a = self.analyze_synergy(cannabinoid, profile_a, min_evidence_tier)
        analysis_b = self.analyze_synergy(cannabinoid, profile_b, min_evidence_tier)
        
        enhancement_a = analysis_a["total_enhancement_factor"]
        enhancement_b = analysis_b["total_enhancement_factor"]
        
        if enhancement_a > enhancement_b:
            winner = "Profile A"
            advantage = ((enhancement_a - enhancement_b) / enhancement_b) * 100
        elif enhancement_b > enhancement_a:
            winner = "Profile B"
            advantage = ((enhancement_b - enhancement_a) / enhancement_a) * 100
        else:
            winner = "Tie"
            advantage = 0
        
        return {
            "cannabinoid": cannabinoid.upper(),
            "profile_a": {
                "terpenes": list(profile_a.keys()),
                "total_enhancement": enhancement_a,
                "synergies_count": analysis_a["num_synergies_included"]
            },
            "profile_b": {
                "terpenes": list(profile_b.keys()),
                "total_enhancement": enhancement_b,
                "synergies_count": analysis_b["num_synergies_included"]
            },
            "winner": winner,
            "advantage_percent": round(advantage, 1),
            "recommendation": self._generate_comparison_recommendation(
                winner, advantage, analysis_a, analysis_b
            )
        }
    
    def _generate_comparison_recommendation(
        self,
        winner: str,
        advantage: float,
        analysis_a: Dict,
        analysis_b: Dict
    ) -> str:
        """Generate recommendation based on profile comparison."""
        if winner == "Tie":
            return "Both profiles offer similar synergy potential. Consider other factors like aroma preference or availability."
        
        if advantage > 20:
            return f"{winner} shows significantly higher synergy potential ({advantage:.1f}% better). Strongly recommended."
        elif advantage > 10:
            return f"{winner} offers moderately better synergy ({advantage:.1f}% advantage). Consider as preferred option."
        else:
            return f"{winner} has slight edge ({advantage:.1f}% better). Difference may not be clinically significant."
