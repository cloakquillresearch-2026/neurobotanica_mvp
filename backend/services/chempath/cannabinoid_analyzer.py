"""
ChemPath Cannabis Module - Cannabinoid Analyzer

Cannabis-specific molecular analysis with entourage effect modeling
and terpene-cannabinoid synergy prediction.

Trade Secret: Synergy algorithms derived from indigenous cannabis knowledge,
entourage effect coefficients, receptor binding synergy matrices.

Schedule III Compliance: Designed for FDA botanical drug pathway.
"""

from enum import Enum
from typing import Dict, List, Optional, Tuple, Set
from pydantic import BaseModel, Field
from datetime import datetime
import uuid
import math


class CannabinoidType(str, Enum):
    """Major cannabinoid classes."""
    # Primary Cannabinoids
    THC = "delta-9-THC"
    THCA = "THCA"
    CBD = "CBD"
    CBDA = "CBDA"
    CBG = "CBG"
    CBGA = "CBGA"
    CBN = "CBN"
    CBC = "CBC"
    CBCA = "CBCA"
    THCV = "THCV"
    CBDV = "CBDV"
    
    # Minor Cannabinoids
    DELTA8_THC = "delta-8-THC"
    THCP = "THCP"
    CBDP = "CBDP"
    CBL = "CBL"
    CBE = "CBE"
    CBT = "CBT"


class TerpeneType(str, Enum):
    """Major cannabis terpenes."""
    # Monoterpenes
    MYRCENE = "myrcene"
    LIMONENE = "limonene"
    PINENE_ALPHA = "alpha-pinene"
    PINENE_BETA = "beta-pinene"
    TERPINOLENE = "terpinolene"
    OCIMENE = "ocimene"
    
    # Sesquiterpenes
    CARYOPHYLLENE = "beta-caryophyllene"
    HUMULENE = "humulene"
    BISABOLOL = "bisabolol"
    GUAIOL = "guaiol"
    NEROLIDOL = "nerolidol"
    
    # Other
    LINALOOL = "linalool"
    TERPINEOL = "terpineol"
    EUCALYPTOL = "eucalyptol"
    GERANIOL = "geraniol"
    BORNEOL = "borneol"


class ReceptorTarget(str, Enum):
    """Endocannabinoid system receptors and related targets."""
    CB1 = "CB1"              # Primary psychoactive receptor
    CB2 = "CB2"              # Immune/inflammatory receptor
    TRPV1 = "TRPV1"          # Vanilloid receptor (pain/heat)
    GPR55 = "GPR55"          # Orphan receptor (bone, cancer)
    GPR18 = "GPR18"          # Microglial receptor
    PPAR_GAMMA = "PPAR-γ"    # Nuclear receptor (metabolic)
    PPAR_ALPHA = "PPAR-α"    # Nuclear receptor (lipid)
    SEROTONIN_5HT1A = "5-HT1A"  # Anxiolytic effects
    ADENOSINE_A2A = "A2A"    # Sleep/sedation
    GABA_A = "GABA-A"        # Sedation/anxiolytic
    TRPA1 = "TRPA1"          # Pain receptor
    TRPM8 = "TRPM8"          # Cold/menthol receptor


class TherapeuticCategory(str, Enum):
    """Therapeutic effect categories."""
    ANALGESIC = "analgesic"
    ANXIOLYTIC = "anxiolytic"
    ANTI_INFLAMMATORY = "anti_inflammatory"
    ANTIEMETIC = "antiemetic"
    ANTISPASMODIC = "antispasmodic"
    NEUROPROTECTIVE = "neuroprotective"
    SEDATIVE = "sedative"
    APPETITE_STIMULANT = "appetite_stimulant"
    ANTICONVULSANT = "anticonvulsant"
    ANTIDEPRESSANT = "antidepressant"
    ANTITUMOR = "antitumor"
    ANTIBACTERIAL = "antibacterial"


class CannabinoidProfile(BaseModel):
    """Chemical profile of cannabinoids in a formulation."""
    profile_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    
    # Cannabinoid concentrations (percentage)
    cannabinoids: Dict[str, float] = Field(default_factory=dict)
    
    # Derived metrics
    total_cannabinoids: float = 0.0
    thc_cbd_ratio: Optional[float] = None
    psychoactive_percentage: float = 0.0
    
    # Analysis metadata
    analysis_date: datetime = Field(default_factory=datetime.utcnow)
    analysis_method: str = "HPLC"


class TerpeneProfile(BaseModel):
    """Chemical profile of terpenes in a formulation."""
    profile_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    
    # Terpene concentrations (percentage)
    terpenes: Dict[str, float] = Field(default_factory=dict)
    
    # Derived metrics
    total_terpenes: float = 0.0
    dominant_terpene: Optional[str] = None
    monoterpene_percentage: float = 0.0
    sesquiterpene_percentage: float = 0.0
    
    # Analysis metadata
    analysis_date: datetime = Field(default_factory=datetime.utcnow)
    analysis_method: str = "GC-MS"


class SynergyInteraction(BaseModel):
    """Individual synergy interaction between compounds."""
    compound_a: str
    compound_b: str
    interaction_type: str  # "synergistic", "additive", "antagonistic"
    synergy_coefficient: float  # >1 synergistic, =1 additive, <1 antagonistic
    affected_receptors: List[str] = Field(default_factory=list)
    therapeutic_impact: List[str] = Field(default_factory=list)
    evidence_strength: str = "moderate"  # "strong", "moderate", "preliminary"
    traditional_knowledge_source: Optional[str] = None


class EntourageProfile(BaseModel):
    """Complete entourage effect analysis."""
    profile_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    formulation_name: str
    
    # Input profiles
    cannabinoid_profile: CannabinoidProfile
    terpene_profile: TerpeneProfile
    
    # Entourage analysis
    overall_entourage_score: float = 0.0  # 0-100 scale
    synergy_interactions: List[SynergyInteraction] = Field(default_factory=list)
    
    # Receptor binding predictions
    receptor_activity: Dict[str, float] = Field(default_factory=dict)  # 0-1 scale
    
    # Therapeutic predictions
    therapeutic_predictions: Dict[str, float] = Field(default_factory=dict)  # 0-1 scale
    primary_indication: Optional[str] = None
    secondary_indications: List[str] = Field(default_factory=list)
    
    # Traditional knowledge attribution
    tk_attributed: bool = False
    tk_communities: List[str] = Field(default_factory=list)
    
    # Metadata
    analysis_timestamp: datetime = Field(default_factory=datetime.utcnow)
    confidence_score: float = 0.0


class FormulationRecommendation(BaseModel):
    """Recommended formulation adjustments."""
    recommendation_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    target_indication: str
    
    # Current vs recommended
    current_entourage_score: float
    projected_entourage_score: float
    improvement_percentage: float
    
    # Specific adjustments
    cannabinoid_adjustments: Dict[str, Dict[str, float]] = Field(default_factory=dict)
    terpene_adjustments: Dict[str, Dict[str, float]] = Field(default_factory=dict)
    
    # Rationale
    rationale: str
    traditional_preparation_reference: Optional[str] = None
    clinical_evidence_citations: List[str] = Field(default_factory=list)


class CannabinoidAnalyzer:
    """
    Cannabis-specific molecular analysis.
    Entourage effect modeling, terpene-cannabinoid synergy prediction.
    
    Trade Secret: Synergy algorithms from indigenous cannabis knowledge,
    receptor binding matrices, entourage coefficients.
    """
    
    # ==========================================================================
    # TRADE SECRET: Cannabinoid-Receptor Binding Affinity Matrix
    # Values represent relative binding affinity (0-1 scale)
    # ==========================================================================
    _CANNABINOID_RECEPTOR_AFFINITY: Dict[str, Dict[str, float]] = {
        CannabinoidType.THC.value: {
            ReceptorTarget.CB1.value: 0.95,
            ReceptorTarget.CB2.value: 0.40,
            ReceptorTarget.GPR55.value: 0.60,
            ReceptorTarget.TRPV1.value: 0.30,
            ReceptorTarget.PPAR_GAMMA.value: 0.25,
        },
        CannabinoidType.CBD.value: {
            ReceptorTarget.CB1.value: 0.05,  # Negative allosteric modulator
            ReceptorTarget.CB2.value: 0.15,
            ReceptorTarget.SEROTONIN_5HT1A.value: 0.70,
            ReceptorTarget.TRPV1.value: 0.65,
            ReceptorTarget.GPR55.value: 0.40,
            ReceptorTarget.PPAR_GAMMA.value: 0.55,
            ReceptorTarget.ADENOSINE_A2A.value: 0.45,
        },
        CannabinoidType.CBG.value: {
            ReceptorTarget.CB1.value: 0.10,
            ReceptorTarget.CB2.value: 0.25,
            ReceptorTarget.TRPV1.value: 0.50,
            ReceptorTarget.PPAR_GAMMA.value: 0.60,
            ReceptorTarget.SEROTONIN_5HT1A.value: 0.40,
        },
        CannabinoidType.CBN.value: {
            ReceptorTarget.CB1.value: 0.25,
            ReceptorTarget.CB2.value: 0.30,
            ReceptorTarget.TRPV1.value: 0.20,
            ReceptorTarget.GABA_A.value: 0.55,  # Sedative effect
        },
        CannabinoidType.CBC.value: {
            ReceptorTarget.CB2.value: 0.20,
            ReceptorTarget.TRPV1.value: 0.70,
            ReceptorTarget.TRPA1.value: 0.65,
        },
        CannabinoidType.THCV.value: {
            ReceptorTarget.CB1.value: 0.40,  # Antagonist at low doses
            ReceptorTarget.CB2.value: 0.35,
            ReceptorTarget.GPR55.value: 0.50,
        },
    }
    
    # ==========================================================================
    # TRADE SECRET: Terpene-Receptor Modulation Matrix
    # ==========================================================================
    _TERPENE_RECEPTOR_MODULATION: Dict[str, Dict[str, float]] = {
        TerpeneType.MYRCENE.value: {
            ReceptorTarget.CB1.value: 0.30,  # Enhances THC absorption
            ReceptorTarget.GABA_A.value: 0.45,
            ReceptorTarget.TRPV1.value: 0.35,
        },
        TerpeneType.LIMONENE.value: {
            ReceptorTarget.SEROTONIN_5HT1A.value: 0.50,
            ReceptorTarget.ADENOSINE_A2A.value: 0.35,
        },
        TerpeneType.PINENE_ALPHA.value: {
            ReceptorTarget.CB1.value: -0.15,  # Counteracts THC memory impairment
            ReceptorTarget.GABA_A.value: 0.25,
        },
        TerpeneType.CARYOPHYLLENE.value: {
            ReceptorTarget.CB2.value: 0.65,  # Direct CB2 agonist
            ReceptorTarget.PPAR_GAMMA.value: 0.40,
        },
        TerpeneType.LINALOOL.value: {
            ReceptorTarget.GABA_A.value: 0.55,
            ReceptorTarget.SEROTONIN_5HT1A.value: 0.40,
        },
        TerpeneType.HUMULENE.value: {
            ReceptorTarget.CB2.value: 0.30,
            ReceptorTarget.PPAR_GAMMA.value: 0.35,
        },
        TerpeneType.BISABOLOL.value: {
            ReceptorTarget.TRPV1.value: 0.40,
            ReceptorTarget.GABA_A.value: 0.30,
        },
    }
    
    # ==========================================================================
    # TRADE SECRET: Synergy Coefficient Matrix (Compound Pairs)
    # >1.0 = synergistic, 1.0 = additive, <1.0 = antagonistic
    # ==========================================================================
    _SYNERGY_COEFFICIENTS: Dict[Tuple[str, str], float] = {
        # CBD-THC interactions
        (CannabinoidType.THC.value, CannabinoidType.CBD.value): 1.35,
        # Myrcene enhances THC
        (CannabinoidType.THC.value, TerpeneType.MYRCENE.value): 1.45,
        # Caryophyllene + CBD for inflammation
        (CannabinoidType.CBD.value, TerpeneType.CARYOPHYLLENE.value): 1.40,
        # Linalool + CBD for anxiety
        (CannabinoidType.CBD.value, TerpeneType.LINALOOL.value): 1.30,
        # Pinene counteracts THC memory effects
        (CannabinoidType.THC.value, TerpeneType.PINENE_ALPHA.value): 0.85,
        # Limonene + CBD for mood
        (CannabinoidType.CBD.value, TerpeneType.LIMONENE.value): 1.25,
        # CBN + Myrcene for sedation
        (CannabinoidType.CBN.value, TerpeneType.MYRCENE.value): 1.50,
        # CBG + Caryophyllene for inflammation
        (CannabinoidType.CBG.value, TerpeneType.CARYOPHYLLENE.value): 1.35,
        # THC + Limonene for mood elevation
        (CannabinoidType.THC.value, TerpeneType.LIMONENE.value): 1.20,
        # CBC + Myrcene for pain
        (CannabinoidType.CBC.value, TerpeneType.MYRCENE.value): 1.30,
    }
    
    # ==========================================================================
    # TRADE SECRET: Therapeutic Prediction Weights
    # ==========================================================================
    _THERAPEUTIC_WEIGHTS: Dict[str, Dict[str, float]] = {
        TherapeuticCategory.ANALGESIC.value: {
            ReceptorTarget.CB1.value: 0.30,
            ReceptorTarget.CB2.value: 0.25,
            ReceptorTarget.TRPV1.value: 0.25,
            ReceptorTarget.TRPA1.value: 0.10,
            ReceptorTarget.GABA_A.value: 0.10,
        },
        TherapeuticCategory.ANXIOLYTIC.value: {
            ReceptorTarget.SEROTONIN_5HT1A.value: 0.35,
            ReceptorTarget.GABA_A.value: 0.30,
            ReceptorTarget.CB1.value: 0.20,
            ReceptorTarget.ADENOSINE_A2A.value: 0.15,
        },
        TherapeuticCategory.ANTI_INFLAMMATORY.value: {
            ReceptorTarget.CB2.value: 0.40,
            ReceptorTarget.PPAR_GAMMA.value: 0.30,
            ReceptorTarget.TRPV1.value: 0.15,
            ReceptorTarget.GPR55.value: 0.15,
        },
        TherapeuticCategory.ANTIEMETIC.value: {
            ReceptorTarget.CB1.value: 0.50,
            ReceptorTarget.SEROTONIN_5HT1A.value: 0.30,
            ReceptorTarget.TRPV1.value: 0.20,
        },
        TherapeuticCategory.SEDATIVE.value: {
            ReceptorTarget.GABA_A.value: 0.40,
            ReceptorTarget.CB1.value: 0.25,
            ReceptorTarget.ADENOSINE_A2A.value: 0.20,
            ReceptorTarget.SEROTONIN_5HT1A.value: 0.15,
        },
        TherapeuticCategory.ANTICONVULSANT.value: {
            ReceptorTarget.GABA_A.value: 0.35,
            ReceptorTarget.TRPV1.value: 0.25,
            ReceptorTarget.GPR55.value: 0.20,
            ReceptorTarget.SEROTONIN_5HT1A.value: 0.20,
        },
        TherapeuticCategory.NEUROPROTECTIVE.value: {
            ReceptorTarget.CB1.value: 0.25,
            ReceptorTarget.CB2.value: 0.25,
            ReceptorTarget.PPAR_GAMMA.value: 0.25,
            ReceptorTarget.SEROTONIN_5HT1A.value: 0.25,
        },
        TherapeuticCategory.APPETITE_STIMULANT.value: {
            ReceptorTarget.CB1.value: 0.60,
            ReceptorTarget.GPR55.value: 0.25,
            ReceptorTarget.GABA_A.value: 0.15,
        },
    }
    
    # Traditional preparation references (from TKPath communities)
    _TRADITIONAL_PREPARATIONS: Dict[str, Dict] = {
        "ayurvedic_bhang": {
            "community": "Indian Ayurvedic Practitioners",
            "optimal_thc_cbd_ratio": (1, 2),  # 1:2 THC:CBD
            "key_terpenes": ["myrcene", "linalool"],
            "primary_use": "anxiolytic",
        },
        "jamaican_ganja": {
            "community": "Jamaican Rastafari Communities",
            "optimal_thc_cbd_ratio": (10, 1),  # High THC
            "key_terpenes": ["myrcene", "limonene", "humulene"],
            "primary_use": "mood_elevation",
        },
        "moroccan_kif": {
            "community": "Moroccan Rif Mountain Cultivators",
            "optimal_thc_cbd_ratio": (5, 1),
            "key_terpenes": ["myrcene", "caryophyllene"],
            "primary_use": "relaxation",
        },
        "himalayan_charras": {
            "community": "Hindu Kush Traditional Cultivators",
            "optimal_thc_cbd_ratio": (15, 1),
            "key_terpenes": ["myrcene", "pinene", "limonene"],
            "primary_use": "analgesic",
        },
    }
    
    def __init__(self):
        self._analysis_cache: Dict[str, EntourageProfile] = {}
    
    def analyze_entourage_effect(
        self,
        cannabinoids: Dict[str, float],
        terpenes: Dict[str, float],
        formulation_name: Optional[str] = None,
    ) -> EntourageProfile:
        """
        Predict synergistic effects based on cannabinoid-terpene interactions.
        
        Trade Secret: Synergy algorithms from indigenous cannabis knowledge,
        receptor binding predictions, therapeutic effect modeling.
        
        Args:
            cannabinoids: Dict of cannabinoid concentrations (%)
            terpenes: Dict of terpene concentrations (%)
            formulation_name: Optional name for the formulation
            
        Returns:
            EntourageProfile with complete synergy analysis
        """
        # Build profiles
        cannabinoid_profile = self._build_cannabinoid_profile(cannabinoids)
        terpene_profile = self._build_terpene_profile(terpenes)
        
        # Create entourage profile
        profile = EntourageProfile(
            formulation_name=formulation_name or f"Formulation_{uuid.uuid4().hex[:8]}",
            cannabinoid_profile=cannabinoid_profile,
            terpene_profile=terpene_profile,
        )
        
        # Calculate receptor activity
        profile.receptor_activity = self._calculate_receptor_activity(
            cannabinoids, terpenes
        )
        
        # Identify synergy interactions
        profile.synergy_interactions = self._identify_synergy_interactions(
            cannabinoids, terpenes
        )
        
        # Calculate overall entourage score
        profile.overall_entourage_score = self._calculate_entourage_score(
            profile.synergy_interactions,
            cannabinoid_profile,
            terpene_profile,
        )
        
        # Predict therapeutic effects
        profile.therapeutic_predictions = self._predict_therapeutic_effects(
            profile.receptor_activity
        )
        
        # Identify primary/secondary indications
        sorted_predictions = sorted(
            profile.therapeutic_predictions.items(),
            key=lambda x: x[1],
            reverse=True
        )
        if sorted_predictions:
            profile.primary_indication = sorted_predictions[0][0]
            profile.secondary_indications = [
                ind for ind, score in sorted_predictions[1:4] if score > 0.3
            ]
        
        # Check for TK attribution
        tk_match = self._match_traditional_preparation(cannabinoids, terpenes)
        if tk_match:
            profile.tk_attributed = True
            profile.tk_communities = [tk_match["community"]]
        
        # Calculate confidence
        profile.confidence_score = self._calculate_confidence(
            cannabinoid_profile, terpene_profile
        )
        
        # Cache result
        self._analysis_cache[profile.profile_id] = profile
        
        return profile
    
    def _build_cannabinoid_profile(
        self, cannabinoids: Dict[str, float]
    ) -> CannabinoidProfile:
        """Build cannabinoid profile from concentrations."""
        profile = CannabinoidProfile(cannabinoids=cannabinoids)
        
        profile.total_cannabinoids = sum(cannabinoids.values())
        
        # Calculate THC:CBD ratio
        thc = cannabinoids.get("THC", cannabinoids.get("delta-9-THC", 0))
        cbd = cannabinoids.get("CBD", 0)
        if cbd > 0:
            profile.thc_cbd_ratio = thc / cbd
        
        # Calculate psychoactive percentage
        psychoactive = ["THC", "delta-9-THC", "delta-8-THC", "THCV", "THCP"]
        profile.psychoactive_percentage = sum(
            cannabinoids.get(c, 0) for c in psychoactive
        )
        
        return profile
    
    def _build_terpene_profile(self, terpenes: Dict[str, float]) -> TerpeneProfile:
        """Build terpene profile from concentrations."""
        profile = TerpeneProfile(terpenes=terpenes)
        
        profile.total_terpenes = sum(terpenes.values())
        
        # Identify dominant terpene
        if terpenes:
            profile.dominant_terpene = max(terpenes, key=terpenes.get)
        
        # Classify terpenes
        monoterpenes = ["myrcene", "limonene", "alpha-pinene", "beta-pinene", 
                       "terpinolene", "ocimene"]
        sesquiterpenes = ["beta-caryophyllene", "humulene", "bisabolol", 
                         "guaiol", "nerolidol"]
        
        profile.monoterpene_percentage = sum(
            terpenes.get(t, 0) for t in monoterpenes
        )
        profile.sesquiterpene_percentage = sum(
            terpenes.get(t, 0) for t in sesquiterpenes
        )
        
        return profile
    
    def _calculate_receptor_activity(
        self,
        cannabinoids: Dict[str, float],
        terpenes: Dict[str, float],
    ) -> Dict[str, float]:
        """
        Calculate predicted receptor activity levels.
        Trade Secret: Receptor binding algorithm.
        """
        receptor_activity: Dict[str, float] = {}
        
        # Initialize all receptors
        for receptor in ReceptorTarget:
            receptor_activity[receptor.value] = 0.0
        
        # Add cannabinoid contributions
        for cannabinoid, concentration in cannabinoids.items():
            normalized_conc = min(concentration / 30.0, 1.0)  # Normalize to 30%
            
            if cannabinoid in self._CANNABINOID_RECEPTOR_AFFINITY:
                for receptor, affinity in self._CANNABINOID_RECEPTOR_AFFINITY[cannabinoid].items():
                    receptor_activity[receptor] += affinity * normalized_conc
        
        # Add terpene modulation
        for terpene, concentration in terpenes.items():
            normalized_conc = min(concentration / 2.0, 1.0)  # Normalize to 2%
            
            if terpene in self._TERPENE_RECEPTOR_MODULATION:
                for receptor, modulation in self._TERPENE_RECEPTOR_MODULATION[terpene].items():
                    receptor_activity[receptor] += modulation * normalized_conc * 0.5
        
        # Normalize to 0-1 scale
        for receptor in receptor_activity:
            receptor_activity[receptor] = min(receptor_activity[receptor], 1.0)
        
        return receptor_activity
    
    def _identify_synergy_interactions(
        self,
        cannabinoids: Dict[str, float],
        terpenes: Dict[str, float],
    ) -> List[SynergyInteraction]:
        """
        Identify synergy interactions between compounds.
        Trade Secret: Synergy detection algorithm.
        """
        interactions = []
        all_compounds = {**cannabinoids, **terpenes}
        compound_list = list(all_compounds.keys())
        
        for i, compound_a in enumerate(compound_list):
            for compound_b in compound_list[i+1:]:
                # Check if we have synergy data for this pair
                pair_key = (compound_a, compound_b)
                reverse_key = (compound_b, compound_a)
                
                coefficient = self._SYNERGY_COEFFICIENTS.get(
                    pair_key, 
                    self._SYNERGY_COEFFICIENTS.get(reverse_key, 1.0)
                )
                
                if coefficient != 1.0:
                    # Determine interaction type
                    if coefficient > 1.0:
                        interaction_type = "synergistic"
                    elif coefficient < 1.0:
                        interaction_type = "antagonistic"
                    else:
                        interaction_type = "additive"
                    
                    interaction = SynergyInteraction(
                        compound_a=compound_a,
                        compound_b=compound_b,
                        interaction_type=interaction_type,
                        synergy_coefficient=coefficient,
                    )
                    
                    # Identify affected receptors
                    receptors_a = set()
                    receptors_b = set()
                    
                    if compound_a in self._CANNABINOID_RECEPTOR_AFFINITY:
                        receptors_a = set(self._CANNABINOID_RECEPTOR_AFFINITY[compound_a].keys())
                    if compound_a in self._TERPENE_RECEPTOR_MODULATION:
                        receptors_a = set(self._TERPENE_RECEPTOR_MODULATION[compound_a].keys())
                    if compound_b in self._CANNABINOID_RECEPTOR_AFFINITY:
                        receptors_b = set(self._CANNABINOID_RECEPTOR_AFFINITY[compound_b].keys())
                    if compound_b in self._TERPENE_RECEPTOR_MODULATION:
                        receptors_b = set(self._TERPENE_RECEPTOR_MODULATION[compound_b].keys())
                    
                    interaction.affected_receptors = list(receptors_a | receptors_b)
                    
                    interactions.append(interaction)
        
        return interactions
    
    def _calculate_entourage_score(
        self,
        interactions: List[SynergyInteraction],
        cannabinoid_profile: CannabinoidProfile,
        terpene_profile: TerpeneProfile,
    ) -> float:
        """
        Calculate overall entourage effect score (0-100).
        Trade Secret: Entourage scoring algorithm.
        """
        score = 50.0  # Baseline
        
        # Synergy bonus
        synergistic_count = sum(
            1 for i in interactions if i.interaction_type == "synergistic"
        )
        antagonistic_count = sum(
            1 for i in interactions if i.interaction_type == "antagonistic"
        )
        
        score += synergistic_count * 5
        score -= antagonistic_count * 3
        
        # Diversity bonus (more compounds = higher entourage potential)
        compound_count = len(cannabinoid_profile.cannabinoids) + len(terpene_profile.terpenes)
        score += min(compound_count * 2, 20)
        
        # Terpene presence bonus
        if terpene_profile.total_terpenes > 0.5:
            score += 10
        if terpene_profile.total_terpenes > 1.5:
            score += 10
        
        # CBD:THC ratio bonus (balanced ratios enhance entourage)
        if cannabinoid_profile.thc_cbd_ratio:
            if 0.5 <= cannabinoid_profile.thc_cbd_ratio <= 2.0:
                score += 15  # Balanced ratio
        
        # High synergy coefficient bonus
        avg_coefficient = 1.0
        if interactions:
            avg_coefficient = sum(i.synergy_coefficient for i in interactions) / len(interactions)
        if avg_coefficient > 1.2:
            score += 10
        
        return max(0, min(100, score))
    
    def _predict_therapeutic_effects(
        self,
        receptor_activity: Dict[str, float],
    ) -> Dict[str, float]:
        """
        Predict therapeutic effects based on receptor activity.
        Trade Secret: Therapeutic prediction algorithm.
        """
        predictions: Dict[str, float] = {}
        
        for category, weights in self._THERAPEUTIC_WEIGHTS.items():
            score = 0.0
            for receptor, weight in weights.items():
                score += receptor_activity.get(receptor, 0) * weight
            predictions[category] = min(score, 1.0)
        
        return predictions
    
    def _match_traditional_preparation(
        self,
        cannabinoids: Dict[str, float],
        terpenes: Dict[str, float],
    ) -> Optional[Dict]:
        """
        Match formulation to traditional preparation for TK attribution.
        Trade Secret: Traditional preparation matching algorithm.
        """
        thc = cannabinoids.get("THC", cannabinoids.get("delta-9-THC", 0))
        cbd = cannabinoids.get("CBD", 0)
        
        for prep_name, prep_info in self._TRADITIONAL_PREPARATIONS.items():
            # Check THC:CBD ratio
            target_ratio = prep_info["optimal_thc_cbd_ratio"]
            if cbd > 0:
                actual_ratio = thc / cbd
                target = target_ratio[0] / target_ratio[1]
                if 0.5 * target <= actual_ratio <= 2.0 * target:
                    # Check key terpenes
                    key_terpenes_present = sum(
                        1 for t in prep_info["key_terpenes"] 
                        if terpenes.get(t, 0) > 0.1
                    )
                    if key_terpenes_present >= len(prep_info["key_terpenes"]) * 0.5:
                        return {
                            "preparation": prep_name,
                            "community": prep_info["community"],
                        }
        
        return None
    
    def _calculate_confidence(
        self,
        cannabinoid_profile: CannabinoidProfile,
        terpene_profile: TerpeneProfile,
    ) -> float:
        """Calculate confidence score for analysis."""
        confidence = 0.5
        
        # More data = higher confidence
        if len(cannabinoid_profile.cannabinoids) >= 3:
            confidence += 0.15
        if len(terpene_profile.terpenes) >= 3:
            confidence += 0.15
        
        # Known compounds = higher confidence
        known_cannabinoids = sum(
            1 for c in cannabinoid_profile.cannabinoids 
            if c in self._CANNABINOID_RECEPTOR_AFFINITY
        )
        known_terpenes = sum(
            1 for t in terpene_profile.terpenes 
            if t in self._TERPENE_RECEPTOR_MODULATION
        )
        
        confidence += min(known_cannabinoids * 0.05, 0.15)
        confidence += min(known_terpenes * 0.05, 0.15)
        
        return min(confidence, 1.0)
    
    def get_formulation_recommendations(
        self,
        current_profile: EntourageProfile,
        target_indication: str,
    ) -> FormulationRecommendation:
        """
        Generate recommendations to optimize formulation for target indication.
        
        Trade Secret: Optimization algorithm based on traditional preparation methods.
        """
        recommendation = FormulationRecommendation(
            target_indication=target_indication,
            current_entourage_score=current_profile.overall_entourage_score,
            projected_entourage_score=0.0,
            improvement_percentage=0.0,
            rationale="",
        )
        
        current_therapeutic = current_profile.therapeutic_predictions.get(target_indication, 0)
        
        # Get optimal receptor profile for indication
        if target_indication in self._THERAPEUTIC_WEIGHTS:
            optimal_receptors = self._THERAPEUTIC_WEIGHTS[target_indication]
            
            # Identify cannabinoids that enhance target receptors
            for cannabinoid, receptors in self._CANNABINOID_RECEPTOR_AFFINITY.items():
                for receptor, weight in optimal_receptors.items():
                    if receptor in receptors and receptors[receptor] > 0.3:
                        current = current_profile.cannabinoid_profile.cannabinoids.get(cannabinoid, 0)
                        if current < 5.0:  # Room for increase
                            recommendation.cannabinoid_adjustments[cannabinoid] = {
                                "current": current,
                                "recommended": min(current + 5.0, 20.0),
                                "reason": f"Enhances {receptor} activity for {target_indication}",
                            }
            
            # Identify terpenes that enhance target receptors
            for terpene, receptors in self._TERPENE_RECEPTOR_MODULATION.items():
                for receptor, weight in optimal_receptors.items():
                    if receptor in receptors and receptors[receptor] > 0.3:
                        current = current_profile.terpene_profile.terpenes.get(terpene, 0)
                        if current < 0.5:  # Room for increase
                            recommendation.terpene_adjustments[terpene] = {
                                "current": current,
                                "recommended": min(current + 0.3, 1.0),
                                "reason": f"Enhances {receptor} activity for {target_indication}",
                            }
        
        # Project improvement
        projected_improvement = 0.15  # Base improvement estimate
        projected_improvement += len(recommendation.cannabinoid_adjustments) * 0.05
        projected_improvement += len(recommendation.terpene_adjustments) * 0.03
        
        recommendation.projected_entourage_score = min(
            current_profile.overall_entourage_score * (1 + projected_improvement),
            100.0
        )
        recommendation.improvement_percentage = projected_improvement * 100
        
        recommendation.rationale = (
            f"Optimizing for {target_indication} by enhancing receptor activity "
            f"through targeted cannabinoid/terpene adjustments. "
            f"Projected {recommendation.improvement_percentage:.1f}% improvement in therapeutic potential."
        )
        
        return recommendation
    
    def get_analysis(self, profile_id: str) -> Optional[EntourageProfile]:
        """Retrieve cached analysis by ID."""
        return self._analysis_cache.get(profile_id)
