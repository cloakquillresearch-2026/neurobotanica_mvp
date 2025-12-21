"""
TKPath Verification Lab

In-house analytical lab verifying TK attribution claims.
Prevents companies from taking TK credit without actual use.

Trade Secret: Matching algorithms for natural product variation,
provenance verification protocols, chemical fingerprinting methods.

Lab Functions:
1. Chemical Fingerprinting - Verify formulations match traditional preparations
2. Provenance Tracking - Blockchain-verified supply chain from indigenous growers
3. Quality Control - GMP compliance while preserving TK methods
"""

from enum import Enum
from typing import Dict, List, Optional, Tuple
from pydantic import BaseModel, Field
from datetime import datetime
import uuid
import hashlib


class VerificationStatus(str, Enum):
    """Status of verification analysis."""
    PENDING = "pending"
    IN_PROGRESS = "in_progress"
    VERIFIED = "verified"
    PARTIAL_MATCH = "partial_match"
    NO_MATCH = "no_match"
    FAILED = "failed"
    DISPUTED = "disputed"


class ProvenanceStatus(str, Enum):
    """Status of provenance verification."""
    UNVERIFIED = "unverified"
    BLOCKCHAIN_VERIFIED = "blockchain_verified"
    THIRD_PARTY_VERIFIED = "third_party_verified"
    SELF_ATTESTED = "self_attested"
    DISPUTED = "disputed"


class ChemicalProfile(BaseModel):
    """Chemical fingerprint of a preparation."""
    profile_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    
    # Cannabinoid profile
    cannabinoids: Dict[str, float] = Field(default_factory=dict)  # e.g., {"THC": 18.5, "CBD": 2.1}
    
    # Terpene profile
    terpenes: Dict[str, float] = Field(default_factory=dict)  # e.g., {"myrcene": 0.45}
    
    # Flavonoid profile
    flavonoids: Dict[str, float] = Field(default_factory=dict)
    
    # Other markers
    other_markers: Dict[str, float] = Field(default_factory=dict)
    
    # Metadata
    analysis_method: str = "HPLC_GC-MS"
    analysis_date: datetime = Field(default_factory=datetime.utcnow)
    lab_id: Optional[str] = None
    
    def get_total_cannabinoids(self) -> float:
        """Calculate total cannabinoid content."""
        return sum(self.cannabinoids.values())
    
    def get_total_terpenes(self) -> float:
        """Calculate total terpene content."""
        return sum(self.terpenes.values())


class TraditionalPreparation(BaseModel):
    """Reference profile for traditional preparation method."""
    preparation_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    preparation_name: str
    community_origin: str
    region: str
    
    # Expected chemical profile ranges
    expected_cannabinoid_ranges: Dict[str, Tuple[float, float]] = Field(default_factory=dict)
    expected_terpene_ranges: Dict[str, Tuple[float, float]] = Field(default_factory=dict)
    
    # Traditional method details
    extraction_method: str
    processing_steps: List[str] = Field(default_factory=list)
    timing_requirements: Optional[str] = None
    
    # Verification parameters
    minimum_match_threshold: float = 0.70  # 70% match required


class ProvenanceRecord(BaseModel):
    """Blockchain-verified supply chain record."""
    record_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    
    # Source information
    source_community: str
    source_region: str
    source_country: str
    gps_coordinates: Optional[Tuple[float, float]] = None
    
    # Cultivation details
    cultivation_method: str
    harvest_date: Optional[datetime] = None
    harvest_timing_notes: Optional[str] = None  # e.g., "full moon harvest"
    
    # Processing chain
    processing_location: Optional[str] = None
    processing_method: Optional[str] = None
    processing_date: Optional[datetime] = None
    
    # Verification
    status: ProvenanceStatus = ProvenanceStatus.UNVERIFIED
    blockchain_tx_hash: Optional[str] = None
    verification_date: Optional[datetime] = None
    verifier_id: Optional[str] = None
    
    def generate_blockchain_hash(self) -> str:
        """Generate hash for blockchain recording."""
        data = (
            f"{self.record_id}:{self.source_community}:{self.harvest_date}:"
            f"{self.processing_method}:{self.gps_coordinates}"
        )
        return hashlib.sha256(data.encode()).hexdigest()


class VerificationReport(BaseModel):
    """Complete verification report for TK attribution claim."""
    report_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    
    # Claim details
    claimed_formulation: str
    claimed_traditional_method: str
    claimed_community: str
    
    # Analysis results
    status: VerificationStatus = VerificationStatus.PENDING
    overall_match_score: float = 0.0
    
    # Component scores
    chemical_match_score: float = 0.0
    method_match_score: float = 0.0
    provenance_score: float = 0.0
    
    # Detailed findings
    chemical_profile_analysis: Optional[ChemicalProfile] = None
    provenance_record: Optional[ProvenanceRecord] = None
    
    # Discrepancies found
    discrepancies: List[str] = Field(default_factory=list)
    recommendations: List[str] = Field(default_factory=list)
    
    # Metadata
    analysis_date: datetime = Field(default_factory=datetime.utcnow)
    analyst_id: Optional[str] = None
    lab_certification: Optional[str] = None
    
    # Verification hash
    verification_hash: Optional[str] = None
    
    def generate_verification_hash(self) -> str:
        """Generate tamper-evident verification hash."""
        data = (
            f"{self.report_id}:{self.claimed_formulation}:"
            f"{self.overall_match_score}:{self.status.value}:"
            f"{self.analysis_date.isoformat()}"
        )
        self.verification_hash = hashlib.sha256(data.encode()).hexdigest()
        return self.verification_hash


class TKVerificationLab:
    """
    In-house analytical lab verifying TK attribution claims.
    Prevents companies from taking TK credit without actual use.
    
    Trade Secret: Matching algorithms for natural product variation,
    provenance verification protocols, chemical fingerprinting methods.
    """
    
    # Match thresholds (Trade Secret)
    _CHEMICAL_MATCH_THRESHOLD = 0.70
    _METHOD_MATCH_THRESHOLD = 0.60
    _PROVENANCE_THRESHOLD = 0.80
    _OVERALL_VERIFICATION_THRESHOLD = 0.75
    
    # Known traditional preparations database
    _TRADITIONAL_PREPARATIONS: Dict[str, TraditionalPreparation] = {}
    
    def __init__(self):
        self._verification_reports: Dict[str, VerificationReport] = {}
        self._provenance_records: Dict[str, ProvenanceRecord] = {}
        self._initialize_traditional_preparations()
    
    def _initialize_traditional_preparations(self):
        """Initialize database of known traditional preparations."""
        self._TRADITIONAL_PREPARATIONS = {
            "charras_himalayan": TraditionalPreparation(
                preparation_name="Himalayan Charras",
                community_origin="Hindu Kush Traditional Cultivators",
                region="Hindu Kush Mountains",
                expected_cannabinoid_ranges={
                    "THC": (15.0, 25.0),
                    "CBD": (0.5, 3.0),
                    "CBN": (0.1, 1.0),
                },
                expected_terpene_ranges={
                    "myrcene": (0.3, 0.8),
                    "limonene": (0.2, 0.5),
                    "pinene": (0.1, 0.4),
                },
                extraction_method="hand_rubbing",
                processing_steps=[
                    "Fresh flower hand rubbing",
                    "Collection of resin",
                    "Manual pressing",
                    "Aging (2-4 weeks)",
                ],
            ),
            "moroccan_hash": TraditionalPreparation(
                preparation_name="Moroccan Traditional Hash",
                community_origin="Moroccan Rif Mountain Cultivators",
                region="Rif Mountains",
                expected_cannabinoid_ranges={
                    "THC": (10.0, 20.0),
                    "CBD": (1.0, 5.0),
                },
                expected_terpene_ranges={
                    "myrcene": (0.2, 0.6),
                    "caryophyllene": (0.3, 0.7),
                },
                extraction_method="dry_sift",
                processing_steps=[
                    "Dry curing of flowers",
                    "Dry sifting through screens",
                    "Collection of kief",
                    "Pressing into blocks",
                ],
            ),
            "bhang_ayurvedic": TraditionalPreparation(
                preparation_name="Ayurvedic Bhang",
                community_origin="Indian Ayurvedic Practitioners",
                region="Northern India",
                expected_cannabinoid_ranges={
                    "THC": (5.0, 15.0),
                    "CBD": (2.0, 8.0),
                },
                expected_terpene_ranges={
                    "linalool": (0.2, 0.5),
                    "myrcene": (0.2, 0.5),
                },
                extraction_method="water_extraction",
                processing_steps=[
                    "Grinding of leaves and flowers",
                    "Soaking in water",
                    "Straining and filtering",
                    "Addition of spices (optional)",
                ],
            ),
            "ganja_jamaican": TraditionalPreparation(
                preparation_name="Jamaican Ganja",
                community_origin="Jamaican Rastafari Communities",
                region="Jamaica",
                expected_cannabinoid_ranges={
                    "THC": (12.0, 22.0),
                    "CBD": (0.5, 2.0),
                },
                expected_terpene_ranges={
                    "myrcene": (0.3, 0.7),
                    "limonene": (0.2, 0.6),
                    "humulene": (0.1, 0.3),
                },
                extraction_method="sun_dried",
                processing_steps=[
                    "Outdoor sun cultivation",
                    "Sun drying of flowers",
                    "Hand trimming",
                    "Curing",
                ],
            ),
        }
    
    def verify_traditional_preparation(
        self,
        claimed_formulation: str,
        traditional_recipe: str,
        chemical_profile: Optional[ChemicalProfile] = None,
        provenance: Optional[ProvenanceRecord] = None,
    ) -> VerificationReport:
        """
        Compare chemical profiles to validate TK attribution.
        
        Trade Secret: Matching algorithms for natural product variation.
        
        Args:
            claimed_formulation: Name of formulation being verified
            traditional_recipe: Key of traditional preparation to match
            chemical_profile: Actual chemical analysis of formulation
            provenance: Supply chain provenance record
            
        Returns:
            VerificationReport with match scores and findings
        """
        report = VerificationReport(
            claimed_formulation=claimed_formulation,
            claimed_traditional_method=traditional_recipe,
            claimed_community="",
        )
        
        # Get traditional preparation reference
        traditional = self._TRADITIONAL_PREPARATIONS.get(traditional_recipe)
        if not traditional:
            report.status = VerificationStatus.FAILED
            report.discrepancies.append(
                f"Unknown traditional preparation: {traditional_recipe}"
            )
            return report
        
        report.claimed_community = traditional.community_origin
        
        # Chemical profile matching
        if chemical_profile:
            report.chemical_profile_analysis = chemical_profile
            report.chemical_match_score = self._calculate_chemical_match(
                chemical_profile, traditional
            )
        else:
            report.discrepancies.append("No chemical profile provided for verification")
            report.recommendations.append("Submit samples for HPLC/GC-MS analysis")
        
        # Method matching (simplified - would analyze preparation steps)
        report.method_match_score = self._calculate_method_match(
            claimed_formulation, traditional
        )
        
        # Provenance verification
        if provenance:
            report.provenance_record = provenance
            report.provenance_score = self._calculate_provenance_score(
                provenance, traditional
            )
            self._provenance_records[provenance.record_id] = provenance
        else:
            report.provenance_score = 0.3  # Partial credit for claim without verification
            report.recommendations.append(
                "Provide blockchain-verified provenance for higher score"
            )
        
        # Calculate overall score
        report.overall_match_score = (
            report.chemical_match_score * 0.50 +
            report.method_match_score * 0.25 +
            report.provenance_score * 0.25
        )
        
        # Determine status
        if report.overall_match_score >= self._OVERALL_VERIFICATION_THRESHOLD:
            report.status = VerificationStatus.VERIFIED
        elif report.overall_match_score >= 0.50:
            report.status = VerificationStatus.PARTIAL_MATCH
            report.recommendations.append(
                f"Score {report.overall_match_score:.2f} below threshold "
                f"{self._OVERALL_VERIFICATION_THRESHOLD}. Address discrepancies."
            )
        else:
            report.status = VerificationStatus.NO_MATCH
            report.discrepancies.append(
                "Formulation does not match claimed traditional preparation"
            )
        
        # Generate verification hash
        report.generate_verification_hash()
        
        # Store report
        self._verification_reports[report.report_id] = report
        
        return report
    
    def _calculate_chemical_match(
        self,
        profile: ChemicalProfile,
        traditional: TraditionalPreparation,
    ) -> float:
        """
        Calculate chemical profile match score.
        Trade Secret: Matching algorithm accounting for natural variation.
        """
        matches = 0
        total_checks = 0
        
        # Check cannabinoids
        for cannabinoid, (min_val, max_val) in traditional.expected_cannabinoid_ranges.items():
            total_checks += 1
            actual = profile.cannabinoids.get(cannabinoid, 0)
            
            if min_val <= actual <= max_val:
                matches += 1.0
            elif actual > 0:
                # Partial credit for presence outside range
                if actual < min_val:
                    matches += max(0, 1 - (min_val - actual) / min_val * 0.5)
                else:
                    matches += max(0, 1 - (actual - max_val) / max_val * 0.5)
        
        # Check terpenes
        for terpene, (min_val, max_val) in traditional.expected_terpene_ranges.items():
            total_checks += 1
            actual = profile.terpenes.get(terpene, 0)
            
            if min_val <= actual <= max_val:
                matches += 1.0
            elif actual > 0:
                matches += 0.5  # Partial credit for presence
        
        if total_checks == 0:
            return 0.5  # No reference data
        
        return matches / total_checks
    
    def _calculate_method_match(
        self,
        formulation_name: str,
        traditional: TraditionalPreparation,
    ) -> float:
        """
        Calculate method match based on preparation description.
        Trade Secret: Method matching algorithm.
        """
        # Simplified matching - real implementation would be more sophisticated
        formulation_lower = formulation_name.lower()
        
        method_keywords = {
            "hand_rubbing": ["charras", "hand-rubbed", "hand rubbed", "fresh"],
            "dry_sift": ["hash", "kief", "sift", "dry"],
            "water_extraction": ["bhang", "water", "extract", "tincture"],
            "sun_dried": ["sun", "outdoor", "natural", "ganja"],
        }
        
        method = traditional.extraction_method
        keywords = method_keywords.get(method, [])
        
        matches = sum(1 for kw in keywords if kw in formulation_lower)
        return min(matches / max(len(keywords), 1), 1.0)
    
    def _calculate_provenance_score(
        self,
        provenance: ProvenanceRecord,
        traditional: TraditionalPreparation,
    ) -> float:
        """
        Calculate provenance verification score.
        Trade Secret: Provenance scoring algorithm.
        """
        score = 0.0
        
        # Community match
        if provenance.source_community.lower() in traditional.community_origin.lower():
            score += 0.3
        
        # Region match
        if provenance.source_region.lower() in traditional.region.lower():
            score += 0.2
        
        # Blockchain verification
        if provenance.status == ProvenanceStatus.BLOCKCHAIN_VERIFIED:
            score += 0.3
        elif provenance.status == ProvenanceStatus.THIRD_PARTY_VERIFIED:
            score += 0.2
        elif provenance.status == ProvenanceStatus.SELF_ATTESTED:
            score += 0.1
        
        # GPS coordinates provided
        if provenance.gps_coordinates:
            score += 0.1
        
        # Harvest timing documented
        if provenance.harvest_timing_notes:
            score += 0.1
        
        return min(score, 1.0)
    
    def get_verification_report(self, report_id: str) -> Optional[VerificationReport]:
        """Retrieve verification report by ID."""
        return self._verification_reports.get(report_id)
    
    def list_traditional_preparations(self) -> List[Dict]:
        """List available traditional preparations for verification."""
        return [
            {
                "key": key,
                "name": prep.preparation_name,
                "community": prep.community_origin,
                "region": prep.region,
                "method": prep.extraction_method,
            }
            for key, prep in self._TRADITIONAL_PREPARATIONS.items()
        ]
    
    def record_provenance(self, provenance: ProvenanceRecord) -> Dict:
        """Record provenance information for future verification."""
        # Generate blockchain hash
        blockchain_hash = provenance.generate_blockchain_hash()
        provenance.blockchain_tx_hash = blockchain_hash
        provenance.verification_date = datetime.utcnow()
        
        self._provenance_records[provenance.record_id] = provenance
        
        return {
            "success": True,
            "record_id": provenance.record_id,
            "blockchain_hash": blockchain_hash,
            "message": "Provenance record stored for verification",
        }
