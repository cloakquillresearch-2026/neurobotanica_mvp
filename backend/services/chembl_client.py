"""
ChEMBL Integration Client - NeuroBotanica Week 5
Client for ChEMBL database API integration.

Provides:
- SMILES-based compound search
- Receptor binding data with full assay context
- Bioactivity data retrieval
- Assay metadata preservation
"""
import requests
from typing import Dict, List, Optional, Any
from dataclasses import dataclass, field
from datetime import datetime
import logging
import time

logger = logging.getLogger(__name__)


@dataclass
class ChEMBLActivity:
    """Single bioactivity measurement from ChEMBL."""
    activity_id: str
    molecule_chembl_id: str
    target_chembl_id: str
    target_name: str
    target_organism: str
    assay_type: str
    assay_description: str
    standard_type: str  # IC50, Ki, EC50, etc.
    standard_value: Optional[float]
    standard_units: str
    standard_relation: str  # =, <, >, etc.
    pchembl_value: Optional[float]  # Standardized -log10 value
    data_validity_comment: Optional[str]
    document_chembl_id: str
    document_year: Optional[int]
    src_id: int
    
    def to_provenance_dict(self) -> Dict:
        """Convert to provenance-rich dictionary format."""
        return {
            "affinity_value": self.standard_value,
            "affinity_unit": self.standard_units,
            "affinity_type": self.standard_type,
            "affinity_relation": self.standard_relation,
            "pchembl_value": self.pchembl_value,
            "assay_type": self.assay_type,
            "assay_description": self.assay_description,
            "target_name": self.target_name,
            "target_organism": self.target_organism,
            "target_chembl_id": self.target_chembl_id,
            "source": f"ChEMBL:{self.molecule_chembl_id}",
            "activity_id": self.activity_id,
            "document_id": self.document_chembl_id,
            "document_year": self.document_year,
            "data_validity": self.data_validity_comment,
            "retrieved_at": datetime.now().isoformat()
        }


@dataclass
class ChEMBLMolecule:
    """Molecule record from ChEMBL database."""
    molecule_chembl_id: str
    pref_name: Optional[str]
    molecule_type: str
    max_phase: int  # Clinical trial phase (0-4)
    therapeutic_flag: bool
    natural_product: bool
    oral: bool
    parenteral: bool
    topical: bool
    molecular_weight: Optional[float]
    alogp: Optional[float]
    psa: Optional[float]
    hba: Optional[int]
    hbd: Optional[int]
    num_ro5_violations: Optional[int]
    canonical_smiles: Optional[str]
    standard_inchi: Optional[str]
    standard_inchi_key: Optional[str]
    
    def to_dict(self) -> Dict:
        return {
            "chembl_id": self.molecule_chembl_id,
            "name": self.pref_name,
            "molecule_type": self.molecule_type,
            "max_phase": self.max_phase,
            "therapeutic_flag": self.therapeutic_flag,
            "natural_product": self.natural_product,
            "routes": {
                "oral": self.oral,
                "parenteral": self.parenteral,
                "topical": self.topical
            },
            "properties": {
                "molecular_weight": self.molecular_weight,
                "alogp": self.alogp,
                "psa": self.psa,
                "hba": self.hba,
                "hbd": self.hbd,
                "ro5_violations": self.num_ro5_violations
            },
            "smiles": self.canonical_smiles,
            "inchi": self.standard_inchi,
            "inchi_key": self.standard_inchi_key
        }


@dataclass
class ChEMBLSearchResult:
    """Result from ChEMBL API search."""
    success: bool
    molecule: Optional[ChEMBLMolecule] = None
    activities: List[ChEMBLActivity] = field(default_factory=list)
    error: Optional[str] = None
    request_time_ms: float = 0.0
    
    def to_dict(self) -> Dict:
        return {
            "success": self.success,
            "molecule": self.molecule.to_dict() if self.molecule else None,
            "activities_count": len(self.activities),
            "error": self.error,
            "request_time_ms": self.request_time_ms
        }


class ChEMBLClient:
    """Client for ChEMBL database REST API.
    
    Provides access to ChEMBL compound and bioactivity data with full
    provenance tracking for scientific rigor.
    
    API Docs: https://www.ebi.ac.uk/chembl/api/data/docs
    """
    
    BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"
    
    # Cannabinoid-relevant targets in ChEMBL
    CANNABINOID_TARGETS = {
        "CB1": "CHEMBL218",   # Cannabinoid CB1 receptor
        "CB2": "CHEMBL253",   # Cannabinoid CB2 receptor
        "GPR55": "CHEMBL2414", # G protein-coupled receptor 55
        "TRPV1": "CHEMBL4105", # Transient receptor potential cation channel subfamily V member 1
        "FAAH": "CHEMBL3421",  # Fatty acid amide hydrolase
        "MAGL": "CHEMBL4840",  # Monoacylglycerol lipase
        "PPARα": "CHEMBL239", # Peroxisome proliferator-activated receptor alpha
        "PPARγ": "CHEMBL235", # Peroxisome proliferator-activated receptor gamma
        "5-HT1A": "CHEMBL214", # 5-hydroxytryptamine receptor 1A
        "TRPA1": "CHEMBL4794", # Transient receptor potential cation channel subfamily A member 1
    }
    
    def __init__(
        self,
        timeout: int = 30,
        max_retries: int = 3,
        retry_delay: float = 1.0
    ):
        """Initialize ChEMBL client.
        
        Args:
            timeout: Request timeout in seconds
            max_retries: Maximum number of retry attempts
            retry_delay: Delay between retries in seconds
        """
        self.session = requests.Session()
        self.session.headers.update({
            "Accept": "application/json",
            "User-Agent": "NeuroBotanica/1.0 (Cannabis Research Platform)"
        })
        self.timeout = timeout
        self.max_retries = max_retries
        self.retry_delay = retry_delay
        logger.info("ChEMBLClient initialized")
    
    def _make_request(
        self,
        endpoint: str,
        params: Optional[Dict] = None
    ) -> Optional[Dict]:
        """Make request to ChEMBL API with retry logic."""
        url = f"{self.BASE_URL}/{endpoint}"
        
        for attempt in range(self.max_retries):
            try:
                response = self.session.get(
                    url,
                    params=params,
                    timeout=self.timeout
                )
                
                if response.status_code == 200:
                    return response.json()
                elif response.status_code == 404:
                    return None
                elif response.status_code == 429:
                    # Rate limited - wait and retry
                    wait_time = self.retry_delay * (attempt + 1) * 2
                    logger.warning(f"Rate limited, waiting {wait_time}s")
                    time.sleep(wait_time)
                else:
                    logger.error(f"ChEMBL API error: {response.status_code}")
                    
            except requests.Timeout:
                logger.warning(f"Request timeout (attempt {attempt + 1})")
                time.sleep(self.retry_delay)
            except requests.RequestException as e:
                logger.error(f"Request error: {e}")
                time.sleep(self.retry_delay)
        
        return None
    
    def search_by_smiles(
        self,
        smiles: str,
        similarity: int = 100
    ) -> Optional[str]:
        """Search ChEMBL by SMILES string.
        
        Args:
            smiles: SMILES string of compound
            similarity: Minimum similarity threshold (40-100)
            
        Returns:
            ChEMBL ID if found, None otherwise
        """
        start_time = time.time()
        
        # First try exact match
        data = self._make_request(
            "molecule",
            params={
                "molecule_structures__canonical_smiles__flexmatch": smiles,
                "limit": 1
            }
        )
        
        if data and data.get("molecules"):
            chembl_id = data["molecules"][0]["molecule_chembl_id"]
            logger.info(f"Found exact SMILES match: {chembl_id}")
            return chembl_id
        
        # Try similarity search if no exact match
        if similarity < 100:
            data = self._make_request(
                f"similarity/{smiles}/{similarity}",
                params={"limit": 1}
            )
            
            if data and data.get("molecules"):
                chembl_id = data["molecules"][0]["molecule_chembl_id"]
                logger.info(f"Found similar compound: {chembl_id}")
                return chembl_id
        
        logger.info(f"No ChEMBL match found for SMILES (took {time.time()-start_time:.2f}s)")
        return None
    
    def search_by_inchi_key(self, inchi_key: str) -> Optional[str]:
        """Search ChEMBL by InChI Key.
        
        Args:
            inchi_key: InChI Key of compound
            
        Returns:
            ChEMBL ID if found, None otherwise
        """
        data = self._make_request(
            "molecule",
            params={
                "molecule_structures__standard_inchi_key": inchi_key,
                "limit": 1
            }
        )
        
        if data and data.get("molecules"):
            return data["molecules"][0]["molecule_chembl_id"]
        return None
    
    def search_by_name(self, name: str) -> Optional[str]:
        """Search ChEMBL by compound name.
        
        Args:
            name: Compound name (e.g., "cannabidiol")
            
        Returns:
            ChEMBL ID if found, None otherwise
        """
        data = self._make_request(
            "molecule",
            params={
                "pref_name__icontains": name,
                "limit": 5
            }
        )
        
        if data and data.get("molecules"):
            # Try to find exact match first
            for mol in data["molecules"]:
                if mol.get("pref_name", "").lower() == name.lower():
                    return mol["molecule_chembl_id"]
            # Return first result if no exact match
            return data["molecules"][0]["molecule_chembl_id"]
        return None
    
    def get_molecule(self, chembl_id: str) -> Optional[ChEMBLMolecule]:
        """Get molecule details from ChEMBL.
        
        Args:
            chembl_id: ChEMBL molecule ID (e.g., "CHEMBL123")
            
        Returns:
            ChEMBLMolecule if found, None otherwise
        """
        data = self._make_request(f"molecule/{chembl_id}")
        
        if not data:
            return None
        
        structs = data.get("molecule_structures", {}) or {}
        props = data.get("molecule_properties", {}) or {}
        
        return ChEMBLMolecule(
            molecule_chembl_id=data["molecule_chembl_id"],
            pref_name=data.get("pref_name"),
            molecule_type=data.get("molecule_type", "Unknown"),
            max_phase=data.get("max_phase", 0) or 0,
            therapeutic_flag=data.get("therapeutic_flag", False),
            natural_product=data.get("natural_product", -1) == 1,
            oral=data.get("oral", False),
            parenteral=data.get("parenteral", False),
            topical=data.get("topical", False),
            molecular_weight=props.get("full_mwt"),
            alogp=props.get("alogp"),
            psa=props.get("psa"),
            hba=props.get("hba"),
            hbd=props.get("hbd"),
            num_ro5_violations=props.get("num_ro5_violations"),
            canonical_smiles=structs.get("canonical_smiles"),
            standard_inchi=structs.get("standard_inchi"),
            standard_inchi_key=structs.get("standard_inchi_key")
        )
    
    def get_receptor_affinities(
        self,
        chembl_id: str,
        target_name: Optional[str] = None,
        assay_type: str = "B",
        limit: int = 500
    ) -> List[ChEMBLActivity]:
        """Get receptor binding data with full assay context.
        
        Args:
            chembl_id: ChEMBL molecule ID
            target_name: Optional target filter (e.g., "CB1", "CB2")
            assay_type: Assay type filter (B=Binding, F=Functional, A=ADMET)
            limit: Maximum number of activities to retrieve
            
        Returns:
            List of ChEMBLActivity records with provenance
        """
        params = {
            "molecule_chembl_id": chembl_id,
            "limit": limit
        }
        
        if assay_type:
            params["assay_type"] = assay_type
        
        # If target specified, filter by target ChEMBL ID
        if target_name and target_name.upper() in self.CANNABINOID_TARGETS:
            params["target_chembl_id"] = self.CANNABINOID_TARGETS[target_name.upper()]
        
        data = self._make_request("activity", params=params)
        
        if not data:
            return []
        
        activities = []
        for act in data.get("activities", []):
            # Filter by target name if specified but not in our known targets
            if target_name and target_name.upper() not in self.CANNABINOID_TARGETS:
                target_pref = act.get("target_pref_name", "")
                if target_name.lower() not in target_pref.lower():
                    continue
            
            try:
                activity = ChEMBLActivity(
                    activity_id=str(act.get("activity_id", "")),
                    molecule_chembl_id=act.get("molecule_chembl_id", chembl_id),
                    target_chembl_id=act.get("target_chembl_id", ""),
                    target_name=act.get("target_pref_name", ""),
                    target_organism=act.get("target_organism", ""),
                    assay_type=act.get("assay_type", ""),
                    assay_description=act.get("assay_description", ""),
                    standard_type=act.get("standard_type", ""),
                    standard_value=self._safe_float(act.get("standard_value")),
                    standard_units=act.get("standard_units", ""),
                    standard_relation=act.get("standard_relation", "="),
                    pchembl_value=self._safe_float(act.get("pchembl_value")),
                    data_validity_comment=act.get("data_validity_comment"),
                    document_chembl_id=act.get("document_chembl_id", ""),
                    document_year=act.get("document_year"),
                    src_id=act.get("src_id", 0)
                )
                activities.append(activity)
            except Exception as e:
                logger.warning(f"Error parsing activity: {e}")
                continue
        
        logger.info(f"Retrieved {len(activities)} activities for {chembl_id}")
        return activities
    
    def get_cannabinoid_target_activities(
        self,
        chembl_id: str
    ) -> Dict[str, List[ChEMBLActivity]]:
        """Get activities for all cannabinoid-relevant targets.
        
        Args:
            chembl_id: ChEMBL molecule ID
            
        Returns:
            Dictionary mapping target names to activity lists
        """
        results = {}
        
        for target_name, target_chembl_id in self.CANNABINOID_TARGETS.items():
            activities = self.get_receptor_affinities(
                chembl_id=chembl_id,
                target_name=target_name,
                limit=100
            )
            if activities:
                results[target_name] = activities
        
        return results
    
    def full_compound_search(
        self,
        smiles: Optional[str] = None,
        name: Optional[str] = None,
        inchi_key: Optional[str] = None,
        include_activities: bool = True
    ) -> ChEMBLSearchResult:
        """Perform full compound search with molecule info and activities.
        
        Args:
            smiles: SMILES string
            name: Compound name
            inchi_key: InChI Key
            include_activities: Whether to fetch bioactivity data
            
        Returns:
            ChEMBLSearchResult with molecule and activity data
        """
        start_time = time.time()
        
        # Search for ChEMBL ID
        chembl_id = None
        
        if inchi_key:
            chembl_id = self.search_by_inchi_key(inchi_key)
        
        if not chembl_id and smiles:
            chembl_id = self.search_by_smiles(smiles)
        
        if not chembl_id and name:
            chembl_id = self.search_by_name(name)
        
        if not chembl_id:
            return ChEMBLSearchResult(
                success=False,
                error="Compound not found in ChEMBL",
                request_time_ms=(time.time() - start_time) * 1000
            )
        
        # Get molecule details
        molecule = self.get_molecule(chembl_id)
        
        # Get activities if requested
        activities = []
        if include_activities:
            activities = self.get_receptor_affinities(chembl_id, limit=200)
        
        return ChEMBLSearchResult(
            success=True,
            molecule=molecule,
            activities=activities,
            request_time_ms=(time.time() - start_time) * 1000
        )
    
    def _safe_float(self, value: Any) -> Optional[float]:
        """Safely convert value to float."""
        if value is None:
            return None
        try:
            return float(value)
        except (ValueError, TypeError):
            return None
    
    def get_known_cannabinoid_ids(self) -> Dict[str, str]:
        """Get ChEMBL IDs for well-known cannabinoids.
        
        Returns cached mapping for common cannabinoids.
        """
        # These are verified ChEMBL IDs for major cannabinoids
        return {
            "THC": "CHEMBL465",            # Δ9-THC
            "CBD": "CHEMBL482897",         # Cannabidiol
            "CBN": "CHEMBL6859",           # Cannabinol
            "CBG": "CHEMBL1214033",        # Cannabigerol
            "CBC": "CHEMBL408306",         # Cannabichromene
            "THCA": "CHEMBL507347",        # THCA
            "CBDA": "CHEMBL1214036",       # CBDA
            "THCV": "CHEMBL407420",        # THCV
            "CBDV": "CHEMBL1234453",       # CBDV
            "Δ8-THC": "CHEMBL560",         # Δ8-THC
            "11-OH-THC": "CHEMBL58612",    # 11-hydroxy-THC
            "Anandamide": "CHEMBL15",      # AEA
            "2-AG": "CHEMBL36036",         # 2-Arachidonoylglycerol
        }
