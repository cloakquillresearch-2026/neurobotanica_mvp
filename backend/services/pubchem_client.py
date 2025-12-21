"""
PubChem Integration Client - NeuroBotanica Week 5
Client for PubChem database API integration.

Provides:
- SMILES/InChI-based compound search
- Compound property retrieval
- Bioassay data with activity outcomes
- Patent information
"""
import requests
from typing import Dict, List, Optional, Any
from dataclasses import dataclass, field
from datetime import datetime
import logging
import time
import urllib.parse

logger = logging.getLogger(__name__)


@dataclass
class PubChemProperty:
    """Compound properties from PubChem."""
    cid: str
    molecular_formula: Optional[str]
    molecular_weight: Optional[float]
    canonical_smiles: Optional[str]
    isomeric_smiles: Optional[str]
    inchi: Optional[str]
    inchi_key: Optional[str]
    iupac_name: Optional[str]
    xlogp: Optional[float]
    exact_mass: Optional[float]
    monoisotopic_mass: Optional[float]
    tpsa: Optional[float]
    complexity: Optional[float]
    charge: Optional[int]
    h_bond_donor_count: Optional[int]
    h_bond_acceptor_count: Optional[int]
    rotatable_bond_count: Optional[int]
    heavy_atom_count: Optional[int]
    atom_stereo_count: Optional[int]
    defined_atom_stereo_count: Optional[int]
    undefined_atom_stereo_count: Optional[int]
    bond_stereo_count: Optional[int]
    covalent_unit_count: Optional[int]
    
    def to_dict(self) -> Dict:
        return {
            "cid": self.cid,
            "molecular_formula": self.molecular_formula,
            "molecular_weight": self.molecular_weight,
            "canonical_smiles": self.canonical_smiles,
            "isomeric_smiles": self.isomeric_smiles,
            "inchi": self.inchi,
            "inchi_key": self.inchi_key,
            "iupac_name": self.iupac_name,
            "xlogp": self.xlogp,
            "exact_mass": self.exact_mass,
            "tpsa": self.tpsa,
            "complexity": self.complexity,
            "h_bond_donors": self.h_bond_donor_count,
            "h_bond_acceptors": self.h_bond_acceptor_count,
            "rotatable_bonds": self.rotatable_bond_count,
            "heavy_atoms": self.heavy_atom_count
        }


@dataclass
class PubChemBioassay:
    """Bioassay activity from PubChem."""
    aid: int
    assay_name: str
    assay_type: str
    activity_outcome: str  # Active, Inactive, Inconclusive, Unspecified
    activity_value: Optional[float]
    activity_unit: Optional[str]
    target_name: Optional[str]
    target_gi: Optional[int]
    pubmed_id: Optional[int]
    
    def to_provenance_dict(self) -> Dict:
        """Convert to provenance-rich dictionary format."""
        return {
            "assay_id": f"PubChem:AID{self.aid}",
            "assay_name": self.assay_name,
            "assay_type": self.assay_type,
            "activity_outcome": self.activity_outcome,
            "activity_value": self.activity_value,
            "activity_unit": self.activity_unit,
            "target_name": self.target_name,
            "target_gi": self.target_gi,
            "pubmed_id": self.pubmed_id,
            "source": "PubChem BioAssay",
            "retrieved_at": datetime.now().isoformat()
        }


@dataclass
class PubChemSynonym:
    """Synonym information from PubChem."""
    cid: str
    synonyms: List[str]
    depositor_names: List[str]
    mesh_headings: List[str]
    
    def get_common_names(self, max_names: int = 10) -> List[str]:
        """Get most common/useful names."""
        all_names = self.synonyms[:max_names]
        return [n for n in all_names if len(n) < 50 and not n.startswith("SCHEMBL")]


@dataclass
class PubChemSearchResult:
    """Result from PubChem API search."""
    success: bool
    cid: Optional[str] = None
    properties: Optional[PubChemProperty] = None
    synonyms: Optional[PubChemSynonym] = None
    bioassays: List[PubChemBioassay] = field(default_factory=list)
    error: Optional[str] = None
    request_time_ms: float = 0.0
    
    def to_dict(self) -> Dict:
        return {
            "success": self.success,
            "cid": self.cid,
            "properties": self.properties.to_dict() if self.properties else None,
            "synonyms": self.synonyms.get_common_names() if self.synonyms else None,
            "bioassays_count": len(self.bioassays),
            "error": self.error,
            "request_time_ms": self.request_time_ms
        }


class PubChemClient:
    """Client for PubChem PUG REST API.
    
    Provides access to PubChem compound data, properties, and bioassay results.
    
    API Docs: https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest
    """
    
    BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    
    # Common property list for compound queries
    DEFAULT_PROPERTIES = [
        "MolecularFormula",
        "MolecularWeight",
        "CanonicalSMILES",
        "IsomericSMILES",
        "InChI",
        "InChIKey",
        "IUPACName",
        "XLogP",
        "ExactMass",
        "MonoisotopicMass",
        "TPSA",
        "Complexity",
        "Charge",
        "HBondDonorCount",
        "HBondAcceptorCount",
        "RotatableBondCount",
        "HeavyAtomCount",
        "AtomStereoCount",
        "DefinedAtomStereoCount",
        "UndefinedAtomStereoCount",
        "BondStereoCount",
        "CovalentUnitCount"
    ]
    
    def __init__(
        self,
        timeout: int = 30,
        max_retries: int = 3,
        retry_delay: float = 0.5
    ):
        """Initialize PubChem client.
        
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
        logger.info("PubChemClient initialized")
    
    def _make_request(
        self,
        url: str,
        params: Optional[Dict] = None
    ) -> Optional[Dict]:
        """Make request to PubChem API with retry logic."""
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
                elif response.status_code == 503:
                    # Service busy - wait and retry
                    wait_time = self.retry_delay * (attempt + 1) * 2
                    logger.warning(f"PubChem busy, waiting {wait_time}s")
                    time.sleep(wait_time)
                else:
                    logger.error(f"PubChem API error: {response.status_code}")
                    
            except requests.Timeout:
                logger.warning(f"Request timeout (attempt {attempt + 1})")
                time.sleep(self.retry_delay)
            except requests.RequestException as e:
                logger.error(f"Request error: {e}")
                time.sleep(self.retry_delay)
        
        return None
    
    def search_by_smiles(self, smiles: str) -> Optional[str]:
        """Search PubChem by SMILES string.
        
        Args:
            smiles: SMILES string of compound
            
        Returns:
            PubChem CID if found, None otherwise
        """
        # URL encode the SMILES
        encoded_smiles = urllib.parse.quote(smiles, safe='')
        url = f"{self.BASE_URL}/compound/smiles/{encoded_smiles}/cids/JSON"
        
        data = self._make_request(url)
        
        if data and data.get("IdentifierList", {}).get("CID"):
            cid = str(data["IdentifierList"]["CID"][0])
            logger.info(f"Found PubChem CID: {cid}")
            return cid
        
        return None
    
    def search_by_inchi_key(self, inchi_key: str) -> Optional[str]:
        """Search PubChem by InChI Key.
        
        Args:
            inchi_key: InChI Key of compound
            
        Returns:
            PubChem CID if found, None otherwise
        """
        url = f"{self.BASE_URL}/compound/inchikey/{inchi_key}/cids/JSON"
        
        data = self._make_request(url)
        
        if data and data.get("IdentifierList", {}).get("CID"):
            return str(data["IdentifierList"]["CID"][0])
        
        return None
    
    def search_by_name(self, name: str) -> Optional[str]:
        """Search PubChem by compound name.
        
        Args:
            name: Compound name
            
        Returns:
            PubChem CID if found, None otherwise
        """
        encoded_name = urllib.parse.quote(name)
        url = f"{self.BASE_URL}/compound/name/{encoded_name}/cids/JSON"
        
        data = self._make_request(url)
        
        if data and data.get("IdentifierList", {}).get("CID"):
            return str(data["IdentifierList"]["CID"][0])
        
        return None
    
    def get_compound_properties(
        self,
        cid: str,
        properties: Optional[List[str]] = None
    ) -> Optional[PubChemProperty]:
        """Get compound properties from PubChem.
        
        Args:
            cid: PubChem Compound ID
            properties: List of property names (uses defaults if None)
            
        Returns:
            PubChemProperty if found, None otherwise
        """
        props = properties or self.DEFAULT_PROPERTIES
        props_str = ",".join(props)
        
        url = f"{self.BASE_URL}/compound/cid/{cid}/property/{props_str}/JSON"
        
        data = self._make_request(url)
        
        if not data or not data.get("PropertyTable", {}).get("Properties"):
            return None
        
        p = data["PropertyTable"]["Properties"][0]
        
        return PubChemProperty(
            cid=cid,
            molecular_formula=p.get("MolecularFormula"),
            molecular_weight=self._safe_float(p.get("MolecularWeight")),
            canonical_smiles=p.get("CanonicalSMILES"),
            isomeric_smiles=p.get("IsomericSMILES"),
            inchi=p.get("InChI"),
            inchi_key=p.get("InChIKey"),
            iupac_name=p.get("IUPACName"),
            xlogp=self._safe_float(p.get("XLogP")),
            exact_mass=self._safe_float(p.get("ExactMass")),
            monoisotopic_mass=self._safe_float(p.get("MonoisotopicMass")),
            tpsa=self._safe_float(p.get("TPSA")),
            complexity=self._safe_float(p.get("Complexity")),
            charge=p.get("Charge"),
            h_bond_donor_count=p.get("HBondDonorCount"),
            h_bond_acceptor_count=p.get("HBondAcceptorCount"),
            rotatable_bond_count=p.get("RotatableBondCount"),
            heavy_atom_count=p.get("HeavyAtomCount"),
            atom_stereo_count=p.get("AtomStereoCount"),
            defined_atom_stereo_count=p.get("DefinedAtomStereoCount"),
            undefined_atom_stereo_count=p.get("UndefinedAtomStereoCount"),
            bond_stereo_count=p.get("BondStereoCount"),
            covalent_unit_count=p.get("CovalentUnitCount")
        )
    
    def get_synonyms(self, cid: str) -> Optional[PubChemSynonym]:
        """Get compound synonyms from PubChem.
        
        Args:
            cid: PubChem Compound ID
            
        Returns:
            PubChemSynonym if found, None otherwise
        """
        url = f"{self.BASE_URL}/compound/cid/{cid}/synonyms/JSON"
        
        data = self._make_request(url)
        
        if not data or not data.get("InformationList", {}).get("Information"):
            return None
        
        info = data["InformationList"]["Information"][0]
        
        return PubChemSynonym(
            cid=cid,
            synonyms=info.get("Synonym", [])[:100],
            depositor_names=[],
            mesh_headings=[]
        )
    
    def get_bioassay_summary(
        self,
        cid: str,
        activity_filter: Optional[str] = None
    ) -> List[PubChemBioassay]:
        """Get bioassay activity summary for compound.
        
        Args:
            cid: PubChem Compound ID
            activity_filter: Filter by outcome ("active", "inactive", "all")
            
        Returns:
            List of PubChemBioassay records
        """
        url = f"{self.BASE_URL}/compound/cid/{cid}/assaysummary/JSON"
        
        data = self._make_request(url)
        
        if not data or not data.get("Table"):
            return []
        
        table = data["Table"]
        columns = table.get("Columns", {}).get("Column", [])
        rows = table.get("Row", [])
        
        # Map column indices
        col_map = {col: idx for idx, col in enumerate(columns)}
        
        bioassays = []
        for row in rows[:200]:  # Limit to 200 assays
            cells = row.get("Cell", [])
            
            try:
                outcome = cells[col_map.get("Activity Outcome", 0)] if "Activity Outcome" in col_map else "Unspecified"
                
                # Apply filter if specified
                if activity_filter:
                    if activity_filter.lower() == "active" and outcome != "Active":
                        continue
                    elif activity_filter.lower() == "inactive" and outcome != "Inactive":
                        continue
                
                bioassay = PubChemBioassay(
                    aid=int(cells[col_map.get("AID", 0)]) if "AID" in col_map else 0,
                    assay_name=cells[col_map.get("Assay Name", 1)] if "Assay Name" in col_map else "",
                    assay_type=cells[col_map.get("Assay Type", 2)] if "Assay Type" in col_map else "",
                    activity_outcome=outcome,
                    activity_value=self._safe_float(cells[col_map.get("Activity Value", 3)]) if "Activity Value" in col_map else None,
                    activity_unit=cells[col_map.get("Activity Unit", 4)] if "Activity Unit" in col_map else None,
                    target_name=cells[col_map.get("Target Name", 5)] if "Target Name" in col_map else None,
                    target_gi=self._safe_int(cells[col_map.get("Target GI", 6)]) if "Target GI" in col_map else None,
                    pubmed_id=self._safe_int(cells[col_map.get("PubMed ID", 7)]) if "PubMed ID" in col_map else None
                )
                bioassays.append(bioassay)
                
            except (IndexError, ValueError) as e:
                logger.debug(f"Error parsing bioassay row: {e}")
                continue
        
        logger.info(f"Retrieved {len(bioassays)} bioassays for CID {cid}")
        return bioassays
    
    def get_cannabinoid_relevant_assays(
        self,
        cid: str
    ) -> List[PubChemBioassay]:
        """Get bioassays relevant to cannabinoid research.
        
        Filters for assays targeting cannabinoid receptors, FAAH, MAGL, etc.
        
        Args:
            cid: PubChem Compound ID
            
        Returns:
            Filtered list of relevant bioassays
        """
        all_assays = self.get_bioassay_summary(cid)
        
        # Keywords for cannabinoid-relevant targets
        cannabinoid_keywords = [
            "cannabinoid", "cb1", "cb2", "cnr1", "cnr2",
            "faah", "magl", "trpv1", "gpr55",
            "endocannabinoid", "anandamide",
            "ppar", "5-ht", "serotonin"
        ]
        
        relevant = []
        for assay in all_assays:
            target = (assay.target_name or "").lower()
            name = (assay.assay_name or "").lower()
            
            if any(kw in target or kw in name for kw in cannabinoid_keywords):
                relevant.append(assay)
        
        return relevant
    
    def full_compound_search(
        self,
        smiles: Optional[str] = None,
        name: Optional[str] = None,
        inchi_key: Optional[str] = None,
        include_bioassays: bool = True
    ) -> PubChemSearchResult:
        """Perform full compound search with properties and bioassays.
        
        Args:
            smiles: SMILES string
            name: Compound name
            inchi_key: InChI Key
            include_bioassays: Whether to fetch bioassay data
            
        Returns:
            PubChemSearchResult with all data
        """
        start_time = time.time()
        
        # Search for CID
        cid = None
        
        if inchi_key:
            cid = self.search_by_inchi_key(inchi_key)
        
        if not cid and smiles:
            cid = self.search_by_smiles(smiles)
        
        if not cid and name:
            cid = self.search_by_name(name)
        
        if not cid:
            return PubChemSearchResult(
                success=False,
                error="Compound not found in PubChem",
                request_time_ms=(time.time() - start_time) * 1000
            )
        
        # Get properties
        properties = self.get_compound_properties(cid)
        
        # Get synonyms
        synonyms = self.get_synonyms(cid)
        
        # Get bioassays if requested
        bioassays = []
        if include_bioassays:
            bioassays = self.get_cannabinoid_relevant_assays(cid)
        
        return PubChemSearchResult(
            success=True,
            cid=cid,
            properties=properties,
            synonyms=synonyms,
            bioassays=bioassays,
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
    
    def _safe_int(self, value: Any) -> Optional[int]:
        """Safely convert value to int."""
        if value is None:
            return None
        try:
            return int(value)
        except (ValueError, TypeError):
            return None
    
    def get_known_cannabinoid_cids(self) -> Dict[str, str]:
        """Get PubChem CIDs for well-known cannabinoids.
        
        Returns cached mapping for common cannabinoids.
        """
        return {
            "THC": "16078",           # Δ9-THC
            "CBD": "644019",          # Cannabidiol
            "CBN": "2543",            # Cannabinol
            "CBG": "5315659",         # Cannabigerol
            "CBC": "30219",           # Cannabichromene
            "THCA": "98523",          # THCA
            "CBDA": "160570",         # CBDA
            "THCV": "100027",         # THCV
            "CBDV": "622251",         # CBDV
            "Δ8-THC": "2977",         # Δ8-THC
            "11-OH-THC": "37482",     # 11-hydroxy-THC
            "Anandamide": "5281969",  # AEA
            "2-AG": "5282280",        # 2-Arachidonoylglycerol
            "CBL": "36494",           # Cannabicyclol
            "CBGV": "67452014",       # CBGV
        }
