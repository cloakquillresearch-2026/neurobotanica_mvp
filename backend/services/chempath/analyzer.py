"""
ChemPath Analyzer - Chemical Characterization Engine

MVP Implementation using RDKit for molecular analysis.
Trade secret: QC rules, scoring weights, and alert thresholds are proprietary.
"""

from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
from typing import Optional, List, Dict, Any, Tuple
import uuid
import re
import math


class QCSeverity(Enum):
    """Severity levels for QC flags."""
    ERROR = "error"
    WARNING = "warning"
    INFO = "info"


@dataclass
class QCFlag:
    """Quality control flag for COA or structure validation."""
    code: str
    severity: QCSeverity
    message: str
    field: Optional[str] = None
    
    def to_dict(self) -> Dict:
        return {
            "code": self.code,
            "severity": self.severity.value,
            "message": self.message,
            "field": self.field
        }


@dataclass
class CompoundInput:
    """Input model for compound structure."""
    name: str
    smiles: Optional[str] = None
    inchi: Optional[str] = None
    inchi_key: Optional[str] = None
    source: str = "customer_private"
    
    def to_dict(self) -> Dict:
        return {
            "name": self.name,
            "smiles": self.smiles,
            "inchi": self.inchi,
            "inchi_key": self.inchi_key,
            "source": self.source
        }


@dataclass
class COAInput:
    """Certificate of Analysis input."""
    lab_name: str
    sample_id: str
    units: str  # "%" or "mg/g"
    cannabinoids: List[Dict]  # [{name, value}]
    terpenes: Optional[List[Dict]] = None
    total_thc: Optional[float] = None
    total_cbd: Optional[float] = None
    moisture: Optional[float] = None
    solvents: Optional[List[Dict]] = None
    microbials: Optional[List[Dict]] = None
    heavy_metals: Optional[List[Dict]] = None
    pesticides: Optional[List[Dict]] = None
    lod_loq: Optional[Dict] = None
    
    def to_dict(self) -> Dict:
        return {
            "lab_name": self.lab_name,
            "sample_id": self.sample_id,
            "units": self.units,
            "cannabinoids": self.cannabinoids,
            "terpenes": self.terpenes,
            "total_thc": self.total_thc,
            "total_cbd": self.total_cbd,
            "moisture": self.moisture,
            "solvents": self.solvents,
            "microbials": self.microbials,
            "heavy_metals": self.heavy_metals,
            "pesticides": self.pesticides,
            "lod_loq": self.lod_loq
        }


@dataclass 
class ChemPathRequest:
    """Request model for ChemPath analysis."""
    compound_input: CompoundInput
    coa_input: Optional[COAInput] = None
    compute_3d: bool = False
    
    def to_dict(self) -> Dict:
        return {
            "compound_input": self.compound_input.to_dict(),
            "coa_input": self.coa_input.to_dict() if self.coa_input else None,
            "compute_3d": self.compute_3d
        }


@dataclass
class NormalizedStructure:
    """Normalized molecular structure output."""
    canonical_smiles: str
    inchi: Optional[str] = None
    inchi_key: Optional[str] = None
    molecular_formula: Optional[str] = None
    normalization_notes: List[str] = field(default_factory=list)
    is_valid: bool = True
    
    def to_dict(self) -> Dict:
        return {
            "canonical_smiles": self.canonical_smiles,
            "inchi": self.inchi,
            "inchi_key": self.inchi_key,
            "molecular_formula": self.molecular_formula,
            "normalization_notes": self.normalization_notes,
            "is_valid": self.is_valid
        }


@dataclass
class DescriptorSet:
    """Computed molecular descriptors."""
    descriptors_2d: Dict[str, float]
    descriptors_3d: Optional[Dict[str, float]] = None
    conformer_meta: Optional[Dict] = None
    
    def to_dict(self) -> Dict:
        return {
            "descriptors_2d": self.descriptors_2d,
            "descriptors_3d": self.descriptors_3d,
            "conformer_meta": self.conformer_meta
        }


@dataclass
class ChemPathResponse:
    """Response model for ChemPath analysis."""
    chempath_job_id: str
    compound_name: str
    normalized_structure: NormalizedStructure
    computed_descriptors: DescriptorSet
    coa_qc_flags: List[QCFlag]
    structure_qc_flags: List[QCFlag]
    data_completeness_score: float
    status: str = "succeeded"
    error: Optional[str] = None
    analyzed_at: str = field(default_factory=lambda: datetime.utcnow().isoformat())
    
    def to_dict(self) -> Dict:
        return {
            "chempath_job_id": self.chempath_job_id,
            "compound_name": self.compound_name,
            "normalized_structure": self.normalized_structure.to_dict(),
            "computed_descriptors": self.computed_descriptors.to_dict(),
            "coa_qc_flags": [f.to_dict() for f in self.coa_qc_flags],
            "structure_qc_flags": [f.to_dict() for f in self.structure_qc_flags],
            "data_completeness_score": self.data_completeness_score,
            "status": self.status,
            "error": self.error,
            "analyzed_at": self.analyzed_at
        }


class ChemPathAnalyzer:
    """Chemical Characterization Engine.
    
    MVP Implementation Notes:
    - Uses RDKit for molecular processing (if available)
    - Falls back to rule-based estimation when RDKit unavailable
    - QC rules and thresholds are trade secrets
    
    Trade Secret Elements (not exposed in API):
    - QC scoring weights
    - Alert threshold values  
    - Decision tree logic
    - Normalization heuristics
    """
    
    # =========================================================================
    # TRADE SECRET: QC Codes and Thresholds
    # These codes are exposed but the exact trigger conditions are proprietary
    # =========================================================================
    
    QC_CODES = {
        # Structure validation
        "STRUCTURE_INVALID": "SMILES/InChI parse failure",
        "STRUCTURE_SALT_NORMALIZED": "Salt form detected and normalized",
        "STRUCTURE_STEREOCHEMISTRY_UNDEFINED": "Undefined stereochemistry detected",
        "STRUCTURE_UNUSUAL_ELEMENTS": "Unusual elements detected in structure",
        "STRUCTURE_TOO_LARGE": "Molecular weight exceeds typical cannabinoid range",
        "STRUCTURE_FRAGMENT_DETECTED": "Multiple fragments detected",
        
        # COA validation
        "COA_TOTAL_EXCEEDS_100": "Sum of cannabinoids + terpenes > 100%",
        "COA_UNIT_MISMATCH": "Mixed units in same COA",
        "COA_MISSING_LOD": "LOD/LOQ not provided",
        "COA_NEGATIVE_VALUE": "Negative concentration value detected",
        "COA_MISSING_CANNABINOIDS": "No cannabinoid data provided",
        "COA_IMPLAUSIBLE_THC": "THC value outside expected range",
        "COA_IMPLAUSIBLE_CBD": "CBD value outside expected range",
        "COA_MICROBIAL_FAIL": "Microbial contamination detected",
        "COA_HEAVY_METAL_FAIL": "Heavy metal contamination detected",
        "COA_PESTICIDE_FAIL": "Pesticide contamination detected",
        "COA_SOLVENT_FAIL": "Residual solvent above limits"
    }
    
    # TRADE SECRET: Threshold values (not exposed)
    _THRESHOLDS = {
        "max_mw": 1500,  # Maximum molecular weight for cannabinoids
        "max_total_percent": 100.5,  # Allow small rounding errors
        "max_thc_percent": 40,  # Max plausible THC
        "max_cbd_percent": 50,  # Max plausible CBD
        "min_value": 0,
    }
    
    # Known cannabinoid SMILES for validation
    KNOWN_CANNABINOIDS = {
        "CBD": "CCCCCC1=CC(=C(C(=C1)O)C2C=C(CCC2C(=C)C)C)O",
        "THC": "CCCCCC1=CC(=C2C3CC(=CCC3C(OC2=C1)(C)C)C)O",
        "CBG": "CCCCCC1=CC(=C(C(=C1)O)CC=C(C)CCC=C(C)C)O",
        "CBC": "CCCCCC1=CC(=C2C=CC(OC2=C1)(C)CCC=C(C)C)O",
        "CBN": "CCCCCC1=CC2=C(C=CC(=C2O)C)C(=C1)O"
    }
    
    def __init__(self):
        """Initialize ChemPath analyzer."""
        self._rdkit_available = self._check_rdkit()
    
    def _check_rdkit(self) -> bool:
        """Check if RDKit is available."""
        try:
            from rdkit import Chem
            return True
        except ImportError:
            return False
    
    def analyze(self, request: ChemPathRequest) -> ChemPathResponse:
        """Run full ChemPath analysis pipeline.
        
        Pipeline steps:
        1. Parse and normalize structure
        2. Compute 2D descriptors
        3. Optionally compute 3D conformers and descriptors
        4. Validate COA if provided
        5. Calculate data completeness score
        """
        job_id = str(uuid.uuid4())
        structure_flags: List[QCFlag] = []
        coa_flags: List[QCFlag] = []
        
        # Step 1: Parse and normalize structure
        normalized = self._normalize_structure(
            request.compound_input,
            structure_flags
        )
        
        # Step 2: Compute descriptors
        descriptors = self._compute_descriptors(
            normalized,
            compute_3d=request.compute_3d
        )
        
        # Step 3: Validate COA if provided
        if request.coa_input:
            coa_flags = self._validate_coa(request.coa_input)
        
        # Step 4: Calculate completeness score
        completeness = self._calculate_completeness(
            request.compound_input,
            request.coa_input,
            normalized,
            descriptors
        )
        
        return ChemPathResponse(
            chempath_job_id=job_id,
            compound_name=request.compound_input.name,
            normalized_structure=normalized,
            computed_descriptors=descriptors,
            coa_qc_flags=coa_flags,
            structure_qc_flags=structure_flags,
            data_completeness_score=completeness,
            status="succeeded" if normalized.is_valid else "failed"
        )
    
    def _normalize_structure(
        self,
        compound: CompoundInput,
        flags: List[QCFlag]
    ) -> NormalizedStructure:
        """Normalize molecular structure.
        
        TRADE SECRET: Specific normalization rules and salt handling logic.
        """
        notes: List[str] = []
        canonical_smiles = ""
        inchi = None
        inchi_key = None
        mol_formula = None
        is_valid = True
        
        # Try to get SMILES from input or InChI
        input_smiles = compound.smiles
        
        if self._rdkit_available:
            canonical_smiles, inchi, inchi_key, mol_formula, is_valid = \
                self._normalize_with_rdkit(compound, flags, notes)
        else:
            # Fallback: basic validation without RDKit
            canonical_smiles, is_valid = self._normalize_fallback(
                compound, flags, notes
            )
        
        return NormalizedStructure(
            canonical_smiles=canonical_smiles,
            inchi=inchi or compound.inchi,
            inchi_key=inchi_key or compound.inchi_key,
            molecular_formula=mol_formula,
            normalization_notes=notes,
            is_valid=is_valid
        )
    
    def _normalize_with_rdkit(
        self,
        compound: CompoundInput,
        flags: List[QCFlag],
        notes: List[str]
    ) -> Tuple[str, Optional[str], Optional[str], Optional[str], bool]:
        """Normalize structure using RDKit."""
        from rdkit import Chem
        from rdkit.Chem import Descriptors, rdMolDescriptors, inchi as rdInchi
        
        mol = None
        canonical_smiles = ""
        inchi_str = None
        inchi_key = None
        mol_formula = None
        
        # Try SMILES first
        if compound.smiles:
            try:
                mol = Chem.MolFromSmiles(compound.smiles)
                if mol:
                    notes.append("Parsed from SMILES")
            except:
                pass
        
        # Try InChI if SMILES failed
        if mol is None and compound.inchi:
            try:
                mol = rdInchi.MolFromInchi(compound.inchi)
                if mol:
                    notes.append("Parsed from InChI")
            except:
                pass
        
        # Check for parse failure
        if mol is None:
            flags.append(QCFlag(
                code="STRUCTURE_INVALID",
                severity=QCSeverity.ERROR,
                message=self.QC_CODES["STRUCTURE_INVALID"]
            ))
            return compound.smiles or "", None, None, None, False
        
        # Check for fragments (salts)
        frags = Chem.GetMolFrags(mol, asMols=True)
        if len(frags) > 1:
            flags.append(QCFlag(
                code="STRUCTURE_SALT_NORMALIZED",
                severity=QCSeverity.INFO,
                message=self.QC_CODES["STRUCTURE_SALT_NORMALIZED"]
            ))
            notes.append(f"Salt form with {len(frags)} fragments - using largest")
            # Keep largest fragment by atom count
            frags_by_size = [(f.GetNumAtoms(), f) for f in frags]
            frags_by_size.sort(reverse=True)
            mol = frags_by_size[0][1]
        
        # Check molecular weight
        mw = Descriptors.MolWt(mol)
        if mw > self._THRESHOLDS["max_mw"]:
            flags.append(QCFlag(
                code="STRUCTURE_TOO_LARGE",
                severity=QCSeverity.WARNING,
                message=f"{self.QC_CODES['STRUCTURE_TOO_LARGE']} (MW: {mw:.1f})"
            ))
        
        # Check stereochemistry
        stereo_info = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        undefined = [s for s in stereo_info if s[1] == '?']
        if undefined:
            flags.append(QCFlag(
                code="STRUCTURE_STEREOCHEMISTRY_UNDEFINED",
                severity=QCSeverity.WARNING,
                message=f"{self.QC_CODES['STRUCTURE_STEREOCHEMISTRY_UNDEFINED']} ({len(undefined)} centers)"
            ))
            notes.append(f"Undefined stereochemistry at {len(undefined)} centers")
        
        # Generate canonical SMILES
        canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
        
        # Generate InChI
        try:
            inchi_str = rdInchi.MolToInchi(mol)
            inchi_key = rdInchi.MolToInchiKey(mol)
        except:
            notes.append("InChI generation failed")
        
        # Get molecular formula
        mol_formula = rdMolDescriptors.CalcMolFormula(mol)
        
        return canonical_smiles, inchi_str, inchi_key, mol_formula, True
    
    def _normalize_fallback(
        self,
        compound: CompoundInput,
        flags: List[QCFlag],
        notes: List[str]
    ) -> Tuple[str, bool]:
        """Basic structure validation without RDKit."""
        notes.append("RDKit not available - using basic validation")
        
        smiles = compound.smiles or ""
        
        # Basic SMILES validation
        if smiles:
            # Check for obviously invalid SMILES
            if not re.match(r'^[A-Za-z0-9@+\-\[\]\(\)\\/#=%\.\*]+$', smiles):
                flags.append(QCFlag(
                    code="STRUCTURE_INVALID",
                    severity=QCSeverity.ERROR,
                    message="Invalid characters in SMILES"
                ))
                return smiles, False
            
            # Check bracket balance
            if smiles.count('[') != smiles.count(']'):
                flags.append(QCFlag(
                    code="STRUCTURE_INVALID",
                    severity=QCSeverity.ERROR,
                    message="Unbalanced brackets in SMILES"
                ))
                return smiles, False
            
            # Check parenthesis balance
            if smiles.count('(') != smiles.count(')'):
                flags.append(QCFlag(
                    code="STRUCTURE_INVALID",
                    severity=QCSeverity.ERROR,
                    message="Unbalanced parentheses in SMILES"
                ))
                return smiles, False
            
            # Check for fragments (dots indicate multiple molecules)
            if '.' in smiles:
                flags.append(QCFlag(
                    code="STRUCTURE_FRAGMENT_DETECTED",
                    severity=QCSeverity.WARNING,
                    message=self.QC_CODES["STRUCTURE_FRAGMENT_DETECTED"]
                ))
                notes.append("Multiple fragments detected - may be salt form")
            
            return smiles, True
        
        # No SMILES provided
        flags.append(QCFlag(
            code="STRUCTURE_INVALID",
            severity=QCSeverity.ERROR,
            message="No valid structure input provided"
        ))
        return "", False
    
    def _compute_descriptors(
        self,
        normalized: NormalizedStructure,
        compute_3d: bool = False
    ) -> DescriptorSet:
        """Compute molecular descriptors.
        
        TRADE SECRET: Descriptor selection and weighting for downstream ML.
        """
        descriptors_2d: Dict[str, float] = {}
        descriptors_3d: Optional[Dict[str, float]] = None
        conformer_meta: Optional[Dict] = None
        
        if not normalized.is_valid:
            return DescriptorSet(
                descriptors_2d={},
                descriptors_3d=None,
                conformer_meta=None
            )
        
        if self._rdkit_available:
            descriptors_2d = self._compute_2d_rdkit(normalized.canonical_smiles)
            
            if compute_3d:
                descriptors_3d, conformer_meta = self._compute_3d_rdkit(
                    normalized.canonical_smiles
                )
        else:
            descriptors_2d = self._estimate_descriptors(normalized.canonical_smiles)
        
        return DescriptorSet(
            descriptors_2d=descriptors_2d,
            descriptors_3d=descriptors_3d,
            conformer_meta=conformer_meta
        )
    
    def _compute_2d_rdkit(self, smiles: str) -> Dict[str, float]:
        """Compute 2D descriptors using RDKit."""
        from rdkit import Chem
        from rdkit.Chem import Descriptors, rdMolDescriptors, Lipinski
        
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return {}
        
        return {
            # Basic properties
            "molecular_weight": round(Descriptors.MolWt(mol), 2),
            "exact_mass": round(Descriptors.ExactMolWt(mol), 4),
            "heavy_atom_count": rdMolDescriptors.CalcNumHeavyAtoms(mol),
            
            # Lipinski descriptors
            "logp": round(Descriptors.MolLogP(mol), 2),
            "tpsa": round(Descriptors.TPSA(mol), 2),
            "hbd": Descriptors.NumHDonors(mol),
            "hba": Descriptors.NumHAcceptors(mol),
            "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
            
            # Ring information
            "num_rings": rdMolDescriptors.CalcNumRings(mol),
            "num_aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
            "num_aliphatic_rings": rdMolDescriptors.CalcNumAliphaticRings(mol),
            
            # Complexity
            "fraction_csp3": round(rdMolDescriptors.CalcFractionCSP3(mol), 3),
            "num_stereocenters": len(Chem.FindMolChiralCenters(mol, includeUnassigned=True)),
            
            # Drug-likeness
            "lipinski_violations": self._count_lipinski_violations(mol),
            "qed": round(Descriptors.qed(mol), 3),
            
            # Surface area
            "labute_asa": round(rdMolDescriptors.CalcLabuteASA(mol), 2)
        }
    
    def _count_lipinski_violations(self, mol) -> int:
        """Count Lipinski Rule of 5 violations."""
        from rdkit.Chem import Descriptors
        
        violations = 0
        if Descriptors.MolWt(mol) > 500:
            violations += 1
        if Descriptors.MolLogP(mol) > 5:
            violations += 1
        if Descriptors.NumHDonors(mol) > 5:
            violations += 1
        if Descriptors.NumHAcceptors(mol) > 10:
            violations += 1
        return violations
    
    def _compute_3d_rdkit(
        self,
        smiles: str
    ) -> Tuple[Optional[Dict[str, float]], Optional[Dict]]:
        """Compute 3D descriptors using RDKit."""
        from rdkit import Chem
        from rdkit.Chem import AllChem, Descriptors3D
        
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return None, None
        
        # Add hydrogens for 3D
        mol = Chem.AddHs(mol)
        
        # Generate conformer
        try:
            result = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
            if result == -1:
                return None, {"error": "Conformer generation failed"}
            
            # Optimize
            AllChem.MMFFOptimizeMolecule(mol)
            
            # Compute 3D descriptors
            descriptors_3d = {
                "radius_of_gyration": round(Descriptors3D.RadiusOfGyration(mol), 3),
                "inertial_shape_factor": round(Descriptors3D.InertialShapeFactor(mol), 3),
                "eccentricity": round(Descriptors3D.Eccentricity(mol), 3),
                "asphericity": round(Descriptors3D.Asphericity(mol), 3),
                "spherocity_index": round(Descriptors3D.SpherocityIndex(mol), 3),
                "pmi1": round(Descriptors3D.PMI1(mol), 3),
                "pmi2": round(Descriptors3D.PMI2(mol), 3),
                "pmi3": round(Descriptors3D.PMI3(mol), 3),
                "npr1": round(Descriptors3D.NPR1(mol), 3),
                "npr2": round(Descriptors3D.NPR2(mol), 3)
            }
            
            conformer_meta = {
                "method": "ETKDG + MMFF",
                "num_conformers": 1,
                "optimization": "MMFF force field"
            }
            
            return descriptors_3d, conformer_meta
            
        except Exception as e:
            return None, {"error": str(e)}
    
    def _estimate_descriptors(self, smiles: str) -> Dict[str, float]:
        """Estimate descriptors without RDKit (rule-based fallback)."""
        if not smiles:
            return {}
        
        # Count atoms (very rough estimation)
        carbon_count = smiles.upper().count('C') - smiles.upper().count('CL')
        oxygen_count = smiles.upper().count('O')
        nitrogen_count = smiles.upper().count('N')
        
        # Estimate molecular weight (very rough)
        est_mw = carbon_count * 12 + oxygen_count * 16 + nitrogen_count * 14 + \
                 (carbon_count * 2) * 1  # Rough H estimation
        
        # Estimate logP (hydrophobicity)
        est_logp = carbon_count * 0.5 - oxygen_count * 1.0 - nitrogen_count * 0.5
        
        return {
            "molecular_weight": round(est_mw, 1),
            "logp": round(max(-2, min(10, est_logp)), 1),
            "hbd": oxygen_count + nitrogen_count,  # Very rough
            "hba": oxygen_count + nitrogen_count,
            "num_rings": smiles.count('1') + smiles.count('2'),
            "estimation_method": "rule_based_fallback"
        }
    
    def _validate_coa(self, coa: COAInput) -> List[QCFlag]:
        """Validate Certificate of Analysis data.
        
        TRADE SECRET: Specific threshold values and validation rules.
        """
        flags: List[QCFlag] = []
        
        # Check for missing cannabinoid data
        if not coa.cannabinoids:
            flags.append(QCFlag(
                code="COA_MISSING_CANNABINOIDS",
                severity=QCSeverity.ERROR,
                message=self.QC_CODES["COA_MISSING_CANNABINOIDS"]
            ))
            return flags
        
        # Validate individual values
        total_cannabinoids = 0.0
        for cb in coa.cannabinoids:
            value = cb.get("value", 0)
            name = cb.get("name", "unknown")
            
            # Check for negative values
            if value < 0:
                flags.append(QCFlag(
                    code="COA_NEGATIVE_VALUE",
                    severity=QCSeverity.ERROR,
                    message=f"Negative value for {name}: {value}",
                    field=f"cannabinoids.{name}"
                ))
            
            total_cannabinoids += max(0, value)
        
        # Add terpenes to total if present
        total_terpenes = 0.0
        if coa.terpenes:
            for terp in coa.terpenes:
                total_terpenes += max(0, terp.get("value", 0))
        
        # Check total doesn't exceed 100%
        total = total_cannabinoids + total_terpenes
        if coa.units == "%" and total > self._THRESHOLDS["max_total_percent"]:
            flags.append(QCFlag(
                code="COA_TOTAL_EXCEEDS_100",
                severity=QCSeverity.ERROR,
                message=f"{self.QC_CODES['COA_TOTAL_EXCEEDS_100']} (Total: {total:.1f}%)"
            ))
        
        # Validate THC range
        if coa.total_thc is not None:
            if coa.total_thc > self._THRESHOLDS["max_thc_percent"]:
                flags.append(QCFlag(
                    code="COA_IMPLAUSIBLE_THC",
                    severity=QCSeverity.WARNING,
                    message=f"THC {coa.total_thc}% exceeds typical maximum"
                ))
        
        # Validate CBD range
        if coa.total_cbd is not None:
            if coa.total_cbd > self._THRESHOLDS["max_cbd_percent"]:
                flags.append(QCFlag(
                    code="COA_IMPLAUSIBLE_CBD",
                    severity=QCSeverity.WARNING,
                    message=f"CBD {coa.total_cbd}% exceeds typical maximum"
                ))
        
        # Check for LOD/LOQ
        if not coa.lod_loq:
            flags.append(QCFlag(
                code="COA_MISSING_LOD",
                severity=QCSeverity.INFO,
                message=self.QC_CODES["COA_MISSING_LOD"]
            ))
        
        # Check contaminants
        flags.extend(self._check_contaminants(coa))
        
        return flags
    
    def _check_contaminants(self, coa: COAInput) -> List[QCFlag]:
        """Check contaminant levels in COA.
        
        TRADE SECRET: Specific limits and state compliance thresholds.
        """
        flags: List[QCFlag] = []
        
        # Check heavy metals
        if coa.heavy_metals:
            for metal in coa.heavy_metals:
                # TRADE SECRET: Actual limits vary by state/jurisdiction
                if metal.get("result", "pass").lower() == "fail":
                    flags.append(QCFlag(
                        code="COA_HEAVY_METAL_FAIL",
                        severity=QCSeverity.ERROR,
                        message=f"Heavy metal {metal.get('name', 'unknown')} above limit"
                    ))
        
        # Check pesticides
        if coa.pesticides:
            for pest in coa.pesticides:
                if pest.get("result", "pass").lower() == "fail":
                    flags.append(QCFlag(
                        code="COA_PESTICIDE_FAIL",
                        severity=QCSeverity.ERROR,
                        message=f"Pesticide {pest.get('name', 'unknown')} detected"
                    ))
        
        # Check residual solvents
        if coa.solvents:
            for solvent in coa.solvents:
                if solvent.get("result", "pass").lower() == "fail":
                    flags.append(QCFlag(
                        code="COA_SOLVENT_FAIL",
                        severity=QCSeverity.ERROR,
                        message=f"Residual solvent {solvent.get('name', 'unknown')} above limit"
                    ))
        
        # Check microbials
        if coa.microbials:
            for micro in coa.microbials:
                if micro.get("result", "pass").lower() == "fail":
                    flags.append(QCFlag(
                        code="COA_MICROBIAL_FAIL",
                        severity=QCSeverity.ERROR,
                        message=f"Microbial {micro.get('name', 'unknown')} above limit"
                    ))
        
        return flags
    
    def _calculate_completeness(
        self,
        compound: CompoundInput,
        coa: Optional[COAInput],
        normalized: NormalizedStructure,
        descriptors: DescriptorSet
    ) -> float:
        """Calculate data completeness score.
        
        TRADE SECRET: Specific weights for each component.
        """
        score = 0.0
        max_score = 0.0
        
        # Structure completeness (40% weight)
        max_score += 40
        if normalized.is_valid:
            score += 20  # Valid structure
            if normalized.canonical_smiles:
                score += 5
            if normalized.inchi_key:
                score += 5
            if normalized.molecular_formula:
                score += 5
            if len(descriptors.descriptors_2d) >= 10:
                score += 5
        
        # 3D descriptors bonus (10% weight)
        max_score += 10
        if descriptors.descriptors_3d:
            score += 10
        
        # COA completeness (50% weight)
        max_score += 50
        if coa:
            score += 10  # COA provided
            if coa.cannabinoids:
                score += 10
            if coa.terpenes:
                score += 5
            if coa.total_thc is not None or coa.total_cbd is not None:
                score += 5
            if coa.heavy_metals:
                score += 5
            if coa.pesticides:
                score += 5
            if coa.solvents:
                score += 5
            if coa.lod_loq:
                score += 5
        
        return round((score / max_score) * 100, 1)
    
    def get_qc_codes(self) -> Dict[str, str]:
        """Get available QC codes and descriptions."""
        return self.QC_CODES.copy()
    
    def get_known_cannabinoids(self) -> List[str]:
        """Get list of known cannabinoid names."""
        return list(self.KNOWN_CANNABINOIDS.keys())
    
    def validate_smiles(self, smiles: str) -> Dict[str, Any]:
        """Quick SMILES validation without full analysis."""
        flags: List[QCFlag] = []
        notes: List[str] = []
        
        if self._rdkit_available:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(smiles)
            is_valid = mol is not None
            if is_valid:
                canonical = Chem.MolToSmiles(mol, canonical=True)
            else:
                canonical = smiles
                flags.append(QCFlag(
                    code="STRUCTURE_INVALID",
                    severity=QCSeverity.ERROR,
                    message="Failed to parse SMILES"
                ))
        else:
            canonical, is_valid = self._normalize_fallback(
                CompoundInput(name="validation", smiles=smiles),
                flags,
                notes
            )
        
        return {
            "is_valid": is_valid,
            "canonical_smiles": canonical,
            "flags": [f.to_dict() for f in flags],
            "notes": notes
        }
