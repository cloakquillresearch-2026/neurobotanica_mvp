"""
ChemPath API Router - Chemical Characterization Endpoints

REST API for ChemPath molecular characterization services.
"""

from fastapi import APIRouter, HTTPException, Query
from pydantic import BaseModel, Field
from typing import Optional, List, Dict, Any
from enum import Enum

from backend.services.chempath import (
    ChemPathAnalyzer,
    ChemPathReportGenerator,
    ReportFormat,
    CompoundInput as CompoundInputModel,
    COAInput as COAInputModel,
    ChemPathRequest as ChemPathRequestModel
)

router = APIRouter(prefix="/api/v1/chempath", tags=["chempath"])


# ============================================================================
# Pydantic Request/Response Models
# ============================================================================

class CompoundInputSchema(BaseModel):
    """Compound structure input."""
    name: str = Field(..., description="Compound name")
    smiles: Optional[str] = Field(None, description="SMILES notation")
    inchi: Optional[str] = Field(None, description="InChI notation")
    inchi_key: Optional[str] = Field(None, description="InChIKey")
    source: str = Field("customer_private", description="Data source")


class CannabinoidEntry(BaseModel):
    """Single cannabinoid entry in COA."""
    name: str
    value: float


class TerpeneEntry(BaseModel):
    """Single terpene entry in COA."""
    name: str
    value: float


class ContaminantEntry(BaseModel):
    """Contaminant test result."""
    name: str
    value: Optional[float] = None
    result: str = "pass"  # pass/fail


class COAInputSchema(BaseModel):
    """Certificate of Analysis input."""
    lab_name: str = Field(..., description="Testing laboratory name")
    sample_id: str = Field(..., description="Sample identifier")
    units: str = Field(..., description="Concentration units: '%' or 'mg/g'")
    cannabinoids: List[CannabinoidEntry] = Field(..., description="Cannabinoid concentrations")
    terpenes: Optional[List[TerpeneEntry]] = Field(None, description="Terpene concentrations")
    total_thc: Optional[float] = Field(None, description="Total THC percentage")
    total_cbd: Optional[float] = Field(None, description="Total CBD percentage")
    moisture: Optional[float] = Field(None, description="Moisture content percentage")
    solvents: Optional[List[ContaminantEntry]] = Field(None, description="Residual solvent results")
    microbials: Optional[List[ContaminantEntry]] = Field(None, description="Microbial test results")
    heavy_metals: Optional[List[ContaminantEntry]] = Field(None, description="Heavy metal results")
    pesticides: Optional[List[ContaminantEntry]] = Field(None, description="Pesticide results")
    lod_loq: Optional[Dict[str, float]] = Field(None, description="LOD/LOQ values")


class ChemPathAnalyzeRequest(BaseModel):
    """Request model for ChemPath analysis."""
    compound: CompoundInputSchema
    coa: Optional[COAInputSchema] = Field(None, description="Certificate of Analysis data")
    compute_3d: bool = Field(False, description="Compute 3D conformer and descriptors")


class SMILESValidationRequest(BaseModel):
    """Request for SMILES validation."""
    smiles: str = Field(..., description="SMILES string to validate")


class ReportFormatEnum(str, Enum):
    """Report format options."""
    markdown = "markdown"
    html = "html"
    json = "json"


# ============================================================================
# Helper Functions
# ============================================================================

def _convert_compound_input(schema: CompoundInputSchema) -> CompoundInputModel:
    """Convert Pydantic schema to dataclass."""
    return CompoundInputModel(
        name=schema.name,
        smiles=schema.smiles,
        inchi=schema.inchi,
        inchi_key=schema.inchi_key,
        source=schema.source
    )


def _convert_coa_input(schema: COAInputSchema) -> COAInputModel:
    """Convert COA Pydantic schema to dataclass."""
    return COAInputModel(
        lab_name=schema.lab_name,
        sample_id=schema.sample_id,
        units=schema.units,
        cannabinoids=[{"name": c.name, "value": c.value} for c in schema.cannabinoids],
        terpenes=[{"name": t.name, "value": t.value} for t in schema.terpenes] if schema.terpenes else None,
        total_thc=schema.total_thc,
        total_cbd=schema.total_cbd,
        moisture=schema.moisture,
        solvents=[{"name": s.name, "value": s.value, "result": s.result} for s in schema.solvents] if schema.solvents else None,
        microbials=[{"name": m.name, "value": m.value, "result": m.result} for m in schema.microbials] if schema.microbials else None,
        heavy_metals=[{"name": h.name, "value": h.value, "result": h.result} for h in schema.heavy_metals] if schema.heavy_metals else None,
        pesticides=[{"name": p.name, "value": p.value, "result": p.result} for p in schema.pesticides] if schema.pesticides else None,
        lod_loq=schema.lod_loq
    )


# ============================================================================
# Analysis Endpoints
# ============================================================================

@router.post("/analyze")
async def analyze_compound(request: ChemPathAnalyzeRequest) -> Dict[str, Any]:
    """Run full ChemPath analysis pipeline.
    
    Performs:
    - Structure normalization and validation
    - 2D descriptor computation
    - Optional 3D conformer generation
    - COA quality check validation
    - Data completeness scoring
    
    Returns complete analysis results including QC flags and descriptors.
    """
    try:
        analyzer = ChemPathAnalyzer()
        
        # Convert to internal models
        compound = _convert_compound_input(request.compound)
        coa = _convert_coa_input(request.coa) if request.coa else None
        
        # Build request
        chempath_request = ChemPathRequestModel(
            compound_input=compound,
            coa_input=coa,
            compute_3d=request.compute_3d
        )
        
        # Run analysis
        response = analyzer.analyze(chempath_request)
        
        return {
            "status": "success",
            "analysis": response.to_dict()
        }
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/validate-smiles")
async def validate_smiles(request: SMILESValidationRequest) -> Dict[str, Any]:
    """Quick SMILES validation without full analysis.
    
    Returns validity status, canonical form, and any structural flags.
    """
    try:
        analyzer = ChemPathAnalyzer()
        result = analyzer.validate_smiles(request.smiles)
        
        return {
            "status": "success",
            "validation": result
        }
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


# ============================================================================
# Report Endpoints
# ============================================================================

@router.post("/report")
async def generate_report(
    request: ChemPathAnalyzeRequest,
    format: ReportFormatEnum = Query(ReportFormatEnum.markdown, description="Report format")
) -> Dict[str, Any]:
    """Analyze compound and generate formatted report.
    
    Runs full analysis pipeline then generates report in requested format.
    """
    try:
        analyzer = ChemPathAnalyzer()
        report_gen = ChemPathReportGenerator()
        
        # Convert to internal models
        compound = _convert_compound_input(request.compound)
        coa = _convert_coa_input(request.coa) if request.coa else None
        
        # Build request
        chempath_request = ChemPathRequestModel(
            compound_input=compound,
            coa_input=coa,
            compute_3d=request.compute_3d
        )
        
        # Run analysis
        response = analyzer.analyze(chempath_request)
        
        # Generate report
        report_format = ReportFormat(format.value)
        report = report_gen.generate_report(response, format=report_format)
        
        return {
            "status": "success",
            "job_id": response.chempath_job_id,
            "format": format.value,
            "report": report
        }
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/summary")
async def generate_summary(request: ChemPathAnalyzeRequest) -> Dict[str, Any]:
    """Generate brief analysis summary.
    
    Returns key metrics and status without full report.
    """
    try:
        analyzer = ChemPathAnalyzer()
        report_gen = ChemPathReportGenerator()
        
        # Convert to internal models
        compound = _convert_compound_input(request.compound)
        coa = _convert_coa_input(request.coa) if request.coa else None
        
        # Build request
        chempath_request = ChemPathRequestModel(
            compound_input=compound,
            coa_input=coa,
            compute_3d=request.compute_3d
        )
        
        # Run analysis
        response = analyzer.analyze(chempath_request)
        
        # Generate summary
        summary = report_gen.generate_summary(response)
        
        return {
            "status": "success",
            "summary": summary
        }
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


# ============================================================================
# Reference Endpoints
# ============================================================================

@router.get("/qc-codes")
async def get_qc_codes() -> Dict[str, Any]:
    """Get available QC flag codes and descriptions.
    
    Returns all possible QC codes that may appear in analysis results.
    """
    analyzer = ChemPathAnalyzer()
    
    return {
        "status": "success",
        "qc_codes": analyzer.get_qc_codes()
    }


@router.get("/known-cannabinoids")
async def get_known_cannabinoids() -> Dict[str, Any]:
    """Get list of known cannabinoid compounds.
    
    Returns cannabinoids with validated SMILES in the system.
    """
    analyzer = ChemPathAnalyzer()
    
    return {
        "status": "success",
        "cannabinoids": analyzer.get_known_cannabinoids()
    }


@router.get("/report-formats")
async def get_report_formats() -> Dict[str, Any]:
    """Get available report output formats."""
    return {
        "status": "success",
        "formats": [
            {
                "format": "markdown",
                "description": "Markdown text format",
                "content_type": "text/markdown"
            },
            {
                "format": "html",
                "description": "HTML document format",
                "content_type": "text/html"
            },
            {
                "format": "json",
                "description": "JSON structured data",
                "content_type": "application/json"
            }
        ]
    }


@router.get("/descriptor-list")
async def get_descriptor_list() -> Dict[str, Any]:
    """Get list of computed molecular descriptors.
    
    Returns descriptors computed during analysis with descriptions.
    """
    return {
        "status": "success",
        "descriptors_2d": {
            "molecular_weight": "Molecular weight in Daltons",
            "exact_mass": "Exact monoisotopic mass",
            "heavy_atom_count": "Number of non-hydrogen atoms",
            "logp": "Calculated partition coefficient (octanol/water)",
            "tpsa": "Topological polar surface area in Å²",
            "hbd": "Number of hydrogen bond donors",
            "hba": "Number of hydrogen bond acceptors",
            "rotatable_bonds": "Number of rotatable bonds",
            "num_rings": "Total ring count",
            "num_aromatic_rings": "Aromatic ring count",
            "num_aliphatic_rings": "Aliphatic ring count",
            "fraction_csp3": "Fraction of sp3 carbons",
            "num_stereocenters": "Number of stereogenic centers",
            "lipinski_violations": "Rule of 5 violations (0-4)",
            "qed": "Quantitative Estimate of Drug-likeness (0-1)",
            "labute_asa": "Labute Approximate Surface Area"
        },
        "descriptors_3d": {
            "radius_of_gyration": "3D size measure",
            "inertial_shape_factor": "Shape anisotropy",
            "eccentricity": "Deviation from spherical shape",
            "asphericity": "Measure of asphericity",
            "spherocity_index": "Sphericity measure",
            "pmi1": "Principal moment of inertia 1",
            "pmi2": "Principal moment of inertia 2",
            "pmi3": "Principal moment of inertia 3",
            "npr1": "Normalized PMI ratio 1",
            "npr2": "Normalized PMI ratio 2"
        }
    }


# ============================================================================
# Health Check
# ============================================================================

@router.get("/health")
async def health_check() -> Dict[str, Any]:
    """ChemPath service health check."""
    analyzer = ChemPathAnalyzer()
    
    return {
        "status": "healthy",
        "service": "chempath",
        "rdkit_available": analyzer._rdkit_available,
        "version": "1.0.0"
    }
