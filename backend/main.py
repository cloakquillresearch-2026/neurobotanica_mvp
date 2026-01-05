"""
NeuroBotanica API - Main Application
Dimeric Cannabinoid Therapeutic Prediction System

FastAPI application with:
- 320-study evidence database access
- FDA compliance documentation generation
- Pharmacology data package APIs
- Comparative efficacy analysis
- 3D Conformer generation (ETKDG method)
- Patient & Treatment management
- OmniPath integration for provenance tracking
- Clinical Evidence API with confidence weighting
- Receptor Affinity with provenance and heterogeneity analysis
"""
from fastapi import FastAPI, Depends, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from sqlalchemy.orm import Session
from sqlalchemy import func
from contextlib import asynccontextmanager
import os

from backend.models.database import get_db, init_db, engine, Base
from backend.api import studies, compounds, fda_compliance, conformers
from backend.api import omnipath
from backend.api import evidence, receptor_affinity
from backend.api import dimers
from backend.routers import patentpath, terpenes, chempath, toxpath, regpath, security, genomepath, biopath, clinpath, dispensary
from backend.middleware.token_validation import TokenValidationMiddleware


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Application lifespan - initialize database on startup."""
    print("🌿 NeuroBotanica API starting...")
    init_db()
    print("✅ Database initialized")
    yield
    print("🛑 NeuroBotanica API shutting down...")


app = FastAPI(
    title="NeuroBotanica API",
    description="""
    ## Dimeric Cannabinoid Therapeutic Prediction System
    
    **Patent Claims Support:**
    - FDA Schedule III compliance documentation (Claim 1j)
    - Pharmacology data packages with mechanism of action
    - Comparative efficacy analysis from 368-study database
    - CMC templates for cannabis pharmaceutical approval
    - 3D conformer generation for molecular modeling
    - Clinical evidence aggregation with confidence weighting
    - Receptor affinity data with full provenance tracking
    - Trade Secret Engines: ChemPath, ToxPath, RegPath, BioPath, ClinPath, GenomePath
    
    **Data Assets:**
    - 368 clinical studies across 22+ conditions
    - 63 cannabinoid compounds with molecular descriptors
    - 10,084 dimer/entourage validated entries
    - 3D conformer ensembles (ETKDG method)
    - Dimeric triangulation predictions
    - FDA-approved drug evidence (Epidiolex, Marinol, Cesamet, Sativex)
    - Receptor affinity database with assay-level provenance
    """,
    version="0.4.0",
    # lifespan=lifespan
)

# CORS middleware for frontend access
app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        "http://localhost:3000",
        "http://localhost:5173",  # Vite dev server
        "https://neurobotanica.cloakandquill.org"
    ],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


# Include routers
app.include_router(studies.router, prefix="/api/v1/studies", tags=["Clinical Studies"])
app.include_router(compounds.router, prefix="/api/v1/compounds", tags=["Cannabinoid Compounds"])
app.include_router(fda_compliance.router, prefix="/api/v1/fda", tags=["FDA Compliance"])
app.include_router(conformers.router, prefix="/api/v1", tags=["3D Conformers"])
app.include_router(omnipath.router, prefix="/api/v1", tags=["OmniPath Integration"])
app.include_router(evidence.router, prefix="/api/v1", tags=["Clinical Evidence"])
app.include_router(receptor_affinity.router, prefix="/api/v1", tags=["Receptor Affinity"])
app.include_router(dimers.router, tags=["Dimers"])
app.include_router(patentpath.router, prefix="/api/v1", tags=["PatentPath Lite"])
app.include_router(terpenes.router, prefix="/api/v1", tags=["Terpene Analysis"])
app.include_router(chempath.router, tags=["ChemPath"])
app.include_router(toxpath.router, tags=["ToxPath"])
app.include_router(regpath.router, tags=["RegPath"])
app.include_router(genomepath.router, tags=["GenomePath"])  # Router already has /api/genomepath prefix
app.include_router(biopath.router, prefix="/api/biopath", tags=["BioPath"])
app.include_router(clinpath.router, prefix="/api/clinpath", tags=["ClinPath"])
app.include_router(dispensary.router, prefix="/api/dispensary", tags=["Dispensary"])
app.include_router(security.router, tags=["Security"])

# Add token validation middleware (disabled by default for development)
# Enable with NEUROBOTANICA_TOKEN_VALIDATION=true environment variable
if os.getenv("NEUROBOTANICA_TOKEN_VALIDATION", "false").lower() == "true":
    app.add_middleware(TokenValidationMiddleware, validate_tokens=True)
else:
    # Add middleware but don't validate (allows testing without auth)
    app.add_middleware(TokenValidationMiddleware, validate_tokens=False)


@app.get("/")
async def root():
    """Root endpoint with API information."""
    try:
        return {
            "message": "🌿 NeuroBotanica API",
            "version": "0.4.0",
            "description": "Dimeric Cannabinoid Therapeutic Prediction System",
            "patent_claims": ["1j - FDA Schedule III Compliance Support"],
            "data_assets": {
                "clinical_studies": 368,
                "conditions": 22,
                "compounds": 63,
                "dimer_entourage_entries": 10084,
                "fda_approved_drugs": ["Epidiolex", "Marinol", "Cesamet", "Sativex"]
            },
            "trade_secret_engines": [
                "ChemPath", "ToxPath", "RegPath", 
                "BioPath", "ClinPath", "GenomePath"
            ],
            "endpoints": {
                "studies": "/api/v1/studies",
                "compounds": "/api/v1/compounds",
                "fda_compliance": "/api/v1/fda",
                "conformers": "/api/v1/conformers",
                "omnipath": "/api/v1/omnipath",
                "dimers": "/api/dimers",
                "chempath": "/api/chempath",
                "toxpath": "/api/toxpath",
                "regpath": "/api/regpath",
                "genomepath": "/api/genomepath",
                "biopath": "/api/biopath",
                "clinpath": "/api/clinpath",
                "docs": "/docs"
            }
        }
    except Exception as e:
        import logging
        logger = logging.getLogger(__name__)
        logger.error(f"Error in root endpoint: {e}")
        import traceback
        traceback.print_exc()
        return {"status": "error", "message": str(e)}


@app.get("/test")
async def test_endpoint():
    """Simple test endpoint without database dependency."""
    try:
        return {"status": "ok", "message": "Test endpoint working"}
    except Exception as e:
        import logging
        logger = logging.getLogger(__name__)
        logger.error(f"Error in test endpoint: {e}")
        import traceback
        traceback.print_exc()
        return {"status": "error", "message": str(e)}


@app.get("/api/v1/stats")
async def get_stats(db: Session = Depends(get_db)):
    """Get database statistics."""
    from backend.models.study import ClinicalStudy
    from backend.models.compound import Cannabinoid, DimericPrediction
    
    try:
        study_count = db.query(ClinicalStudy).count()
        compound_count = db.query(Cannabinoid).count()
        dimer_count = db.query(DimericPrediction).count()
        
        # Get condition breakdown
        conditions = db.query(
            ClinicalStudy.condition,
            func.count(ClinicalStudy.id)
        ).group_by(ClinicalStudy.condition).all()
        
        return {
            "total_studies": study_count,
            "total_compounds": compound_count,
            "dimeric_predictions": dimer_count,
            "conditions": {c[0]: c[1] for c in conditions},
            "fda_approved_coverage": ["Epidiolex", "Marinol", "Cesamet", "Sativex"]
        }
    except Exception as e:
        return {
            "total_studies": 0,
            "total_compounds": 0,
            "dimeric_predictions": 0,
            "message": "Database not yet populated",
            "note": "Run data loading scripts to populate"
        }


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000, reload=True)
