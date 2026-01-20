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
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware


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


@app.get("/")
async def root():
    """Root endpoint returning API information."""
    return {
        "message": "🌿 NeuroBotanica API - Dimeric Cannabinoid Therapeutic Prediction System",
        "version": "0.4.0",
        "status": "operational",
        "endpoints": {
            "health": "/health",
            "api_docs": "/docs",
            "api_redoc": "/redoc"
        }
    }


@app.get("/health")
async def health_check():
    """Health check endpoint."""
    return {"status": "healthy", "timestamp": "2024-01-01T00:00:00Z"}