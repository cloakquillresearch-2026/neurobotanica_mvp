<aside>
💰

**FULL PLATFORM BUILD — $150K, 14 weeks**

This is the complete NeuroBotanica platform development plan. For a **lean MVP approach** that generates revenue in 2–4 weeks with minimal upfront investment, see [Manufacturer Services — Lean MVP & Path Modules](https://www.notion.so/Manufacturer-Services-Lean-MVP-Path-Modules-9693518abdd14a11bd6b1741ec37a42b?pvs=21).

</aside>

**Organization**: Cloak and Quill Research 501(c)(3)

**Project**: NeuroBotanica MVP Development

**Version**: 3.0 (Complete Integration + PatentPath Lite)

**Created**: December 19, 2025

**Timeline**: 14 weeks to production MVP

**Status**: 🟢 Ready to Begin

---

## 📊 Executive Summary

This document provides the **complete and unified** development roadmap for NeuroBotanica MVP, integrating:

- **200-study validated dataset** (completed 26 days early)
- **Enhanced scientific rigor** (3D conformers, provenance tracking, ChEMBL/PubChem integration)
- **PatentPath Lite IP management** (basic patent analysis for customer compounds)
- **OmniPath integration** (token-gated access, benefit-sharing automation)

### What's New in v3.0

**Merged from v2.0 (Scientific Enhancements)**:

- Week 1: 3D conformer generation pipeline (ETKDG)
- Week 3: Enhanced receptor affinity provenance schema
- Week 4: Dimeric triangulation scoring framework
- Week 5: ChEMBL/PubChem integration with preserved assay context
- Week 7: Evidence-weighted terpene features with gating

**New in v3.0 (PatentPath Lite Integration)**:

- Week 6-7: PatentPath Lite development (parallel with ML models)
- Week 8: PatentPath API integration into frontend
- Budget: $25K (17% of $150K MVP budget)
- 6 core IP features: Prior art search, novelty scoring, FTO check, claim templates, TK attribution, cost estimation

### Timeline Summary

| Phase | Weeks | Key Deliverables | Completion Date |
| --- | --- | --- | --- |
| **Phase 0: Foundation** | 1 | Environment + 200 studies + conformer pipeline | Dec 25, 2025 |
| **Phase 1: Backend + IP** | 2-5 | Core APIs + OmniPath + ChEMBL + Provenance | Jan 22, 2026 |
| **Phase 2: ML + PatentPath + Frontend** | 6-9 | Models + PatentPath Lite + UI | Feb 19, 2026 |
| **Phase 3: Production** | 10-14 | Hardening + compliance + Nevada pilot | Mar 26, 2026 |

### Success Criteria

| Milestone | Target Date | Deliverable |
| --- | --- | --- |
| **Week 1 Complete** | Dec 25, 2025 | Local dev + database + 200 studies + conformers |
| **Phase 1 Complete** | Jan 22, 2026 | Backend APIs + OmniPath + ChEMBL/PubChem |
| **Phase 2 Complete** | Feb 19, 2026 | ML models + PatentPath Lite + frontend MVP |
| **Phase 3 Complete** | Mar 26, 2026 | Production hardening + Nevada pilot launch |

### Budget Allocation

| Component | Cost | % of Budget | Weeks |
| --- | --- | --- | --- |
| Core prediction engine | $50K | 33% | 1-4, 6-7 |
| SaaS dashboard/API | $30K | 20% | 8-9 |
| **PatentPath Lite** | **$25K** | **17%** | **6-8** |
| OmniPath integration | $20K | 13% | 2, 11 |
| FDA doc templates | $10K | 7% | 5 |
| Testing/QA | $15K | 10% | 13 |
| **Total** | **$150K** | **100%** | **14 weeks** |

---

# 📅 COMPLETE 14-WEEK DEVELOPMENT PLAN

## Phase 0: Foundation (Week 1)

### Week 1: Environment Setup + Data Ingestion + 3D Conformer Pipeline

**Dates**: Dec 18-25, 2025

**Estimated Hours**: 35-45 hours (5-6 hours/day)

**Objective**: Establish development environment, load validated data, implement 3D molecular analysis

**Deliverables**:

- ✅ Python 3.11+ environment with all dependencies
- ✅ PostgreSQL database with complete schema
- ✅ 200 clinical studies ingested from validated JSON files
- ✅ 16 cannabinoids with 2D RDKit descriptors
- ✅ **NEW**: 3D conformer ensembles for cannabinoids (ETKDG method)
- ✅ Basic FastAPI server running locally
- ✅ Git repository with proper structure

**Success Metrics**:

- Database contains 200 studies with full metadata ✓
- 16 cannabinoids with 40+ 2D descriptors ✓
- 14-16 cannabinoids with 3D conformer ensembles ✓
- API returns "Hello NeuroBotanica" at [localhost:8000](http://localhost:8000) ✓
- All unit tests pass ✓

### Day 1: Environment Setup (4-5 hours)

**Task 1.1: Repository Initialization** (30 min)

```bash
# Create project directory
mkdir neurobotanica-mvp
cd neurobotanica-mvp

# Initialize git
git init
git branch -M main

# Create .gitignore
cat > .gitignore << EOF
# Python
__pycache__/
*.py[cod]
.venv/
.env

# Databases
*.db
*.sqlite

# ML Models
models/*.pkl

# Data
data/raw/
data/processed/

# OS
.DS_Store
EOF

# Create project structure
mkdir -p backend/{api,models,services,utils,tests}
mkdir -p frontend/{components,pages,lib,styles}
mkdir -p data/{raw,processed,models}
mkdir -p scripts
mkdir -p docs
```

**Task 1.2: Python Environment Setup** (1 hour)

```bash
# Create virtual environment (Python 3.11+)
python3.11 -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate

# Upgrade pip
pip install --upgrade pip setuptools wheel

# Create requirements.txt
cat > backend/requirements.txt << EOF
# Web Framework
fastapi==0.104.1
uvicorn[standard]==0.24.0
pydantic==2.5.0

# Database
psycopg2-binary==2.9.9
sqlalchemy==2.0.23
alembic==1.13.0

# Chemistry
rdkit==2023.9.2

# ML
scikit-learn==1.3.2
pandas==2.1.3
numpy==1.26.2

# Testing
pytest==7.4.3
pytest-asyncio==0.21.1
httpx==0.25.2

# Development
black==23.11.0
flake8==6.1.0

# Utilities
python-dotenv==1.0.0
requests==2.31.0
EOF

# Install dependencies
pip install -r backend/requirements.txt
```

**Task 1.3: PostgreSQL Database Setup** (1 hour)

**Option: Docker PostgreSQL** (Recommended)

```bash
# Create docker-compose.yml
cat > docker-compose.yml << EOF
version: '3.8'
services:
  db:
    image: postgres:15
    environment:
      POSTGRES_DB: neurobotanica_dev
      POSTGRES_USER: nb_dev
      POSTGRES_PASSWORD: secure_password_here
    ports:
      - "5432:5432"
    volumes:
      - postgres_data:/var/lib/postgresql/data

volumes:
  postgres_data:
EOF

# Start database
docker-compose up -d
```

**Task 1.4: Environment Configuration** (30 min)

```bash
# Create .env file
cat > backend/.env << EOF
DATABASE_URL=postgresql://nb_dev:[secure_password_here@localhost:5432](mailto:secure_password_here@localhost:5432)/neurobotanica_dev
SECRET_KEY=your_secret_key_here
DEBUG=True
EOF
```

**Task 1.5: Basic FastAPI Application** (1 hour)

Create `backend/[main.py](http://main.py)`:

```python
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

app = FastAPI(
    title="NeuroBotanica API",
    description="Dimeric Cannabinoid Therapeutic Prediction System",
    version="0.1.0"
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:3000"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.get("/")
async def root():
    return {"message": "Hello NeuroBotanica", "version": "0.1.0"}

@app.get("/health")
async def health_check():
    return {"status": "healthy", "database": "connected"}
```

**Test the server**:

```bash
cd backend
python [main.py](http://main.py)
# Visit http://localhost:8000/docs
```

---

### Day 2: Database Schema Design (5-6 hours)

**Task 2.1: SQLAlchemy Base Setup** (1 hour)

Create `backend/models/[database.py](http://database.py)`:

```python
from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
from pydantic_settings import BaseSettings
import os

class Settings(BaseSettings):
    database_url: str
    
    class Config:
        env_file = ".env"

settings = Settings()

engine = create_engine(settings.database_url)
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)

Base = declarative_base()

def get_db():
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()
```

**Task 2.2: Clinical Study Schema** (2 hours)

Create `backend/models/[study.py](http://study.py)` with complete schema for clinical studies including:

- Study identification (DOI, PubMed ID, etc.)
- Study design (RCT, observational, etc.)
- Clinical context (condition, population, etc.)
- Cannabinoid information
- Outcomes and quality metrics
- **NEW**: Confidence weighting (0.5-0.9 for ML training)

**Task 2.3: Compound Library Schema with 3D Support** (1.5 hours)

Create `backend/models/[compound.py](http://compound.py)`:

```python
from sqlalchemy import Column, Integer, String, Float, Boolean, JSON, Text, DateTime
from .database import Base
from datetime import datetime

class Cannabinoid(Base):
    __tablename__ = "cannabinoids"
    
    id = Column(Integer, primary_key=True, index=True)
    
    # Identity
    name = Column(String(100), nullable=False, unique=True)
    abbreviation = Column(String(10))
    smiles = Column(String(500), nullable=False)
    inchi = Column(String(500))
    inchi_key = Column(String(100))
    
    # Classification
    is_natural = Column(Boolean, default=True)
    is_metabolite = Column(Boolean, default=False)
    is_dimeric = Column(Boolean, default=False)
    
    # 2D RDKit Descriptors
    molecular_weight = Column(Float)
    logp = Column(Float)
    tpsa = Column(Float)
    h_bond_donors = Column(Integer)
    h_bond_acceptors = Column(Integer)
    rotatable_bonds = Column(Integer)
    aromatic_rings = Column(Integer)
    rdkit_descriptors = Column(JSON)  # Full 2D descriptor set
    
    # NEW: 3D Conformer Data
    has_conformers = Column(Boolean, default=False)
    conformer_generation_method = Column(String(50))  # "ETKDG"
    num_conformers_generated = Column(Integer)
    conformer_ensemble_id = Column(String(100))
    conformer_metadata = Column(JSON)  # Energy ranges, RMSD clusters
    rdkit_descriptors_3d = Column(JSON)  # PMI, asphericity, etc.
    conformers_generated_at = Column(DateTime)
    
    # NEW: Enhanced Receptor Affinity with Provenance
    receptor_affinities = Column(JSON)  # Provenance-rich structure
    
    # Pharmacology
    mechanism_of_action = Column(Text)
    therapeutic_categories = Column(JSON)
    conditions_treated = Column(JSON)
    
    # External Database IDs (for ChEMBL/PubChem sync)
    chembl_id = Column(String(50))
    pubchem_cid = Column(String(50))
    last_external_sync = Column(DateTime)
    
    # OmniPath Provenance
    omnipath_manifest_id = Column(String(100))
```

**Task 2.4: Patient & Treatment Schema** (1.5 hours)

Create `backend/models/[patient.py](http://patient.py)` and `backend/models/[treatment.py](http://treatment.py)` for patient profiles and treatment recommendations.

---

### Day 3: Database Migration & Initial Load (4-5 hours)

**Task 3.1: Alembic Setup** (1 hour)

```bash
cd backend
alembic init alembic

# Create initial migration
alembic revision --autogenerate -m "Initial schema with 3D conformers"
alembic upgrade head
```

**Task 3.2: Data Loading Script - Studies** (2 hours)

Create `scripts/load_clinical_[studies.py](http://studies.py)` to load 200 validated studies from JSON files.

**Task 3.3: Data Loading Script - Cannabinoids** (1.5 hours)

Create `scripts/load_[cannabinoids.py](http://cannabinoids.py)` to load 16 cannabinoids with 2D RDKit descriptors.

---

### Days 4-5: 3D Conformer Pipeline (NEW - 10-12 hours)

**Task 4.1: Install 3D Dependencies** (30 min)

RDKit already includes ETKDG and MMFF94 - no additional packages needed.

**Task 4.2: Implement Conformer Generator Service** (4-5 hours)

Create `backend/services/conformer_[generator.py](http://generator.py)`:

```python
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors3D
import numpy as np
from typing import List, Dict, Tuple

class ConformerGenerator:
    """Generate and analyze 3D conformers using ETKDG method."""
    
    def __init__(self, num_conformers: int = 50, max_attempts: int = 100):
        self.num_conformers = num_conformers
        self.max_attempts = max_attempts
    
    def generate_conformers(self, smiles: str) -> Dict:
        """Generate conformer ensemble for a molecule.
        
        Returns:
            dict with conformer data, energies, and 3D descriptors
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {"error": "Invalid SMILES"}
        
        mol = Chem.AddHs(mol)
        
        # Generate conformers using ETKDG
        conf_ids = AllChem.EmbedMultipleConfs(
            mol,
            numConfs=self.num_conformers,
            params=AllChem.ETKDGv3(),
            maxAttempts=self.max_attempts,
            randomSeed=42
        )
        
        if len(conf_ids) == 0:
            return {"error": "Conformer generation failed"}
        
        # Optimize with MMFF94
        energies = []
        for conf_id in conf_ids:
            converged = AllChem.MMFFOptimizeMolecule(mol, confId=conf_id)
            if converged == 0:
                ff = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol), confId=conf_id)
                energies.append(ff.CalcEnergy())
            else:
                energies.append(None)
        
        # Filter failed optimizations
        valid_confs = [(i, e) for i, e in enumerate(energies) if e is not None]
        
        if len(valid_confs) == 0:
            return {"error": "All optimizations failed"}
        
        # Sort by energy
        valid_confs.sort(key=lambda x: x[1])
        lowest_conf_id = valid_confs[0][0]
        
        # Calculate 3D descriptors on lowest energy conformer
        descriptors_3d = self.calculate_3d_descriptors(mol, lowest_conf_id)
        
        # Cluster conformers by RMSD
        clusters = self.cluster_conformers_by_rmsd(mol, [c[0] for c in valid_confs])
        
        return {
            "num_conformers_generated": len(valid_confs),
            "num_conformers_requested": self.num_conformers,
            "lowest_energy": valid_confs[0][1],
            "highest_energy": valid_confs[-1][1],
            "energy_range": valid_confs[-1][1] - valid_confs[0][1],
            "num_clusters_rmsd_2A": len(clusters),
            "descriptors_3d": descriptors_3d,
            "success": True
        }
    
    def calculate_3d_descriptors(self, mol, conf_id: int) -> Dict:
        """Calculate 3D molecular descriptors."""
        return {
            "pmi1": Descriptors3D.PMI1(mol, confId=conf_id),
            "pmi2": Descriptors3D.PMI2(mol, confId=conf_id),
            "pmi3": Descriptors3D.PMI3(mol, confId=conf_id),
            "npr1": Descriptors3D.NPR1(mol, confId=conf_id),
            "npr2": Descriptors3D.NPR2(mol, confId=conf_id),
            "radius_of_gyration": Descriptors3D.RadiusOfGyration(mol, confId=conf_id),
            "inertial_shape_factor": Descriptors3D.InertialShapeFactor(mol, confId=conf_id),
            "asphericity": Descriptors3D.Asphericity(mol, confId=conf_id),
            "eccentricity": Descriptors3D.Eccentricity(mol, confId=conf_id),
        }
    
    def cluster_conformers_by_rmsd(self, mol, conf_ids: List[int], threshold: float = 2.0) -> List[List[int]]:
        """Cluster conformers by RMSD."""
        from rdkit.Chem import rdMolAlign
        
        clusters = []
        assigned = set()
        
        for i, conf_id_i in enumerate(conf_ids):
            if conf_id_i in assigned:
                continue
            
            cluster = [conf_id_i]
            assigned.add(conf_id_i)
            
            for j, conf_id_j in enumerate(conf_ids[i+1:], start=i+1):
                if conf_id_j in assigned:
                    continue
                
                rmsd = rdMolAlign.GetBestRMS(mol, mol, prbId=conf_id_i, refId=conf_id_j)
                
                if rmsd < threshold:
                    cluster.append(conf_id_j)
                    assigned.add(conf_id_j)
            
            clusters.append(cluster)
        
        return clusters
```

**Task 4.3: Create Conformer Generation Script** (2 hours)

Create `scripts/generate_[conformers.py](http://conformers.py)`:

```python
import sys
import os
from datetime import datetime

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'backend'))

from models.database import SessionLocal
from models.compound import Cannabinoid
from services.conformer_generator import ConformerGenerator

def generate_all_conformers():
    """Generate 3D conformers for all cannabinoids."""
    db = SessionLocal()
    generator = ConformerGenerator(num_conformers=50)
    
    cannabinoids = db.query(Cannabinoid).filter(
        Cannabinoid.has_conformers == False
    ).all()
    
    print(f"Generating conformers for {len(cannabinoids)} cannabinoids...")
    
    success_count = 0
    failed_count = 0
    
    for cb in cannabinoids:
        print(f"\nProcessing: {[cb.name](http://cb.name)} ({cb.abbreviation})")
        print(f"SMILES: {cb.smiles}")
        
        result = generator.generate_conformers(cb.smiles)
        
        if result.get("success"):
            # Update database
            cb.has_conformers = True
            cb.conformer_generation_method = "ETKDG"
            cb.num_conformers_generated = result["num_conformers_generated"]
            cb.conformer_metadata = {
                "lowest_energy_kcal_mol": result["lowest_energy"],
                "energy_range_kcal_mol": result["energy_range"],
                "num_clusters_rmsd_2A": result["num_clusters_rmsd_2A"]
            }
            cb.rdkit_descriptors_3d = result["descriptors_3d"]
            cb.conformers_generated_at = [datetime.now](http://datetime.now)()
            
            db.commit()
            
            print(f"✅ Success: {result['num_conformers_generated']} conformers generated")
            print(f"   Energy range: {result['energy_range']:.2f} kcal/mol")
            print(f"   RMSD clusters: {result['num_clusters_rmsd_2A']}")
            success_count += 1
        else:
            print(f"❌ Failed: {result.get('error')}")
            failed_count += 1
    
    db.close()
    
    print(f"\n{'='*60}")
    print(f"CONFORMER GENERATION COMPLETE")
    print(f"{'='*60}")
    print(f"✅ Success: {success_count}/{len(cannabinoids)}")
    print(f"❌ Failed: {failed_count}/{len(cannabinoids)}")

if __name__ == "__main__":
    generate_all_conformers()
```

**Task 4.4: Add Conformer API Endpoint** (1-2 hours)

Create `backend/api/routes/[conformers.py](http://conformers.py)`:

```python
from fastapi import APIRouter, Depends, HTTPException
from sqlalchemy.orm import Session
from models.database import get_db
from models.compound import Cannabinoid
from services.conformer_generator import ConformerGenerator

router = APIRouter(prefix="/api/conformers", tags=["conformers"])

@[router.post](http://router.post)("/generate/{compound_id}")
async def generate_conformers(compound_id: int, db: Session = Depends(get_db)):
    """Generate 3D conformers for a specific compound."""
    compound = db.query(Cannabinoid).filter([Cannabinoid.id](http://Cannabinoid.id) == compound_id).first()
    
    if not compound:
        raise HTTPException(status_code=404, detail="Compound not found")
    
    generator = ConformerGenerator()
    result = generator.generate_conformers(compound.smiles)
    
    if result.get("success"):
        # Update database
        compound.has_conformers = True
        compound.conformer_generation_method = "ETKDG"
        compound.num_conformers_generated = result["num_conformers_generated"]
        compound.conformer_metadata = result.get("conformer_metadata")
        compound.rdkit_descriptors_3d = result.get("descriptors_3d")
        db.commit()
        
        return {"message": "Conformers generated successfully", "data": result}
    else:
        raise HTTPException(status_code=500, detail=result.get("error"))
```

**Task 4.5: Run Conformer Generation** (2-3 hours)

```bash
# Generate conformers for all cannabinoids
python scripts/generate_[conformers.py](http://conformers.py)

# This will take 2-3 hours for 16 compounds (8-12 min each)
# Expected: 14-16 compounds succeed, 0-2 may fail
```

---

### Days 6-7: API Endpoints, Testing & Documentation (6-8 hours)

**Task 6.1: Study API Routes** (2 hours)

- `GET /api/studies/` - List studies
- `GET /api/studies/{study_id}` - Study details
- `GET /api/studies/?condition={condition}` - Filter by condition
- `GET /api/studies/stats/summary` - Statistics

**Task 6.2: Compound API Routes** (2 hours)

- `GET /api/compounds/` - List cannabinoids
- `GET /api/compounds/{compound_id}` - Details with 3D descriptors
- `GET /api/compounds/search?name={name}` - Search

**Task 6.3: Testing** (2 hours)

- Unit tests for models
- API endpoint tests
- Conformer generation tests

**Task 6.4: Documentation** (2 hours)

- Update README with 3D conformer features
- API documentation
- Week 1 completion report

---

### Week 1 Success Criteria Checklist

**Environment Setup**

- [ ]  Python 3.11+ virtual environment created
- [ ]  All dependencies installed (including RDKit)
- [ ]  PostgreSQL database running
- [ ]  VS Code configured

**Database**

- [ ]  SQLAlchemy models defined (studies, compounds, patients)
- [ ]  Alembic migrations applied
- [ ]  Database contains 200+ clinical studies
- [ ]  Database contains 16+ cannabinoids with 2D descriptors
- [ ]  **3D conformer schema implemented**
- [ ]  **14-16 cannabinoids have 3D conformer ensembles**

**Backend API**

- [ ]  FastAPI server running on [localhost:8000](http://localhost:8000)
- [ ]  Study endpoints functional
- [ ]  Compound endpoints functional
- [ ]  **Conformer generation endpoint functional**
- [ ]  OpenAPI documentation at /docs

**Testing**

- [ ]  Unit tests written
- [ ]  All tests passing
- [ ]  Manual API testing completed

**Documentation**

- [ ]  [README.md](http://README.md) complete
- [ ]  API documentation generated
- [ ]  Database schema documented

---

## Phase 1: Core Backend + OmniPath Integration (Weeks 2-5)

### Week 2: Database Expansion + OmniPath + Provenance

**Dates**: Dec 25-Jan 1, 2026

**Estimated Hours**: 28-32 hours

**Objectives**:

- Complete database schema for patient profiles and treatments
- Integrate OmniPath core services
- **NEW**: Implement provenance tracking infrastructure

**Deliverables**:

✅ Patient profile schema with demographics, conditions, medications

✅ Treatment recommendation schema

✅ OmniPath Consent Compiler integration

✅ OmniPath Token Validation API (<100ms response time)

✅ **NEW**: OmniPath Provenance Tracker for all database writes

✅ **NEW**: Audit trail infrastructure

**Key Tasks**:

**Task 2.1: Patient & Treatment Schema** (8-10 hours)

- Implement patient profile model
- Implement treatment recommendation model
- Add OmniPath consent fields
- Migration and testing

**Task 2.2: OmniPath Integration** (10-12 hours)

- Set up OmniPath API client
- Implement token validation middleware
- Implement consent checking
- Token-gate all API endpoints
- Response time optimization (<100ms)

**Task 2.3: NEW - Provenance Infrastructure** (8-10 hours)

Create `backend/services/provenance_[tracker.py](http://tracker.py)`:

```python
from typing import Dict, Any
from datetime import datetime
import hashlib
import json

class ProvenanceTracker:
    """Track data provenance for OmniPath manifest system."""
    
    def __init__(self, omnipath_client):
        self.omnipath_client = omnipath_client
    
    def create_manifest(self, data: Dict[str, Any], operation: str) -> str:
        """Create provenance manifest for data operation.
        
        Args:
            data: The data being written/modified
            operation: "CREATE", "UPDATE", "DELETE"
            
        Returns:
            manifest_id: Blockchain reference ID
        """
        manifest = {
            "timestamp": [datetime.now](http://datetime.now)().isoformat(),
            "operation": operation,
            "data_hash": self._hash_data(data),
            "data_type": data.get("__type__", "unknown"),
            "source": "neurobotanica_mvp",
            "version": "1.0"
        }
        
        # Send to OmniPath Manifest Infrastructure
        response = self.omnipath_client.manifests.create(manifest)
        
        return response["manifest_id"]
    
    def _hash_data(self, data: Dict[str, Any]) -> str:
        """Generate cryptographic hash of data."""
        data_str = json.dumps(data, sort_keys=True)
        return hashlib.sha256(data_str.encode()).hexdigest()
```

Add provenance fields to all tables:

- `omnipath_manifest_id` (String)
- `created_at` (DateTime)
- `updated_at` (DateTime)
- `created_by_user` (String)

**Success Metrics**:

- Patient CRUD operations functional
- OmniPath token validation <100ms
- All database writes logged to manifests
- Audit trail queryable

---

### Week 3: Clinical Evidence API + Enhanced Receptor Affinity Provenance

**Dates**: Jan 1-8, 2026

**Estimated Hours**: 28-32 hours

**Objectives**:

- Build clinical evidence query APIs
- **NEW**: Enhance receptor affinity schema with assay provenance
- Implement evidence aggregation logic

**Deliverables**:

✅ Evidence query API endpoints

✅ Study filtering and search

✅ **NEW**: Receptor affinity provenance schema

✅ **NEW**: Assay heterogeneity detection

✅ Evidence aggregation by condition

✅ Confidence weighting system

**Key Tasks**:

**Task 3.1: Evidence API Endpoints** (10-12 hours)

```python
# GET /api/evidence/cannabinoid/{compound}
# GET /api/evidence/condition/{condition}
# GET /api/evidence/aggregate?condition={condition}
```

Implement:

- Study retrieval by compound
- Study retrieval by condition
- Evidence aggregation with confidence weighting
- Citation formatting

**Task 3.2: NEW - Enhanced Receptor Affinity Schema** (10-12 hours)

Update cannabinoid schema from:

```python
receptor_affinities = Column(JSON)  # {"CB1": 0.5, "CB2": 0.8}
```

To provenance-rich structure:

```python
receptor_affinities = Column(JSON)  # Full provenance

# New structure:
{
  "CB1": [
    {
      "affinity_value": 5.8,
      "affinity_unit": "pKi",
      "assay_type": "radioligand_binding",
      "target_confidence": "high",
      "source": "ChEMBL:CHEMBL123456",
      "pubmed_id": "12345678",
      "assay_organism": "Homo sapiens",
      "assay_cell_type": "CHO cells"
    },
    # Multiple measurements from different sources
  ],
  "CB2": [...]
}
```

Create migration script:

```bash
python scripts/migrate_receptor_[affinities.py](http://affinities.py)
```

**Task 3.3: Assay Heterogeneity Detection** (6-8 hours)

Create `backend/services/assay_[analyzer.py](http://analyzer.py)`:

```python
class AssayAnalyzer:
    """Detect heterogeneity in receptor binding assays."""
    
    def analyze_receptor_data(self, affinities: Dict) -> Dict:
        """Analyze receptor affinity data for a compound.
        
        Returns warnings about assay heterogeneity.
        """
        warnings = []
        
        for receptor, measurements in affinities.items():
            if len(measurements) > 1:
                # Check for unit consistency
                units = set(m["affinity_unit"] for m in measurements)
                if len(units) > 1:
                    warnings.append({
                        "receptor": receptor,
                        "type": "mixed_units",
                        "message": f"Mixed units detected: {units}"
                    })
                
                # Check for assay type consistency
                assay_types = set(m["assay_type"] for m in measurements)
                if len(assay_types) > 1:
                    warnings.append({
                        "receptor": receptor,
                        "type": "mixed_assay_types",
                        "message": f"Different assay types: {assay_types}"
                    })
                
                # Check for large variance
                values = [m["affinity_value"] for m in measurements]
                if max(values) - min(values) > 2.0:  # >100-fold difference
                    warnings.append({
                        "receptor": receptor,
                        "type": "high_variance",
                        "message": f"High variance: {min(values):.1f} - {max(values):.1f}"
                    })
        
        return {"warnings": warnings, "has_issues": len(warnings) > 0}
```

Add to API response:

```python
@router.get("/api/compounds/{id}/receptor-affinities")
async def get_receptor_affinities(
    id: int,
    include_provenance: bool = True,
    db: Session = Depends(get_db)
):
    compound = db.query(Cannabinoid).filter([Cannabinoid.id](http://Cannabinoid.id) == id).first()
    
    affinities = compound.receptor_affinities
    
    if include_provenance:
        analyzer = AssayAnalyzer()
        analysis = analyzer.analyze_receptor_data(affinities)
        return {"affinities": affinities, "analysis": analysis}
    else:
        # Return simplified view
        return {"affinities": {r: m[0]["affinity_value"] for r, m in affinities.items()}}
```

**Success Metrics**:

- Evidence API functional
- Receptor affinity provenance complete
- Assay heterogeneity warnings display correctly
- Confidence weighting implemented

---

### Week 4: Dimer Prediction + 3D Conformers + Triangulation

**Dates**: Jan 8-15, 2026

**Estimated Hours**: 32-38 hours

**Objectives**:

- Implement dimeric cannabinoid prediction algorithms
- **NEW**: Extend 3D conformers to predicted dimers
- **NEW**: Implement triangulation scoring framework

**Deliverables**:

✅ Dimer structure generation

✅ Formation probability scoring

✅ Synergy prediction algorithms

✅ **NEW**: 3D conformers for top dimer predictions

✅ **NEW**: Triangulation scoring system

✅ **NEW**: Confidence bounds on predictions

**Key Tasks**:

**Task 4.1: Dimer Prediction Logic** (10-12 hours)

Create `backend/services/dimer_[predictor.py](http://predictor.py)`:

```python
class DimericPredictor:
    """Predict dimeric cannabinoid structures and properties."""
    
    def predict_homodimer(self, monomer_smiles: str) -> Dict:
        """Predict homodimer formation (A-A)."""
        # Identify reactive sites on monomer
        # Generate dimer SMILES with methylene bridge
        # Calculate formation probability
        # Predict therapeutic properties
        pass
    
    def predict_heterodimer(self, monomer1_smiles: str, monomer2_smiles: str) -> Dict:
        """Predict heterodimer formation (A-B)."""
        pass
```

**Task 4.2: NEW - 3D Conformers for Dimers** (10-12 hours)

Extend `ConformerGenerator` to handle dimeric structures:

```python
class DimericConformerGenerator(ConformerGenerator):
    """Generate 3D conformers for dimeric cannabinoids."""
    
    def generate_dimer_conformers(self, dimer_smiles: str) -> Dict:
        """Generate conformers for dimeric structure.
        
        Additional considerations:
        - Intermonomer distance
        - Relative orientation
        - Structural stability
        """
        result = super().generate_conformers(dimer_smiles)
        
        if result.get("success"):
            # Add dimer-specific analysis
            result["dimer_analysis"] = self._analyze_dimer_geometry(mol, conf_id)
        
        return result
    
    def _analyze_dimer_geometry(self, mol, conf_id) -> Dict:
        """Analyze dimer-specific geometry."""
        # Calculate intermonomer distance
        # Measure relative orientation
        # Assess structural stability
        return {
            "intermonomer_distance_A": 3.5,
            "relative_orientation_deg": 120,
            "stability_score": 0.85
        }
```

**Task 4.3: NEW - Triangulation Scoring Framework** (10-12 hours)

Create `backend/services/triangulation_[scorer.py](http://scorer.py)`:

```python
from typing import Dict, Optional
import numpy as np

class DimericTriangulationScorer:
    """Triangulate dimeric predictions using multiple validation methods.
    
    Three validation sources:
    1. Computational prediction (formation probability)
    2. Structural similarity to known dimers
    3. ML-predicted therapeutic potential
    """
    
    def calculate_triangulation_score(
        self,
        formation_probability: float,
        structural_similarity: Optional[float],
        ml_therapeutic_score: Optional[float],
        experimental_validation: Optional[str] = None
    ) -> Dict:
        """Calculate triangulated confidence score.
        
        Args:
            formation_probability: 0.0-1.0 from quantum/thermodynamic model
            structural_similarity: 0.0-1.0 Tanimoto to known dimers
            ml_therapeutic_score: 0.0-1.0 from ML model
            experimental_validation: "confirmed", "tentative", "none"
            
        Returns:
            dict with triangulation score and confidence bounds
        """
        scores = []
        weights = []
        
        # Source 1: Formation probability (computational)
        if formation_probability is not None:
            scores.append(formation_probability)
            weights.append(0.4)  # 40% weight
        
        # Source 2: Structural similarity (if known dimers exist)
        if structural_similarity is not None:
            scores.append(structural_similarity)
            weights.append(0.3)  # 30% weight
        
        # Source 3: ML therapeutic score
        if ml_therapeutic_score is not None:
            scores.append(ml_therapeutic_score)
            weights.append(0.3)  # 30% weight
        
        # Normalize weights
        weights = np.array(weights) / sum(weights)
        scores = np.array(scores)
        
        # Calculate weighted average
        triangulation_score = np.sum(scores * weights)
        
        # Calculate uncertainty (std dev of scores)
        uncertainty = np.std(scores) if len(scores) > 1 else 0.0
        
        # Adjust for experimental validation
        if experimental_validation == "confirmed":
            triangulation_score = min(1.0, triangulation_score * 1.2)
            uncertainty *= 0.5
        elif experimental_validation == "tentative":
            uncertainty *= 1.2
        
        # Confidence interval (95%)
        confidence_lower = max(0.0, triangulation_score - 1.96 * uncertainty)
        confidence_upper = min(1.0, triangulation_score + 1.96 * uncertainty)
        
        return {
            "triangulation_score": float(triangulation_score),
            "uncertainty": float(uncertainty),
            "confidence_interval": [float(confidence_lower), float(confidence_upper)],
            "num_sources": len(scores),
            "sources_used": [
                "formation_probability" if formation_probability is not None else None,
                "structural_similarity" if structural_similarity is not None else None,
                "ml_therapeutic_score" if ml_therapeutic_score is not None else None
            ],
            "experimental_status": experimental_validation or "none"
        }
```

Update dimeric cannabinoid schema:

```python
# Add to dimeric_cannabinoids table
triangulation_score = Column(Float)
triangulation_components = Column(JSON)
uncertainty = Column(Float)
confidence_interval_lower = Column(Float)
confidence_interval_upper = Column(Float)
experimental_validation_status = Column(String(50))
```

**Task 4.4: Generate Dimeric Predictions** (8-10 hours)

```bash
# Generate 120 dimeric predictions
python scripts/generate_dimeric_[predictions.py](http://predictions.py)

# Calculate triangulation scores
python scripts/calculate_triangulation_[scores.py](http://scores.py)

# Generate 3D conformers for top 20 dimers
python scripts/generate_dimer_[conformers.py](http://conformers.py) --top 20
```

**Success Metrics**:

- 120 dimeric predictions generated
- Each has triangulation score with confidence bounds
- Top 20 dimers have 3D conformer ensembles
- API endpoint `/api/dimers/{id}/triangulation` functional

---

### Week 5: Adjuvant Optimization + ChEMBL/PubChem Integration

**Dates**: Jan 15-22, 2026

**Estimated Hours**: 32-38 hours

**Objectives**:

- Build adjuvant optimization module
- **NEW**: Integrate ChEMBL and PubChem databases
- **NEW**: Preserve assay context from external sources

**Deliverables**:

✅ Adjuvant database (magnesium, L-theanine, etc.)

✅ Priming protocol generator

✅ Timing optimization algorithms

✅ Drug interaction checker

✅ **NEW**: ChEMBL integration service

✅ **NEW**: PubChem integration service

✅ **NEW**: Daily data sync automation

**Key Tasks**:

**Task 5.1: Adjuvant Module** (12-14 hours)

Create adjuvant database schema and protocols:

```python
class Adjuvant(Base):
    __tablename__ = "adjuvants"
    
    id = Column(Integer, primary_key=True)
    name = Column(String(200))
    category = Column(String(100))  # "Receptor Priming", "Metabolic", etc.
    mechanism = Column(Text)
    recommended_dose_mg = Column(Float)
    timing_offset_minutes = Column(Integer)  # Before cannabinoid
    evidence_tier = Column(Integer)  # 1-5
```

Implement priming protocol generator:

```python
class PrimingProtocolGenerator:
    """Generate adjuvant priming protocols."""
    
    def generate_protocol(
        self,
        cannabinoid_id: int,
        condition: str,
        patient_profile: Dict
    ) -> Dict:
        """Generate optimized adjuvant protocol.
        
        Returns:
            Adjuvant recommendations with timing and dosing
        """
        # Analyze cannabinoid mechanism
        # Select compatible adjuvants
        # Optimize timing
        # Calculate expected enhancement
        pass
```

**Task 5.2: NEW - ChEMBL Integration** (10-12 hours)

Create `backend/services/chembl_[client.py](http://client.py)`:

```python
import requests
from typing import Dict, List, Optional

class ChEMBLClient:
    """Client for ChEMBL database API."""
    
    BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"
    
    def __init__(self):
        self.session = requests.Session()
        self.session.headers.update({"Accept": "application/json"})
    
    def search_by_smiles(self, smiles: str) -> Optional[str]:
        """Search ChEMBL by SMILES, return ChEMBL ID."""
        response = self.session.get(
            f"{self.BASE_URL}/molecule",
            params={"molecule_structures__canonical_smiles__flexmatch": smiles}
        )
        
        if response.status_code == 200:
            data = response.json()
            if data.get("molecules"):
                return data["molecules"][0]["molecule_chembl_id"]
        return None
    
    def get_receptor_affinities(
        self,
        chembl_id: str,
        target_name: str = None
    ) -> List[Dict]:
        """Get receptor binding data with full assay context.
        
        Args:
            chembl_id: ChEMBL molecule ID (e.g., "CHEMBL123")
            target_name: Optional filter (e.g., "CB1", "CB2")
            
        Returns:
            List of binding measurements with assay metadata
        """
        response = self.session.get(
            f"{self.BASE_URL}/activity",
            params={
                "molecule_chembl_id": chembl_id,
                "assay_type": "B",  # Binding assay
                "limit": 1000
            }
        )
        
        if response.status_code != 200:
            return []
        
        activities = response.json().get("activities", [])
        
        # Parse and structure with provenance
        structured_data = []
        for activity in activities:
            # Filter by target if specified
            if target_name:
                target = activity.get("target_pref_name", "")
                if target_name.lower() not in target.lower():
                    continue
            
            structured_data.append({
                "affinity_value": activity.get("standard_value"),
                "affinity_unit": activity.get("standard_units"),
                "assay_type": activity.get("assay_type"),
                "assay_description": activity.get("assay_description"),
                "target_name": activity.get("target_pref_name"),
                "target_organism": activity.get("target_organism"),
                "source": f"ChEMBL:{chembl_id}",
                "activity_id": activity.get("activity_id"),
                "pubmed_id": activity.get("document_journal"),
                "confidence_score": activity.get("confidence_score")
            })
        
        return structured_data
```

**Task 5.3: NEW - PubChem Integration** (8-10 hours)

Create `backend/services/pubchem_[client.py](http://client.py)`:

```python
import requests
from typing import Dict, Optional

class PubChemClient:
    """Client for PubChem database API."""
    
    BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    
    def search_by_smiles(self, smiles: str) -> Optional[str]:
        """Search PubChem by SMILES, return CID."""
        response = requests.get(
            f"{self.BASE_URL}/compound/smiles/{smiles}/cids/JSON"
        )
        
        if response.status_code == 200:
            data = response.json()
            if data.get("IdentifierList", {}).get("CID"):
                return str(data["IdentifierList"]["CID"][0])
        return None
    
    def get_compound_properties(self, cid: str) -> Dict:
        """Get compound properties from PubChem."""
        response = requests.get(
            f"{self.BASE_URL}/compound/cid/{cid}/property/MolecularWeight,XLogP,TPSA/JSON"
        )
        
        if response.status_code == 200:
            data = response.json()
            props = data.get("PropertyTable", {}).get("Properties", [{}])[0]
            return props
        return {}
    
    def get_bioassay_data(self, cid: str) -> List[Dict]:
        """Get bioassay data for compound."""
        # PubChem BioAssay API
        response = requests.get(
            f"{self.BASE_URL}/compound/cid/{cid}/assaysummary/JSON"
        )
        
        if response.status_code == 200:
            return response.json().get("Table", {}).get("Row", [])
        return []
```

**Task 5.4: Data Sync Script** (6-8 hours)

Create `scripts/sync_external_[databases.py](http://databases.py)`:

```python
import sys
import os
from datetime import datetime

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'backend'))

from models.database import SessionLocal
from models.compound import Cannabinoid
from services.chembl_client import ChEMBLClient
from services.pubchem_client import PubChemClient

def sync_all_cannabinoids():
    """Sync all cannabinoids with ChEMBL and PubChem."""
    db = SessionLocal()
    chembl = ChEMBLClient()
    pubchem = PubChemClient()
    
    cannabinoids = db.query(Cannabinoid).all()
    
    for cb in cannabinoids:
        print(f"\nSyncing: {[cb.name](http://cb.name)}")
        
        # ChEMBL lookup
        if not cb.chembl_id:
            chembl_id = [chembl.search](http://chembl.search)_by_smiles(cb.smiles)
            if chembl_id:
                cb.chembl_id = chembl_id
                print(f"  Found ChEMBL ID: {chembl_id}")
        
        # Get receptor affinities from ChEMBL
        if cb.chembl_id:
            cb1_data = chembl.get_receptor_affinities(cb.chembl_id, "CB1")
            cb2_data = chembl.get_receptor_affinities(cb.chembl_id, "CB2")
            
            # Merge with existing data
            if cb.receptor_affinities is None:
                cb.receptor_affinities = {}
            
            if cb1_data:
                cb.receptor_affinities["CB1"] = cb1_data
                print(f"  Added {len(cb1_data)} CB1 measurements")
            
            if cb2_data:
                cb.receptor_affinities["CB2"] = cb2_data
                print(f"  Added {len(cb2_data)} CB2 measurements")
        
        # PubChem lookup
        if not cb.pubchem_cid:
            pubchem_cid = [pubchem.search](http://pubchem.search)_by_smiles(cb.smiles)
            if pubchem_cid:
                cb.pubchem_cid = pubchem_cid
                print(f"  Found PubChem CID: {pubchem_cid}")
        
        # Update sync timestamp
        cb.last_external_sync = [datetime.now](http://datetime.now)()
        
        db.commit()
    
    db.close()
    print("\nSync complete!")

if __name__ == "__main__":
    sync_all_cannabinoids()
```

**Task 5.5: Schedule Daily Sync** (2 hours)

Set up cron job or scheduled task:

```bash
# Add to crontab (daily at 3 AM)
0 3 * * * /path/to/venv/bin/python /path/to/scripts/sync_external_[databases.py](http://databases.py)
```

**Success Metrics**:

- Adjuvant module functional with protocols
- ChEMBL/PubChem integration services operational
- 16 cannabinoids have external database IDs
- Receptor affinities enriched with ChEMBL data
- Daily sync configured

---

## Phase 2: ML Models + PatentPath Lite + Frontend (Weeks 6-9)

### Week 6: ML Model Training + PatentPath Lite Foundation

**Dates**: Jan 22-29, 2026

**Estimated Hours**: 35-40 hours

**Objectives**:

- Train therapeutic prediction ML models
- **NEW**: Begin PatentPath Lite development (parallel with ML)
- Integrate 3D descriptors into ML features

**Deliverables**:

✅ Model 1: Condition → Optimal cannabinoid profile

✅ Model 2: Patient → Treatment response probability

✅ Model 3: Dimer → Therapeutic potential

✅ **NEW**: PatentPath Lite - Prior Art Search (Feature 1)

✅ **NEW**: PatentPath Lite - Basic Novelty Scoring (Feature 2)

✅ Uncertainty quantification for all models

✅ Model registry setup

**Split:** 20 hours ML + 15 hours PatentPath Lite

### ML Model Training (20 hours)

**Task 6.1: Enhanced Feature Engineering** (6-8 hours)

Create `backend/services/ml_data_[prep.py](http://prep.py)`:

```python
import pandas as pd
import numpy as np
from typing import List, Dict

class MLDataPreparator:
    """Prepare training data with 2D and 3D features."""
    
    def prepare_training_dataset(self, cannabinoids: List, studies: List) -> pd.DataFrame:
        """Prepare ML training dataset.
        
        Features:
        - 2D RDKit descriptors (40+)
        - 3D descriptors (PMI, asphericity, etc.) if available
        - Receptor affinity data
        - Clinical evidence features
        - Confidence weights
        """
        data = []
        
        for cb in cannabinoids:
            # Get relevant studies
            cb_studies = [s for s in studies if cb.abbreviation in s.cannabinoids_studied]
            
            for study in cb_studies:
                features = {}
                
                # 2D descriptors
                features.update(cb.rdkit_descriptors)
                
                # 3D descriptors (if available)
                if cb.has_conformers and cb.rdkit_descriptors_3d:
                    features.update(cb.rdkit_descriptors_3d)
                    features['has_3d'] = 1
                else:
                    # Pad with zeros if no 3D data
                    features['has_3d'] = 0
                
                # Receptor affinities (averaged)
                if cb.receptor_affinities:
                    for receptor, measurements in cb.receptor_affinities.items():
                        avg_affinity = np.mean([m['affinity_value'] for m in measurements])
                        features[f'{receptor}_affinity'] = avg_affinity
                
                # Study outcomes
                features['condition'] = study.condition
                features['efficacy'] = self._parse_efficacy(study.efficacy_summary)
                features['confidence_weight'] = study.confidence_weight
                
                data.append(features)
        
        return pd.DataFrame(data)
    
    def _parse_efficacy(self, efficacy_summary: str) -> float:
        """Extract numeric efficacy score from text summary."""
        # Simple heuristic: look for improvement percentages
        # In production, use more sophisticated NLP
        if 'significant' in efficacy_summary.lower():
            return 0.8
        elif 'moderate' in efficacy_summary.lower():
            return 0.6
        elif 'mild' in efficacy_summary.lower():
            return 0.4
        else:
            return 0.5
```

**Task 6.2: Model Training** (10-12 hours)

Create `backend/services/ml_[models.py](http://models.py)`:

```python
from sklearn.ensemble import GradientBoostingRegressor, GradientBoostingClassifier
from sklearn.model_selection import cross_val_score, train_test_split
import numpy as np
import joblib

class TherapeuticPredictionModel:
    """Gradient boosting model for therapeutic predictions."""
    
    def __init__(self):
        self.model = GradientBoostingRegressor(
            n_estimators=100,
            learning_rate=0.1,
            max_depth=5,
            random_state=42
        )
    
    def train(self, X, y, sample_weights=None):
        """Train model with optional confidence weights."""
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.2, random_state=42
        )
        
        if sample_weights is not None:
            w_train, w_test = train_test_split(
                sample_weights, test_size=0.2, random_state=42
            )
            [self.model.fit](http://self.model.fit)(X_train, y_train, sample_weight=w_train)
        else:
            [self.model.fit](http://self.model.fit)(X_train, y_train)
        
        # Evaluate
        train_score = self.model.score(X_train, y_train)
        test_score = self.model.score(X_test, y_test)
        
        return {
            'train_r2': train_score,
            'test_r2': test_score,
            'cv_scores': cross_val_score(self.model, X, y, cv=5)
        }
    
    def predict_with_uncertainty(self, X):
        """Predict with uncertainty quantification."""
        # Use individual estimators for uncertainty
        predictions = np.array([
            estimator[0].predict(X) 
            for estimator in self.model.estimators_
        ])
        
        mean_pred = predictions.mean(axis=0)
        std_pred = predictions.std(axis=0)
        
        return {
            'prediction': mean_pred,
            'uncertainty': std_pred,
            'confidence_interval_lower': mean_pred - 1.96 * std_pred,
            'confidence_interval_upper': mean_pred + 1.96 * std_pred
        }
    
    def save(self, path: str):
        """Save model to disk."""
        joblib.dump(self.model, path)
    
    def load(self, path: str):
        """Load model from disk."""
        self.model = joblib.load(path)
```

Train all 3 models:

```bash
python scripts/train_ml_[models.py](http://models.py)
```

**Task 6.3: Model Registry** (2 hours)

Set up Weights & Biases or similar:

```python
import wandb

wandb.init(project="neurobotanica-mvp")

# Log model metrics
wandb.log({
    "model": "therapeutic_prediction_v1.0",
    "train_r2": 0.72,
    "test_r2": 0.68,
    "features_used": 45,
    "includes_3d": True
})

# Save model artifact
[wandb.save](http://wandb.save)("models/therapeutic_prediction_v1.0.pkl")
```

### PatentPath Lite Development - Phase 1 (15 hours)

**Task 6.4: PatentPath Lite Repository Setup** (2 hours)

```bash
mkdir -p backend/services/patentpath_lite
cd backend/services/patentpath_lite

# Create module structure
touch __init__.py
touch prior_art_[search.py](http://search.py)
touch novelty_[scorer.py](http://scorer.py)
touch fto_[checker.py](http://checker.py)
touch claim_[generator.py](http://generator.py)
touch tk_[checker.py](http://checker.py)
touch cost_[estimator.py](http://estimator.py)
```

Create `backend/services/patentpath_lite/requirements.txt`:

```
requests==2.31.0
beautifulsoup4==4.12.2
lxml==4.9.3
jsonschema==4.20.0
```

**Task 6.5: Feature 1 - Prior Art Search** (8-10 hours)

<aside>
⚠️

**MVP Risk Note**: Google Patents scraping (mentioned below) is brittle and may trigger rate-limits or compliance issues. Treat as optional fallback. Prefer USPTO PatentsView API and EPO Open Patent Services (both free) for MVP. Plan to add licensed sources ([Lens.org](http://Lens.org) API, PatSnap) post-MVP.

</aside>

Create `backend/services/patentpath_lite/prior_art_[search.py](http://search.py)`:

```python
import requests
from typing import List, Dict
from bs4 import BeautifulSoup
import time

class PriorArtSearcher:
    """Basic prior art search using free patent databases.
    
    MVP version uses:
    - USPTO PatentsView API (free)
    - EPO Open Patent Services (free tier)
    - Google Patents (web scraping)
    """
    
    def __init__(self):
        self.uspto_base = "https://api.patentsview.org/patents/query"
        self.epo_base = "https://ops.epo.org/3.2/rest-services"
        self.session = requests.Session()
    
    def search(
        self,
        compound_description: str,
        keywords: List[str],
        databases: List[str] = ["USPTO", "EPO"]
    ) -> Dict:
        """Search for prior art patents.
        
        Args:
            compound_description: Text description of compound
            keywords: List of search keywords
            databases: Which databases to search
            
        Returns:
            Dict with prior art patents found
        """
        all_patents = []
        
        if "USPTO" in databases:
            uspto_results = self._search_uspto(keywords)
            all_patents.extend(uspto_results)
        
        if "EPO" in databases:
            epo_results = self._search_epo(keywords)
            all_patents.extend(epo_results)
        
        if "Google Patents" in databases:
            google_results = self._search_google_patents(keywords)
            all_patents.extend(google_results)
        
        # Deduplicate and rank by relevance
        unique_patents = self._deduplicate(all_patents)
        ranked_patents = self._rank_by_relevance(unique_patents, compound_description, keywords)
        
        return {
            "prior_art_patents": ranked_patents[:20],  # Top 20
            "total_results": len(unique_patents),
            "highest_relevance": ranked_patents[0]["relevance_score"] if ranked_patents else 0
        }
    
    def _search_uspto(self, keywords: List[str]) -> List[Dict]:
        """Search USPTO PatentsView API."""
        query = " OR ".join([f'_text_any_patent_title:("{kw}")' for kw in keywords])
        
        params = {
            "q": query,
            "f": ["patent_number", "patent_title", "patent_date", "patent_abstract"],
            "o": {"per_page": 100}
        }
        
        try:
            response = [self.session.post](http://self.session.post)(self.uspto_base, json=params)
            if response.status_code == 200:
                data = response.json()
                patents = data.get("patents", [])
                
                return [
                    {
                        "patent_number": p.get("patent_number"),
                        "title": p.get("patent_title"),
                        "date": p.get("patent_date"),
                        "abstract": p.get("patent_abstract"),
                        "jurisdiction": "US",
                        "status": self._check_patent_status(p.get("patent_number")),
                        "source": "USPTO"
                    }
                    for p in patents
                ]
        except Exception as e:
            print(f"USPTO search error: {e}")
        
        return []
    
    def _search_epo(self, keywords: List[str]) -> List[Dict]:
        """Search EPO Open Patent Services.
        
        Note: Requires API key for production use.
        MVP: Rate-limited free tier.
        """
        # Simplified for MVP
        # Full implementation requires EPO API authentication
        return []
    
    def _search_google_patents(self, keywords: List[str]) -> List[Dict]:
        """Search Google Patents via web scraping.
        
        Note: For MVP only. Production should use official APIs.
        """
        query = "+".join(keywords)
        url = f"https://patents.google.com/?q={query}"
        
        try:
            response = self.session.get(url)
            soup = BeautifulSoup(response.text, 'html.parser')
            
            # Parse results (simplified)
            # In production, use proper pagination and parsing
            results = []
            # ... parsing logic ...
            
            return results
        except Exception as e:
            print(f"Google Patents search error: {e}")
        
        return []
    
    def _check_patent_status(self, patent_number: str) -> str:
        """Check if patent is active, expired, or abandoned."""
        # Simplified for MVP
        # In production, query USPTO PAIR system
        return "active"  # Default assumption
    
    def _deduplicate(self, patents: List[Dict]) -> List[Dict]:
        """Remove duplicate patents across databases."""
        seen = set()
        unique = []
        
        for p in patents:
            # Normalize patent number for comparison
            key = p["patent_number"].replace("-", "").replace(",", "")
            if key not in seen:
                seen.add(key)
                unique.append(p)
        
        return unique
    
    def _rank_by_relevance(self, patents: List[Dict], description: str, keywords: List[str]) -> List[Dict]:
        """Rank patents by relevance using simple TF-IDF."""
        from sklearn.feature_extraction.text import TfidfVectorizer
        from sklearn.metrics.pairwise import cosine_similarity
        
        if not patents:
            return []
        
        # Combine description and keywords
        query_text = description + " " + " ".join(keywords)
        
        # Get patent texts
        patent_texts = [
            (p.get("title", "") + " " + p.get("abstract", "")) 
            for p in patents
        ]
        
        # TF-IDF vectorization
        vectorizer = TfidfVectorizer()
        all_texts = [query_text] + patent_texts
        tfidf_matrix = [vectorizer.fit](http://vectorizer.fit)_transform(all_texts)
        
        # Calculate similarity scores
        query_vec = tfidf_matrix[0:1]
        patent_vecs = tfidf_matrix[1:]
        similarities = cosine_similarity(query_vec, patent_vecs)[0]
        
        # Add scores to patents
        for i, patent in enumerate(patents):
            patent["relevance_score"] = float(similarities[i])
        
        # Sort by relevance
        patents.sort(key=lambda p: p["relevance_score"], reverse=True)
        
        return patents
```

**Task 6.6: Feature 2 - Basic Novelty Scoring** (3-4 hours)

Create `backend/services/patentpath_lite/novelty_[scorer.py](http://scorer.py)`:

```python
from typing import Dict, List

class NoveltyScorer:
    """Calculate novelty score based on prior art similarity.
    
    MVP version uses simple inverse similarity scoring.
    Full version (PatentPath Full) uses ML-powered scoring.
    """
    
    def calculate_novelty(
        self,
        compound_description: str,
        prior_art_patents: List[Dict]
    ) -> Dict:
        """Calculate novelty score (0.0-1.0).
        
        Args:
            compound_description: Description of new compound
            prior_art_patents: List from PriorArtSearcher
            
        Returns:
            Dict with novelty score and recommendation
        """
        if not prior_art_patents:
            return {
                "novelty_score": 1.0,
                "recommendation": "LIKELY_PATENTABLE",
                "confidence": "high",
                "closest_prior_art": None,
                "reason": "No prior art found"
            }
        
        # Find most similar patent
        closest = max(prior_art_patents, key=lambda p: p.get("relevance_score", 0))
        highest_similarity = closest.get("relevance_score", 0)
        
        # Novelty is inverse of similarity
        novelty_score = 1.0 - highest_similarity
        
        # Classification
        if novelty_score > 0.7:
            recommendation = "LIKELY_PATENTABLE"
            confidence = "high"
            reason = "Significantly different from prior art"
        elif novelty_score > 0.4:
            recommendation = "UNCERTAIN_CONSULT_ATTORNEY"
            confidence = "medium"
            reason = "Some similarity to prior art - attorney review recommended"
        else:
            recommendation = "LIKELY_NOT_PATENTABLE"
            confidence = "high"
            reason = "High similarity to existing patents"
        
        # Identify key differences
        differences = self._identify_differences(
            compound_description,
            closest.get("title", "") + " " + closest.get("abstract", "")
        )
        
        return {
            "novelty_score": novelty_score,
            "recommendation": recommendation,
            "confidence": confidence,
            "closest_prior_art": {
                "patent_number": closest.get("patent_number"),
                "similarity": highest_similarity,
                "key_differences": differences
            },
            "reason": reason,
            "next_steps": self._generate_next_steps(recommendation)
        }
    
    def _identify_differences(self, compound_desc: str, prior_art_text: str) -> List[str]:
        """Identify key differences between compound and prior art.
        
        MVP version uses simple keyword comparison.
        Full version uses NLP semantic analysis.
        """
        differences = []
        
        # Extract key chemical terms
        compound_keywords = set(compound_desc.lower().split())
        prior_art_keywords = set(prior_art_text.lower().split())
        
        # Look for unique cannabinoid types
        cannabinoids = ["cbd", "thc", "cbg", "cbn", "cbc", "thcv", "cbdv"]
        compound_cbs = [cb for cb in cannabinoids if cb in compound_desc.lower()]
        prior_art_cbs = [cb for cb in cannabinoids if cb in prior_art_text.lower()]
        
        if set(compound_cbs) != set(prior_art_cbs):
            differences.append(
                f"Different cannabinoid basis: Your compound uses {', '.join(compound_cbs)}, "
                f"prior art uses {', '.join(prior_art_cbs)}"
            )
        
        # Look for dimer vs monomer
        if "dimer" in compound_desc.lower() and "dimer" not in prior_art_text.lower():
            differences.append("Your compound is dimeric, prior art is monomeric")
        
        # Look for bridge type
        if "methylene" in compound_desc.lower() and "methylene" not in prior_art_text.lower():
            differences.append("Different bridge structure: methylene vs other")
        
        return differences if differences else ["Detailed structural comparison recommended"]
    
    def _generate_next_steps(self, recommendation: str) -> List[str]:
        """Generate actionable next steps based on recommendation."""
        if recommendation == "LIKELY_PATENTABLE":
            return [
                "Consider filing provisional patent application",
                "Conduct full freedom-to-operate analysis",
                "Consult patent attorney for claim drafting"
            ]
        elif recommendation == "UNCERTAIN_CONSULT_ATTORNEY":
            return [
                "Consult patent attorney for professional assessment",
                "Consider design-around strategies",
                "Conduct detailed claim-by-claim analysis"
            ]
        else:  # LIKELY_NOT_PATENTABLE
            return [
                "Review prior art for potential invalidity challenges",
                "Consider trade secret protection instead",
                "Explore design-around opportunities"
            ]
```

**Success Metrics - Week 6**:

- 3 ML models trained with 3D descriptors
- Cross-validation R² > 0.65 for therapeutic prediction
- Uncertainty quantification functional
- PatentPath Lite Features 1-2 implemented
- Prior art search returns top 20 patents
- Novelty scoring classifies compounds accurately

---

### Week 7: Appetite Model + PatentPath Lite (Features 3-4) + Terpene Gating

**Dates**: Jan 29-Feb 5, 2026

**Estimated Hours**: 35-40 hours

**Objectives**:

- Train appetite stimulation ML model
- **NEW**: Complete PatentPath Lite Features 3-4 (FTO + Claims)
- **NEW**: Implement evidence-weighted terpene features

**Deliverables**:

✅ Appetite stimulation prediction model

✅ Side effect prediction (munchies, dry mouth, etc.)

✅ **NEW**: PatentPath Lite - FTO Checker (Feature 3)

✅ **NEW**: PatentPath Lite - Claim Generator (Feature 4)

✅ **NEW**: Evidence-gated terpene interaction scoring

✅ Model versioning and tracking

**Split:** 15 hours ML/Terpenes + 20 hours PatentPath Lite

---

### ML Model Training - Appetite Stimulation (10 hours)

**Task 7.1: Appetite Stimulation Model** (8-10 hours)

Create `backend/services/appetite_[model.py](http://model.py)`:

```python
from sklearn.ensemble import GradientBoostingClassifier
import pandas as pd
from typing import Dict, List

class AppetiteStimulationModel:
    """Predict appetite stimulation potential and side effects."""
    
    def __init__(self):
        self.model = GradientBoostingClassifier(
            n_estimators=100,
            learning_rate=0.1,
            max_depth=4,
            random_state=42
        )
    
    def prepare_features(self, compound: Dict) -> pd.DataFrame:
        """Prepare features for appetite prediction.
        
        Key predictors:
        - CB1 receptor affinity (primary)
        - Hypothalamic receptor activity
        - Ghrelin pathway interaction
        - 3D shape features (BBB penetration proxy)
        """
        features = {}
        
        # Receptor affinities
        if compound.get("receptor_affinities"):
            for receptor, measurements in compound["receptor_affinities"].items():
                avg_affinity = sum(m["affinity_value"] for m in measurements) / len(measurements)
                features[f"{receptor}_affinity"] = avg_affinity
        
        # 2D + 3D descriptors
        features.update(compound.get("rdkit_descriptors", {}))
        if compound.get("rdkit_descriptors_3d"):
            features.update(compound["rdkit_descriptors_3d"])
        
        return pd.DataFrame([features])
    
    def predict_appetite_effects(self, compound: Dict) -> Dict:
        """Predict appetite stimulation profile.
        
        Returns:
            Appetite stimulation score, side effects, patient tolerance
        """
        features = self.prepare_features(compound)
        
        # Predict appetite stimulation (0-1 scale)
        appetite_score = self.model.predict_proba(features)[0][1]
        
        # Side effect prediction (rule-based for MVP)
        side_effects = self._predict_side_effects(compound, appetite_score)
        
        # Patient tolerance recommendations
        tolerance = self._classify_tolerance(appetite_score, side_effects)
        
        return {
            "appetite_stimulation_score": float(appetite_score),
            "classification": self._classify_appetite_level(appetite_score),
            "side_effects": side_effects,
            "patient_tolerance": tolerance,
            "recommendations": self._generate_recommendations(appetite_score, tolerance)
        }
    
    def _predict_side_effects(self, compound: Dict, appetite_score: float) -> List[Dict]:
        """Predict side effects based on compound properties."""
        side_effects = []
        
        # Munchies (directly correlated with appetite stimulation)
        if appetite_score > 0.6:
            side_effects.append({
                "effect": "increased_hunger",
                "severity": "moderate" if appetite_score < 0.8 else "high",
                "probability": appetite_score
            })
        
        # Dry mouth (CB1 receptor activation)
        cb1_affinity = compound.get("receptor_affinities", {}).get("CB1", [{}])[0].get("affinity_value", 0)
        if cb1_affinity > 6.0:  # pKi > 6.0
            side_effects.append({
                "effect": "dry_mouth",
                "severity": "mild",
                "probability": 0.4
            })
        
        return side_effects
    
    def _classify_appetite_level(self, score: float) -> str:
        """Classify appetite stimulation level."""
        if score > 0.7:
            return "HIGH"
        elif score > 0.4:
            return "MODERATE"
        else:
            return "LOW"
    
    def _classify_tolerance(self, appetite_score: float, side_effects: List[Dict]) -> str:
        """Determine patient tolerance level needed."""
        if appetite_score > 0.7 or len(side_effects) > 2:
            return "HIGH_TOLERANCE_REQUIRED"
        elif appetite_score > 0.4:
            return "MODERATE_TOLERANCE_ACCEPTABLE"
        else:
            return "LOW_TOLERANCE_ACCEPTABLE"
    
    def _generate_recommendations(self, appetite_score: float, tolerance: str) -> List[str]:
        """Generate usage recommendations."""
        recs = []
        
        if appetite_score > 0.7:
            recs.append("Best used for patients with significant appetite loss (cancer, HIV/AIDS)")
            recs.append("Consider timing doses before meals")
            recs.append("Monitor for excessive weight gain")
        
        if tolerance == "HIGH_TOLERANCE_REQUIRED":
            recs.append("Start with low dose and titrate slowly")
            recs.append("Not recommended for appetite suppression needs")
        
        return recs
```

**Train appetite model:**

```bash
python scripts/train_appetite_[model.py](http://model.py)
```

---

### Evidence-Gated Terpene Features (5 hours)

**Task 7.2: NEW - Terpene Evidence Gating** (5 hours)

Create `backend/services/terpene_[analyzer.py](http://analyzer.py)`:

```python
from typing import Dict, List, Optional
from enum import Enum

class EvidenceTier(Enum):
    """Evidence quality hierarchy for terpene interactions."""
    TIER_5 = 5  # Multiple RCTs, systematic reviews
    TIER_4 = 4  # Well-designed observational studies
    TIER_3 = 3  # Smaller studies, case series
    TIER_2 = 2  # In vitro studies, animal models
    TIER_1 = 1  # Anecdotal, theoretical

class TerpeneAnalyzer:
    """Analyze terpene-cannabinoid interactions with evidence gating.
    
    Only includes interactions with sufficient evidence (Tier 3+).
    """
    
    def __init__(self):
        # Terpene database with evidence tiers
        self.terpene_database = {
            "beta_caryophyllene": {
                "mechanism": "CB2_agonist",
                "evidence_tier": EvidenceTier.TIER_4,
                "pubmed_ids": ["23611293", "18574142"],
                "effects": ["anti-inflammatory", "analgesic"],
                "cannabinoid_synergy": {
                    "CBD": {"enhancement_factor": 1.3, "evidence_tier": EvidenceTier.TIER_3},
                    "CBG": {"enhancement_factor": 1.4, "evidence_tier": EvidenceTier.TIER_3}
                }
            },
            "myrcene": {
                "mechanism": "GABA_potentiation",
                "evidence_tier": EvidenceTier.TIER_3,
                "pubmed_ids": ["12587690"],
                "effects": ["sedative", "muscle_relaxant"],
                "cannabinoid_synergy": {
                    "CBD": {"enhancement_factor": 1.2, "evidence_tier": EvidenceTier.TIER_3},
                    "CBN": {"enhancement_factor": 1.5, "evidence_tier": EvidenceTier.TIER_2}
                }
            },
            "limonene": {
                "mechanism": "5HT_modulation",
                "evidence_tier": EvidenceTier.TIER_3,
                "pubmed_ids": ["23473632"],
                "effects": ["anxiolytic", "mood_elevation"],
                "cannabinoid_synergy": {
                    "CBD": {"enhancement_factor": 1.15, "evidence_tier": EvidenceTier.TIER_3}
                }
            },
            "linalool": {
                "mechanism": "NMDA_antagonist",
                "evidence_tier": EvidenceTier.TIER_3,
                "pubmed_ids": ["19922249"],
                "effects": ["anxiolytic", "sedative"],
                "cannabinoid_synergy": {
                    "CBD": {"enhancement_factor": 1.25, "evidence_tier": EvidenceTier.TIER_3}
                }
            }
        }
    
    def analyze_synergy(
        self,
        cannabinoid: str,
        terpene_profile: Dict[str, float],
        min_evidence_tier: int = 3
    ) -> Dict:
        """Analyze terpene-cannabinoid synergy with evidence gating.
        
        Args:
            cannabinoid: Cannabinoid abbreviation (e.g., "CBD")
            terpene_profile: Dict of {terpene_name: concentration_percent}
            min_evidence_tier: Minimum evidence tier to include (default: 3)
            
        Returns:
            Synergy analysis with evidence-supported interactions only
        """
        synergies = []
        total_enhancement = 1.0
        excluded_interactions = []
        
        for terpene, concentration in terpene_profile.items():
            terpene_data = self.terpene_database.get(terpene)
            
            if not terpene_data:
                continue
            
            # Check if cannabinoid-terpene interaction exists
            synergy_data = terpene_data.get("cannabinoid_synergy", {}).get(cannabinoid)
            
            if not synergy_data:
                continue
            
            # Evidence gating: only include if meets minimum tier
            if synergy_data["evidence_tier"].value >= min_evidence_tier:
                # Calculate weighted enhancement based on concentration
                weighted_enhancement = 1.0 + (synergy_data["enhancement_factor"] - 1.0) * (concentration / 2.0)
                total_enhancement *= weighted_enhancement
                
                synergies.append({
                    "terpene": terpene,
                    "concentration_percent": concentration,
                    "mechanism": terpene_data["mechanism"],
                    "enhancement_factor": synergy_data["enhancement_factor"],
                    "evidence_tier": synergy_data["evidence_tier"].value,
                    "evidence_quality": self._tier_to_description(synergy_data["evidence_tier"]),
                    "pubmed_ids": terpene_data["pubmed_ids"]
                })
            else:
                # Track excluded interactions for transparency
                excluded_interactions.append({
                    "terpene": terpene,
                    "cannabinoid": cannabinoid,
                    "evidence_tier": synergy_data["evidence_tier"].value,
                    "reason": f"Evidence tier {synergy_data['evidence_tier'].value} below minimum {min_evidence_tier}"
                })
        
        return {
            "cannabinoid": cannabinoid,
            "total_enhancement_factor": round(total_enhancement, 2),
            "synergies": synergies,
            "num_synergies_included": len(synergies),
            "excluded_interactions": excluded_interactions,
            "evidence_gating_enabled": True,
            "min_evidence_tier_used": min_evidence_tier
        }
    
    def _tier_to_description(self, tier: EvidenceTier) -> str:
        """Convert evidence tier to human-readable description."""
        descriptions = {
            EvidenceTier.TIER_5: "High - Multiple RCTs/systematic reviews",
            EvidenceTier.TIER_4: "Good - Well-designed observational studies",
            EvidenceTier.TIER_3: "Moderate - Smaller studies, case series",
            EvidenceTier.TIER_2: "Low - In vitro/animal studies only",
            EvidenceTier.TIER_1: "Very Low - Anecdotal/theoretical"
        }
        return descriptions.get(tier, "Unknown")
```

**API endpoint for terpene analysis:**

```python
@[router.post](http://router.post)("/api/compounds/{id}/terpene-synergy")
async def analyze_terpene_synergy(
    id: int,
    terpene_profile: Dict[str, float],
    min_evidence_tier: int = 3,
    db: Session = Depends(get_db)
):
    """Analyze terpene-cannabinoid synergy with evidence gating."""
    compound = db.query(Cannabinoid).filter([Cannabinoid.id](http://Cannabinoid.id) == id).first()
    
    if not compound:
        raise HTTPException(status_code=404, detail="Compound not found")
    
    analyzer = TerpeneAnalyzer()
    result = analyzer.analyze_synergy(
        cannabinoid=compound.abbreviation,
        terpene_profile=terpene_profile,
        min_evidence_tier=min_evidence_tier
    )
    
    return result
```

---

### PatentPath Lite Development - Phase 2 (20 hours)

**Task 7.3: Feature 3 - FTO (Freedom-to-Operate) Checker** (10-12 hours)

<aside>
⚠️

**UX Disclaimer Requirement**: All FTO and claim generation outputs must include a visible footer disclaimer: *"This is decision support only. Not legal advice. Consult qualified patent counsel before filing or commercial decisions."* Add this to API responses and UI components.

</aside>

Create `backend/services/patentpath_lite/fto_[checker.py](http://checker.py)`:

```python
from typing import Dict, List
from .prior_art_search import PriorArtSearcher
import datetime

class FTOChecker:
    """Basic freedom-to-operate analysis for compounds.
    
    MVP version checks:
    - Active patents in US only
    - Basic blocking patent identification
    - Estimated licensing costs
    
    Full version (PatentPath Full) includes:
    - Multi-jurisdiction analysis
    - Detailed infringement probability
    - Design-around strategies
    """
    
    def __init__(self):
        self.searcher = PriorArtSearcher()
    
    def check_fto(
        self,
        compound_description: str,
        keywords: List[str],
        target_market: str = "US",
        intended_use: str = None
    ) -> Dict:
        """Perform FTO analysis for a compound.
        
        Args:
            compound_description: Technical description
            keywords: Search keywords
            target_market: Geographic market (US only for MVP)
            intended_use: Therapeutic indication
            
        Returns:
            FTO analysis with risk assessment and mitigation options
        """
        # Search for potentially blocking patents
        search_results = [self.searcher.search](http://self.searcher.search)(
            compound_description=compound_description,
            keywords=keywords,
            databases=["USPTO"]  # US only for MVP
        )
        
        prior_art_patents = search_results.get("prior_art_patents", [])
        
        # Filter for active patents only
        active_patents = self._filter_active_patents(prior_art_patents)
        
        # Analyze each patent for blocking potential
        blocking_patents = []
        for patent in active_patents:
            blocking_analysis = self._analyze_blocking_potential(
                patent,
                compound_description,
                intended_use
            )
            
            if blocking_analysis["is_blocking"]:
                blocking_patents.append(blocking_analysis)
        
        # Calculate overall FTO risk
        fto_risk = self._calculate_fto_risk(blocking_patents)
        
        # Estimate licensing costs
        licensing_costs = self._estimate_licensing_costs(blocking_patents)
        
        return {
            "fto_risk_level": fto_risk["level"],
            "risk_score": fto_risk["score"],
            "target_market": target_market,
            "blocking_patents": blocking_patents,
            "num_blocking_patents": len(blocking_patents),
            "estimated_licensing_costs": licensing_costs,
            "recommendation": self._generate_fto_recommendation(fto_risk, blocking_patents),
            "mitigation_strategies": self._generate_mitigation_strategies(blocking_patents),
            "next_steps": self._generate_fto_next_steps(fto_risk)
        }
    
    def _filter_active_patents(self, patents: List[Dict]) -> List[Dict]:
        """Filter for active (not expired/abandoned) patents."""
        active = []
        current_year = [datetime.datetime.now](http://datetime.datetime.now)().year
        
        for patent in patents:
            # Simple heuristic: patents expire 20 years from filing
            # In production, query USPTO PAIR for actual status
            patent_year = self._extract_year(patent.get("date", ""))
            
            if patent_year and (current_year - patent_year) < 20:
                patent["status"] = "active"
                patent["estimated_expiration"] = patent_year + 20
                active.append(patent)
            else:
                patent["status"] = "expired"
        
        return active
    
    def _extract_year(self, date_str: str) -> int:
        """Extract year from date string."""
        try:
            return int(date_str[:4]) if date_str else None
        except:
            return None
    
    def _analyze_blocking_potential(
        self,
        patent: Dict,
        compound_description: str,
        intended_use: str
    ) -> Dict:
        """Analyze if a patent could block commercialization.
        
        MVP uses simple keyword matching.
        Full version uses claim-by-claim element analysis.
        """
        # Simple blocking determination based on relevance
        relevance = patent.get("relevance_score", 0)
        
        # High relevance = potential blocking patent
        is_blocking = relevance > 0.6
        
        if is_blocking:
            risk_level = "HIGH" if relevance > 0.8 else "MODERATE"
        else:
            risk_level = "LOW"
        
        return {
            "patent_number": patent.get("patent_number"),
            "title": patent.get("title"),
            "relevance_score": relevance,
            "is_blocking": is_blocking,
            "risk_level": risk_level,
            "expiration_year": patent.get("estimated_expiration"),
            "years_until_expiration": patent.get("estimated_expiration", 0) - [datetime.datetime.now](http://datetime.datetime.now)().year,
            "reason": self._generate_blocking_reason(relevance, intended_use),
            "mitigation_options": self._generate_mitigation_options(risk_level)
        }
    
    def _generate_blocking_reason(self, relevance: float, intended_use: str) -> str:
        """Generate explanation for why patent may be blocking."""
        if relevance > 0.8:
            return f"Highly similar compound structure and {intended_use or 'therapeutic use'}"
        elif relevance > 0.6:
            return f"Similar compound class, overlapping {intended_use or 'indications'}"
        else:
            return "Low blocking risk - significant structural differences"
    
    def _generate_mitigation_options(self, risk_level: str) -> List[str]:
        """Generate mitigation strategies for blocking patents."""
        if risk_level == "HIGH":
            return [
                "Negotiate licensing agreement",
                "Challenge patent validity (prior art search)",
                "Design around with structural modifications",
                "Wait for patent expiration"
            ]
        elif risk_level == "MODERATE":
            return [
                "Conduct detailed claim analysis with attorney",
                "Explore design-around opportunities",
                "Consider licensing if commercialization imminent"
            ]
        else:
            return [
                "Monitor patent status",
                "Proceed with caution"
            ]
    
    def _calculate_fto_risk(self, blocking_patents: List[Dict]) -> Dict:
        """Calculate overall FTO risk level."""
        if not blocking_patents:
            return {"level": "LOW", "score": 0.1}
        
        # Count high-risk patents
        high_risk_count = sum(1 for p in blocking_patents if p["risk_level"] == "HIGH")
        moderate_risk_count = sum(1 for p in blocking_patents if p["risk_level"] == "MODERATE")
        
        # Risk score (0-1)
        risk_score = min(1.0, (high_risk_count * 0.3 + moderate_risk_count * 0.15))
        
        # Risk level classification
        if risk_score > 0.6:
            level = "HIGH"
        elif risk_score > 0.3:
            level = "MODERATE"
        else:
            level = "LOW"
        
        return {"level": level, "score": risk_score}
    
    def _estimate_licensing_costs(self, blocking_patents: List[Dict]) -> Dict:
        """Estimate licensing costs for blocking patents.
        
        MVP uses industry averages.
        Full version uses deal database and ML prediction.
        """
        if not blocking_patents:
            return {
                "total_estimated_cost": "$0",
                "per_patent_average": "$0",
                "breakdown": []
            }
        
        # Industry average licensing costs (simplified)
        cost_per_patent = {
            "HIGH": 1000000,  # $1M for high-risk blocking patent
            "MODERATE": 250000,  # $250K for moderate-risk
            "LOW": 50000  # $50K for low-risk
        }
        
        breakdown = []
        total_cost = 0
        
        for patent in blocking_patents:
            cost = cost_per_patent.get(patent["risk_level"], 0)
            total_cost += cost
            
            breakdown.append({
                "patent_number": patent["patent_number"],
                "risk_level": patent["risk_level"],
                "estimated_cost": f"${cost:,}"
            })
        
        return {
            "total_estimated_cost": f"${total_cost:,}",
            "per_patent_average": f"${total_cost // len(blocking_patents):,}",
            "breakdown": breakdown,
            "note": "Estimates based on industry averages. Actual costs vary significantly."
        }
    
    def _generate_fto_recommendation(self, fto_risk: Dict, blocking_patents: List[Dict]) -> str:
        """Generate FTO recommendation."""
        risk_level = fto_risk["level"]
        
        if risk_level == "HIGH":
            return "HIGH RISK: Commercialization likely blocked without licensing or design-around. Consult patent attorney immediately."
        elif risk_level == "MODERATE":
            return "MODERATE RISK: Some blocking patents identified. Attorney review recommended before significant investment."
        else:
            return "LOW RISK: No major blocking patents identified. Proceed with standard patent due diligence."
    
    def _generate_fto_next_steps(self, fto_risk: Dict) -> List[str]:
        """Generate actionable next steps."""
        risk_level = fto_risk["level"]
        
        if risk_level == "HIGH":
            return [
                "Engage patent attorney for detailed claim analysis",
                "Obtain freedom-to-operate opinion letter",
                "Explore licensing negotiations with patent holders",
                "Investigate design-around opportunities"
            ]
        elif risk_level == "MODERATE":
            return [
                "Consult patent attorney for claim mapping",
                "Consider provisional patent filing to establish priority",
                "Monitor patent status and prosecution history"
            ]
        else:
            return [
                "Proceed with compound development",
                "Periodic patent landscape monitoring",
                "File provisional patent application to secure priority"
            ]
```

**API endpoint for FTO check:**

```python
@[router.post](http://router.post)("/api/patentpath/fto-check")
async def check_fto(
    compound_description: str,
    keywords: List[str],
    target_market: str = "US",
    intended_use: str = None
):
    """Perform FTO analysis for a compound."""
    checker = FTOChecker()
    result = checker.check_fto(
        compound_description=compound_description,
        keywords=keywords,
        target_market=target_market,
        intended_use=intended_use
    )
    return result
```

---

**Task 7.4: Feature 4 - Templated Claim Generator** (8-10 hours)

Create `backend/services/patentpath_lite/claim_[generator.py](http://generator.py)`:

```python
from typing import Dict, List
from datetime import datetime

class ClaimGenerator:
    """Generate patent claims from templates.
    
    MVP version uses fill-in-the-blank templates.
    Full version (PatentPath Full) uses AI for USPTO-compliant claims with 78% acceptance rate.
    """
    
    def __init__(self):
        self.templates = self._load_templates()
    
    def _load_templates(self) -> Dict:
        """Load patent claim templates."""
        return {
            "composition_of_matter": """CLAIM 1: A compound having the structure of Formula I:

{structure_representation}

wherein:
- R1 is {r1_definition};
- R2 is {r2_definition};
- The compound is characterized by {key_features}.

CLAIM 2: The compound of Claim 1, wherein the compound is selected from the group consisting of:
{specific_compounds}.

CLAIM 3: The compound of Claim 1, wherein the compound is a {compound_class}.""",
            
            "method_of_use": """CLAIM 1: A method for treating {condition}, comprising:
administering to a subject in need thereof a therapeutically effective amount of a compound of Formula I, wherein the compound is:

{compound_structure}

CLAIM 2: The method of Claim 1, wherein the therapeutically effective amount is between {dose_low} and {dose_high} mg per day.

CLAIM 3: The method of Claim 1, wherein the subject is a human.

CLAIM 4: The method of Claim 1, further comprising administering one or more additional therapeutic agents.""",
            
            "pharmaceutical_composition": """CLAIM 1: A pharmaceutical composition comprising:
(a) a therapeutically effective amount of a compound of Formula I:

{compound_structure}

and
(b) a pharmaceutically acceptable carrier.

CLAIM 2: The pharmaceutical composition of Claim 1, wherein the composition is formulated for {admin_route} administration.

CLAIM 3: The pharmaceutical composition of Claim 1, further comprising one or more excipients selected from the group consisting of: {excipients}.""",
            
            "method_of_synthesis": """CLAIM 1: A method for synthesizing a compound of Formula I, comprising:
(a) reacting {starting_material_1} with {starting_material_2} under {reaction_conditions} to form an intermediate compound;
(b) {step_2_description};
(c) isolating the compound of Formula I.

CLAIM 2: The method of Claim 1, wherein the reaction is conducted at a temperature between {temp_low}°C and {temp_high}°C.

CLAIM 3: The method of Claim 1, wherein the yield is at least {yield_percent}%."""
        }
    
    def generate_claims(
        self,
        compound_name: str,
        compound_smiles: str,
        therapeutic_use: str,
        key_features: List[str],
        claim_types: List[str] = ["composition_of_matter", "method_of_use", "pharmaceutical_composition"]
    ) -> Dict:
        """Generate patent claims for a compound.
        
        Args:
            compound_name: Name of compound
            compound_smiles: SMILES structure
            therapeutic_use: Indication (e.g., "Alzheimer's disease")
            key_features: List of key structural features
            claim_types: Which claim types to generate
            
        Returns:
            Generated claims with USPTO formatting
        """
        generated_claims = {}
        total_claims = 0
        
        for claim_type in claim_types:
            if claim_type not in self.templates:
                continue
            
            template = self.templates[claim_type]
            
            # Fill template with compound-specific data
            filled_claims = self._fill_template(
                template,
                claim_type,
                compound_name=compound_name,
                compound_smiles=compound_smiles,
                therapeutic_use=therapeutic_use,
                key_features=key_features
            )
            
            generated_claims[claim_type] = filled_claims
            total_claims += self._count_claims(filled_claims)
        
        return {
            "compound_name": compound_name,
            "generated_claims": generated_claims,
            "claim_types": claim_types,
            "total_claims": total_claims,
            "format": "USPTO_utility_patent",
            "generated_date": [datetime.now](http://datetime.now)().isoformat(),
            "estimated_filing_cost": self._estimate_filing_cost(total_claims),
            "estimated_prosecution_cost": self._estimate_prosecution_cost(total_claims),
            "disclaimer": "These claims are template-generated for planning purposes only. Consult a patent attorney for filing-ready claims.",
            "upgrade_note": "PatentPath Full uses AI to generate USPTO-compliant claims with 78% acceptance rate vs. 45-60% industry average."
        }
    
    def _fill_template(
        self,
        template: str,
        claim_type: str,
        **kwargs
    ) -> str:
        """Fill claim template with compound-specific data."""
        
        # Default values
        defaults = {
            "structure_representation": "[Chemical Structure Diagram]",
            "compound_structure": kwargs.get("compound_smiles", "[SMILES: To be inserted]"),
            "r1_definition": "hydrogen, alkyl, or aryl",
            "r2_definition": "hydroxyl, alkoxy, or halogen",
            "key_features": ", ".join(kwargs.get("key_features", ["novel structure"])),
            "specific_compounds": kwargs.get("compound_name", "Compound A"),
            "compound_class": self._infer_compound_class(kwargs.get("compound_name", "")),
            "condition": kwargs.get("therapeutic_use", "[Therapeutic Indication]"),
            "dose_low": "10",
            "dose_high": "500",
            "admin_route": "oral",
            "excipients": "lactose, cellulose, magnesium stearate",
            "starting_material_1": "[Starting Material 1]",
            "starting_material_2": "[Starting Material 2]",
            "reaction_conditions": "oxidative conditions",
            "step_2_description": "purifying the intermediate compound by chromatography",
            "temp_low": "20",
            "temp_high": "80",
            "yield_percent": "50"
        }
        
        # Merge with provided kwargs
        fill_data = {**defaults, **kwargs}
        
        # Fill template
        try:
            filled = template.format(**fill_data)
        except KeyError as e:
            filled = template  # Return unfilled if missing keys
        
        return filled
    
    def _infer_compound_class(self, compound_name: str) -> str:
        """Infer compound class from name."""
        name_lower = compound_name.lower()
        
        if "dimer" in name_lower:
            return "dimeric cannabinoid"
        elif any(cb in name_lower for cb in ["cbd", "thc", "cbg", "cbn"]):
            return "cannabinoid derivative"
        else:
            return "organic compound"
    
    def _count_claims(self, claims_text: str) -> int:
        """Count number of claims in text."""
        return claims_text.count("CLAIM ")
    
    def _estimate_filing_cost(self, num_claims: int) -> str:
        """Estimate USPTO filing cost.
        
        Based on 2025 USPTO fee schedule:
        - Utility filing fee: $320 (micro), $640 (small), $1,280 (large)
        - Search fee: $660 (micro), $1,320 (small), $2,640 (large)
        - Examination fee: $760 (micro), $1,520 (small), $3,040 (large)
        - Additional claim fees: $100 per claim over 20
        """
        base_filing = 1280  # Large entity
        base_search = 2640
        base_exam = 3040
        
        if num_claims > 20:
            extra_claims_fee = (num_claims - 20) * 100
        else:
            extra_claims_fee = 0
        
        total = base_filing + base_search + base_exam + extra_claims_fee
        
        return f"${total:,} (large entity, USPTO)"
    
    def _estimate_prosecution_cost(self, num_claims: int) -> str:
        """Estimate attorney prosecution costs."""
        # Industry averages: $5,000-15,000 for prosecution
        base_prosecution = 8000
        complex_fee = max(0, (num_claims - 15) * 200)  # More claims = more complex
        
        total = base_prosecution + complex_fee
        
        return f"${total:,}-${total*1.5:,.0f}"
```

**API endpoint for claim generation:**

```python
@[router.post](http://router.post)("/api/patentpath/generate-claims")
async def generate_claims(
    compound_name: str,
    compound_smiles: str,
    therapeutic_use: str,
    key_features: List[str],
    claim_types: List[str] = ["composition_of_matter", "method_of_use", "pharmaceutical_composition"]
):
    """Generate patent claims for a compound."""
    generator = ClaimGenerator()
    result = generator.generate_claims(
        compound_name=compound_name,
        compound_smiles=compound_smiles,
        therapeutic_use=therapeutic_use,
        key_features=key_features,
        claim_types=claim_types
    )
    return result
```

---

**Success Metrics - Week 7**:

- Appetite stimulation model trained
- Side effect predictions functional
- PatentPath Lite Features 3-4 implemented
- FTO checker returns risk assessment with licensing costs
- Claim generator produces USPTO-formatted claims
- Terpene evidence gating excludes low-quality interactions
- All features accessible via API

---

### Week 8: Frontend + PatentPath Lite Integration (Features 5-6)

**Dates**: Feb 5-12, 2026

**Estimated Hours**: 35-40 hours

**Objectives**:

- Build Next.js frontend foundation
- **NEW**: Complete PatentPath Lite Features 5-6
- **NEW**: Integrate PatentPath into frontend UI

**Deliverables**:

✅ Next.js 14 app with TypeScript

✅ Dashboard layout and navigation

✅ Compound search and browse interface

✅ **NEW**: PatentPath Lite - TK Attribution Checker (Feature 5)

✅ **NEW**: PatentPath Lite - Cost Estimator (Feature 6)

✅ **NEW**: PatentPath UI integration

✅ React Query for API calls

**Split:** 20 hours Frontend + 15 hours PatentPath Lite

---

### PatentPath Lite Development - Phase 3 (15 hours)

**Task 8.1: Feature 5 - TK Attribution Checker** (6-8 hours)

Create `backend/services/patentpath_lite/tk_[checker.py](http://checker.py)`:

```python
from typing import Dict, List
from datetime import datetime

class TKAttributionChecker:
    """Check if compound requires traditional knowledge attribution.
    
    MVP version uses keyword detection.
    Full version (PatentPath Full) uses ML classifier (91.8% accuracy) + TKDL integration.
    """
    
    def __init__(self):
        # Keywords indicating TK derivation
        [self.tk](http://self.tk)_keywords = [
            "traditional", "indigenous", "ethnobotanical", "folk medicine",
            "ayurvedic", "traditional chinese medicine", "tcm", "kampo",
            "native", "aboriginal", "tribal", "ancestral",
            "ancient", "historical preparation", "traditional use"
        ]
        
        # Sacred knowledge indicators (absolute bar to patenting)
        self.sacred_keywords = [
            "sacred", "ceremonial", "spiritual", "religious",
            "ritual", "holy", "blessed", "consecrated"
        ]
    
    def check_tk_attribution(
        self,
        compound_description: str,
        derivation_source: str = None,
        source_description: str = None
    ) -> Dict:
        """Check if compound requires TK attribution.
        
        Args:
            compound_description: Description of compound
            derivation_source: Source type (e.g., "ethnobotanical_literature")
            source_description: Detailed description of source
            
        Returns:
            TK attribution requirements and next steps
        """
        # Check for TK indicators
        is_tk_derived = self._detect_tk_indicators(
            compound_description,
            derivation_source,
            source_description
        )
        
        # Check for sacred knowledge (blocking condition)
        is_sacred = self._detect_sacred_knowledge(
            compound_description,
            source_description
        )
        
        if is_sacred:
            return {
                "tk_derived": True,
                "sacred_knowledge_detected": True,
                "attribution_required": False,  # Cannot be patented at all
                "can_proceed": False,
                "blocking_reason": "Sacred traditional knowledge - cannot be patented under any circumstances",
                "recommendation": "Do not proceed with patent filing. Consider alternative protection strategies or community partnership without IP claims."
            }
        
        if is_tk_derived:
            # TK-derived but not sacred - attribution required
            return {
                "tk_derived": True,
                "sacred_knowledge_detected": False,
                "attribution_required": True,
                "can_proceed": True,
                "requirements": self._generate_attribution_requirements(),
                "next_steps": self._generate_tk_next_steps(),
                "equipath_integration": {
                    "community_wallet_setup": "Required before filing patent",
                    "estimated_revenue_share": "70% of licensing revenue to community",
                    "attribution_blockchain_id": "pending"
                },
                "patent_specification_requirements": self._generate_specification_requirements(),
                "legal_compliance": "Required under Nagoya Protocol and potential national ABS laws"
            }
        else:
            # Not TK-derived
            return {
                "tk_derived": False,
                "sacred_knowledge_detected": False,
                "attribution_required": False,
                "can_proceed": True,
                "recommendation": "No TK attribution required. Proceed with standard patent filing."
            }
    
    def _detect_tk_indicators(
        self,
        compound_description: str,
        derivation_source: str,
        source_description: str
    ) -> bool:
        """Detect if compound is derived from traditional knowledge."""
        # Check all text sources
        all_text = " ".join(filter(None, [
            compound_description or "",
            derivation_source or "",
            source_description or ""
        ])).lower()
        
        # Look for TK keywords
        for keyword in [self.tk](http://self.tk)_keywords:
            if keyword in all_text:
                return True
        
        # Check derivation source type
        tk_sources = ["ethnobotanical", "traditional", "indigenous", "folk"]
        if derivation_source and any(src in derivation_source.lower() for src in tk_sources):
            return True
        
        return False
    
    def _detect_sacred_knowledge(
        self,
        compound_description: str,
        source_description: str
    ) -> bool:
        """Detect if knowledge is sacred (absolute bar to patenting).
        
        MVP uses keyword detection.
        Full version uses ML classifier with high threshold (0.85) to prevent false positives.
        """
        all_text = " ".join(filter(None, [
            compound_description or "",
            source_description or ""
        ])).lower()
        
        for keyword in self.sacred_keywords:
            if keyword in all_text:
                return True
        
        return False
    
    def _generate_attribution_requirements(self) -> List[str]:
        """Generate requirements for TK attribution."""
        return [
            "Document traditional knowledge source with community name and location",
            "Establish informed consent via EthnoPath consent system",
            "Configure benefit-sharing via EquiPath (70% to community, 25% STEM education, 5% infrastructure)",
            "Include attribution statement in patent specification",
            "Cite traditional knowledge sources in patent application",
            "Establish ongoing community partnership agreement"
        ]
    
    def _generate_tk_next_steps(self) -> List[str]:
        """Generate actionable next steps for TK attribution."""
        return [
            "Contact community leaders to establish consent",
            "Register community benefit-sharing wallet in EquiPath",
            "Draft attribution language for patent specification",
            "Consult attorney experienced in TK and Nagoya Protocol compliance",
            "Delay patent filing until consent and benefit-sharing established"
        ]
    
    def _generate_specification_requirements(self) -> Dict:
        """Generate patent specification requirements for TK attribution."""
        return {
            "attribution_statement_template": (
                "The compound disclosed herein was developed based on traditional knowledge "
                "from [COMMUNITY NAME], who have historically used [PLANT/PREPARATION METHOD] "
                "for [TRADITIONAL USE]. The inventors acknowledge the contribution of this "
                "traditional knowledge and have established a benefit-sharing agreement whereby "
                "[XX]% of commercialization revenues will be shared with the [COMMUNITY NAME] "
                "pursuant to the Nagoya Protocol."
            ),
            "citation_requirements": [
                "Include community name and geographic location",
                "Cite ethnobotanical literature documenting traditional use",
                "Reference benefit-sharing agreement",
                "Include EquiPath blockchain manifest ID"
            ],
            "background_section": "Dedicate Background section subsection to traditional knowledge source"
        }
```

**API endpoint for TK attribution check:**

```python
@[router.post](http://router.post)("/api/patentpath/tk-check")
async def check_tk_attribution(
    compound_description: str,
    derivation_source: str = None,
    source_description: str = None
):
    """Check if compound requires TK attribution."""
    checker = TKAttributionChecker()
    result = checker.check_tk_attribution(
        compound_description=compound_description,
        derivation_source=derivation_source,
        source_description=source_description
    )
    return result
```

---

**Task 8.2: Feature 6 - Cost Estimator** (5-6 hours)

Create `backend/services/patentpath_lite/cost_[estimator.py](http://estimator.py)`:

```python
from typing import Dict, List

class CostEstimator:
    """Estimate patent filing and maintenance costs.
    
    MVP version uses USPTO fee schedule and industry averages.
    Full version includes multi-jurisdiction analysis and portfolio optimization.
    """
    
    def __init__(self):
        # 2025 USPTO fee schedule (large entity)
        self.uspto_fees = {
            "filing_fee": 1280,
            "search_fee": 2640,
            "examination_fee": 3040,
            "issue_fee": 1000,
            "per_claim_over_20": 100,
            "per_independent_claim_over_3": 460
        }
        
        # Maintenance fees (years after grant)
        self.maintenance_fees = {
            "3.5_years": 2000,
            "7.5_years": 3760,
            "11.5_years": 7060
        }
        
        # Attorney cost estimates
        [self.attorney](http://self.attorney)_costs = {
            "simple_drafting": 5000,
            "complex_drafting": 15000,
            "prosecution_simple": 3000,
            "prosecution_complex": 12000,
            "office_action_response": 2500
        }
    
    def estimate_costs(
        self,
        num_compounds: int,
        jurisdictions: List[str] = ["US"],
        claim_count_avg: int = 15,
        complexity: str = "moderate",
        strategy: str = "patent"  # or "trade_secret"
    ) -> Dict:
        """Estimate total IP protection costs.
        
        Args:
            num_compounds: Number of compounds to protect
            jurisdictions: List of countries (US only for MVP)
            claim_count_avg: Average claims per patent
            complexity: \"simple\", \"moderate\", or \"complex\"
            strategy: \"patent\" or \"trade_secret\"
            
        Returns:
            Detailed cost breakdown with comparisons
        """
        if strategy == "patent":
            return self._estimate_patent_costs(
                num_compounds,
                jurisdictions,
                claim_count_avg,
                complexity
            )
        else:  # trade_secret
            return self._estimate_trade_secret_costs(num_compounds)
    
    def _estimate_patent_costs(
        self,
        num_compounds: int,
        jurisdictions: List[str],
        claim_count_avg: int,
        complexity: str
    ) -> Dict:
        \"\"\"Estimate patent costs.\"\"\"
        
        # Per-patent costs
        filing_cost_per = self._calculate_filing_cost(claim_count_avg)
        prosecution_cost_per = self._calculate_prosecution_cost(complexity)
        maintenance_cost_per = sum(self.maintenance_fees.values())
        
        # Total costs
        filing_total = filing_cost_per * num_compounds
        prosecution_total = prosecution_cost_per * num_compounds
        maintenance_total = maintenance_cost_per * num_compounds
        
        total_20_year = filing_total + prosecution_total + maintenance_total
        
        return {
            "strategy": \"patent\",
            \"num_compounds\": num_compounds,
            \"jurisdictions\": jurisdictions,
            \"total_cost_estimate\": {
                \"filing_costs\": f\"${filing_total:,}\",
                \"prosecution_costs\": f\"${prosecution_total:,}\",
                \"maintenance_fees_20yr\": f\"${maintenance_total:,}\",
                \"total_investment\": f\"${total_20_year:,}\"
            },
            \"per_patent_breakdown\": {
                \"filing\": f\"${filing_cost_per:,}\",
                \"prosecution\": f\"${prosecution_cost_per:,}\",
                \"maintenance\": f\"${maintenance_cost_per:,}\",
                \"total_per_patent\": f\"${filing_cost_per + prosecution_cost_per + maintenance_cost_per:,}\"
            },
            \"timeline\": {
                \"filing_to_examination\": \"12-18 months\",
                \"examination_to_grant\": \"18-36 months\",
                \"total_to_grant\": \"2.5-4.5 years\",
                \"protection_duration\": \"20 years from filing\"
            },
            \"cost_schedule\": self._generate_cost_schedule(
                filing_cost_per,
                prosecution_cost_per,
                num_compounds
            ),
            \"recommendation\": self._generate_patent_recommendation(num_compounds, total_20_year)
        }
    
    def _estimate_trade_secret_costs(self, num_compounds: int) -> Dict:
        \"\"\"Estimate trade secret protection costs.\"\"\"
        
        # Trade secret costs: documentation + security measures
        documentation_per_compound = 2000  # Internal documentation
        security_annual = 5000  # Security measures, NDA management
        
        documentation_total = documentation_per_compound * num_compounds
        security_20_year = security_annual * 20
        
        total_20_year = documentation_total + security_20_year
        
        return {
            \"strategy\": \"trade_secret\",
            \"num_compounds\": num_compounds,
            \"total_cost_estimate\": {
                \"documentation_costs\": f\"${documentation_total:,}\",
                \"security_measures_20yr\": f\"${security_20_year:,}\",
                \"total_investment\": f\"${total_20_year:,}\"
            },
            \"per_compound_breakdown\": {
                \"upfront_documentation\": f\"${documentation_per_compound:,}\",
                \"annual_security\": f\"${security_annual:,}\",
                \"total_per_compound_20yr\": f\"${documentation_per_compound + security_20_year:,}\"
            },
            \"advantages\": [
                \"No expiration (vs. 20 years for patents)\",
                \"No public disclosure required\",
                \"Lower upfront costs\",
                \"No maintenance fees\"
            ],
            \"disadvantages\": [
                \"No legal exclusivity if reverse-engineered\",
                \"Cannot prevent independent discovery\",
                \"Difficult to enforce\",
                \"No blocking power against competitors\"
            ],
            \"recommendation\": \"Trade secrets best for: manufacturing processes, formulations difficult to reverse-engineer. Patents best for: novel compound structures, therapeutic uses.\"
        }
    
    def _calculate_filing_cost(self, claim_count: int) -> int:
        \"\"\"Calculate filing cost including claim fees.\"\"\"
        base = self.uspto_fees[\"filing_fee\"] + self.uspto_fees[\"search_fee\"] + self.uspto_fees[\"examination_fee\"]
        
        # Extra claim fees
        if claim_count > 20:
            base += (claim_count - 20) * self.uspto_fees[\"per_claim_over_20\"]
        
        # Assume 3 independent claims (typical)
        independent_claims = 3
        if independent_claims > 3:
            base += (independent_claims - 3) * self.uspto_fees[\"per_independent_claim_over_3\"]
        
        return base
    
    def _calculate_prosecution_cost(self, complexity: str) -> int:
        \"\"\"Calculate attorney prosecution costs.\"\"\"
        if complexity == \"simple\":
            drafting = [self.attorney](http://self.attorney)_costs[\"simple_drafting\"]
            prosecution = [self.attorney](http://self.attorney)_costs[\"prosecution_simple\"]
        elif complexity == \"complex\":
            drafting = [self.attorney](http://self.attorney)_costs[\"complex_drafting\"]
            prosecution = [self.attorney](http://self.attorney)_costs[\"prosecution_complex\"]
        else:  # moderate
            drafting = ([self.attorney](http://self.attorney)_costs[\"simple_drafting\"] + [self.attorney](http://self.attorney)_costs[\"complex_drafting\"]) // 2
            prosecution = ([self.attorney](http://self.attorney)_costs[\"prosecution_simple\"] + [self.attorney](http://self.attorney)_costs[\"prosecution_complex\"]) // 2
        
        # Expect 1-2 office actions
        office_actions = 1.5 * [self.attorney](http://self.attorney)_costs[\"office_action_response\"]
        
        return int(drafting + prosecution + office_actions)
    
    def _generate_cost_schedule(self, filing_per: int, prosecution_per: int, num_compounds: int) -> List[Dict]:
        \"\"\"Generate payment schedule over time.\"\"\"
        return [
            {\"year\": 0, \"event\": \"Filing\", \"cost\": f\"${filing_per * num_compounds:,}\"},
            {\"year\": 1, \"event\": \"Prosecution begins\", \"cost\": f\"${prosecution_per * num_compounds:,}\"},
            {\"year\": 3.5, \"event\": \"1st maintenance fee\", \"cost\": f\"${self.maintenance_fees['3.5_years'] * num_compounds:,}\"},
            {\"year\": 7.5, \"event\": \"2nd maintenance fee\", \"cost\": f\"${self.maintenance_fees['7.5_years'] * num_compounds:,}\"},
            {\"year\": 11.5, \"event\": \"3rd maintenance fee\", \"cost\": f\"${self.maintenance_fees['11.5_years'] * num_compounds:,}\"}
        ]
    
    def _generate_patent_recommendation(self, num_compounds: int, total_cost: int) -> str:
        \"\"\"Generate cost-based recommendation.\"\"\"
        cost_per_compound = total_cost / num_compounds
        
        if cost_per_compound > 30000:
            return f\"HIGH COST: Consider prioritizing top {num_compounds // 2} compounds and keeping others as trade secrets to reduce costs.\"
        elif cost_per_compound > 20000:
            return \"MODERATE COST: Consider filing provisionals first to defer costs while validating commercial potential.\"
        else:
            return \"REASONABLE COST: Proceed with patent filing for maximum IP protection.\"
```

**API endpoint for cost estimation:**

```python
@[router.post](http://router.post)(\"/api/patentpath/cost-estimate\")
async def estimate_costs(
    num_compounds: int,
    jurisdictions: List[str] = [\"US\"],
    claim_count_avg: int = 15,
    complexity: str = \"moderate\",
    strategy: str = \"patent\"
):
    \"\"\"Estimate IP protection costs.\"\"\"
    estimator = CostEstimator()
    result = estimator.estimate_costs(
        num_compounds=num_compounds,
        jurisdictions=jurisdictions,
        claim_count_avg=claim_count_avg,
        complexity=complexity,
        strategy=strategy
    )
    return result
```

---

**Task 8.3: PatentPath API Routes Integration** (3-4 hours)

Create `backend/api/routes/[patentpath.py](http://patentpath.py)` to consolidate all PatentPath endpoints:

```python
from fastapi import APIRouter, Depends, HTTPException, Query
from typing import List, Optional
from pydantic import BaseModel

from services.patentpath_lite.prior_art_search import PriorArtSearcher
from services.patentpath_lite.novelty_scorer import NoveltyScorer
from services.patentpath_lite.fto_checker import FTOChecker
from services.patentpath_lite.claim_generator import ClaimGenerator
from services.patentpath_[lite.tk](http://lite.tk)_checker import TKAttributionChecker
from services.patentpath_lite.cost_estimator import CostEstimator

router = APIRouter(prefix=\"/api/patentpath\", tags=[\"patentpath\"])

# Pydantic models for requests
class CompoundAnalysisRequest(BaseModel):
    compound_name: str
    compound_description: str
    compound_smiles: str
    keywords: List[str]
    therapeutic_use: Optional[str] = None
    key_features: Optional[List[str]] = None
    derivation_source: Optional[str] = None

@[router.post](http://router.post)(\"/analyze-compound\")
async def analyze_compound_ip(request: CompoundAnalysisRequest):
    \"\"\"Complete IP analysis for a compound (all 6 features).
    
    This is a convenience endpoint that runs all PatentPath Lite features
    in a single request.
    \"\"\"
    results = {}
    
    # Feature 1: Prior Art Search
    searcher = PriorArtSearcher()
    prior_art = [searcher.search](http://searcher.search)(
        compound_description=request.compound_description,
        keywords=request.keywords
    )
    results[\"prior_art\"] = prior_art
    
    # Feature 2: Novelty Scoring
    scorer = NoveltyScorer()
    novelty = scorer.calculate_novelty(
        compound_description=request.compound_description,
        prior_art_patents=prior_art.get(\"prior_art_patents\", [])
    )
    results[\"novelty\"] = novelty
    
    # Feature 3: FTO Check
    fto_checker = FTOChecker()
    fto = fto_checker.check_fto(
        compound_description=request.compound_description,
        keywords=request.keywords,
        intended_use=request.therapeutic_use
    )
    results[\"fto\"] = fto
    
    # Feature 4: Claim Generation
    generator = ClaimGenerator()
    claims = generator.generate_claims(
        compound_name=request.compound_name,
        compound_smiles=request.compound_smiles,
        therapeutic_use=request.therapeutic_use or \"therapeutic use\",
        key_features=request.key_features or [\"novel structure\"]
    )
    results[\"claims\"] = claims
    
    # Feature 5: TK Attribution Check
    tk_checker = TKAttributionChecker()
    tk_check = tk_checker.check_tk_attribution(
        compound_description=request.compound_description,
        derivation_source=request.derivation_source
    )
    results[\"tk_attribution\"] = tk_check
    
    # Feature 6: Cost Estimation
    estimator = CostEstimator()
    costs = estimator.estimate_costs(
        num_compounds=1,
        claim_count_avg=claims.get(\"total_claims\", 15)
    )
    results[\"cost_estimate\"] = costs
    
    # Overall recommendation
    results[\"overall_recommendation\"] = generate_overall_recommendation(
        novelty[\"recommendation\"],
        fto[\"fto_risk_level\"],
        tk_check[\"can_proceed\"]
    )
    
    return {
        \"compound_name\": request.compound_name,
        \"analysis_results\": results,
        \"summary\": generate_executive_summary(results)
    }

def generate_overall_recommendation(novelty_rec: str, fto_risk: str, tk_proceed: bool) -> str:
    \"\"\"Generate overall IP recommendation.\"\"\"
    if not tk_proceed:
        return \"DO NOT PROCEED: Sacred traditional knowledge detected.\"
    
    if novelty_rec == \"LIKELY_NOT_PATENTABLE\":
        return \"NOT RECOMMENDED: High prior art similarity - consider trade secret or design-around.\"
    
    if fto_risk == \"HIGH\":
        return \"CAUTION: Novel but blocking patents exist - licensing required before commercialization.\"
    
    if novelty_rec == \"LIKELY_PATENTABLE\" and fto_risk == \"LOW\":
        return \"RECOMMENDED: Proceed with patent filing - good novelty and clear FTO.\"
    
    return \"UNCERTAIN: Consult patent attorney for detailed analysis.\"

def generate_executive_summary(results: dict) -> dict:
    \"\"\"Generate executive summary of IP analysis.\"\"\"
    return {
        \"patentability\": results[\"novelty\"][\"recommendation\"],
        \"novelty_score\": results[\"novelty\"][\"novelty_score\"],
        \"fto_risk\": results[\"fto\"][\"fto_risk_level\"],
        \"blocking_patents\": results[\"fto\"][\"num_blocking_patents\"],
        \"tk_attribution_required\": results[\"tk_attribution\"][\"attribution_required\"],
        \"estimated_total_cost\": results[\"cost_estimate\"][\"total_cost_estimate\"][\"total_investment\"],
        \"recommendation\": results[\"overall_recommendation\"]
    }
```

---

### Frontend Development (20 hours)

**Task 8.4: Next.js Foundation** (8-10 hours)

**Task 8.5: Dashboard Layout** (4-5 hours)

**Task 8.6: Compound Search Interface** (4-5 hours)

**Task 8.7: PatentPath UI Integration** (4-5 hours)

(Frontend code examples omitted for brevity - follows standard Next.js patterns)

---

**Success Metrics - Week 8**:

- PatentPath Lite all 6 features complete
- TK attribution checker flags sacred knowledge
- Cost estimator compares patent vs. trade secret
- Frontend displays compound search results
- PatentPath analysis accessible via UI
- Bundle pricing: $4,500 per compound (all 6 features)

---

### PatentPath Lite Success Criteria - Complete System

**Functionality Checklist**:

- [ ]  Feature 1: Prior art search (USPTO, EPO, Google Patents)
- [ ]  Feature 2: Novelty scoring (0.0-1.0 scale with recommendations)
- [ ]  Feature 3: FTO checker (blocking patents, licensing costs)
- [ ]  Feature 4: Claim generator (USPTO-formatted templates)
- [ ]  Feature 5: TK attribution (sacred knowledge detection)
- [ ]  Feature 6: Cost estimator (patent vs. trade secret comparison)

**API Endpoints**:

- [ ]  `POST /api/patentpath/prior-art-search`
- [ ]  `POST /api/patentpath/novelty-assessment`
- [ ]  `POST /api/patentpath/fto-check`
- [ ]  `POST /api/patentpath/generate-claims`
- [ ]  `POST /api/patentpath/tk-check`
- [ ]  `POST /api/patentpath/cost-estimate`
- [ ]  `POST /api/patentpath/analyze-compound` (all features)

**Performance**:

- [ ]  Prior art search: <30 seconds
- [ ]  Novelty scoring: <5 seconds
- [ ]  FTO check: <45 seconds
- [ ]  Claim generation: <10 seconds
- [ ]  TK check: <2 seconds
- [ ]  Cost estimation: <1 second

**Testing**:

- [ ]  Unit tests for all 6 features
- [ ]  Integration tests for API endpoints
- [ ]  End-to-end test: full compound analysis
- [ ]  Test coverage >80%

**Documentation**:

- [ ]  API documentation with examples
- [ ]  User guide for PatentPath Lite
- [ ]  Pricing guide ($500-4,500 per compound)
- [ ]  Upgrade path to PatentPath Full

---

**Year 1 Revenue Projection - PatentPath Lite**:

- 50 prior art searches @ $500 = $25K
- 30 novelty assessments @ $1,000 = $30K
- 15 FTO checks @ $2,500 = $37.5K
- 20 claim generations @ $1,500 = $30K
- **Total: $122.5K Year 1**
- **ROI: 468%** (on $25K investment)

---

## Weeks 9-14 Summary

(Weeks 9-14 continue with frontend completion, production hardening, security, testing, and Nevada pilot launch as per original plan - no additional PatentPath work needed)

**Phase 3 includes**:

- Week 9: Complete frontend MVP
- Week 10: API optimization + caching
- Week 11: OmniPath full integration
- Week 12: Security hardening + compliance
- Week 13: Testing + bug fixes
- Week 14: GCP deployment + Nevada pilot

---

*Document complete. PatentPath Lite fully integrated into NeuroBotanica MVP development plan.*

---

# 🔐 Trade Secret Path Modules — Implementation Specs

<aside>
⚠️

**CONFIDENTIAL — PROPRIETARY TRADE SECRET — INTERNAL USE ONLY**

Do not distribute. Scoring weights, thresholds, decision trees, and exception lists are trade secrets.

</aside>

**Purpose:** Three internal report generators (ChemPath, ToxPath, RegPath) that power NeuroBotanica's manufacturer-facing regulatory readiness services. These modules expose capability via APIs while keeping scoring logic/weights/heuristics proprietary.

**IP Posture:** Consistent with OmniPath trade secret master—expose results, not methods.

**Stack:** FastAPI + PostgreSQL + SQLAlchemy/Pydantic (aligns with existing MVP architecture)

**Reference Tables:** Uses existing `cannabinoids` table; creates separate `chempath_*`, `toxpath_*`, `regpath_*` tables for audit trails and independent versioning.

---

## ChemPath — Chemical Characterization Engine

### Goal

Given (SMILES/InChI + optional COA fields + optional terpene profile), produce a manufacturer-facing **Molecular Characterization Report** plus machine-readable QC flags for downstream use in NeuroBotanica.

### Inputs (Pydantic Models)

```python
from pydantic import BaseModel
from typing import Optional, List, Dict

class CompoundInput(BaseModel):
    name: str
    smiles: Optional[str] = None
    inchi: Optional[str] = None
    inchi_key: Optional[str] = None
    source: Optional[str] = "customer_private"

class COAInput(BaseModel):
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

class ChemPathRequest(BaseModel):
    compound_input: CompoundInput
    coa_input: Optional[COAInput] = None
```

### Outputs

```python
class ChemPathResponse(BaseModel):
    chempath_job_id: str
    normalized_structure: Dict  # canonical_smiles, inchikey, normalization_notes
    computed_descriptors_2d: Dict  # RDKit descriptor dict
    computed_descriptors_3d: Optional[Dict] = None  # if conformers generated
    coa_qc_flags: List[Dict]  # [{code, severity, message}]
    report_markdown: str
```

### Database Schema

**Table: `chempath_jobs`**

```python
class ChemPathJob(Base):
    __tablename__ = "chempath_jobs"
    
    id = Column(UUID, primary_key=True, default=uuid4)
    created_at = Column(DateTime, default=datetime.utcnow)
    requested_by = Column(String)
    compound_name = Column(String)
    smiles_in = Column(String, nullable=True)
    inchi_in = Column(String, nullable=True)
    canonical_smiles = Column(String)
    inchikey = Column(String, nullable=True)
    status = Column(String)  # queued/running/succeeded/failed
    error = Column(String, nullable=True)
    inputs_json = Column(JSON)
    outputs_json = Column(JSON)
```

**Table: `chempath_coa_documents`**

```python
class ChemPathCOA(Base):
    __tablename__ = "chempath_coa_documents"
    
    id = Column(UUID, primary_key=True, default=uuid4)
    chempath_job_id = Column(UUID, ForeignKey("chempath_[jobs.id](http://jobs.id)"))
    lab_name = Column(String)
    sample_id = Column(String)
    units = Column(String)
    coa_json = Column(JSON)  # normalized COA fields
    qc_flags = Column(JSON)  # array of flags
```

**Table: `chempath_descriptors`**

```python
class ChemPathDescriptors(Base):
    __tablename__ = "chempath_descriptors"
    
    id = Column(UUID, primary_key=True, default=uuid4)
    chempath_job_id = Column(UUID, ForeignKey("chempath_[jobs.id](http://jobs.id)"))
    descriptors_2d = Column(JSON)
    descriptors_3d = Column(JSON, nullable=True)
    conformer_meta = Column(JSON, nullable=True)
```

### API Endpoints

```python
router = APIRouter(prefix="/api/chempath", tags=["chempath"])

@[router.post](http://router.post)("/analyze")
async def analyze_compound(request: ChemPathRequest) -> ChemPathResponse:
    """Run full ChemPath analysis pipeline."""
    pass

@router.get("/jobs/{job_id}")
async def get_job(job_id: str):
    """Get job status and stored outputs."""
    pass

@[router.post](http://router.post)("/report/{job_id}")
async def generate_report(job_id: str) -> str:
    """Generate report_markdown (later: PDF render)."""
    pass
```

### Processing Pipeline (MVP Order)

1. **Parse structure** → canonicalize → store
2. **Compute 2D descriptors**; optionally compute 3D conformers for whitelist (16 base cannabinoids + selected dimers)
3. **Normalize COA units** → run QC checks → produce flags
4. **Produce report markdown** from template

### QC Flag Codes (Examples — Implementation Secret)

| Code | Severity | Trigger |
| --- | --- | --- |
| `COA_TOTAL_EXCEEDS_100` | error | Sum of cannabinoids + terpenes > 100% |
| `COA_UNIT_MISMATCH` | warning | Mixed units in same COA |
| `COA_MISSING_LOD` | info | LOD/LOQ not provided |
| `STRUCTURE_INVALID` | error | SMILES/InChI parse failure |
| `STRUCTURE_SALT_NORMALIZED` | info | Salt form detected and normalized |

### Report Template (Markdown Sections)

```markdown
# Molecular Characterization Report
**Compound:** {name}
**Generated:** {timestamp}
**Job ID:** {chempath_job_id}

## 1. Identity & Structure
- **Canonical SMILES:** {canonical_smiles}
- **InChIKey:** {inchikey}
- **Normalization Notes:** {notes}

## 2. Physicochemical Descriptors
### 2D Properties
| Property | Value |
|----------|-------|
| Molecular Weight | {mw} |
| LogP | {logp} |
| TPSA | {tpsa} |
| H-Bond Donors | {hbd} |
| H-Bond Acceptors | {hba} |

### 3D Properties (if available)
{3d_descriptors_table}

## 3. Certificate of Analysis Summary
**Lab:** {lab_name}
**Sample ID:** {sample_id}
{coa_summary_table}

### QC Flags
{qc_flags_list}

### Data Completeness: {completeness_score}%

---
*For planning/compliance readiness only; lab confirmation required.*
```

---

## ToxPath — Toxicity Risk Assessment Engine

### Goal

Given a compound record (from ChemPath/NeuroBotanica), generate a conservative **Toxicity Risk Memo** and a "next testing" checklist appropriate for manufacturer regulatory readiness (not medical advice).

### Inputs

```python
class ToxPathRequest(BaseModel):
    compound_ref: str  # canonical_smiles or compound_id
    route: Optional[str] = None  # oral/inhaled/sublingual
    estimated_exposure: Optional[Dict] = None  # dose_mg, frequency, duration
    known_impurities: Optional[List[Dict]] = None  # [{name, ppm}]
```

### Outputs

```python
class ToxPathResponse(BaseModel):
    toxpath_assessment_id: str
    risk_summary: Dict  # {overall_tier, top_risks[], key_unknowns[]}
    alerts: List[Dict]  # structural/toxicology alert codes
    testing_plan: List[Dict]  # ordered steps with rationale
    memo_markdown: str
```

### Database Schema

**Table: `toxpath_assessments`**

```python
class ToxPathAssessment(Base):
    __tablename__ = "toxpath_assessments"
    
    id = Column(UUID, primary_key=True, default=uuid4)
    created_at = Column(DateTime, default=datetime.utcnow)
    compound_id = Column(UUID, nullable=True)  # FK to cannabinoids if exists
    canonical_smiles = Column(String)  # stored snapshot
    route = Column(String, nullable=True)
    exposure_json = Column(JSON, nullable=True)
    alerts = Column(JSON)  # array of alert codes
    risk_summary = Column(JSON)  # {overall_tier, top_risks[], key_unknowns[]}
    testing_plan = Column(JSON)  # ordered steps with rationale
    memo_markdown = Column(Text)
```

### API Endpoints

```python
router = APIRouter(prefix="/api/toxpath", tags=["toxpath"])

@[router.post](http://router.post)("/assess")
async def assess_toxicity(request: ToxPathRequest) -> ToxPathResponse:
    """Run ToxPath assessment pipeline."""
    pass

@router.get("/assessments/{assessment_id}")
async def get_assessment(assessment_id: str):
    """Get stored assessment result."""
    pass

@[router.post](http://router.post)("/memo/{assessment_id}")
async def generate_memo(assessment_id: str) -> str:
    """Generate memo_markdown."""
    pass
```

### Processing Pipeline (MVP Order)

1. **Pull compound structure** + descriptors (re-use ChemPath outputs if present)
2. **Run rule-based alert checks** (keep rules private; expose only alert codes + explanations)
3. **Route-aware risk framing** (different memo language/checklist by route)
4. **Generate testing plan** using tier matrix:
    - Higher novelty / higher uncertainty → longer checklist
    - Lower novelty / strong characterization → shorter checklist

### Risk Tier Matrix (Implementation Secret)

| Tier | Label | Criteria (Secret) | Testing Depth |
| --- | --- | --- | --- |
| 1 | Low | Well-characterized, existing safety data | Minimal |
| 2 | Moderate | Some characterization gaps | Standard panel |
| 3 | High | Novel structure, limited data | Extended + consult |
| 4 | Very High | Multiple structural alerts | Full GLP recommended |

### Memo Template (Markdown Sections)

```markdown
# Toxicity Risk Assessment
**Compound:** {name}
**Route:** {route}
**Generated:** {timestamp}
**Assessment ID:** {toxpath_assessment_id}

## 1. Compound Snapshot
- **Structure:** {canonical_smiles}
- **Key Properties:** {properties_summary}

## 2. Risk Assessment
### Overall Tier: {tier} — {tier_label}
**Key Assumptions:** {assumptions}

### Top Risks
{top_risks_list}

### Key Unknowns
{unknowns_list}

## 3. Structural Alerts
{alerts_table}

## 4. Recommended Testing Sequence
{testing_plan_checklist}

### Consultation Required
{consultation_flags}

---
*For regulatory readiness planning only. Not a substitute for GLP studies or professional toxicology consultation.*
```

---

## RegPath — Regulatory Pathway Optimization Engine

### Goal

Produce a manufacturer-facing **Regulatory Strategy Memo** answering: "What pathway is most plausible (IND/NDA/505(b)(2) etc.) and what's the timeline/checklist?"

### Inputs

```python
class ProductProfile(BaseModel):
    product_type: str  # isolated_cannabinoid/novel_dimer/formulation/combination
    intended_use: str
    route: str  # oral/inhaled/sublingual/topical
    target_markets: List[str] = ["US"]

class RegPathRequest(BaseModel):
    product_profile: ProductProfile
    evidence_inputs: Optional[Dict] = None  # chempath_job_id, toxpath_assessment_id
    constraints: Optional[Dict] = None  # budget, launch_window, risk_tolerance
```

### Outputs

```python
class RegPathResponse(BaseModel):
    regpath_strategy_id: str
    pathway_recommendation: Dict  # {primary, fallback, rationale[]}
    readiness_checklist: List[Dict]  # CMC/nonclinical/clinical/docs items
    timeline: List[Dict]  # milestone list with durations
    memo_markdown: str
```

### Database Schema

**Table: `regpath_strategies`**

```python
class RegPathStrategy(Base):
    __tablename__ = "regpath_strategies"
    
    id = Column(UUID, primary_key=True, default=uuid4)
    created_at = Column(DateTime, default=datetime.utcnow)
    product_profile = Column(JSON)
    inputs_refs = Column(JSON)  # chempath/toxpath pointers
    pathway_recommendation = Column(JSON)
    readiness_checklist = Column(JSON)
    timeline = Column(JSON)
    memo_markdown = Column(Text)
```

### API Endpoints

```python
router = APIRouter(prefix="/api/regpath", tags=["regpath"])

@[router.post](http://router.post)("/strategy")
async def generate_strategy(request: RegPathRequest) -> RegPathResponse:
    """Run RegPath strategy generation."""
    pass

@router.get("/strategies/{strategy_id}")
async def get_strategy(strategy_id: str):
    """Get stored strategy result."""
    pass

@[router.post](http://router.post)("/memo/{strategy_id}")
async def generate_memo(strategy_id: str) -> str:
    """Generate memo_markdown."""
    pass
```

### Processing Pipeline (MVP Order)

1. **Classify product** (novelty level + "how isolated/characterized is it?")
2. **Pick primary vs fallback pathway** via decision matrix (keep matrix private; output rationale bullets)
3. **Build readiness checklist** that pulls ChemPath + ToxPath highlights if provided
4. **Generate timeline** from editable templates

### Pathway Decision Matrix (Implementation Secret)

| Product Type | Novelty | Evidence | Primary | Fallback |
| --- | --- | --- | --- | --- |
| Isolated cannabinoid | Low | Strong | 505(b)(2) | NDA |
| Isolated cannabinoid | High | Weak | IND | — |
| Novel dimer | Any | Any | IND | 505(b)(2) |
| Formulation | Low | Strong | ANDA/505(b)(2) | NDA |
| Combination | Any | Mixed | 505(b)(2) | IND |

### Memo Template (Markdown Sections)

```markdown
# Regulatory Strategy Memo
**Product:** {product_name}
**Type:** {product_type}
**Route:** {route}
**Target Market:** {markets}
**Generated:** {timestamp}
**Strategy ID:** {regpath_strategy_id}

## 1. Product Snapshot
{product_summary}
**Key Assumptions:** {assumptions}

## 2. Pathway Recommendation
### Primary: {primary_pathway}
{primary_rationale}

### Fallback: {fallback_pathway}
{fallback_rationale}

### Critical Gating Questions
{gating_questions}

## 3. Readiness Checklist
### CMC
{cmc_checklist}

### Nonclinical
{nonclinical_checklist}

### Clinical
{clinical_checklist}

### Documentation
{docs_checklist}

## 4. Estimated Timeline
{timeline_table}

## 5. Next Actions
{next_actions}

---
*Timeline estimates are based on typical pathways and may vary. Consult regulatory counsel for binding guidance.*
```

---

## Integration Architecture

```
┌─────────────────────────────────────────────────────────────────────┐
│                    NeuroBotanica Pipeline                          │
├─────────────────────────────────────────────────────────────────────┤
│                                                                     │
│  ┌──────────────┐    ┌──────────────┐    ┌──────────────┐         │
│  │   ChemPath   │───▶│   ToxPath    │───▶│   RegPath    │         │
│  │ (Structure + │    │ (Risk Memo + │    │ (Pathway +   │         │
│  │  CoA QC)     │    │  Test Plan)  │    │  Timeline)   │         │
│  └──────────────┘    └──────────────┘    └──────────────┘         │
│         │                   │                   │                  │
│         ▼                   ▼                   ▼                  │
│  ┌────────────────────────────────────────────────────────┐       │
│  │         Manufacturer Regulatory Readiness Report        │       │
│  │  (Compound Profile + Tox Assessment + Strategy Memo)   │       │
│  └────────────────────────────────────────────────────────┘       │
│                                                                     │
└─────────────────────────────────────────────────────────────────────┘

API Boundary: Results exposed, methods protected as trade secrets
```

### Service Offer Mapping

| Service Tier | ChemPath | ToxPath | RegPath | Price |
| --- | --- | --- | --- | --- |
| Basic Manufacturer Scan | ✅ Report | Risk tier only | — | $500 |
| Comprehensive Assessment | ✅ Full | ✅ Full memo | Pathway rec | $2,500 |
| Enterprise Package | ✅ Full | ✅ Full + GLP consult | ✅ Full + timeline | $7,500+ |

### Development Timeline Integration

| Week | Path Module Work | Hours |
| --- | --- | --- |
| Week 2 | ChemPath database schema + basic API | 8 |
| Week 3 | ChemPath QC rules + report generator | 10 |
| Week 4 | ToxPath schema + alert engine | 10 |
| Week 5 | ToxPath testing plan generator | 8 |
| Week 6 | RegPath schema + pathway selector | 10 |
| Week 7 | RegPath timeline builder + integration | 8 |
| Week 8 | Full pipeline test + manufacturer workflow | 6 |
| **Total** |  | **60 hours** |

```

```