# Phase 4 Data Expansion & Model Roadmap

## ✅ What Has Been Accomplished

- **NORML Dataset Extraction**
  - 320+ studies extracted, covering 17+ indications (chronic pain, cancer, anxiety, epilepsy, etc.)
  - Data quality: 45% clinical trials, 68 high-quality RCTs, 12 meta-analyses, 24 systematic reviews
  - Individual and consolidated JSON files generated and validated

- **Data Validation & Integration**
  - All NORML studies robustly ingested and validated (359 valid, 9 invalid skipped)
  - Merged into main training dataset (`neurobotanica_enriched_dataset.json`)
  - No PubMed/EuropePMC automation used due to technical issues—NORML is the current foundation

- **Model Retraining**
  - Retrained Therapeutic Prediction and Patient Response models using only validated NORML data
  - Dimer Potential Model skipped (insufficient data)
  - All models and training reports saved

- **Trade Secret & Patent Protection**
  - Trade secret documentation and legal framework established
  - Patent claims updated for chemical-genomic integration

---

## 🔜 What Needs To Be Done Next

### 1. Data Expansion (Highest Priority)
- Fill gaps for high-priority conditions (Opioid Use Disorder, Fibromyalgia, Autism, Migraine, Menopause, ALS, Spinal Cord Injury, HIV/AIDS, Endometriosis)
- Reinforce thin areas (PTSD, Epilepsy, Tourette, Glaucoma)
- Target: 500+ studies, 24+ conditions, 50%+ RCTs
- Continue manual/curated evidence extraction if automation is unreliable

### 2. Adjuvant & Entourage Data Integration
- Build and validate adjuvant (terpene, flavonoid, lipid, alkaloid, polyphenol) mechanism database
- Integrate adjuvant evidence into training data (bioavailability, potentiation, side effect mitigation)

### 3. Incremental Model Layering
- After NORML+adjuvant base is robust, incrementally add:
  - GenomePath (genomic correlations)
  - ChemPath (chemical similarity/structure)
  - BioPath (bias-corrected efficacy validation)
  - ClinPath (adaptive clinical trial optimization)
- Each layer must be validated and benchmarked before full integration

### 4. Documentation & Reporting
- Maintain clear documentation of data sources, validation steps, and model changes
- Update roadmap and progress reports after each major milestone

---

## 📈 Immediate Next Steps

1. Gap analysis: Identify missing/underrepresented conditions in current dataset
2. Curate or extract additional studies for priority conditions
3. Validate and merge new evidence into the training set
4. Retrain and benchmark models after each data expansion
5. Begin adjuvant integration planning and schema design

---

## 🧬 Dimer Potential: Gap & Technical Integration Plan

### Current Status
---

## 🛠️ Dimer Schema & Integration: Implementation Steps
---

## 🧩 Dimer Integration: Python Code Templates

### 1. Dimer Schema Definition (Python Dataclass)
```python
from dataclasses import dataclass
from typing import Optional, Tuple

@dataclass
class DimerEvidence:
  compound_1: str
  compound_2: str
  dimer_type: str  # e.g., 'cannabinoid-cannabinoid', 'cannabinoid-terpene'
  interaction_effect: str  # 'synergistic', 'antagonistic', 'additive', 'unknown'
  evidence_type: str  # 'clinical', 'preclinical', 'in silico', 'prediction'
  effect_size: Optional[str] = None
  study_type: Optional[str] = None
  source: Optional[str] = None
  confidence: Optional[str] = None
  notes: Optional[str] = None
```

### 2. Dimer Evidence Validation Function
```python
def validate_dimer_entry(entry: dict) -> bool:
  required = ['compound_1', 'compound_2', 'dimer_type', 'interaction_effect', 'evidence_type']
  for field in required:
    if field not in entry or not entry[field]:
      return False
  return True
```

### 3. Dimer Extraction Example (from JSON)
```python
import json
from typing import List

def load_dimer_evidence(json_path: str) -> List[DimerEvidence]:
  with open(json_path, 'r') as f:
    data = json.load(f)
  dimers = []
  for entry in data:
    if validate_dimer_entry(entry):
      dimers.append(DimerEvidence(**entry))
  return dimers
```

### 4. Integration into Training Dataset
```python
def merge_dimeric_evidence(training_data: list, dimeric_data: list) -> list:
  # Assumes both are lists of dicts or dataclasses
  # Optionally deduplicate by (compound_1, compound_2, condition, source)
  merged = training_data.copy()
  for dimer in dimeric_data:
    if dimer not in merged:
      merged.append(dimer)
  return merged
```

### 5. Model Pipeline Update (Feature Example)
```python
def dimer_to_features(dimer: DimerEvidence) -> dict:
  # Example: one-hot encode dimer_type, interaction_effect, etc.
  features = {
    'compound_1': dimer.compound_1,
    'compound_2': dimer.compound_2,
    'dimer_type': dimer.dimer_type,
    'interaction_effect': dimer.interaction_effect,
    # Add more feature engineering as needed
  }
  return features
```

---

### Step 1: Schema Extension
- Define new fields in the evidence/training schema:
  - `compound_1`, `compound_2` (or `dimer_pair`)
  - `dimer_type` (e.g., cannabinoid-cannabinoid, cannabinoid-terpene)
  - `interaction_effect` (synergistic, antagonistic, additive, unknown)
  - `evidence_type` (clinical, preclinical, in silico, prediction)
  - Retain `source`, `confidence`, `notes`, `effect_size`, `study_type`
- Update data validation scripts to enforce new schema for dimeric entries

### Step 2: Data Extraction & Curation
- Systematically review literature and existing datasets for studies on compound pairings
- Extract, structure, and validate dimeric evidence (clinical, preclinical, in silico)
- Cross-reference and validate dimeric predictions with available evidence
- Tag and log all dimeric entries for traceability

### Step 3: Data Integration
- Merge validated dimeric evidence into the main training dataset
- Ensure deduplication and clear separation of monomeric vs dimeric entries
- Update data loaders and preprocessing scripts to handle dimeric features

### Step 4: Model Pipeline Update
- Extend model input pipeline to accept dimeric features
- Retrain or fine-tune TherapeuticPredictionModel to learn pairwise effects
- Benchmark dimeric model performance against monomeric baseline
- If dimeric data is sparse, develop ensemble or hybrid models

### Step 5: Documentation & Trade Secret Management
- Document all schema changes, extraction logic, and integration steps
- Maintain dimer/entourage logic as internal trade secret until public disclosure is strategic
- Prepare technical documentation for future patent filings and IP defense

### Step 6: Roadmap & Milestones
- Q1 2026: Complete schema extension, begin extraction/curation, validate initial dimeric entries
- Q2 2026: Integrate dimeric data, retrain/benchmark models, document trade secret logic
- Ongoing: Expand dimer/entourage evidence, update models, and maintain documentation

---

- Schema and validation pipelines must be extended to support dimeric (and higher-order) relationships.

### Technical Plan for Dimer Schema & Integration

**1. Schema Design**
  - Add new fields to evidence and training data:
    - `compound_1`, `compound_2` (or `dimer_pair` as a tuple/list)
    - `dimer_type` (e.g., cannabinoid-cannabinoid, cannabinoid-terpene)
    - `interaction_effect` (synergistic, antagonistic, additive, unknown)
    - `evidence_type` (clinical, preclinical, in silico, prediction)
    - `source`, `confidence`, `notes`, `effect_size`, `study_type` (as for monomers)
  - Update model input pipeline to accept dimeric features.

**2. Data Extraction & Validation**
  - Systematically extract and validate studies reporting on specific compound pairings (dimers) from literature and existing datasets.
  - Cross-reference dimeric predictions with available clinical/preclinical evidence for validation.
  - Tag and log all dimeric evidence for traceability and future audit.

**3. Model Integration**
  - Extend the TherapeuticPredictionModel to accept dimeric features and learn pairwise effects.
  - Benchmark dimeric model predictions against monomeric baselines and available evidence.
  - Develop ensemble or hybrid models if dimeric data is sparse.

**4. Trade Secret & Patent Alignment**
  - Keep dimer/entourage schema and extraction logic as internal trade secret until public disclosure is strategically advantageous.
  - Document all schema changes and integration steps for future patent filings and IP defense.

**5. Roadmap**
  - Q1 2026: Schema design, initial extraction, and validation of dimeric evidence
  - Q2 2026: Model integration, benchmarking, and trade secret documentation
  - Ongoing: Continuous expansion and validation of dimer/entourage evidence

---
