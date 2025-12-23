# Phase 4 Data Expansion & Model Roadmap

> **Note:** This summary has been superseded by the consolidated [PHASE_4_MASTER_BENCHMARK_REPORT.md](PHASE_4_MASTER_BENCHMARK_REPORT.md). Please refer to the master report for the latest benchmarks, integration status, and technical accomplishments.


## ✅ What Has Been Accomplished (as of 2025-12-23)

- **NORML Dataset Extraction & Validation**
  - 368 studies extracted, covering 22+ conditions (chronic pain, cancer, anxiety, epilepsy, ALS, autism, etc.)
  - 100% validation rate for all priority conditions; all studies curated and validated
  - Individual and consolidated JSON files generated, validated, and merged

- **Adjuvant/Entourage Data Structuring**
  - Terpene, flavonoid, and adjuvant synergy data integrated (9,884 synergies, evidence tiers)
  - Output: `neurobotanica_enriched_synergies.json` (353 compounds with structured entourage data)

- **Evidence Integration & Model Retraining**
  - All new evidence and synergy data merged into the main training set
  - Retrained all models: Therapeutic Prediction, Dimer Potential, Patient Response
  - Training report and model artifacts updated

- **Benchmarking & Validation of New Model Layers**
  - GenomePath, ChemPath, BioPath, ClinPath, and TKPath layers fully benchmarked (100% test pass)
  - MetaPath orchestration and ethical pricing logic validated (100% test pass)

- **Documentation & Roadmap Updates**
  - PHASE_4_DATA_EXPANSION_SUMMARY.md and technical roadmap updated after each expansion
  - All schema, validation, and integration steps documented for IP/trade secret protection

- **Trade Secret & Patent Protection**
  - Trade secret documentation and legal framework established
  - Patent claims updated for chemical-genomic and dimer/entourage integration

---


## 🚀 Current Status & Next Steps

- All evidence curation, validation, and integration steps for Phase 4 are complete
- All model layers (including GenomePath, ChemPath, BioPath, ClinPath, TKPath, MetaPath) are benchmarked and validated
- Documentation and technical roadmap are up to date

**Next Milestones:**
- Continue periodic evidence expansion as new studies emerge
- Monitor model performance and retrain as needed
- Prepare for Q1 2026: dimer/entourage schema extension, new evidence curation, and next-generation model integration

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
