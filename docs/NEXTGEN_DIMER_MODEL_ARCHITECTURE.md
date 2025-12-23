# Next-Gen Dimer/Entourage Model Architecture Plan

## 1. Schema Foundation
- Uses `NextGenDimerEvidence` dataclass (see nextgen_dimer_schema.py)
- Supports omics, pathway, AI, patient stratification, time series, and regulatory fields

## 2. Data Ingestion Pipeline
- Ingests from unified next-gen JSON (e.g., data/processed/nextgen_dimers.json)
- Validates and normalizes all new fields (omics, pathway, population, time series)
- Supports batch and streaming ingestion for real-world/longitudinal data

## 3. Model Architecture (Draft)

### a. Input Modules
- **Compound Embeddings:**
  - SMILES-based molecular embeddings (RDKit, ECFP4, ChemBERTa)
  - Pathway/omics signature embeddings (if available)
- **Contextual Features:**
  - Patient stratification (age, sex, genotype, comorbidities)
  - Regulatory status, benefit-risk ratio
  - Time series (for real-world evidence)

### b. Core Model
- **Multi-Modal Fusion:**
  - Concatenate/fuse compound, omics, pathway, and patient features
  - Attention layers for context weighting (e.g., transformer or cross-attention)
- **Prediction Heads:**
  - Synergy/antagonism classification
  - Effect size regression
  - Adverse effect prediction
  - Population-specific efficacy

### c. Training & Validation
- **Loss Functions:**
  - Multi-task: classification (synergy/antagonism), regression (effect size, benefit-risk)
  - Weighted by evidence confidence, population size, and regulatory status
- **Evaluation:**
  - Stratified by population, evidence type, and time series
  - Real-world/omics validation sets

### d. Output
- **Explainable AI:**
  - Feature attribution for regulatory/clinical transparency
  - Pathway and omics signature reporting
- **API/Integration:**
  - Ready for edge deployment and clinical/dispensary integration

---

*This plan is the foundation for Q1 2026 next-gen dimer/entourage model development.*
