# Clinical Evidence Expansion Guide
**Goal:** Expand from 3 positive samples (1.5%) to 200+ positive samples (40-50%)

---

## Priority: Systematic Reviews & Meta-Analyses

### Why Start Here?
- **Highest quality evidence** - synthesize multiple studies
- **Effect sizes clearly stated** - easier to extract large/medium/small ratings
- **Fewer studies to review** - one meta-analysis covers many RCTs
- **Higher confidence** - established clinical consensus

### Where to Search

**1. PubMed (FREE)**
- Go to: https://pubmed.ncbi.nlm.nih.gov/
- Use advanced search with filters:
  - Publication Type: Meta-Analysis, Systematic Review
  - Publication Date: Last 10 years (2015-2025)

**2. Cochrane Database (FREE)**
- Go to: https://www.cochranelibrary.com/
- Highest quality systematic reviews
- Cannabis therapeutics section: https://www.cochranelibrary.com/search?q=cannabis

**3. Google Scholar (FREE)**
- Go to: https://scholar.google.com/
- Search: "[compound] [condition] meta-analysis"
- Filter by date: Since 2015

---

## Search Strategy by Compound

### THC (Delta-9-Tetrahydrocannabinol)
**Current:** 246 studies in NORML dataset
**Target:** Extract effect sizes for 30+ conditions

**Priority Searches:**
```
1. "THC chronic pain meta-analysis"
2. "dronabinol cancer pain systematic review"
3. "THC PTSD randomized controlled trial"
4. "THC nausea chemotherapy meta-analysis"
5. "THC multiple sclerosis spasticity Cochrane"
6. "THC sleep disorders systematic review"
7. "THC Alzheimer's cognitive function RCT"
8. "THC Parkinson's tremor clinical trial"
```

**Expected Yield:** 20-30 high-quality studies → 15-20 positive samples

---

### CBD (Cannabidiol)
**Current:** 282 studies in NORML dataset  
**Target:** Extract effect sizes for 30+ conditions

**Priority Searches:**
```
1. "CBD epilepsy meta-analysis" (GOLD STANDARD - Epidiolex FDA approved)
2. "cannabidiol anxiety systematic review"
3. "CBD Dravet syndrome clinical trial"
4. "CBD Lennox-Gastaut syndrome RCT"
5. "CBD social anxiety disorder meta-analysis"
6. "CBD schizophrenia psychosis systematic review"
7. "CBD inflammatory bowel disease RCT"
8. "CBD addiction treatment meta-analysis"
9. "CBD sleep quality systematic review"
10. "CBD autism spectrum disorder clinical trial"
```

**Expected Yield:** 25-35 high-quality studies → 20-25 positive samples

---

### CBN (Cannabinol)
**Current:** Minimal studies
**Target:** 5-10 positive samples

**Priority Searches:**
```
1. "CBN sleep insomnia clinical trial"
2. "cannabinol sedative systematic review"
3. "CBN antibacterial activity study"
```

**Expected Yield:** 5-10 studies → 3-5 positive samples

---

### CBG (Cannabigerol)
**Current:** Minimal studies
**Target:** 5-10 positive samples

**Priority Searches:**
```
1. "CBG inflammatory bowel disease study"
2. "cannabigerol glaucoma intraocular pressure"
3. "CBG antibacterial MRSA study"
4. "CBG neuroprotection study"
```

**Expected Yield:** 5-10 studies → 3-5 positive samples

---

## How to Extract Effect Sizes

### From Meta-Analyses
Look for these terms in abstract/conclusions:

**LARGE effect:**
- "Strong evidence for efficacy"
- "Clinically significant improvement"
- "Number needed to treat (NNT) < 5"
- "Effect size d > 0.8" or "SMD > 0.8"
- ">50% symptom reduction"
- "FDA approval based on..."

**MEDIUM effect:**
- "Moderate evidence"
- "Significant clinical benefit"
- "NNT 5-10"
- "Effect size d 0.5-0.8"
- "30-50% improvement"

**SMALL effect:**
- "Modest benefit"
- "Statistically significant but..."
- "NNT > 10"
- "Effect size d 0.2-0.5"
- "10-30% improvement"

**MINIMAL/NONE:**
- "No significant difference"
- "Insufficient evidence"
- "Not clinically meaningful"
- "p > 0.05"

---

## Evidence Entry Workflow

### Step 1: Find Study
- Search PubMed/Cochrane
- Filter for systematic reviews/meta-analyses first
- Download PDF or save PubMed ID

### Step 2: Extract Key Data
- Compound name
- Condition/indication
- Study type (meta-analysis/systematic review/RCT)
- Year published
- Number of participants (total across all studies)
- Effect size rating (large/medium/small/minimal/none)
- Confidence level (high/medium/low)

### Step 3: Record Evidence
```python
# Run this to create a template
python scripts/expand_clinical_evidence.py template THC "chronic pain"
```

Then fill in the template and save to:
`data/clinical_evidence/thc_chronic_pain_20251223.md`

### Step 4: Import to Dataset
```python
# Batch import all new evidence files
python scripts/import_evidence.py
```

This will:
- Parse all .md files in `data/clinical_evidence/`
- Add to `neurobotanica_fully_enriched_fixed.json`
- Update compound `clinical_studies` arrays
- Regenerate training dataset with new samples

---

## Quality Control Checklist

For each evidence entry, verify:

- [ ] **Study type identified** (prefer meta-analysis > systematic review > RCT)
- [ ] **Effect size clearly categorized** (large/medium/small/minimal/none)
- [ ] **Source documented** (PubMed ID or DOI)
- [ ] **Participant count noted** (for confidence assessment)
- [ ] **Confidence level assigned** (based on study quality)
- [ ] **Compound name standardized** (THC not Delta-9-THC)
- [ ] **Condition name matches targets** (check therapeutic_targets in dataset)

---

## Target Milestones

### Week 1 (Dec 23-29, 2025)
- [ ] THC: 20 high-quality studies extracted → 15 positive samples
- [ ] CBD: 25 high-quality studies extracted → 20 positive samples
- **Target:** 35 positive samples (vs current 3)

### Week 2 (Dec 30 - Jan 5, 2026)
- [ ] THC: Additional 15 studies → 10 positive samples
- [ ] CBD: Additional 15 studies → 12 positive samples
- [ ] CBN: 5 studies → 3 positive samples
- [ ] CBG: 5 studies → 3 positive samples
- **Target:** 63 total positive samples

### Week 3 (Jan 6-12, 2026)
- [ ] Expand to THCV, CBDA, THCA, CBC
- [ ] Fill gaps in existing conditions
- [ ] Target rare conditions with strong evidence
- **Target:** 100+ total positive samples

### Week 4 (Jan 13-19, 2026)
- [ ] Final push to 200 positive samples
- [ ] Quality review and validation
- [ ] Retrain binary classifier
- **Target:** 200+ positive samples, 40-50% positive rate

---

## Resources

**Free Databases:**
- PubMed: https://pubmed.ncbi.nlm.nih.gov/
- Cochrane Library: https://www.cochranelibrary.com/
- Google Scholar: https://scholar.google.com/
- ClinicalTrials.gov: https://clinicaltrials.gov/

**Cannabis Research Databases:**
- Project CBD: https://www.projectcbd.org/science/cannabis-pharmacology
- NORML: https://norml.org/marijuana/library/recent-medical-marijuana-research/
- International Association for Cannabinoid Medicines: https://www.cannabis-med.org/

**Effect Size Calculators:**
- Effect Size Calculator: https://www.psychometrica.de/effect_size.html
- NNT Calculator: https://www.thennt.com/

---

## Tips for Efficient Extraction

1. **Start with Cochrane reviews** - highest quality, clear conclusions
2. **Focus on FDA-approved indications** - guaranteed large effect sizes (e.g., CBD for epilepsy)
3. **Look for "Table 1" or "Summary of Findings"** - quick effect size extraction
4. **Abstract is usually enough** - effect sizes stated in conclusions
5. **Skip negative studies initially** - focus on positive evidence first
6. **Batch similar conditions** - search all pain studies together
7. **Use citation networks** - one good review cites others
8. **Check supplementary materials** - often have detailed effect sizes

---

## Common Pitfalls to Avoid

❌ **Don't:**
- Mix preclinical (animal) studies with clinical evidence
- Include case reports as "observational studies" (need n>10)
- Rate statistically significant as "large" - must be clinically meaningful
- Duplicate evidence from same meta-analysis
- Include in-vitro studies

✅ **Do:**
- Focus on human clinical trials only
- Require clinical significance for large/medium ratings
- Cross-reference meta-analyses to avoid duplication
- Prioritize recent studies (2015-2025)
- Document source clearly for verification

---

## Questions? Issues?

1. **Unclear effect size rating?** → Default to "small" and add note
2. **Conflicting meta-analyses?** → Use most recent with largest sample
3. **Missing participant count?** → Estimate from study type (meta-analysis ~1000, RCT ~100)
4. **Study paywall?** → Search PubMed Central for free full text or use abstract only

---

**Next Step:** Run `python scripts/expand_clinical_evidence.py template [COMPOUND] [CONDITION]` to start!
