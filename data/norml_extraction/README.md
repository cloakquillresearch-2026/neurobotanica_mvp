# NORML Manual Extraction Workflow

**Status**: ✅ ACTIVE  
**Start Date**: December 17, 2025  
**Target Completion**: January 15, 2026  
**Progress**: 10/200 studies (5%)

---

## 🎯 Mission

Extract **200+ clinical studies** from NORML's medical cannabis research library to build an FDA-compliant therapeutic validation dataset for **NeuroBotanica AI platform**.

**Market Catalyst**: Trump Executive Order Cannabis Rescheduling (Dec 17, 2025) creates immediate demand for clinical evidence.

---

## 📋 Quick Start Guide

### Step 1: Choose Your Condition

**Week 1** (Dec 17-23): Chronic Pain  
**Week 2** (Dec 24-30): Anxiety  
**Week 3** (Dec 31-Jan 6): Depression  
**Week 4** (Jan 7-13): Arthritis

### Step 2: Access NORML Page

**Chronic Pain**: https://norml.org/marijuana/library/recent-medical-marijuana-research/chronic-pain/  
**Anxiety**: https://norml.org/marijuana/library/recent-medical-marijuana-research/anxiety/  
**Depression**: https://norml.org/marijuana/library/recent-medical-marijuana-research/depression/  
**Arthritis**: https://norml.org/marijuana/library/recent-medical-marijuana-research/rheumatoid-arthritis/

### Step 3: Extract Studies

For each reference on the NORML page, fill out this JSON structure:

```json
{
  "study_id": "[CONDITION]_[TYPE]_[NUMBER]",
  "study_type": "RCT",
  "condition": "CHRONIC_PAIN",
  "study_title": "[exact title from citation]",
  "citation": "[author et al. year. journal vol: pages]",
  "publication_year": 2023,
  "journal": "[Journal Name]",
  
  "intervention": {
    "cannabis_type": "vaporized cannabis",
    "cannabinoid_profile": "THC-dominant",
    "delivery_method": "vaporized",
    "dosing_information": "low dose protocol",
    "treatment_duration": "4 weeks"
  },
  
  "outcomes": {
    "key_findings": [
      "Primary finding from study",
      "Secondary finding",
      "Additional result"
    ],
    "outcome_measures": ["VAS pain scale", "PHQ-9"],
    "quantitative_results": [],
    "adverse_effects": ["mild drowsiness", "dry mouth"],
    "effect_size_category": "medium"
  },
  
  "study_quality": {
    "sample_size": 120,
    "institution": "University Name",
    "country": "United States",
    "randomized": true,
    "placebo_controlled": true
  },
  
  "notes": "Extracted from NORML [condition] research library"
}
```

### Step 4: Add to JSON File

Add your study to the appropriate file:
- `data/norml_extraction/chronic_pain_studies.json`
- `data/norml_extraction/anxiety_studies.json`
- `data/norml_extraction/depression_studies.json`
- `data/norml_extraction/arthritis_studies.json`

### Step 5: Validate

Run the validator:
```bash
python scripts/validate_norml_extraction.py
```

---

## 🔍 What to Extract from NORML

### From the Citation
- **Study title** (exact)
- **Authors** (first author et al.)
- **Year**
- **Journal name**
- **Volume and pages**

### From the Text
- **Study type**: RCT, observational, case series
- **Sample size**: Look for "N = XX patients"
- **Cannabis type**: CBD, THC, whole plant, specific cannabinoid
- **Delivery method**: oral, inhaled, vaporized, sublingual, topical
- **Key findings**: Main results paragraph
- **Outcomes**: Pain scales, symptom measures
- **Adverse effects**: Side effects mentioned

---

## 📊 Quality Standards

### Required Fields (Must Have)
- study_id
- study_type
- condition
- study_title
- citation
- publication_year
- journal

### Highly Recommended
- sample_size
- key_findings (at least 1)
- delivery_method
- country

### Optional (If Available)
- cannabinoid_profile
- dosing_information
- quantitative_results
- institution

---

## 🎓 Study Type Classification

**RCT** (Randomized Controlled Trial):
- Keywords: "randomized", "placebo-controlled", "double-blind"
- Gold standard for efficacy

**Observational**:
- Keywords: "observational", "cohort", "registry", "real-world"
- Shows real-world effectiveness

**Case Series**:
- Keywords: "case series", "case report", "small sample"
- Pilot data, exploratory

---

## ✅ Validation Checklist

Before moving to next study, verify:

- [ ] Study ID follows format: `[CONDITION]_[TYPE]_[NUMBER]`
- [ ] All required fields present
- [ ] Citation format matches NORML reference
- [ ] Year is 4 digits (not 20XX)
- [ ] Key findings are substantive (not just "showed results")
- [ ] JSON syntax is valid (check commas, brackets)

---

## 🚀 Productivity Tips

### Time Management
- **Target**: 8 studies/day
- **Time per study**: ~30 minutes
- **Daily quota**: 4 hours of focused extraction

### Batch Processing
1. **Morning**: Read NORML page, identify 10 studies
2. **Midday**: Extract 5 studies
3. **Afternoon**: Extract 5 studies
4. **Evening**: Validate and commit

### Quality Over Speed
- **Don't rush**: Accurate extraction > quantity
- **When stuck**: Leave field blank, add note
- **When uncertain**: Mark study for review

---

## 📈 Progress Tracking

### Daily Log Template

```markdown
## December 17, 2025

**Condition**: Chronic Pain  
**Studies Extracted**: 10  
**Total Time**: 5 hours  
**Challenges**: Finding full citations  
**Notes**: Focused on HIV neuropathy studies

**Next Session**:
- Continue with opioid-sparing studies
- Extract references #15-25
```

---

## 🛠️ Tools & Resources

### Validation Script
```bash
# Validate current extraction
python scripts/validate_norml_extraction.py

# Merge all phases
python scripts/validate_norml_extraction.py --merge
```

### File Locations
- **Extractions**: `data/norml_extraction/`
- **Scripts**: `scripts/`
- **Status**: `data/norml_extraction/EXTRACTION_STATUS.md`

---

## 🎯 Success Metrics

### Week 1 Goals (Chronic Pain)
- [ ] 40-50 studies extracted
- [ ] 100% validation pass rate
- [ ] All major subconditions covered (neuropathy, cancer, fibromyalgia)

### Week 2 Goals (Anxiety)
- [ ] 30-40 studies extracted
- [ ] Diverse study types (RCT, observational)
- [ ] Multiple cannabinoid profiles

### Week 3 Goals (Depression)
- [ ] 30-40 studies extracted
- [ ] Comorbidity studies included
- [ ] Adolescent and adult populations

### Week 4 Goals (Arthritis)
- [ ] 20-30 studies extracted
- [ ] Topical delivery methods highlighted
- [ ] Inflammatory marker data captured

---

## 💪 Staying Motivated

**Remember**:
- 200+ studies = **FDA-ready dataset**
- Each study = **Commercial value** for dispensaries
- Trump EO = **Market timing** is perfect
- Beta partners = **Waiting for this data**

**Impact**:
- VeriTrad → $5,000 MRR target
- NeuroBotanica → $6,000 MRR target
- Combined → **$132,000 annual revenue** potential

---

## 🆘 Troubleshooting

### "Can't find full citation"
→ Use what's available in NORML reference, mark for review

### "Unclear cannabinoid profile"
→ Use "not specified" or "varied"

### "No quantitative results"
→ Leave array empty `[]`, focus on qualitative findings

### "Duplicate study across conditions"
→ Extract once, note cross-relevance in notes field

---

## 📞 Questions?

**Status Updates**: Check `EXTRACTION_STATUS.md`  
**Validation Issues**: Run `validate_norml_extraction.py`  
**JSON Syntax Help**: Use online JSON validator

---

**Let's extract! The cannabis industry needs this data NOW. 🚀**

---

**Last Updated**: December 17, 2025 - Start of extraction campaign
