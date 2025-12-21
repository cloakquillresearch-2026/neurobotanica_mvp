# Quick Reference: Study Extraction Template

**Print this for reference while extracting studies**

---

## 📋 JSON Template (Copy/Paste)

```json
{
  "study_id": "CHRONIC_PAIN_RCT_011",
  "study_type": "RCT",
  "condition": "CHRONIC_PAIN",
  "study_title": "",
  "citation": "",
  "publication_year": 2023,
  "journal": "",
  
  "intervention": {
    "cannabis_type": "",
    "cannabinoid_profile": "",
    "delivery_method": "",
    "dosing_information": "",
    "treatment_duration": ""
  },
  
  "outcomes": {
    "key_findings": [
      "",
      "",
      ""
    ],
    "outcome_measures": [""],
    "quantitative_results": [],
    "adverse_effects": [""],
    "effect_size_category": "medium"
  },
  
  "study_quality": {
    "sample_size": null,
    "institution": "",
    "country": "",
    "randomized": true,
    "placebo_controlled": true
  },
  
  "notes": "Extracted from NORML chronic pain research library"
}
```

---

## 🔑 Field Definitions

### Study ID Format
```
[CONDITION]_[TYPE]_[NUMBER]

Examples:
CHRONIC_PAIN_RCT_001
ANXIETY_OBSERVATIONAL_001
DEPRESSION_CASE_SERIES_001
```

### Study Type Options
- `"RCT"` - Randomized Controlled Trial
- `"observational"` - Cohort, registry, real-world
- `"case_series"` - Small pilot, case reports

### Condition Values
- `"CHRONIC_PAIN"`
- `"ANXIETY"`
- `"DEPRESSION"`
- `"ARTHRITIS"`

### Delivery Method Options
- `"inhaled"` - Smoked cannabis
- `"vaporized"` - Vaporizer
- `"oral"` - Capsules, oils, tinctures
- `"sublingual"` - Under tongue
- `"topical"` - Creams, lotions
- `"metered-dose"` - Inhaler device

### Cannabinoid Profile Examples
- `"THC-dominant"`
- `"CBD-dominant"`
- `"balanced THC:CBD"`
- `"THC:CBD 1:1"`
- `"whole plant extract"`
- `"not specified"`

### Effect Size Category
- `"small"` - Modest effect
- `"medium"` - Moderate effect
- `"large"` - Strong effect
- `"unknown"` - Not reported

### Country Options (common)
- `"United States"`
- `"Canada"`
- `"United Kingdom"`
- `"Israel"`
- `"Netherlands"`
- `"Australia"`

---

## ⚡ Quick Extraction Steps

### 1. Read NORML Reference
Example from page:
```
[7] Abrams et al. 2007. Cannabis in painful HIV-associated 
sensory neuropathy: a randomized placebo-controlled trial. 
Neurology 68: 515-521.
```

### 2. Fill Template

**Study ID**: 
```
CHRONIC_PAIN_RCT_007
(increment number from previous)
```

**Citation**:
```
"Abrams et al. 2007. Neurology 68: 515-521."
```

**Study Title**:
```
"Cannabis in painful HIV-associated sensory neuropathy: 
a randomized placebo-controlled trial"
```

**Year**: `2007`

**Journal**: `"Neurology"`

**Study Type**: `"RCT"` (says "randomized placebo-controlled")

### 3. Read Paragraph Description

Look for:
- Sample size: "38 patients enrolled..."
- Intervention: "smoked cannabis 3x daily..."
- Results: "30% reduction in pain..."
- Side effects: "mild dizziness reported..."

### 4. Fill Intervention Block

```json
"intervention": {
  "cannabis_type": "smoked medicinal cannabis",
  "cannabinoid_profile": "not specified",
  "delivery_method": "inhaled",
  "dosing_information": "3 times daily",
  "treatment_duration": "5 days"
}
```

### 5. Fill Outcomes Block

```json
"outcomes": {
  "key_findings": [
    "Cannabis significantly reduced HIV neuropathic pain",
    "30% reduction in pain scores vs placebo",
    "Well-tolerated in HIV population"
  ],
  "outcome_measures": [
    "Neuropathic pain scale",
    "Pain intensity VAS"
  ],
  "quantitative_results": [],
  "adverse_effects": ["mild dizziness", "dry mouth"],
  "effect_size_category": "large"
}
```

### 6. Fill Study Quality Block

```json
"study_quality": {
  "sample_size": 38,
  "institution": "University of California San Francisco",
  "country": "United States",
  "randomized": true,
  "placebo_controlled": true
}
```

---

## 🚨 Common Mistakes to Avoid

### ❌ DON'T
```json
"study_id": "001"  // Missing condition and type
"publication_year": "2007"  // Should be number not string
"delivery_method": "smoking"  // Use "inhaled"
"randomized": "yes"  // Should be true/false
```

### ✅ DO
```json
"study_id": "CHRONIC_PAIN_RCT_001"  // Complete format
"publication_year": 2007  // Number
"delivery_method": "inhaled"  // Standard term
"randomized": true  // Boolean
```

---

## 📝 When Information is Missing

### Not in Text?
```json
"cannabinoid_profile": "not specified"
"dosing_information": "not specified"
"sample_size": null  // Use null, not 0
"quantitative_results": []  // Empty array
```

### Uncertain?
Add to notes:
```json
"notes": "Extracted from NORML chronic pain library. 
Sample size unclear from reference text."
```

---

## ✅ Validation Checklist

Before moving to next study:

- [ ] Study ID incremented correctly
- [ ] All required fields filled
- [ ] Year is number (not string)
- [ ] Boolean fields are true/false (not "yes"/"no")
- [ ] Arrays use brackets `[]`
- [ ] Objects use braces `{}`
- [ ] Commas between fields (but not after last field)
- [ ] JSON syntax valid

---

## 🎯 Quality Targets

**Per Study Time**: 20-30 minutes
- 5 min: Read reference
- 10 min: Fill template
- 5 min: Quality check
- 5 min: Add to file

**Per Day Target**: 8 studies
- 4 morning
- 4 afternoon
- 1 validation run

---

## 💾 Saving Your Work

### After Each Study
1. Add study to JSON array
2. Save file
3. Run validation every 5 studies

### Validation Command
```bash
python scripts/validate_norml_extraction.py
```

### Expected Output
```
✓ Study 11: CHRONIC_PAIN_RCT_011 - Complete
✓ Study 12: CHRONIC_PAIN_RCT_012 - Complete
...
Complete studies: 12/12
Completion rate: 100.0%
```

---

## 📱 Quick Reference Card Ends Here

**Save this file for easy access while extracting!**

**File**: `EXTRACTION_QUICKREF.md`  
**Location**: `data/norml_extraction/`

---

*Last updated: December 17, 2025*
