#!/usr/bin/env python3
"""
Clinical Evidence Expansion Workflow
Target: 200+ compound-target pairs with evidence ratings

Current state: 368 NORML studies → 3 positive samples (1.5%)
Target state: 200+ positive samples (40-50% positive rate)

Strategy:
1. Systematic reviews and meta-analyses (highest quality)
2. PubMed searches for well-studied compounds (THC, CBD, CBN)
3. ClinicalTrials.gov completed studies
4. Extract effect sizes from study conclusions
"""

import json
from pathlib import Path
from typing import Dict, List, Optional
from dataclasses import dataclass, asdict
from datetime import datetime


@dataclass
class ClinicalEvidence:
    """Structured clinical evidence for compound-condition pair."""
    compound: str
    condition: str
    effect_size: str  # 'large', 'medium', 'small', 'minimal', 'none'
    study_type: str  # 'meta_analysis', 'systematic_review', 'rct', 'observational'
    source: str  # PubMed ID, DOI, or NORML study ID
    participants: Optional[int] = None
    year: Optional[int] = None
    notes: str = ""
    confidence: str = "medium"  # 'high', 'medium', 'low'
    
    def to_binary_label(self) -> int:
        """Convert to binary classification label (1=effective, 0=not effective)."""
        return 1 if self.effect_size in ['large', 'medium'] else 0


class EvidenceExpander:
    """Workflow manager for clinical evidence expansion."""
    
    def __init__(self):
        self.project_root = Path(__file__).parent.parent
        self.data_dir = self.project_root / 'data' / 'clinical_evidence'
        self.data_dir.mkdir(parents=True, exist_ok=True)
        
        # Load existing NORML evidence
        self.existing_evidence = self._load_norml_evidence()
        
        print(f"Loaded {len(self.existing_evidence)} existing NORML studies")
        print(f"Evidence directory: {self.data_dir}")
    
    def _load_norml_evidence(self) -> List[ClinicalEvidence]:
        """Load existing NORML studies as baseline."""
        norml_file = self.project_root / 'data' / 'training' / 'neurobotanica_fully_enriched_fixed.json'
        
        if not norml_file.exists():
            return []
        
        with open(norml_file, 'r', encoding='utf-8') as f:
            data = json.load(f)
        
        evidence_list = []
        for compound_data in data['compounds']:
            compound = compound_data['compound_name']
            studies = compound_data.get('clinical_studies', [])
            
            for study in studies:
                evidence = ClinicalEvidence(
                    compound=compound,
                    condition=study.get('condition', 'unknown'),
                    effect_size=study.get('effect_size') or 'none',
                    study_type=study.get('study_type', 'observational'),
                    source=f"NORML_{study.get('year', 'unknown')}",
                    participants=study.get('participants'),
                    year=study.get('year'),
                    notes=study.get('conclusion', ''),
                    confidence='medium'
                )
                evidence_list.append(evidence)
        
        return evidence_list
    
    def generate_search_queries(self) -> Dict[str, List[str]]:
        """Generate PubMed search queries for priority compounds."""
        
        # Priority compounds with most studies
        priority_compounds = {
            'THC': ['delta-9-tetrahydrocannabinol', 'THC', 'dronabinol'],
            'CBD': ['cannabidiol', 'CBD', 'epidiolex'],
            'CBN': ['cannabinol', 'CBN'],
            'CBG': ['cannabigerol', 'CBG'],
            'THCV': ['tetrahydrocannabivarin', 'THCV']
        }
        
        # Priority conditions with strong evidence potential
        conditions = [
            'chronic pain', 'neuropathic pain', 'cancer pain',
            'PTSD', 'anxiety', 'depression',
            'epilepsy', 'seizures',
            'multiple sclerosis', 'spasticity',
            'nausea', 'chemotherapy-induced nausea',
            'sleep disorders', 'insomnia',
            'Alzheimer\'s disease', 'dementia',
            'Parkinson\'s disease',
            'inflammatory bowel disease', 'Crohn\'s disease'
        ]
        
        queries = {}
        
        for compound, terms in priority_compounds.items():
            queries[compound] = []
            
            for condition in conditions:
                for term in terms:
                    # Systematic review query
                    query = f'("{term}"[Title/Abstract] AND "{condition}"[Title/Abstract] AND ("systematic review"[Publication Type] OR "meta-analysis"[Publication Type]))'
                    queries[compound].append({
                        'query': query,
                        'condition': condition,
                        'type': 'systematic_review',
                        'priority': 'high'
                    })
                    
                    # RCT query
                    query = f'("{term}"[Title/Abstract] AND "{condition}"[Title/Abstract] AND "randomized controlled trial"[Publication Type])'
                    queries[compound].append({
                        'query': query,
                        'condition': condition,
                        'type': 'rct',
                        'priority': 'medium'
                    })
        
        return queries
    
    def create_evidence_template(self, compound: str, condition: str) -> str:
        """Create template for manual evidence entry."""
        template = f"""
# Clinical Evidence Entry
**Compound:** {compound}
**Condition:** {condition}
**Date:** {datetime.now().strftime('%Y-%m-%d')}

---

## Study Information
- **Study Type:** [meta_analysis / systematic_review / rct / observational]
- **Source:** [PubMed ID or DOI]
- **Year:** [YYYY]
- **Participants:** [number]

## Evidence Quality
- **Effect Size:** [large / medium / small / minimal / none]
  - **large**: Clear, clinically significant benefit (e.g., >50% symptom reduction, NNT<5)
  - **medium**: Moderate benefit (e.g., 30-50% improvement, NNT 5-10)
  - **small**: Minor benefit (e.g., 10-30% improvement, NNT >10)
  - **minimal**: Statistically significant but not clinically meaningful
  - **none**: No significant effect or negative result

- **Confidence:** [high / medium / low]
  - **high**: Meta-analysis or multiple large RCTs (n>500)
  - **medium**: Single large RCT or systematic review
  - **low**: Small RCT (n<100) or observational study

## Study Findings
**Key Results:**
[Summarize main efficacy findings]

**Safety Profile:**
[Note any safety concerns or adverse events]

**Conclusion:**
[Study authors' conclusion about efficacy]

## Notes
[Additional context, limitations, special populations, etc.]

---
## Save Instructions
1. Fill in all bracketed fields above
2. Save as: `{compound.lower()}_{condition.lower().replace(' ', '_')}_YYYYMMDD.md`
3. Run: `python scripts/import_evidence.py` to add to dataset
"""
        return template
    
    def save_search_guide(self):
        """Save comprehensive search guide for evidence gathering."""
        guide_path = self.data_dir / 'EVIDENCE_EXPANSION_GUIDE.md'
        
        guide = """# Clinical Evidence Expansion Guide
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
"""
        
        with open(guide_path, 'w', encoding='utf-8') as f:
            f.write(guide)
        
        print(f"\n✅ Evidence expansion guide saved: {guide_path}")
        return guide_path
    
    def save_search_queries(self):
        """Save PubMed search queries for easy copy-paste."""
        queries_path = self.data_dir / 'pubmed_search_queries.txt'
        
        queries = self.generate_search_queries()
        
        with open(queries_path, 'w', encoding='utf-8') as f:
            f.write("# PubMed Search Queries - Clinical Evidence Expansion\n")
            f.write(f"# Generated: {datetime.now().strftime('%Y-%m-%d')}\n")
            f.write("# Copy and paste these directly into PubMed Advanced Search\n\n")
            
            for compound, query_list in queries.items():
                f.write(f"\n{'='*80}\n")
                f.write(f"# {compound}\n")
                f.write(f"{'='*80}\n\n")
                
                # Group by priority
                high_priority = [q for q in query_list if q['priority'] == 'high']
                medium_priority = [q for q in query_list if q['priority'] == 'medium']
                
                f.write("## HIGH PRIORITY - Systematic Reviews & Meta-Analyses\n\n")
                for i, q in enumerate(high_priority[:10], 1):  # Top 10
                    f.write(f"{i}. {q['condition']}\n")
                    f.write(f"   {q['query']}\n\n")
                
                f.write("\n## MEDIUM PRIORITY - Randomized Controlled Trials\n\n")
                for i, q in enumerate(medium_priority[:5], 1):  # Top 5
                    f.write(f"{i}. {q['condition']}\n")
                    f.write(f"   {q['query']}\n\n")
        
        print(f"✅ Search queries saved: {queries_path}")
        return queries_path
    
    def analyze_current_state(self):
        """Analyze current evidence to identify gaps."""
        print("\n" + "="*80)
        print("CURRENT EVIDENCE ANALYSIS")
        print("="*80)
        
        # Count by compound
        compound_counts = {}
        for evidence in self.existing_evidence:
            compound = evidence.compound
            if compound not in compound_counts:
                compound_counts[compound] = {'total': 0, 'positive': 0, 'negative': 0}
            
            compound_counts[compound]['total'] += 1
            if evidence.to_binary_label() == 1:
                compound_counts[compound]['positive'] += 1
            else:
                compound_counts[compound]['negative'] += 1
        
        print("\nEvidence by Compound:")
        print(f"{'Compound':<15} {'Total':<10} {'Positive':<10} {'Negative':<10} {'Pos %':<10}")
        print("-"*65)
        
        total_samples = sum(c['total'] for c in compound_counts.values())
        total_positive = sum(c['positive'] for c in compound_counts.values())
        
        for compound, counts in sorted(compound_counts.items(), key=lambda x: x[1]['total'], reverse=True):
            pos_pct = (counts['positive'] / counts['total'] * 100) if counts['total'] > 0 else 0
            print(f"{compound:<15} {counts['total']:<10} {counts['positive']:<10} {counts['negative']:<10} {pos_pct:<10.1f}%")
        
        print("-"*65)
        print(f"{'TOTAL':<15} {total_samples:<10} {total_positive:<10} {total_samples - total_positive:<10} {total_positive/total_samples*100:<10.1f}%")
        
        print(f"\n⚠️  Current positive sample rate: {total_positive/total_samples*100:.1f}% ({total_positive}/{total_samples})")
        print(f"🎯 Target positive sample rate: 40-50% (200+ positive samples)")
        print(f"📈 Need to add: {200 - total_positive} positive samples")
        
        # Identify top compounds for expansion
        print("\n" + "="*80)
        print("EXPANSION PRIORITIES")
        print("="*80)
        
        priorities = [
            ('THC', 246, 'Most studied, highest yield potential'),
            ('CBD', 282, 'FDA-approved uses, strong evidence base'),
            ('CBN', 4, 'Sleep/sedation niche, emerging evidence'),
            ('CBG', 0, 'Anti-inflammatory, antibacterial potential')
        ]
        
        print("\n1. HIGH PRIORITY - Evidence-Rich Compounds:")
        for compound, study_count, rationale in priorities[:2]:
            print(f"   • {compound}: {study_count} NORML studies")
            print(f"     Rationale: {rationale}")
            current = compound_counts.get(compound, {}).get('positive', 0)
            print(f"     Current: {current} positive samples")
            print(f"     Target: +{30 if compound == 'THC' else 25} positive samples")
            print()
        
        print("2. MEDIUM PRIORITY - Emerging Compounds:")
        for compound, study_count, rationale in priorities[2:]:
            print(f"   • {compound}: {study_count} NORML studies")
            print(f"     Rationale: {rationale}")
            print(f"     Target: +5-10 positive samples")
            print()


def main():
    """Main workflow."""
    import sys
    
    expander = EvidenceExpander()
    
    if len(sys.argv) > 1:
        command = sys.argv[1].lower()
        
        if command == 'analyze':
            expander.analyze_current_state()
        
        elif command == 'template':
            if len(sys.argv) < 4:
                print("Usage: python scripts/expand_clinical_evidence.py template [COMPOUND] [CONDITION]")
                print('Example: python scripts/expand_clinical_evidence.py template THC "chronic pain"')
                return
            
            compound = sys.argv[2]
            condition = ' '.join(sys.argv[3:])
            
            template = expander.create_evidence_template(compound, condition)
            
            filename = f"{compound.lower()}_{condition.lower().replace(' ', '_')}_{datetime.now().strftime('%Y%m%d')}.md"
            filepath = expander.data_dir / filename
            
            with open(filepath, 'w', encoding='utf-8') as f:
                f.write(template)
            
            print(f"\n✅ Evidence template created: {filepath}")
            print(f"\nNext steps:")
            print(f"1. Open the file and fill in study details")
            print(f"2. Search PubMed for: {compound} {condition} meta-analysis")
            print(f"3. Extract effect size from study conclusions")
            print(f"4. Run: python scripts/import_evidence.py")
        
        elif command == 'setup':
            print("\n🚀 Setting up evidence expansion workflow...\n")
            expander.analyze_current_state()
            expander.save_search_guide()
            expander.save_search_queries()
            
            print("\n" + "="*80)
            print("✅ SETUP COMPLETE")
            print("="*80)
            print(f"\nYour evidence expansion workspace is ready in:")
            print(f"  {expander.data_dir}")
            print(f"\nKey files created:")
            print(f"  • EVIDENCE_EXPANSION_GUIDE.md - Complete search strategy")
            print(f"  • pubmed_search_queries.txt - Ready-to-use PubMed queries")
            print(f"\nNext steps:")
            print(f"  1. Read: {expander.data_dir / 'EVIDENCE_EXPANSION_GUIDE.md'}")
            print(f"  2. Start with THC chronic pain (highest yield)")
            print(f"  3. Create template: python scripts/expand_clinical_evidence.py template THC \"chronic pain\"")
            print(f"  4. Search PubMed and fill in evidence")
            print(f"  5. Import: python scripts/import_evidence.py")
            print()
        
        else:
            print(f"Unknown command: {command}")
            print("\nAvailable commands:")
            print("  analyze  - Analyze current evidence and identify gaps")
            print("  setup    - Create evidence expansion workspace")
            print("  template [COMPOUND] [CONDITION] - Create evidence entry template")
    
    else:
        print("\nClinical Evidence Expansion Workflow")
        print("="*80)
        print("\nCommands:")
        print("  python scripts/expand_clinical_evidence.py setup")
        print("    → Create evidence expansion workspace and guides")
        print()
        print("  python scripts/expand_clinical_evidence.py analyze")
        print("    → Analyze current evidence and identify gaps")
        print()
        print("  python scripts/expand_clinical_evidence.py template THC \"chronic pain\"")
        print("    → Create evidence entry template for manual data entry")
        print()
        print("Start with: python scripts/expand_clinical_evidence.py setup")
        print()


if __name__ == '__main__':
    main()
