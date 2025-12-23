#!/usr/bin/env python3
"""
Automated Clinical Evidence Expansion
- Uses PubMed Entrez API to search for systematic reviews/meta-analyses
- Extracts effect size from abstracts using keyword heuristics
- Auto-generates evidence entries for review

Requires: Biopython (for Entrez)
  pip install biopython
"""

import os
import sys
import time
import json
from pathlib import Path
from datetime import datetime
from typing import List, Dict
from Bio import Entrez

# Set your email for Entrez API
Entrez.email = os.environ.get('ENTREZ_EMAIL', 'your_email@example.com')

# Output directory
PROJECT_ROOT = Path(__file__).parent.parent
EVIDENCE_DIR = PROJECT_ROOT / 'data' / 'clinical_evidence'
EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)

# Compound-condition search map (priority pairs)
COMPOUND_CONDITIONS = {
    'THC': [
        'chronic pain', 'cancer pain', 'PTSD', 'nausea', 'multiple sclerosis', 'sleep disorders', "Alzheimer's disease", "Parkinson's disease"
    ],
    'CBD': [
        'epilepsy', 'anxiety', 'depression', 'inflammatory bowel disease', 'sleep disorders', 'autism', 'addiction', 'schizophrenia'
    ],
    'CBN': ['sleep disorders', 'insomnia'],
    'CBG': ['inflammatory bowel disease', 'glaucoma', 'antibacterial', 'neuroprotection']
}

# Effect size keyword map
EFFECT_SIZE_MAP = [
    ('large', [
        'strong evidence', 'clinically significant', 'NNT < 5', 'effect size d > 0.8', 'SMD > 0.8', '>50% reduction', 'FDA approval', 'robust efficacy', 'marked improvement'
    ]),
    ('medium', [
        'moderate evidence', 'significant clinical benefit', 'NNT 5-10', 'effect size d 0.5', 'effect size d 0.6', 'effect size d 0.7', 'effect size d 0.8', 'SMD 0.5', 'SMD 0.6', 'SMD 0.7', 'SMD 0.8', '30% improvement', '40% improvement', '50% improvement', 'moderate improvement'
    ]),
    ('small', [
        'modest benefit', 'statistically significant', 'NNT > 10', 'effect size d 0.2', 'effect size d 0.3', 'effect size d 0.4', 'effect size d 0.5', 'SMD 0.2', 'SMD 0.3', 'SMD 0.4', 'SMD 0.5', '10% improvement', '20% improvement', 'minor improvement'
    ]),
    ('minimal', [
        'not clinically meaningful', 'minimal effect', 'borderline significance', 'weak evidence', 'marginal benefit'
    ]),
    ('none', [
        'no significant difference', 'insufficient evidence', 'not effective', 'no effect', 'p > 0.05', 'failed to demonstrate', 'negative result'
    ])
]

# PubMed search utility
def pubmed_search(term, retmax=10):
    handle = Entrez.esearch(db="pubmed", term=term, retmax=retmax, sort="relevance")
    record = Entrez.read(handle)
    handle.close()
    return record.get('IdList', [])

# PubMed fetch utility
def pubmed_fetch_details(id_list):
    ids = ','.join(id_list)
    handle = Entrez.efetch(db="pubmed", id=ids, rettype="medline", retmode="text")
    records = handle.read()
    handle.close()
    return records

# Abstract fetch utility (returns list of dicts)
def fetch_abstracts(id_list):
    handle = Entrez.efetch(db="pubmed", id=','.join(id_list), rettype="abstract", retmode="xml")
    records = Entrez.read(handle)
    handle.close()
    results = []
    for article in records['PubmedArticle']:
        try:
            pmid = article['MedlineCitation']['PMID']
            title = article['MedlineCitation']['Article']['ArticleTitle']
            abstract = article['MedlineCitation']['Article']['Abstract']['AbstractText'][0]
            year = int(article['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate'].get('Year', '0'))
            results.append({'pmid': str(pmid), 'title': title, 'abstract': abstract, 'year': year})
        except Exception:
            continue
    return results

# Heuristic effect size extraction
def extract_effect_size(text):
    text = text.lower()
    for label, keywords in EFFECT_SIZE_MAP:
        for kw in keywords:
            if kw in text:
                return label
    return 'none'

# Main automation routine
RETMAX = int(os.environ.get('PUBMED_RETMAX', '20'))


def auto_expand():
    total_saved = 0
    for compound, conditions in COMPOUND_CONDITIONS.items():
        for condition in conditions:
            # Build PubMed query for systematic reviews/meta-analyses
            query = f'("{compound}"[Title/Abstract] AND "{condition}"[Title/Abstract] AND ("systematic review"[Publication Type] OR "meta-analysis"[Publication Type]))'
            print(f"Searching: {compound} + {condition} (retmax={RETMAX})")
            try:
                ids = pubmed_search(query, retmax=RETMAX)
                print(f"  Found {len(ids)} IDs")
                if not ids:
                    print(f"  No results for query: {query}")
                    continue
                abstracts = fetch_abstracts(ids)
                if not abstracts:
                    print("  No abstracts returned by fetch.")
                for ab in abstracts:
                    # ensure abstract text exists
                    abstract_text = ab.get('abstract', '') if isinstance(ab, dict) else ''
                    effect_size = extract_effect_size(abstract_text)
                    entry = {
                        'compound': compound,
                        'condition': condition,
                        'effect_size': effect_size,
                        'study_type': 'systematic_review/meta-analysis',
                        'source': f"PubMed:{ab.get('pmid', '')}",
                        'participants': None,
                        'year': ab.get('year'),
                        'notes': ab.get('title'),
                        'confidence': 'high' if effect_size in ['large', 'medium'] else 'medium',
                        'abstract': abstract_text
                    }
                    # Save as JSON for review
                    fname = f"{compound.lower()}_{condition.lower().replace(' ', '_')}_{ab.get('pmid', 'unknown')}.json"
                    with open(EVIDENCE_DIR / fname, 'w', encoding='utf-8') as f:
                        json.dump(entry, f, indent=2)
                    total_saved += 1
                    print(f"  Saved: {fname} (effect_size: {effect_size})")
                time.sleep(1)  # Be polite to NCBI
            except Exception as e:
                print(f"  Error: {e}")
                continue
    print(f"Total evidence files saved: {total_saved}")

if __name__ == '__main__':
    print("Automated Clinical Evidence Expansion - PubMed Systematic Reviews")
    print("="*70)
    auto_expand()
    print("\n✅ Automation complete. Review new evidence files in data/clinical_evidence/")
