#!/usr/bin/env python3
"""
Automated Clinical Study Expansion Script
Expands NORML clinical studies from 368 to 500+ by pulling from:
- ClinicalTrials.gov API (completed cannabis trials)
- PubMed Entrez API (systematic reviews and meta-analyses)

Usage:
    python scripts/automate_clinical_study_expansion.py [--dry-run] [--limit N]

Requirements:
    pip install requests biopython
    Set ENTREZ_EMAIL environment variable for PubMed API
"""

import os
import sys
import json
import time
import requests
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Any
from dataclasses import dataclass, asdict
import logging
import argparse

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

try:
    from Bio import Entrez
    ENTREZ_AVAILABLE = True
    # Set Entrez email
    Entrez.email = os.environ.get('ENTREZ_EMAIL', 'neurobotanica@cqr.org')
except ImportError:
    ENTREZ_AVAILABLE = False
    Entrez = None
    print("Warning: BioPython not available. PubMed search disabled.")

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

@dataclass
class ClinicalTrial:
    """ClinicalTrials.gov study data."""
    nct_id: str
    title: str
    condition: str
    status: str
    phase: str
    enrollment: Optional[int]
    completion_date: Optional[str]
    interventions: List[str]
    outcomes: List[str]
    url: str

@dataclass
class PubMedStudy:
    """PubMed study data."""
    pmid: str
    title: str
    abstract: str
    authors: List[str]
    journal: str
    year: int
    condition: str
    intervention: str
    study_type: str
    key_findings: List[str]

class ClinicalStudyExpander:
    """Automated clinical study expansion from external sources."""

    def __init__(self, dry_run: bool = False, limit: int = 100):
        self.dry_run = dry_run
        self.limit = limit
        self.project_root = Path(__file__).parent.parent
        self.norml_dir = self.project_root / 'data' / 'norml_extraction'
        self.output_dir = self.project_root / 'data' / 'expanded_studies'
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # ClinicalTrials.gov API base
        self.ct_api_base = "https://clinicaltrials.gov/api/v2/studies"

        # Load existing studies to avoid duplicates
        self.existing_studies = self._load_existing_studies()

    def _load_existing_studies(self) -> set:
        """Load existing study IDs to avoid duplicates."""
        existing = set()
        for json_file in self.norml_dir.glob('*_studies.json'):
            try:
                with open(json_file, 'r', encoding='utf-8') as f:
                    data = json.load(f)
                    for study in data.get('studies', []):
                        existing.add(study.get('study_id', ''))
            except Exception as e:
                print(f"Warning: Error loading {json_file}: {e}")
        return existing

    def search_clinicaltrials_gov(self, condition: str, max_results: int = 50) -> List[ClinicalTrial]:
        """Search ClinicalTrials.gov for cannabis-related studies."""
        trials = []

        # Condition mapping for broader search terms
        condition_mapping = {
            'chronic pain': ['pain', 'chronic pain', 'neuropathic pain', 'nociceptive pain'],
            'cancer pain': ['cancer pain', 'oncology pain', 'malignant pain'],
            'anxiety': ['anxiety', 'anxiety disorders', 'generalized anxiety'],
            'depression': ['depression', 'major depressive disorder', 'depressive disorder'],
            'epilepsy': ['epilepsy', 'seizure disorder', 'epileptic seizures'],
            'multiple sclerosis': ['multiple sclerosis', 'ms', 'disseminated sclerosis'],
            'ptsd': ['ptsd', 'post-traumatic stress disorder', 'post traumatic stress'],
            'insomnia': ['insomnia', 'sleep disorder', 'sleep disturbance'],
            'nausea': ['nausea', 'chemotherapy induced nausea', 'radiation induced nausea'],
            'appetite': ['appetite', 'cachexia', 'anorexia', 'weight loss'],
            'inflammation': ['inflammation', 'inflammatory disease', 'arthritis'],
            'alzheimer': ["alzheimer's disease", 'alzheimer', 'dementia'],
            'parkinson': ["parkinson's disease", 'parkinson', 'parkinsonian disorders'],
            'glaucoma': ['glaucoma', 'ocular hypertension'],
            'ibd': ['inflammatory bowel disease', 'crohn disease', 'ulcerative colitis'],
            'autism': ['autism', 'autism spectrum disorder', 'asd'],
            'schizophrenia': ['schizophrenia', 'psychotic disorder'],
            'addiction': ['addiction', 'substance use disorder', 'dependence'],
            'tourette': ['tourette syndrome', 'tic disorder'],
            'covid': ['covid-19', 'sars-cov-2', 'coronavirus']
        }

        # Get search terms for this condition
        search_conditions = condition_mapping.get(condition.lower(), [condition.lower()])

        # Expanded cannabis search terms
        cannabis_terms = [
            'cannabis', 'cannabinoid', 'THC', 'CBD', 'CBN', 'CBG', 'cannabidiol',
            'dronabinol', 'marinol', 'sativex', 'nabiximols', 'epidolex', 'epidiolex',
            'medical marijuana', 'medical cannabis', 'cannabis sativa', 'hemp extract'
        ]

        # Strategy 1: Search by condition first, then filter for cannabis interventions
        condition_query = ' OR '.join([f'"{cond}"' for cond in search_conditions])
        params = {
            'query.cond': condition_query,
            'filter.overallStatus': 'COMPLETED,TERMINATED,WITHDRAWN',  # Include more statuses
            'countTotal': 'true',
            'pageSize': min(max_results * 3, 100)  # Get more results to filter
        }

        logger.info(f"Searching ClinicalTrials.gov for {condition} with condition query: {condition_query[:100]}...")

        try:
            response = requests.get(self.ct_api_base, params=params, timeout=30)
            response.raise_for_status()
            data = response.json()

            logger.info(f"ClinicalTrials.gov API returned {len(data.get('studies', []))} studies")

            for study in data.get('studies', []):
                protocol = study.get('protocolSection', {})

                # Check if this study is relevant to our condition
                conditions_list = protocol.get('conditionsModule', {}).get('conditions', [])
                conditions_text = ' '.join([
                    cond if isinstance(cond, str) else cond.get('condition', '')
                    for cond in conditions_list
                ]).lower()

                # Check if any of our search conditions match
                condition_match = any(
                    search_cond in conditions_text
                    for search_cond in search_conditions
                )

                # For cannabis studies, be more lenient - include if it has pain-related terms
                # or if it's a general pain study that might be relevant
                pain_related = any(term in conditions_text for term in [
                    'pain', 'chronic', 'neuropathic', 'nociceptive', 'fibromyalgia',
                    'arthritis', 'rheumatoid', 'osteoarthritis', 'headache', 'migraine',
                    'back pain', 'neck pain', 'musculoskeletal', 'inflammation'
                ])

                if not (condition_match or pain_related):
                    logger.debug(f"Skipping study {protocol.get('identificationModule', {}).get('nctId', '')} - conditions: {conditions_text[:100]}...")
                    continue

                # Extract interventions
                interventions = []
                arms = protocol.get('armsInterventionsModule', {}).get('interventions', [])
                logger.debug(f"Study {protocol.get('identificationModule', {}).get('nctId', '')} has {len(arms)} interventions")
                for i, arm in enumerate(arms):
                    arm_name = arm.get('name', '').lower()
                    logger.debug(f"  Intervention {i+1}: '{arm_name}'")
                    if any(term in arm_name for term in cannabis_terms):
                        interventions.append(arm.get('name', ''))
                        logger.debug(f"  -> MATCHED cannabis term")
                    else:
                        logger.debug(f"  -> No match")

                # Check for cannabis in title or description as well
                title = protocol.get('identificationModule', {}).get('briefTitle', '').lower()
                description = protocol.get('descriptionModule', {}).get('briefSummary', '').lower() if protocol.get('descriptionModule') else ''

                has_cannabis_in_text = any(term.lower() in title or term.lower() in description for term in cannabis_terms)

                if not interventions and not has_cannabis_in_text:
                    logger.debug(f"No cannabis interventions or mentions found for study {protocol.get('identificationModule', {}).get('nctId', '')}")
                    continue  # Skip if no cannabis interventions or mentions

                if interventions:
                    logger.info(f"Including study {protocol.get('identificationModule', {}).get('nctId', '')} - conditions: {conditions_text}...")
                elif has_cannabis_in_text:
                    logger.info(f"Including study {protocol.get('identificationModule', {}).get('nctId', '')} - cannabis mentioned in title/description")
                    interventions = ['Cannabis mentioned in study text']  # Placeholder intervention

                logger.info(f"Study {protocol.get('identificationModule', {}).get('nctId', '')} passed all filters - adding to results")

                # Extract outcomes
                outcomes = []
                outcomes_module = protocol.get('outcomesModule', {})
                primary_outcomes = outcomes_module.get('primaryOutcomes', [])
                for outcome in primary_outcomes:
                    if isinstance(outcome, dict):
                        outcomes.append(outcome.get('measure', ''))
                    else:
                        outcomes.append(str(outcome))

                trial = ClinicalTrial(
                    nct_id=protocol.get('identificationModule', {}).get('nctId', ''),
                    title=protocol.get('identificationModule', {}).get('briefTitle', ''),
                    condition=condition,
                    status=protocol.get('statusModule', {}).get('overallStatus', ''),
                    phase=protocol.get('designModule', {}).get('phases', ['N/A'])[0] if protocol.get('designModule', {}).get('phases') else 'N/A',
                    enrollment=protocol.get('designModule', {}).get('enrollmentInfo', {}).get('count'),
                    completion_date=protocol.get('statusModule', {}).get('completionDateStruct', {}).get('date'),
                    interventions=interventions,
                    outcomes=outcomes,
                    url=f"https://clinicaltrials.gov/study/{protocol.get('identificationModule', {}).get('nctId', '')}"
                )

                trials.append(trial)

        except Exception as e:
            logger.error(f"Error searching ClinicalTrials.gov for {condition}: {e}")

        # Strategy 2: If few results, try condition + intervention search
        if len(trials) < 5:
            try:
                condition_query = ' OR '.join([f'"{cond}"' for cond in search_conditions])
                params = {
                    'query.cond': condition_query,
                    'query.intr': ' OR '.join([f'"{term}"' for term in cannabis_terms[:5]]),  # Limit to avoid URL too long
                    'filter.overallStatus': 'COMPLETED',
                    'countTotal': 'true',
                    'pageSize': min(max_results - len(trials), 100)
                }

                logger.info(f"Trying secondary search for {condition} with combined condition + intervention query")

                response = requests.get(self.ct_api_base, params=params, timeout=30)
                response.raise_for_status()
                data = response.json()

                existing_ncts = {t.nct_id for t in trials}

                for study in data.get('studies', []):
                    nct_id = study.get('protocolSection', {}).get('identificationModule', {}).get('nctId', '')
                    if nct_id in existing_ncts:
                        continue

                    protocol = study.get('protocolSection', {})

                    # Extract interventions
                    interventions = []
                    arms = protocol.get('armsInterventionsModule', {}).get('interventions', [])
                    for arm in arms:
                        arm_name = arm.get('name', '').lower()
                        if any(term in arm_name for term in cannabis_terms):
                            interventions.append(arm.get('name', ''))

                    if not interventions:
                        continue

                    # Extract outcomes
                    outcomes = []
                    outcomes_module = protocol.get('outcomesModule', {})
                    primary_outcomes = outcomes_module.get('primaryOutcomes', [])
                    for outcome in primary_outcomes:
                        if isinstance(outcome, dict):
                            outcomes.append(outcome.get('measure', ''))
                        else:
                            outcomes.append(str(outcome))

                    trial = ClinicalTrial(
                        nct_id=nct_id,
                        title=protocol.get('identificationModule', {}).get('briefTitle', ''),
                        condition=condition,
                        status=protocol.get('statusModule', {}).get('overallStatus', ''),
                        phase=protocol.get('designModule', {}).get('phases', ['N/A'])[0] if protocol.get('designModule', {}).get('phases') else 'N/A',
                        enrollment=protocol.get('designModule', {}).get('enrollmentInfo', {}).get('count'),
                        completion_date=protocol.get('statusModule', {}).get('completionDateStruct', {}).get('date'),
                        interventions=interventions,
                        outcomes=outcomes,
                        url=f"https://clinicaltrials.gov/study/{nct_id}"
                    )

                    trials.append(trial)

            except Exception as e:
                logger.warning(f"Error in secondary ClinicalTrials.gov search for {condition}: {e}")

        # Strategy 3: If still few results, try cannabis intervention first, then filter conditions
        if len(trials) < 3:
            try:
                intervention_query = ' OR '.join([f'"{term}"' for term in cannabis_terms[:8]])  # Use more terms
                params = {
                    'query.intr': intervention_query,
                    'filter.overallStatus': 'COMPLETED',
                    'countTotal': 'true',
                    'pageSize': min(max_results * 3, 100)  # Get more results to filter
                }

                logger.info(f"Trying tertiary search for {condition} with cannabis intervention first")

                response = requests.get(self.ct_api_base, params=params, timeout=30)
                response.raise_for_status()
                data = response.json()

                existing_ncts = {t.nct_id for t in trials}

                for study in data.get('studies', []):
                    nct_id = study.get('protocolSection', {}).get('identificationModule', {}).get('nctId', '')
                    if nct_id in existing_ncts:
                        continue

                    protocol = study.get('protocolSection', {})

                    # Check if this study is relevant to our condition
                    conditions_list = protocol.get('conditionsModule', {}).get('conditions', [])
                    conditions_text = ' '.join([
                        cond if isinstance(cond, str) else cond.get('condition', '')
                        for cond in conditions_list
                    ]).lower()

                    # Check if any of our search conditions match
                    condition_match = any(
                        search_cond in conditions_text
                        for search_cond in search_conditions
                    )

                    # For cannabis studies, be more lenient - include if it has pain-related terms
                    pain_related = any(term in conditions_text for term in [
                        'pain', 'chronic', 'neuropathic', 'nociceptive', 'fibromyalgia',
                        'arthritis', 'rheumatoid', 'osteoarthritis', 'headache', 'migraine',
                        'back pain', 'neck pain', 'musculoskeletal', 'inflammation'
                    ])

                    if not (condition_match or pain_related):
                        continue

                    # Extract interventions
                    interventions = []
                    arms = protocol.get('armsInterventionsModule', {}).get('interventions', [])
                    for arm in arms:
                        arm_name = arm.get('name', '').lower()
                        if any(term in arm_name for term in cannabis_terms):
                            interventions.append(arm.get('name', ''))

                    if not interventions:
                        continue

                    # Extract outcomes
                    outcomes = []
                    outcomes_module = protocol.get('outcomesModule', {})
                    primary_outcomes = outcomes_module.get('primaryOutcomes', [])
                    for outcome in primary_outcomes:
                        if isinstance(outcome, dict):
                            outcomes.append(outcome.get('measure', ''))
                        else:
                            outcomes.append(str(outcome))

                    trial = ClinicalTrial(
                        nct_id=nct_id,
                        title=protocol.get('identificationModule', {}).get('briefTitle', ''),
                        condition=condition,
                        status=protocol.get('statusModule', {}).get('overallStatus', ''),
                        phase=protocol.get('designModule', {}).get('phases', ['N/A'])[0] if protocol.get('designModule', {}).get('phases') else 'N/A',
                        enrollment=protocol.get('designModule', {}).get('enrollmentInfo', {}).get('count'),
                        completion_date=protocol.get('statusModule', {}).get('completionDateStruct', {}).get('date'),
                        interventions=interventions,
                        outcomes=outcomes,
                        url=f"https://clinicaltrials.gov/study/{nct_id}"
                    )

                    trials.append(trial)

            except Exception as e:
                logger.warning(f"Error in tertiary ClinicalTrials.gov search for {condition}: {e}")

        return trials[:max_results]

    def search_pubmed(self, condition: str, max_results: int = 50) -> List[PubMedStudy]:
        """Search PubMed for systematic reviews and meta-analyses."""
        if not ENTREZ_AVAILABLE:
            logger.warning("BioPython not available, skipping PubMed search")
            return []

        studies = []

        # Search query for systematic reviews/meta-analyses
        query = f'({condition}[Title/Abstract]) AND cannabis[Title/Abstract] AND (systematic review[Publication Type] OR meta-analysis[Publication Type])'

        try:
            # Search for PMIDs
            handle = Entrez.esearch(db='pubmed', term=query, retmax=max_results)
            record = Entrez.read(handle)
            pmids = record.get('IdList', [])
            handle.close()

            if not pmids:
                return studies

            # Fetch details for first few to avoid rate limits
            pmids_subset = pmids[:10]  # Limit to 10 for testing

            for pmid in pmids_subset:
                try:
                    # Get summary for metadata
                    handle = Entrez.esummary(db='pubmed', id=pmid)
                    summary = Entrez.read(handle)
                    handle.close()

                    if summary and len(summary) > 0:
                        article = summary[0]
                        study = PubMedStudy(
                            pmid=pmid,
                            title=article.get('Title', ''),
                            abstract='',  # Skip abstract for now to avoid rate limits
                            authors=[auth.get('Name', '') for auth in article.get('AuthorList', [])],
                            journal=article.get('Source', ''),
                            year=int(article.get('PubDate', '').split()[0]) if article.get('PubDate') else 2025,
                            condition=condition,
                            intervention='cannabis/cannabinoids',
                            study_type='systematic review/meta-analysis',
                            key_findings=['Automated extraction - review full paper for details']
                        )
                        studies.append(study)

                except Exception as e:
                    logger.warning(f"Error processing PubMed ID {pmid}: {e}")

                time.sleep(0.5)  # Rate limiting

        except Exception as e:
            logger.error(f"Error searching PubMed for {condition}: {e}")

        return studies[:max_results]

    def convert_trial_to_norml_format(self, trial: ClinicalTrial) -> Dict[str, Any]:
        """Convert ClinicalTrial to NORML study format."""
        return {
            "study_id": f"CT_{trial.nct_id}",
            "study_type": "RCT" if trial.phase != "N/A" else "Clinical Trial",
            "condition": trial.condition.upper().replace(' ', '_'),
            "study_title": trial.title,
            "citation": f"ClinicalTrials.gov NCT{trial.nct_id}",
            "publication_year": int(trial.completion_date.split('-')[0]) if trial.completion_date else 2025,
            "journal": "ClinicalTrials.gov",
            "intervention": {
                "cannabis_type": "cannabis/cannabinoids",
                "cannabinoid_profile": ", ".join(trial.interventions),
                "delivery_method": "various",
                "dosing_information": f"Phase {trial.phase}",
                "treatment_duration": "clinical trial duration"
            },
            "outcomes": {
                "key_findings": trial.outcomes,
                "outcome_measures": trial.outcomes,
                "adverse_events": [],
                "evidence_rating": "medium",
                "confidence_level": "medium"
            },
            "metadata": {
                "nct_id": trial.nct_id,
                "enrollment": trial.enrollment,
                "phase": trial.phase,
                "source": "ClinicalTrials.gov",
                "url": trial.url
            }
        }

    def convert_pubmed_to_norml_format(self, study: PubMedStudy) -> Dict[str, Any]:
        """Convert PubMedStudy to NORML study format."""
        return {
            "study_id": f"PM_{study.pmid}",
            "study_type": study.study_type.upper().replace(' ', '_'),
            "condition": study.condition.upper().replace(' ', '_'),
            "study_title": study.title,
            "citation": f"{study.authors[0] if study.authors else 'Unknown'} et al. {study.year}. {study.journal}",
            "publication_year": study.year,
            "journal": study.journal,
            "intervention": {
                "cannabis_type": "cannabis/cannabinoids",
                "cannabinoid_profile": study.intervention,
                "delivery_method": "various",
                "dosing_information": "systematic review",
                "treatment_duration": "various"
            },
            "outcomes": {
                "key_findings": study.key_findings,
                "outcome_measures": study.key_findings,
                "adverse_events": [],
                "evidence_rating": "high",
                "confidence_level": "high"
            },
            "metadata": {
                "pmid": study.pmid,
                "authors": study.authors,
                "abstract": study.abstract[:500] + "..." if len(study.abstract) > 500 else study.abstract,
                "source": "PubMed",
                "url": f"https://pubmed.ncbi.nlm.nih.gov/{study.pmid}/"
            }
        }

    def expand_condition(self, condition: str) -> Dict[str, Any]:
        """Expand studies for a specific condition."""
        logger.info(f"Expanding studies for condition: {condition}")

        # Search ClinicalTrials.gov
        ct_trials = self.search_clinicaltrials_gov(condition, max_results=self.limit // 2)
        logger.info(f"Found {len(ct_trials)} ClinicalTrials.gov studies for {condition}")

        # Search PubMed (temporarily disabled to focus on ClinicalTrials.gov)
        pm_studies = []  # self.search_pubmed(condition, max_results=self.limit // 2)
        logger.info(f"Found {len(pm_studies)} PubMed studies for {condition} (disabled)")

        # Convert to NORML format
        new_studies = []

        for trial in ct_trials:
            study_id = f"CT_{trial.nct_id}"
            if study_id not in self.existing_studies:
                study = self.convert_trial_to_norml_format(trial)
                new_studies.append(study)

        for study in pm_studies:
            study_id = f"PM_{study.pmid}"
            if study_id not in self.existing_studies:
                study = self.convert_pubmed_to_norml_format(study)
                new_studies.append(study)

        # Load existing studies for this condition
        existing_file = self.norml_dir / f"{condition.lower().replace(' ', '_')}_studies.json"
        existing_data = {"studies": []}

        if existing_file.exists():
            try:
                with open(existing_file, 'r', encoding='utf-8') as f:
                    existing_data = json.load(f)
            except Exception as e:
                logger.warning(f"Error loading existing data for {condition}: {e}")

        # Combine existing and new studies
        all_studies = existing_data.get('studies', []) + new_studies

        # Update metadata
        expanded_data = {
            "extraction_date": datetime.now().strftime("%Y-%m-%d"),
            "condition": condition.upper().replace(' ', '_'),
            "phase": existing_data.get('phase', '3'),
            "priority": existing_data.get('priority', 'Expansion'),
            "total_studies": len(all_studies),
            "studies": all_studies,
            "expansion_metadata": {
                "clinicaltrials_gov_added": len([s for s in new_studies if s['study_id'].startswith('CT_')]),
                "pubmed_added": len([s for s in new_studies if s['study_id'].startswith('PM_')]),
                "total_added": len(new_studies),
                "expansion_date": datetime.now().isoformat()
            }
        }

        if not self.dry_run:
            output_file = self.output_dir / f"{condition.lower().replace(' ', '_')}_expanded_studies.json"
            with open(output_file, 'w') as f:
                json.dump(expanded_data, f, indent=2, default=str)
            logger.info(f"Saved expanded data to {output_file}")

        return expanded_data

    def run_expansion(self, conditions: List[str] = None):
        """Run expansion for specified conditions or all available."""
        if conditions is None:
            # Get all conditions from existing files
            conditions = []
            for json_file in self.norml_dir.glob('*_studies.json'):
                condition_name = json_file.stem.replace('_studies', '').replace('_', ' ').title()
                conditions.append(condition_name)

        logger.info(f"Starting expansion for {len(conditions)} conditions")
        logger.info(f"Dry run: {self.dry_run}, Limit per source: {self.limit}")

        total_new_studies = 0
        results = {}

        for condition in conditions:
            try:
                result = self.expand_condition(condition)
                new_count = result.get('expansion_metadata', {}).get('total_added', 0)
                total_new_studies += new_count
                results[condition] = {
                    'total_studies': result['total_studies'],
                    'new_studies': new_count
                }
                logger.info(f"{condition}: {result['total_studies']} total, {new_count} new")
            except Exception as e:
                logger.error(f"Error expanding {condition}: {e}")

        logger.info(f"Expansion complete. Total new studies: {total_new_studies}")
        return results


def main():
    parser = argparse.ArgumentParser(description='Automate clinical study expansion')
    parser.add_argument('--dry-run', action='store_true', help='Run without saving changes')
    parser.add_argument('--limit', type=int, default=50, help='Max studies per source per condition')
    parser.add_argument('--conditions', nargs='*', help='Specific conditions to expand')

    args = parser.parse_args()

    expander = ClinicalStudyExpander(dry_run=args.dry_run, limit=args.limit)
    results = expander.run_expansion(args.conditions)

    print("\nExpansion Results:")
    print("=" * 50)
    for condition, stats in results.items():
        print(f"{condition}: {stats['total_studies']} total, {stats['new_studies']} new")


if __name__ == '__main__':
    main()