#!/usr/bin/env python3
"""
Merge Expanded Studies Back into NORML Extraction Files
=======================================================

Merges the expanded studies from data/expanded_studies/ back into the
original NORML extraction files in data/norml_extraction/.

This ensures the integration script picks up the new studies.
"""

import json
from pathlib import Path
from typing import Dict, List
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def merge_expanded_studies():
    """Merge expanded studies back into original NORML files."""

    project_root = Path(__file__).parent.parent
    norml_dir = project_root / 'data' / 'norml_extraction'
    expanded_dir = project_root / 'data' / 'expanded_studies'

    logger.info("Merging expanded studies back into NORML extraction files...")

    total_merged = 0

    # Manual mapping of expanded filenames to original filenames
    filename_mapping = {
        'alzheimer_expanded_studies.json': 'alzheimers_studies.json',
        'alzheimers_expanded_studies.json': 'alzheimers_studies.json',
        'parkinson_expanded_studies.json': 'parkinsons_studies.json',
        'autism_expanded_studies.json': 'autism_spectrum_disorder_studies.json',
        'back_pain_expanded_studies.json': None,  # No original file
        'cancer_expanded_studies.json': None,  # No original file
        'fibromyalgia_expanded_studies.json': None,  # No original file
        'hiv_expanded_studies.json': None,  # No original file
        'migraine_expanded_studies.json': None,  # No original file
        'neck_pain_expanded_studies.json': None,  # No original file
        'neuropathic_pain_expanded_studies.json': None,  # No original file
        'schizophrenia_expanded_studies.json': None,  # No original file
        'ibd_crohn_expanded_studies.json': 'ibd_crohns_studies.json',
    }

    # Process each expanded file
    for expanded_file in expanded_dir.glob('*_expanded_studies.json'):
        # Extract condition name from filename
        expanded_filename = expanded_file.name

        # Find corresponding original file
        original_filename = filename_mapping.get(expanded_filename)
        if original_filename is None:
            # Try to find by pattern matching
            condition_base = expanded_filename.replace('_expanded_studies.json', '')
            original_file = None
            for orig_file in norml_dir.glob('*_studies.json'):
                orig_base = orig_file.name.replace('_studies.json', '')
                if condition_base in orig_base or orig_base in condition_base:
                    original_file = orig_file
                    break
            if not original_file:
                logger.warning(f"No original file found for condition: {condition_base}")
                continue
        else:
            original_file = norml_dir / original_filename
            if not original_file.exists():
                logger.warning(f"Mapped original file does not exist: {original_filename}")
                continue

        # Load original data
        with open(original_file, 'r', encoding='utf-8') as f:
            original_data = json.load(f)

        # Load expanded data
        with open(expanded_file, 'r', encoding='utf-8') as f:
            expanded_data = json.load(f)

        # Get existing study IDs (check multiple ID fields)
        existing_ids = set()
        for study in original_data.get('studies', []):
            # Check various ID fields that might exist
            study_id = (study.get('nct_id') or study.get('pmid') or study.get('study_id') or '')
            if study_id:
                existing_ids.add(study_id)

        # Add new studies
        new_studies = []
        for study in expanded_data.get('studies', []):
            study_id = study.get('nct_id', study.get('pmid', study.get('study_id', '')))
            if study_id and study_id not in existing_ids:
                new_studies.append(study)
                existing_ids.add(study_id)

        if new_studies:
            original_data['studies'].extend(new_studies)
            original_data['total_studies'] = len(original_data['studies'])

            # Update metadata
            if 'expansion_metadata' not in original_data:
                original_data['expansion_metadata'] = {}
            original_data['expansion_metadata']['last_expanded'] = expanded_data.get('expansion_metadata', {}).get('timestamp')
            original_data['expansion_metadata']['total_expanded'] = len(new_studies)

            # Save updated original file
            with open(original_file, 'w', encoding='utf-8') as f:
                json.dump(original_data, f, indent=2, ensure_ascii=False)

            logger.info(f"{expanded_filename}: Added {len(new_studies)} new studies")
            total_merged += len(new_studies)
        else:
            logger.info(f"{expanded_filename}: No new studies to add")

    logger.info(f"Merge complete. Total studies added: {total_merged}")

if __name__ == '__main__':
    merge_expanded_studies()