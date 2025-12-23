import json
import os
import glob
import argparse
from collections import Counter, defaultdict
from datetime import date


def _print_exposure_ascertainment_summary(all_studies):
    """Print counts for exposure_ascertainment across studies.

    This is intended as a lightweight QA helper to reduce self-report bias by making
    exposure capture methods visible (e.g., urine toxicology vs self-report).
    """
    if not all_studies:
        return

    overall = Counter()
    by_condition = defaultdict(Counter)

    for study in all_studies:
        condition = study.get('condition', 'UNKNOWN')
        ascertainment = study.get('exposure_ascertainment', 'missing')
        overall[ascertainment] += 1
        by_condition[condition][ascertainment] += 1

    total = sum(overall.values())
    missing_count = overall.get('missing', 0)
    set_count = total - missing_count

    # If nothing has this field yet (all missing), don't spam output.
    if set_count == 0:
        return

    print("\n" + "=" * 60)
    print("EXPOSURE ASCERTAINMENT SUMMARY")
    print("=" * 60)
    print("Overall (all conditions):")
    print(f"  - total_studies: {total}")
    print(f"  - exposure_ascertainment_set: {set_count}")
    print(f"  - exposure_ascertainment_missing: {missing_count}")
    print("  - distribution (set values only):")

    for key, count in sorted(
        ((k, v) for k, v in overall.items() if k != 'missing'),
        key=lambda kv: (-kv[1], kv[0]),
    ):
        print(f"    - {key}: {count}")

    print("\nBy condition:")
    for condition in sorted(by_condition.keys()):
        counter = by_condition[condition]
        # Skip conditions that never set the field (all missing)
        if len(counter) == 1 and 'missing' in counter:
            continue
        print(f"  {condition}:")
        for key, count in sorted(counter.items(), key=lambda kv: (-kv[1], kv[0])):
            print(f"    - {key}: {count}")

    print("=" * 60 + "\n")


def _infer_exposure_ascertainment_minimal(study):
    """Infer exposure_ascertainment conservatively.

    Minimal ruleset: only tag clear RCTs as trial_assigned.
    """
    study_type = (study.get('study_type') or '').strip().lower()
    study_id = (study.get('study_id') or '').strip().upper()

    # Tag only when the study is clearly randomized/assigned.
    if study_type in {"rct", "randomized_controlled_trial"}:
        return "trial_assigned"
    if "RANDOM" in study_type:
        return "trial_assigned"
    if "_RCT_" in study_id or study_id.endswith("_RCT"):
        return "trial_assigned"

    return None


def apply_min_exposure_tags(
    condition_files,
    tag_conditions,
    dry_run=True,
):
    """Apply minimal exposure_ascertainment tags to source JSON files.

    This is intended to seed a small set of conditions with consistent tags
    before retraining, without risking widespread mislabeling.
    """
    total_updates = 0
    per_condition_updates = Counter()

    for file_path in condition_files:
        if not os.path.exists(file_path):
            continue

        with open(file_path, encoding='utf-8') as f:
            data = json.load(f)

        condition = data.get('condition')
        if condition not in tag_conditions:
            continue

        studies = data.get('studies', [])
        updated = 0

        for study in studies:
            if 'exposure_ascertainment' in study:
                continue
            inferred = _infer_exposure_ascertainment_minimal(study)
            if inferred:
                study['exposure_ascertainment'] = inferred
                updated += 1

        if updated and not dry_run:
            with open(file_path, 'w', encoding='utf-8') as f:
                json.dump(data, f, indent=2)
                f.write("\n")

        total_updates += updated
        per_condition_updates[condition] += updated

    mode = "DRY RUN" if dry_run else "WRITE"
    print("\n" + "=" * 60)
    print(f"MINIMUM EXPOSURE TAGGING ({mode})")
    print("=" * 60)
    print(f"Total studies tagged: {total_updates}")
    for condition in sorted(per_condition_updates.keys()):
        print(f"  - {condition}: {per_condition_updates[condition]}")
    print("=" * 60 + "\n")

    return total_updates

def validate_study_extraction(filename):
    """Quick validation of extracted studies"""
    with open(filename, encoding='utf-8') as f:
        data = json.load(f)
    
    studies = data.get('studies', [])
    
    condition = data.get('condition', os.path.basename(filename))
    phase = data.get('phase', 'unknown')
    extraction_date = data.get('extraction_date', 'unknown')

    print(f"\n{'='*60}")
    print(f"VALIDATION: {condition}")
    print(f"{'='*60}")
    print(f"Total studies: {len(studies)}")
    print(f"Phase: {phase}")
    print(f"Extraction date: {extraction_date}")
    
    # Check required fields
    required_fields = ['study_id', 'study_type', 'condition', 'study_title', 'citation']
    
    complete_count = 0
    for i, study in enumerate(studies, 1):
        missing = [f for f in required_fields if f not in study]
        if missing:
            print(f"[WARN] Study {i}: Missing fields: {missing}")
        else:
            print(f"[OK] Study {i}: {study['study_id']} - Complete")
            complete_count += 1
    
    print(f"{'='*60}")
    print(f"Complete studies: {complete_count}/{len(studies)}")
    print(f"Completion rate: {(complete_count/len(studies)*100):.1f}%")
    print(f"{'='*60}\n")

def merge_phase_data(output_filename='norml_complete_200plus_studies.json', exposure_summary=False):
    """Merge Phase 1, Phase 2, and Phase 3 data"""

    # Load existing Phase 1 legacy data (if exists)
    phase1_legacy_file = 'data/norml_extraction/norml_training_dataset_FINAL.json'

    # Auto-discover all condition extraction files
    condition_files = sorted(
        glob.glob('data/norml_extraction/*_studies.json')
    )
    
    all_studies = []
    conditions_loaded = []
    parsed_extraction_dates = []
    
    # Load Phase 1 legacy file if exists (skip if individual files present)
    if os.path.exists(phase1_legacy_file) and not condition_files:
        with open(phase1_legacy_file, encoding='utf-8') as f:
            phase1 = json.load(f)
            all_studies.extend(phase1.get('studies', []))
            print(f"[OK] Phase 1 (legacy): {len(phase1.get('studies', []))} studies loaded")

    # Load all discovered condition files
    for file_path in condition_files:
        if os.path.exists(file_path):
            with open(file_path, encoding='utf-8') as f:
                data = json.load(f)
                all_studies.extend(data.get('studies', []))
                condition = data.get('condition')
                if condition and condition not in conditions_loaded:
                    conditions_loaded.append(condition)

                extracted_on = data.get('extraction_date')
                if isinstance(extracted_on, str):
                    try:
                        parsed_extraction_dates.append(date.fromisoformat(extracted_on))
                    except ValueError:
                        pass

                print(f"[OK] {data.get('condition', os.path.basename(file_path))}: {len(data.get('studies', []))} studies loaded")
    
    # Create merged dataset
    merged_extraction_date = (
        max(parsed_extraction_dates).isoformat()
        if parsed_extraction_dates
        else date.today().isoformat()
    )

    merged_data = {
        "extraction_date": merged_extraction_date,
        "phases_complete": ["Phase 1", "Phase 2", "Phase 3"],
        "total_conditions": len(conditions_loaded),
        "total_studies": len(all_studies),
        "conditions": conditions_loaded,
        "studies": all_studies
    }
    
    # Save merged dataset
    with open(output_filename, 'w', encoding='utf-8') as f:
        json.dump(merged_data, f, indent=2)
    
    print(f"\n{'='*60}")
    print(f"[OK] MERGED DATASET CREATED: {output_filename}")
    print(f"Total studies: {merged_data['total_studies']}")
    print(f"Conditions: {len(merged_data['conditions'])}")
    print(f"{'='*60}\n")

    if exposure_summary:
        _print_exposure_ascertainment_summary(all_studies)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Validate and merge NORML extraction JSON files")
    parser.add_argument(
        "--exposure-summary",
        action="store_true",
        help="Print counts of exposure_ascertainment values (overall and per condition)",
    )
    parser.add_argument(
        "--apply-min-exposure-tags",
        action="store_true",
        help="Write minimal exposure_ascertainment tags (trial_assigned for clear RCTs) into a small default set of conditions",
    )
    parser.add_argument(
        "--write",
        action="store_true",
        help="When used with --apply-min-exposure-tags, actually writes changes to JSON files (otherwise dry-run)",
    )
    args = parser.parse_args()

    print("\nNORML Study Extraction Validator")
    print("=" * 60)

    # Validate all available condition files
    condition_files = sorted(glob.glob('data/norml_extraction/*_studies.json'))
    if not condition_files:
        print("[WARN] No condition files found in data/norml_extraction/")

    # Optional: apply minimal tags before validation/merge so summary reflects them.
    if args.apply_min_exposure_tags:
        # Minimum tagging set recommendation: seed a couple high-volume, clean conditions.
        default_tag_conditions = {"ANXIETY", "CHRONIC_PAIN"}
        apply_min_exposure_tags(
            condition_files=condition_files,
            tag_conditions=default_tag_conditions,
            dry_run=(not args.write),
        )

    for file_path in condition_files:
        if os.path.exists(file_path):
            validate_study_extraction(file_path)

    # Try to merge if other files exist
    print("\nAttempting to merge all available data...")
    merge_phase_data(exposure_summary=args.exposure_summary)
