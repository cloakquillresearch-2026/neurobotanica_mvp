import argparse
import json
import os
from typing import Any, Dict, List, Tuple

REQUIRED_FIELDS = ["study_id", "study_type", "condition", "study_title", "citation"]


def _infer_study_type(study: Dict[str, Any]) -> str:
    existing = (study.get("study_type") or "").strip()
    if existing:
        return existing

    study_id = (study.get("study_id") or "").upper()
    if "_RCT_" in study_id or study_id.endswith("_RCT"):
        return "RCT"
    if "_OBSERVATIONAL_" in study_id or study_id.endswith("_OBS"):
        return "OBSERVATIONAL"
    if "_SYSTEMATIC_REVIEW_" in study_id or "_META_ANALYSIS_" in study_id:
        return "SYSTEMATIC_REVIEW"
    if "_GUIDELINE_" in study_id:
        return "GUIDELINE"
    if "_PRECLINICAL_" in study_id:
        return "PRECLINICAL"
    if "_MECHANISTIC_" in study_id:
        return "MECHANISTIC"
    if "_CASE_" in study_id:
        return "CASE_SERIES"

    return "UNKNOWN"


def _compose_citation(study: Dict[str, Any]) -> str:
    # Prefer an existing non-empty citation.
    existing = (study.get("citation") or "").strip()
    if existing:
        return existing

    # Build a conservative citation string from available fields.
    title = (study.get("study_title") or study.get("title") or "").strip()
    authors = (study.get("authors") or "").strip()
    year = study.get("year")
    journal = (study.get("journal") or "").strip()

    parts: List[str] = []
    if authors:
        parts.append(authors.rstrip("."))
    if title:
        parts.append(title.rstrip("."))
    if journal:
        parts.append(journal.rstrip("."))
    if year is not None and str(year).strip():
        parts.append(str(year).strip())

    if parts:
        return ". ".join(parts) + "."

    # Absolute fallback: required field must be present, but we avoid inventing details.
    return "Citation missing (legacy record; needs curation)."


def _ensure_required_fields(data: Dict[str, Any]) -> Tuple[int, int]:
    condition = (data.get("condition") or "").strip() or "UNKNOWN"
    studies: List[Dict[str, Any]] = data.get("studies", [])

    updated = 0
    problems = 0

    for index, study in enumerate(studies, 1):
        changed = False

        # condition
        if not (study.get("condition") or "").strip():
            study["condition"] = condition
            changed = True

        # study_title
        if not (study.get("study_title") or "").strip():
            fallback_title = (study.get("title") or "").strip()
            if fallback_title:
                study["study_title"] = fallback_title
            else:
                study["study_title"] = f"Untitled study ({study.get('study_id', f'{condition}_{index:03d}')})"
            changed = True

        # study_type
        if not (study.get("study_type") or "").strip():
            study["study_type"] = _infer_study_type(study)
            changed = True

        # citation
        if not (study.get("citation") or "").strip():
            study["citation"] = _compose_citation(study)
            changed = True

        # verify required fields exist (study_id is assumed present; if missing, flag)
        missing = [f for f in REQUIRED_FIELDS if not (study.get(f) or "").strip()]
        if missing:
            problems += 1

        if changed:
            updated += 1

    return updated, problems


def main() -> int:
    parser = argparse.ArgumentParser(description="Fix missing required fields in *_studies.json files")
    parser.add_argument("files", nargs="+", help="One or more JSON files to fix")
    args = parser.parse_args()

    total_updated = 0
    total_problems = 0

    for path in args.files:
        if not os.path.exists(path):
            raise FileNotFoundError(path)

        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)

        updated, problems = _ensure_required_fields(data)

        with open(path, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2, ensure_ascii=False)
            f.write("\n")

        print(f"[OK] {os.path.basename(path)}: updated {updated} studies; remaining-problem-studies {problems}")
        total_updated += updated
        total_problems += problems

    print(f"[OK] Total updated studies: {total_updated}")
    if total_problems:
        print(f"[WARN] Studies still missing required fields: {total_problems}")
        return 2

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
