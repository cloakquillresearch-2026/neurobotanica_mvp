"""Post-harden training correlations by attaching high-confidence citations.

This script links each TK practice in the correlation dataset to NORML extraction
studies so the re-training pipeline can rely on grounded, traceable evidence.
It enriches `evidence_sources` with concrete PMIDs/DOIs and stores structured
matches under `hardened_sources` for downstream auditing.
"""

from __future__ import annotations

import argparse
import glob
import json
import re
from collections import defaultdict
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

PMID_RE = re.compile(r"PMID[:\s]*([0-9]{5,9})", re.IGNORECASE)
DOI_RE = re.compile(r"\b10\.\d{4,9}/[-._;()/:A-Z0-9]+\b", re.IGNORECASE)


CONDITION_ALIASES = {
    # Appetite and metabolic
    "appetite_loss": "appetite_cachexia",
    "appetite_stimulation": "appetite_cachexia",
    "cachexia": "appetite_cachexia",
    "wasting_syndrome": "appetite_cachexia",
    "blood_sugar_regulation": "appetite_cachexia",
    "diabetes": "appetite_cachexia",
    "metabolic_syndrome": "appetite_cachexia",
    # Pain and musculoskeletal
    "joint_pain": "arthritis",
    "inflammation": "arthritis",
    "swelling": "arthritis",
    "arthritis": "arthritis",
    "bone_healing": "arthritis",
    "radiculopathy": "chronic_pain",
    "sciatica": "chronic_pain",
    "migraine": "chronic_pain",
    "headache": "chronic_pain",
    "chronic_pain": "chronic_pain",
    "chronic_widespread_pain": "chronic_pain",
    "neuropathic_pain": "chronic_pain",
    "nerve_pain": "chronic_pain",
    "neuralgia": "chronic_pain",
    "bone_pain": "chronic_pain",
    "fibromyalgia": "chronic_pain",
    "lower_back_pain": "chronic_pain",
    "athletic_performance": "chronic_pain",
    "cancer_pain": "cancer_palliative",
    # Neurology
    "multiple_sclerosis": "multiple_sclerosis",
    "ms": "multiple_sclerosis",
    "muscle_recovery": "multiple_sclerosis",
    "muscle_spasticity": "multiple_sclerosis",
    "muscle_spasms": "multiple_sclerosis",
    "weakness": "multiple_sclerosis",
    "neurological_disorders": "multiple_sclerosis",
    "neuroprotection": "alzheimers",
    "tremor": "parkinsons",
    "tremors": "parkinsons",
    "rigidity": "parkinsons",
    "bradykinesia": "parkinsons",
    "bone_fractures": "arthritis",
    "fractures": "arthritis",
    "parkinsons": "parkinsons",
    "parkinsonism": "parkinsons",
    # Mental health
    "anxiety": "anxiety",
    "nervousness": "anxiety",
    "mood_disorders": "depression",
    "anhedonia": "depression",
    "depression": "depression",
    "fatigue": "post_acute_sequelae_of_sars_cov_2",
    "high_blood_pressure": "post_acute_sequelae_of_sars_cov_2",
    "hypertension": "post_acute_sequelae_of_sars_cov_2",
    "cognitive_decline": "alzheimers",
    "dementia": "alzheimers",
    "ptsd": "ptsd",
    "trauma": "ptsd",
    "nightmares": "ptsd",
    "hypervigilance": "ptsd",
    # Sleep
    "insomnia": "insomnia",
    "sleep_disturbance": "insomnia",
    # Respiratory & inflammatory
    "asthma": "covid_19",
    "bronchitis": "covid_19",
    "bronchospasm": "covid_19",
    "cough": "covid_19",
    "fever": "covid_19",
    "infection": "covid_19",
    "respiratory_distress": "covid_19",
    "respiratory_inflammation": "covid_19",
    # GI
    "ibd": "ibd_crohns",
    "crohns": "ibd_crohns",
    "ulcerative_colitis": "ibd_crohns",
    "constipation": "ibd_crohns",
    "gastrointestinal_health": "ibd_crohns",
    "digestive_health": "ibd_crohns",
    # Nausea
    "nausea": "nausea_chemotherapy",
    "vomiting": "nausea_chemotherapy",
    "chemotherapy_nausea": "nausea_chemotherapy",
    # Ophthalmology
    "glaucoma": "glaucoma",
    "intraocular_pressure": "glaucoma",
    # Dermatology & wound care
    "eczema": "dermatology",
    "psoriasis": "dermatology",
    "dermatitis": "dermatology",
    "dry_skin": "dermatology",
    "skin_damage": "dermatology",
    "wound_healing": "dermatology",
    "burns": "dermatology",
    "bruises": "dermatology",
    # Women’s health
    "dysmenorrhea": "womens_healthcare",
    "menstrual_cramps": "womens_healthcare",
    "pelvic_pain": "womens_healthcare",
    "lactation_support": "womens_healthcare",
    "milk_production": "womens_healthcare",
    "postpartum_recovery": "womens_healthcare",
    # Seizure
    "epilepsy": "epilepsy",
    "seizures": "epilepsy",
    # Tourette
    "tourette": "tourette_syndrome",
    "tourette_syndrome": "tourette_syndrome",
}


def normalize_condition(name: Optional[str]) -> str:
    """Normalize condition tokens so TK indications match extraction files."""

    if not name:
        return ""
    return name.strip().lower().replace(" ", "_")


def resolve_condition_key(name: Optional[str]) -> str:
    normalized = normalize_condition(name)
    return CONDITION_ALIASES.get(normalized, normalized)


def classify_study_tier(study_type: Optional[str]) -> str:
    label = (study_type or "").lower()
    if any(term in label for term in ("systematic", "meta")):
        return "systematic"
    if any(term in label for term in ("rct", "random", "placebo", "double")):
        return "rct"
    if any(term in label for term in ("observational", "cohort", "registry", "case", "real-world")):
        return "observational"
    return "other"


def quality_from_confidence(confidence: float) -> str:
    if confidence >= 0.85:
        return "excellent"
    if confidence >= 0.75:
        return "good"
    if confidence >= 0.60:
        return "moderate"
    return "poor"


def normalize_practice_name(name: Optional[str]) -> str:
    """Lowercase helper used to map predicted practice names back to TK records."""

    if not name:
        return ""
    return re.sub(r"\s+", " ", name.strip()).lower()


def coerce_int(value: Optional[object]) -> Optional[int]:
    """Convert mixed-type year fields into integers when possible."""

    if value is None:
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def extract_pmid(text: str) -> Optional[str]:
    match = PMID_RE.search(text or "")
    if match:
        return match.group(1)
    return None


def extract_doi(text: str) -> Optional[str]:
    match = DOI_RE.search(text or "")
    if match:
        return match.group(0).rstrip(".")
    return None


def study_priority(study_type: Optional[str]) -> int:
    """Give randomized/systematic evidence priority during selection."""

    label = (study_type or "").lower()
    if "rct" in label or "random" in label:
        return 0
    if "systematic" in label or "meta" in label:
        return 1
    if "observational" in label or "cohort" in label or "registry" in label:
        return 2
    if "case" in label:
        return 3
    return 4


@dataclass
class HardenedCitation:
    condition: str
    study_id: Optional[str]
    study_title: Optional[str]
    study_type: Optional[str]
    publication_year: Optional[int]
    citation: str
    pmid: Optional[str]
    doi: Optional[str]

    @property
    def dedupe_key(self) -> str:
        if self.pmid:
            return f"pmid:{self.pmid}"
        if self.doi:
            return f"doi:{self.doi.lower()}"
        return self.citation

    def as_dict(self) -> Dict[str, object]:
        payload = {
            "condition": self.condition,
            "study_id": self.study_id,
            "study_title": self.study_title,
            "study_type": self.study_type,
            "publication_year": self.publication_year,
            "citation": self.citation,
        }
        if self.pmid:
            payload["pmid"] = self.pmid
        if self.doi:
            payload["doi"] = self.doi
        return payload


def load_condition_index(pattern: str) -> Dict[str, List[HardenedCitation]]:
    """Load NORML extraction files and build a citation index per condition."""

    condition_index: Dict[str, List[HardenedCitation]] = {}

    for file_path in glob.glob(pattern):
        path_obj = Path(file_path)
        with path_obj.open(encoding="utf-8") as handle:
            payload = json.load(handle)

        condition_name = normalize_condition(
            payload.get("condition") or path_obj.stem.replace("_studies", "")
        )

        studies: List[HardenedCitation] = []
        for study in payload.get("studies", []):
            raw_citation = (study.get("citation") or "").strip()
            if not raw_citation:
                continue

            hc = HardenedCitation(
                condition=condition_name,
                study_id=study.get("study_id"),
                study_title=study.get("study_title") or study.get("title"),
                study_type=study.get("study_type"),
                publication_year=coerce_int(
                    study.get("publication_year") or study.get("year")
                ),
                citation=raw_citation,
                pmid=extract_pmid(raw_citation),
                doi=extract_doi(raw_citation),
            )
            studies.append(hc)

        studies.sort(
            key=lambda entry: (
                study_priority(entry.study_type),
                -1 * (entry.publication_year or 0)
            )
        )
        condition_index[condition_name] = studies

    return condition_index


def gather_citations(
    indications: Sequence[str],
    condition_index: Dict[str, List[HardenedCitation]],
    max_total: int,
    per_condition: int,
) -> List[HardenedCitation]:
    """Collect top studies per indication until reaching the configured cap."""

    selected: List[HardenedCitation] = []
    seen: set[str] = set()

    for indication in indications:
        condition_key = resolve_condition_key(indication)
        if not condition_key:
            continue

        pulled = 0
        for study in condition_index.get(condition_key, []):
            dedupe = study.dedupe_key
            if dedupe in seen:
                continue

            selected.append(study)
            seen.add(dedupe)
            pulled += 1

            if len(selected) >= max_total:
                return selected
            if pulled >= per_condition:
                break

    return selected


def update_evidence_sources(
    existing: Iterable[str],
    citations: Sequence[HardenedCitation]
) -> List[str]:
    """Merge in new PMID/DOI tokens so legacy references are preserved."""

    updated = set(existing or [])
    for citation in citations:
        if citation.pmid:
            updated.add(f"PMID:{citation.pmid}")
        elif citation.doi:
            updated.add(f"DOI:{citation.doi}")
        else:
            updated.add(citation.citation)
    return sorted(updated)


def apply_confidence_boost(correlation: Dict[str, object]) -> float:
    sources = correlation.get("hardened_sources") or []
    if not sources:
        return 0.0

    tier_counts = defaultdict(int)
    for source in sources:
        tier = classify_study_tier(source.get("study_type"))
        tier_counts[tier] += 1

    trial_like = tier_counts.get("rct", 0)
    systematic_like = tier_counts.get("systematic", 0)
    observational_like = tier_counts.get("observational", 0)

    if trial_like == 0 and systematic_like == 0 and observational_like < 2:
        return 0.0

    base_conf = float(correlation.get("confidence") or 0.0)
    boost = 0.0

    if trial_like:
        boost += min(0.04 * trial_like, 0.08)
    if systematic_like:
        boost += min(0.03 * systematic_like, 0.06)
    if boost == 0.0 and observational_like >= 2:
        boost += 0.02

    if boost <= 0:
        return 0.0

    new_conf = max(base_conf, min(base_conf + boost, 0.85))
    delta = new_conf - base_conf
    if delta <= 0:
        return 0.0

    correlation["confidence"] = new_conf
    correlation["quality"] = quality_from_confidence(new_conf)
    correlation["evidence_confidence_boost"] = round(delta, 4)
    return delta


def build_practice_indexes(tk_practices_path: str) -> tuple[Dict[str, dict], Dict[str, dict]]:
    with open(tk_practices_path, encoding="utf-8") as handle:
        payload = json.load(handle)

    by_id = {practice.get("practice_id"): practice for practice in payload.get("practices", [])}
    by_name = {
        normalize_practice_name(practice.get("practice_name")): practice
        for practice in payload.get("practices", [])
    }
    return by_id, by_name


def post_harden(
    correlations_path: str,
    tk_practices_path: str,
    condition_pattern: str,
    max_total: int,
    per_condition: int,
    output_path: Optional[str],
    dry_run: bool,
) -> None:
    with open(correlations_path, encoding="utf-8") as handle:
        corr_payload = json.load(handle)

    correlations = corr_payload.get("correlations", [])
    metadata = corr_payload.get("metadata", {})

    practice_by_id, practice_by_name = build_practice_indexes(tk_practices_path)
    condition_index = load_condition_index(condition_pattern)

    print(f"Loaded {len(condition_index)} condition files into citation index")

    updated = 0
    skipped = 0
    total_citations = 0
    missing_reasons = defaultdict(int)
    confidence_boosts = 0
    confidence_delta_total = 0.0

    for correlation in correlations:
        direction = correlation.get("direction")
        practice = None

        if direction == "tk_to_genomic":
            practice = practice_by_id.get(correlation.get("tk_practice_id"))
        elif direction == "genomic_to_tk":
            practice = practice_by_name.get(
                normalize_practice_name(correlation.get("tk_practice_predicted"))
            )
            if not practice:
                predicted = correlation.get("tk_practice_predicted", "")
                match = re.search(r"(TK_\d+)", predicted)
                if match:
                    practice = practice_by_id.get(match.group(1))

        if not practice:
            skipped += 1
            missing_reasons["no_practice_match"] += 1
            continue

        indications = practice.get("indications", [])
        if not indications:
            skipped += 1
            missing_reasons["no_indications"] += 1
            continue

        citations = gather_citations(indications, condition_index, max_total, per_condition)
        if not citations:
            skipped += 1
            missing_reasons["no_condition_citations"] += 1
            continue

        correlation["hardened_sources"] = [item.as_dict() for item in citations]
        correlation["evidence_sources"] = update_evidence_sources(
            correlation.get("evidence_sources", []),
            citations,
        )

        delta = apply_confidence_boost(correlation)
        if delta > 0:
            confidence_boosts += 1
            confidence_delta_total += delta

        updated += 1
        total_citations += len(citations)

    metadata["post_hardened_on"] = datetime.now(UTC).isoformat()
    metadata["post_hardened_citations"] = total_citations
    metadata["post_hardened_correlations"] = updated
    metadata["post_hardened_confidence_boosts"] = confidence_boosts
    if confidence_boosts:
        metadata["post_hardened_avg_confidence_delta"] = round(
            confidence_delta_total / confidence_boosts, 4
        )
    corr_payload["metadata"] = metadata

    print("\nPost-hardening summary")
    print("=" * 60)
    print(f"Correlations processed: {len(correlations)}")
    print(f"Correlations hardened: {updated}")
    print(f"Citations attached: {total_citations}")
    print(f"Skipped correlations: {skipped}")
    if missing_reasons:
        print("Skipped breakdown:")
        for reason, count in sorted(missing_reasons.items()):
            print(f"  - {reason}: {count}")
    print("=" * 60)

    if dry_run:
        print("Dry run complete – no files written")
        return

    destination = output_path or correlations_path
    with open(destination, "w", encoding="utf-8") as handle:
        json.dump(corr_payload, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    print(f"✅ Hardened dataset written to {destination}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Attach NORML evidence citations to GenomePath training correlations"
    )
    parser.add_argument(
        "--correlations",
        default="data/processed/training_correlations.json",
        help="Path to training_correlations.json",
    )
    parser.add_argument(
        "--tk-practices",
        default="data/processed/tk_practices.json",
        help="Path to tk_practices.json for indication lookup",
    )
    parser.add_argument(
        "--condition-pattern",
        default="data/norml_extraction/*_studies.json",
        help="Glob for NORML extraction files",
    )
    parser.add_argument(
        "--max-total",
        type=int,
        default=5,
        help="Maximum citations to attach per correlation",
    )
    parser.add_argument(
        "--per-condition",
        type=int,
        default=2,
        help="Maximum citations to pull per indication/condition",
    )
    parser.add_argument(
        "--output",
        help="Optional output path (defaults to in-place overwrite)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Preview changes without writing output",
    )
    args = parser.parse_args()

    post_harden(
        correlations_path=args.correlations,
        tk_practices_path=args.tk_practices,
        condition_pattern=args.condition_pattern,
        max_total=args.max_total,
        per_condition=args.per_condition,
        output_path=args.output,
        dry_run=args.dry_run,
    )


if __name__ == "__main__":
    main()
