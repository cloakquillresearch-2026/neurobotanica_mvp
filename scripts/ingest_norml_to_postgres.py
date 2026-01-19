"""Ingest NORML JSON dataset into a PostgreSQL database.

Usage:
  Set environment variable `DATABASE_URL` to a valid SQLAlchemy URL
  (e.g. postgres://user:pass@host:5432/dbname) and run:

    python scripts/ingest_norml_to_postgres.py --file norml_complete_200plus_studies.json

If `DATABASE_URL` is not set, the script will print instructions and exit.
"""
import os
import json
import argparse
from typing import Any, Dict
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.exc import IntegrityError

PROJECT_ROOT = os.path.dirname(os.path.dirname(__file__))
DEFAULT_JSON = os.path.join(PROJECT_ROOT, "norml_complete_200plus_studies.json")


def normalize_record(rec: Dict[str, Any]) -> Dict[str, Any]:
    # Normalize keys, ensure required fields exist, coerce types
    out = {}
    out["study_id"] = str(rec.get("study_id") or rec.get("id") or os.urandom(6).hex())
    out["condition"] = (rec.get("condition") or rec.get("indication") or "UNKNOWN").upper()
    out["study_type"] = (rec.get("study_type") or rec.get("design") or "OBSERVATIONAL")
    out["study_title"] = rec.get("title") or rec.get("study_title")
    out["year"] = int(rec.get("year")) if rec.get("year") else None
    out["cannabinoid"] = rec.get("cannabinoid") or rec.get("intervention") or None
    out["sample_size"] = rec.get("sample_size") or rec.get("n") or None
    out["results_summary"] = rec.get("results") or rec.get("summary") or None
    out["evidence_grade"] = rec.get("evidence_grade") or None
    out["is_pivotal_trial"] = bool(rec.get("is_pivotal_trial") or rec.get("pivotal"))
    out["confidence_weight"] = rec.get("confidence_weight") or None
    return out


def ingest(json_file: str, database_url: str):
    engine = create_engine(database_url)
    Session = sessionmaker(bind=engine)

    # Import here to ensure models are bound to engine metadata
    from backend.models.study import ClinicalStudy

    session = Session()
    added = 0
    skipped = 0
    with open(json_file, "r", encoding="utf-8") as fh:
        data = json.load(fh)
        for rec in data:
            norm = normalize_record(rec)
            cs = ClinicalStudy(
                study_id=norm["study_id"],
                condition=norm["condition"],
                study_type=norm["study_type"],
                study_title=norm.get("study_title"),
                year=norm.get("year"),
                cannabinoid=norm.get("cannabinoid"),
                sample_size=norm.get("sample_size"),
                results_summary=norm.get("results_summary"),
                evidence_grade=norm.get("evidence_grade"),
                is_pivotal_trial=norm.get("is_pivotal_trial", False),
                confidence_weight=norm.get("confidence_weight")
            )
            try:
                session.add(cs)
                session.commit()
                added += 1
            except IntegrityError:
                session.rollback()
                skipped += 1
            except Exception:
                session.rollback()
                skipped += 1

    session.close()
    print(f"Ingestion complete. Added={added}, Skipped={skipped}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", default=DEFAULT_JSON, help="Path to NORML JSON file")
    args = parser.parse_args()

    database_url = os.environ.get("DATABASE_URL")
    if not database_url:
        print("Error: Please set DATABASE_URL environment variable to your Postgres connection string.")
        return

    if not os.path.exists(args.file):
        print(f"JSON file not found: {args.file}")
        return

    ingest(args.file, database_url)


if __name__ == "__main__":
    main()
