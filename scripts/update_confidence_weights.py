"""Compute and persist confidence weights for ClinicalStudy records.

Usage:
  Set DATABASE_URL to point to your DB (Postgres), then run:
    python scripts/update_confidence_weights.py

This computes `compute_confidence_for_study` for each record and writes
the value into the `confidence_weight` column for downstream use.
"""
import os
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from backend.services.confidence import compute_confidence_for_study


def main():
    database_url = os.environ.get("DATABASE_URL")
    if not database_url:
        print("Please set DATABASE_URL to your Postgres database URL.")
        return

    engine = create_engine(database_url)
    Session = sessionmaker(bind=engine)

    from backend.models.study import ClinicalStudy

    session = Session()
    updated = 0
    for s in session.query(ClinicalStudy).all():
        try:
            conf = compute_confidence_for_study(s)
            s.confidence_weight = float(conf)
            session.add(s)
            updated += 1
        except Exception:
            session.rollback()
            continue

    session.commit()
    session.close()
    print(f"Updated confidence_weight for {updated} studies")


if __name__ == "__main__":
    main()
