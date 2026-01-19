"""
Load 320 Clinical Studies from JSON to Database
Populates the evidence database for FDA Compliance Module

Source: norml_complete_200plus_studies.json (320 studies across 16 conditions)
"""
import json
import sys
import os
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from sqlalchemy.orm import Session
from backend.models.database import engine, SessionLocal, Base
from backend.models.study import ClinicalStudy


def parse_effect_size(effect_size_str: str) -> float:
    """Extract numeric effect size from string like 'Large (d = 0.91)'."""
    if not effect_size_str:
        return None
    
    # Look for patterns like "d = 0.91", "d=1.9", "SMD=0.82"
    import re
    patterns = [
        r'd\s*=\s*([\d.]+)',
        r'SMD\s*=\s*([\d.]+)',
        r'ES\s*=\s*([\d.]+)',
        r'\((\d+\.?\d*)\)'
    ]
    
    for pattern in patterns:
        match = re.search(pattern, effect_size_str)
        if match:
            try:
                return float(match.group(1))
            except ValueError:
                continue
    
    return None


def determine_evidence_grade(study_data: dict) -> str:
    """Determine evidence grade based on study characteristics."""
    study_type = study_data.get("study_type", "").upper()
    quality = study_data.get("quality_metrics", {})
    blinding = quality.get("blinding", "")
    
    if study_type == "RCT":
        if "double-blind" in blinding.lower():
            return "Level 1"  # High-quality RCT
        else:
            return "Level 2"  # Lower-quality RCT
    elif study_type in ["SYSTEMATIC REVIEW", "META-ANALYSIS"]:
        return "Level 1"
    elif study_type == "OBSERVATIONAL":
        return "Level 3"
    elif study_type in ["CASE SERIES", "CASE REPORT"]:
        return "Level 4"
    else:
        return "Level 3"  # Default


def is_fda_pivotal(study_data: dict) -> bool:
    """Check if study is a pivotal FDA trial."""
    relevance = study_data.get("regulatory_relevance", "").lower()
    title = study_data.get("study_title", "").lower()
    
    keywords = ["pivotal", "phase 3", "phase iii", "fda approval", "breakthrough", "epidiolex", "marinol"]
    
    return any(kw in relevance or kw in title for kw in keywords)


def extract_fda_drug(study_data: dict) -> str:
    """Extract FDA drug reference from study data."""
    intervention = study_data.get("intervention", {})
    if isinstance(intervention, str):
        cannabinoid = intervention
    else:
        cannabinoid = intervention.get("cannabinoid", "") or ""
    title = study_data.get("study_title", "")
    relevance = study_data.get("regulatory_relevance", "")
    
    all_text = f"{cannabinoid} {title} {relevance}".lower()
    
    if "epidiolex" in all_text or "cbd" in cannabinoid.lower():
        return "Epidiolex"
    elif "marinol" in all_text or "dronabinol" in all_text:
        return "Marinol"
    elif "cesamet" in all_text or "nabilone" in all_text:
        return "Cesamet"
    elif "sativex" in all_text or "nabiximols" in all_text:
        return "Sativex"
    
    return None


def load_studies(json_path: str, db: Session) -> int:
    """Load studies from JSON file into database."""
    print(f"📂 Loading studies from: {json_path}")
    
    with open(json_path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    studies = data.get("studies", [])
    print(f"📊 Found {len(studies)} studies in JSON")
    
    loaded_count = 0
    skipped_count = 0
    seen_ids = set()  # Track IDs within this batch
    
    def to_json_str(val):
        """Convert list/dict to JSON string for Text columns."""
        if val is None:
            return None
        if isinstance(val, (list, dict)):
            return json.dumps(val)
        return str(val)
    
    for study_data in studies:
        study_id = study_data.get("study_id")
        
        # Skip if already seen in this batch
        if study_id in seen_ids:
            skipped_count += 1
            continue
        seen_ids.add(study_id)
        
        # Check if already exists in database
        existing = db.query(ClinicalStudy).filter(
            ClinicalStudy.study_id == study_id
        ).first()
        
        if existing:
            skipped_count += 1
            continue
        
        # Extract intervention details
        intervention = study_data.get("intervention", {})
        # Normalize intervention when it's a simple string in some records
        if isinstance(intervention, str):
            intervention = {"cannabinoid": intervention}
        outcomes = study_data.get("outcomes", {})
        # Normalize possible list values to dict for consistency
        if not isinstance(outcomes, dict):
            outcomes = {"results": outcomes} if outcomes else {}
        safety = study_data.get("safety", {})
        if not isinstance(safety, dict):
            safety = {}
        quality = study_data.get("quality_metrics", {})
        if not isinstance(quality, dict):
            quality = {}
        
        # Extract condition from study_id if not present (e.g., DEPRESSION_RCT_001 -> DEPRESSION)
        condition = study_data.get("condition")
        if not condition and study_id:
            # Study ID format: CONDITION_TYPE_NUMBER
            parts = study_id.split('_')
            if len(parts) >= 1:
                condition = parts[0].replace('_', ' ')
        # Extract study_type from study_id if not present (e.g., DEPRESSION_RCT_001 -> RCT)
        study_type = study_data.get("study_type")
        if not study_type and study_id:
            parts = study_id.split('_')
            if len(parts) >= 2:
                type_map = {
                    'RCT': 'randomized_controlled_trial',
                    'OBSERVATIONAL': 'observational',
                    'SYSTEMATIC': 'systematic_review',
                    'META': 'meta_analysis',
                    'CASE': 'case_series',
                    'GUIDELINE': 'guideline',
                    'MECHANISTIC': 'mechanistic'
                }
                study_type = type_map.get(parts[1], parts[1])
        
        # Create study record
        study = ClinicalStudy(
            study_id=study_id,
            condition=condition or "UNKNOWN",
            study_type=study_type or "unknown",
            study_title=study_data.get("study_title") or study_data.get("title") or f"Study {study_id}",
            
            # Citation
            citation=study_data.get("citation"),
            authors=to_json_str(study_data.get("authors")),
            year=study_data.get("year"),
            journal=study_data.get("journal"),
            
            # Population
            sample_size=study_data.get("sample_size"),
            
            # Intervention
            cannabinoid=intervention.get("cannabinoid"),
            dosage=intervention.get("dosage"),
            duration=intervention.get("duration"),
            delivery_method=intervention.get("delivery_method"),
            intervention_details=intervention,
            
            # Outcomes
            primary_measure=to_json_str(outcomes.get("primary_measure")),
            results_summary=to_json_str(outcomes.get("results")),
            effect_size=outcomes.get("effect_size"),
            effect_size_numeric=parse_effect_size(outcomes.get("effect_size")),
            secondary_outcomes=to_json_str(outcomes.get("secondary_outcomes")),
            outcomes_json=outcomes,
            
            # Safety
            adverse_events=to_json_str(safety.get("adverse_events")),
            serious_adverse_events=to_json_str(safety.get("serious_adverse_events")),
            dropout_rate=safety.get("dropout_rate"),
            safety_json=safety,
            
            # Quality
            randomization_method=quality.get("randomization"),
            blinding=quality.get("blinding"),
            funding_source=quality.get("funding_source"),
            conflicts_of_interest=to_json_str(quality.get("conflicts_of_interest")),
            quality_metrics_json=quality,
            
            # FDA Relevance
            regulatory_relevance=to_json_str(study_data.get("regulatory_relevance")),
            priority_rating=study_data.get("priority"),
            is_pivotal_trial=is_fda_pivotal(study_data),
            fda_drug_reference=extract_fda_drug(study_data),
            
            # Evidence grading
            evidence_grade=determine_evidence_grade(study_data),
            confidence_weight=0.7  # Default, can be adjusted
        )
        
        db.add(study)
        loaded_count += 1
        
        if loaded_count % 50 == 0:
            print(f"  ✅ Loaded {loaded_count} studies...")
            db.commit()
    
    db.commit()
    print(f"✅ Loaded {loaded_count} new studies (skipped {skipped_count} existing)")
    
    return loaded_count


def main():
    """Main entry point."""
    print("🌿 NeuroBotanica Clinical Study Loader")
    print("=" * 50)
    
    # Initialize database
    print("📦 Creating database tables...")
    Base.metadata.create_all(bind=engine)
    
    # Find JSON file
    project_root = Path(__file__).parent.parent
    json_paths = [
        project_root / "norml_complete_200plus_studies.json",
        project_root / "data" / "norml_extraction" / "norml_complete_200plus_studies.json"
    ]
    
    json_path = None
    for path in json_paths:
        if path.exists():
            json_path = path
            break
    
    if not json_path:
        print("❌ Could not find norml_complete_200plus_studies.json")
        print("   Looked in:")
        for path in json_paths:
            print(f"   - {path}")
        sys.exit(1)
    
    # Load studies
    db = SessionLocal()
    try:
        loaded = load_studies(str(json_path), db)
        
        # Print summary
        total = db.query(ClinicalStudy).count()
        conditions = db.query(ClinicalStudy.condition).distinct().count()
        
        print("\n" + "=" * 50)
        print("📊 DATABASE SUMMARY")
        print("=" * 50)
        print(f"Total studies: {total}")
        print(f"Conditions: {conditions}")
        
        # Condition breakdown
        from sqlalchemy import func
        condition_counts = db.query(
            ClinicalStudy.condition,
            func.count(ClinicalStudy.id)
        ).group_by(ClinicalStudy.condition).all()
        
        print("\nStudies by condition:")
        for condition, count in sorted(condition_counts, key=lambda x: -x[1]):
            print(f"  {condition}: {count}")
        
    finally:
        db.close()
    
    print("\n✅ Data loading complete!")


if __name__ == "__main__":
    main()
