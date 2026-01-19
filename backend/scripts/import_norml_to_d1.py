"""
Import NORML clinical studies dataset into Cloudflare D1 database.
Generates SQL INSERT statements and executes via wrangler.
"""

import json
import subprocess
from pathlib import Path
from typing import Dict, Any

# Study quality scoring (used for confidence_score)
STUDY_TYPE_SCORES = {
    'randomized_controlled_trial': 1.0,
    'rct': 1.0,
    'cohort_study': 0.8,
    'cohort': 0.8,
    'case_control': 0.7,
    'case_series': 0.6,
    'case_report': 0.5,
    'observational': 0.6,
    'cross_sectional': 0.6,
    'meta_analysis': 0.9,
    'systematic_review': 0.85
}


def escape_sql_string(value: str) -> str:
    """Escape single quotes for SQL."""
    if value is None:
        return 'NULL'
    return "'" + str(value).replace("'", "''") + "'"


def calculate_confidence_score(study: Dict[str, Any]) -> float:
    """
    Calculate confidence score based on study type and quality.
    
    Returns:
        float: Score between 0.0 and 1.0
    """
    study_type = (study.get('study_type') or '').lower()
    
    # Base score from study type
    base_score = STUDY_TYPE_SCORES.get(study_type, 0.5)
    
    # Adjust based on study quality indicators
    quality = study.get('study_quality', {}) or {}
    adjustments = 0.0
    
    # Peer-reviewed journal bonus
    if quality.get('journal') and quality.get('journal') != 'Not specified':
        adjustments += 0.05
    
    # Recent publication bonus (2020+)
    pub_year = quality.get('publication_year')
    try:
        if pub_year and int(pub_year) >= 2020:
            adjustments += 0.05
    except Exception:
        pass
    
    # Sample size consideration
    design = study.get('study_design', {}) or {}
    sample_size = design.get('sample_size', 0)
    try:
        if isinstance(sample_size, int) and sample_size > 100:
            adjustments += 0.10
        elif isinstance(sample_size, int) and sample_size > 50:
            adjustments += 0.05
    except Exception:
        pass
    
    # Cap at 1.0
    final_score = min(1.0, base_score + adjustments)
    return final_score


def json_to_sql(obj: Any) -> str:
    """Convert Python object to SQL-safe JSON string."""
    if obj is None:
        return 'NULL'
    return escape_sql_string(json.dumps(obj, ensure_ascii=False))


def generate_insert_sql(studies: list) -> str:
    """
    Generate SQL INSERT statements for all studies.
    
    Args:
        studies: List of study dictionaries from NORML JSON
        
    Returns:
        str: SQL statements ready to execute
    """
    sql_statements = []
    
    for study in studies:
        study_id = study.get('study_id') or study.get('id') or ''
        if not study_id:
            continue  # Skip studies without ID
        
        # Calculate confidence score
        confidence = calculate_confidence_score(study)
        
        # Build INSERT statement (use OR IGNORE to avoid UNIQUE constraint failures)
        sql = f"""
    INSERT OR IGNORE INTO clinical_studies (
    study_id,
    study_type,
    condition,
    study_title,
    patient_characteristics,
    intervention,
    outcomes,
    study_design,
    study_quality,
    clinical_significance,
    key_findings,
    citation,
    confidence_score
) VALUES (
    {escape_sql_string(study_id)},
    {escape_sql_string(study.get('study_type'))},
    {escape_sql_string(study.get('condition'))},
    {escape_sql_string(study.get('study_title'))},
    {json_to_sql(study.get('patient_characteristics'))},
    {json_to_sql(study.get('intervention'))},
    {json_to_sql(study.get('outcomes'))},
    {json_to_sql(study.get('study_design'))},
    {json_to_sql(study.get('study_quality'))},
    {escape_sql_string(study.get('clinical_significance'))},
    {json_to_sql(study.get('key_findings'))},
    {escape_sql_string(study.get('citation'))},
    {confidence}
);
""".strip()
        
        sql_statements.append(sql)
    
    return '\n\n'.join(sql_statements)


def generate_conditions_sql(studies: list) -> str:
    """
    Generate conditions table from study data.
    
    Args:
        studies: List of study dictionaries
        
    Returns:
        str: SQL INSERT statements for conditions
    """
    # Aggregate conditions from studies
    conditions_map = {}
    
    for study in studies:
        condition = (study.get('condition') or '').strip()
        if not condition or condition.lower() == 'not specified':
            continue
        
        condition_upper = condition.upper()
        
        if condition_upper not in conditions_map:
            conditions_map[condition_upper] = {
                'name': condition,
                'studies': [],
                'cannabinoids': set(),
                'effects': set()
            }
        
        conditions_map[condition_upper]['studies'].append(study.get('study_id') or study.get('id'))
        
        # Extract cannabinoid info from intervention (handle dict or string)
        intervention = study.get('intervention') or {}
        cannabinoid_profile = ''
        if isinstance(intervention, dict):
            cannabinoid_profile = intervention.get('cannabinoid_profile') or intervention.get('cannabinoid') or ''
        else:
            # Some source records store intervention as a free-form string
            try:
                cannabinoid_profile = str(intervention)
            except Exception:
                cannabinoid_profile = ''

        try:
            profile_up = (cannabinoid_profile or '').upper()
            if 'CBD' in profile_up:
                conditions_map[condition_upper]['cannabinoids'].add('CBD')
            if 'THC' in profile_up:
                conditions_map[condition_upper]['cannabinoids'].add('THC')
        except Exception:
            pass
    
    # Generate INSERT statements
    sql_statements = []
    
    for condition_id, data in conditions_map.items():
        # Determine category from condition name
        condition_lower = condition_id.lower()
        if 'anxiety' in condition_lower or 'ptsd' in condition_lower:
            category = 'anxiety'
        elif 'pain' in condition_lower or 'arthritis' in condition_lower:
            category = 'pain'
        elif 'sleep' in condition_lower or 'insomnia' in condition_lower:
            category = 'sleep'
        elif 'epilepsy' in condition_lower or 'seizure' in condition_lower:
            category = 'epilepsy'
        elif 'depression' in condition_lower:
            category = 'depression'
        else:
            category = 'other'
        
        cannabinoids = list(data['cannabinoids']) if data['cannabinoids'] else []
        
        sql = f"""
INSERT OR IGNORE INTO conditions (
    condition_id,
    condition_name,
    category,
    recommended_cannabinoids,
    evidence_count
) VALUES (
    {escape_sql_string(condition_id)},
    {escape_sql_string(data['name'])},
    {escape_sql_string(category)},
    {json_to_sql(cannabinoids)},
    {len(data['studies'])}
);
""".strip()
        
        sql_statements.append(sql)
    
    return '\n\n'.join(sql_statements)


def main():
    """Main import function."""
    print("🔄 Starting NORML dataset import to D1...")
    
    # Find NORML JSON file (try expected filename, fallback to known file)
    norml_file = Path('norml_training_dataset_FINAL.json')
    if not norml_file.exists():
        alt = Path('..') / 'norml_complete_200plus_studies.json'
        if alt.exists():
            norml_file = alt
        else:
            print(f"❌ Error: Could not find norml_training_dataset_FINAL.json or {alt}")
            return
    
    # Load NORML data (support multiple top-level shapes)
    print(f"📖 Loading {norml_file}...")
    with open(norml_file, 'r', encoding='utf-8') as f:
        data = json.load(f)

    # Normalize to a list of study dicts
    if isinstance(data, dict):
        if 'studies' in data and isinstance(data['studies'], list):
            studies = data['studies']
        elif 'data' in data and isinstance(data['data'], list):
            studies = data['data']
        else:
            # If it's a single study object, wrap it
            # Or if dict maps ids->study objects, take values
            maybe_list = list(data.values())
            if all(isinstance(x, dict) for x in maybe_list):
                studies = maybe_list
            else:
                # Unknown shape: wrap the dict as single study
                studies = [data]
    elif isinstance(data, list):
        studies = data
    else:
        print("❌ Unsupported NORML JSON structure")
        return

    # Filter out non-dict entries (some source files may contain strings/newlines)
    initial_count = len(studies)
    studies = [s for s in studies if isinstance(s, dict)]
    skipped = initial_count - len(studies)
    print(f"✅ Loaded {len(studies)} study objects (skipped {skipped} non-object entries)")
    
    # Generate SQL for studies
    print("🔨 Generating SQL for clinical studies...")
    studies_sql = generate_insert_sql(studies)
    
    # Generate SQL for conditions
    print("🔨 Generating SQL for conditions...")
    conditions_sql = generate_conditions_sql(studies)
    
    # Combine SQL
    full_sql = studies_sql + '\n\n' + conditions_sql
    
    # Write to file
    output_file = Path('backend/db/norml_import.sql')
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    print(f"💾 Writing SQL to {output_file}...")
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(full_sql)
    
    print(f"✅ Generated {output_file.stat().st_size:,} bytes of SQL")
    
    # Execute on D1
    print("☁️  Executing on Cloudflare D1...")
    result = subprocess.run(
        ['npx', 'wrangler', 'd1', 'execute', 'neurobotanica-clinical-evidence', 
         '--file', str(output_file), '--remote', '--config', '../workers/api-proxy/wrangler.toml'],
        capture_output=True,
        text=True
    )
    
    if result.returncode == 0:
        print("✅ Import completed successfully!")
        print(result.stdout)
    else:
        print("❌ Import failed:")
        print(result.stderr)
        return
    
    print("\n🎉 NORML dataset imported to D1!")
    print(f"   - {len(studies)} clinical studies")
    print(f"   - Conditions table populated")
    print(f"   - Ready for budtender recommendations")


if __name__ == '__main__':
    main()
