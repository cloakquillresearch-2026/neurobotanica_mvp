-- Migration: create clinical_studies + conditions tables populated from legacy data

-- Stub table to satisfy foreign key dependencies in later migrations.
CREATE TABLE IF NOT EXISTS omnipath_consent_artifacts (
  consent_id TEXT PRIMARY KEY,
  community_id TEXT,
  knowledge_type TEXT,
  consent_status TEXT,
  consent_scope TEXT,
  tk_labels TEXT,
  bc_labels TEXT,
  signature_rsa4096 TEXT,
  signature_ecdsa TEXT,
  manifest_hash TEXT,
  community_vote_pct REAL,
  revocation_authority TEXT,
  granted_at TEXT,
  expires_at TEXT,
  revoked_at TEXT,
  last_verified_at TEXT,
  created_at TEXT
);

-- Ensure the legacy neurobotanica_clinical_studies table exists so the copy
-- statements below do not fail even on empty/new databases.
CREATE TABLE IF NOT EXISTS neurobotanica_clinical_studies (
  study_id TEXT,
  study_type TEXT,
  condition TEXT,
  intervention TEXT,
  outcomes TEXT,
  key_findings TEXT,
  citation TEXT,
  pubmed_id TEXT,
  doi TEXT,
  confidence_score REAL,
  sample_size INTEGER,
  study_duration_days TEXT,
  publication_year INTEGER,
  consent_id TEXT,
  traditional_knowledge_source INTEGER,
  community_attribution TEXT,
  created_at TEXT
);

-- Clinical studies table mirrors the schema expected by the worker.
CREATE TABLE IF NOT EXISTS clinical_studies (
  study_id TEXT PRIMARY KEY,
  study_type TEXT,
  condition TEXT,
  intervention TEXT,
  outcomes TEXT,
  key_findings TEXT,
  citation TEXT,
  confidence_score REAL,
  sample_size INTEGER,
  study_duration_days TEXT,
  publication_year INTEGER,
  created_at TEXT DEFAULT (datetime('now'))
);

-- Populate from the neurobotanica_clinical_studies dataset if present.
INSERT OR IGNORE INTO clinical_studies (
  study_id,
  study_type,
  condition,
  intervention,
  outcomes,
  key_findings,
  citation,
  confidence_score,
  sample_size,
  study_duration_days,
  publication_year,
  created_at
)
SELECT
  study_id,
  study_type,
  condition,
  intervention,
  outcomes,
  key_findings,
  citation,
  confidence_score,
  sample_size,
  study_duration_days,
  publication_year,
  COALESCE(created_at, datetime('now'))
FROM neurobotanica_clinical_studies;

CREATE INDEX IF NOT EXISTS idx_clinical_studies_condition ON clinical_studies(condition);
CREATE INDEX IF NOT EXISTS idx_clinical_studies_type ON clinical_studies(study_type);

-- Conditions reference table used by the frontend + worker.
CREATE TABLE IF NOT EXISTS conditions (
  condition_id TEXT PRIMARY KEY,
  condition_name TEXT NOT NULL,
  category TEXT,
  recommended_cannabinoids TEXT,
  evidence_count INTEGER DEFAULT 0,
  description TEXT,
  created_at TEXT DEFAULT (datetime('now')),
  updated_at TEXT DEFAULT (datetime('now'))
);

-- Refresh condition metadata from the clinical_studies table.
INSERT INTO conditions (
  condition_id,
  condition_name,
  category,
  recommended_cannabinoids,
  evidence_count,
  description,
  updated_at
)
SELECT
  LOWER(REPLACE(TRIM(condition), ' ', '_')) AS condition_id,
  condition AS condition_name,
  CASE
    WHEN LOWER(condition) LIKE '%pain%' OR LOWER(condition) LIKE '%arthritis%' THEN 'pain'
    WHEN LOWER(condition) LIKE '%anxiety%' OR LOWER(condition) LIKE '%stress%' THEN 'anxiety'
    WHEN LOWER(condition) LIKE '%sleep%' OR LOWER(condition) LIKE '%insomnia%' THEN 'sleep'
    WHEN LOWER(condition) LIKE '%seizure%' OR LOWER(condition) LIKE '%epilepsy%' THEN 'neurological'
    WHEN LOWER(condition) LIKE '%ptsd%' THEN 'trauma'
    WHEN LOWER(condition) LIKE '%glaucoma%' THEN 'ophthalmic'
    WHEN LOWER(condition) LIKE '%inflammation%' OR LOWER(condition) LIKE '%ibd%' THEN 'inflammation'
    ELSE 'general'
  END AS category,
  'CBD,THC' AS recommended_cannabinoids,
  COUNT(*) AS evidence_count,
  'Auto-generated from clinical evidence corpus' AS description,
  datetime('now') AS updated_at
FROM clinical_studies
WHERE TRIM(condition) != ''
GROUP BY condition;

CREATE INDEX IF NOT EXISTS idx_conditions_name ON conditions(condition_name);
CREATE INDEX IF NOT EXISTS idx_conditions_category ON conditions(category);
