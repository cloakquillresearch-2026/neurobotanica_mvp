-- Clinical Studies Table (505+ NORML studies)
CREATE TABLE IF NOT EXISTS clinical_studies (
    study_id TEXT PRIMARY KEY,
    study_type TEXT NOT NULL,
    condition TEXT NOT NULL,
    study_title TEXT,
    patient_characteristics TEXT, -- JSON blob
    intervention TEXT, -- JSON blob: cannabinoid_profile, delivery_method, dosing
    outcomes TEXT, -- JSON blob: clinical_observations, symptom_changes
    study_design TEXT, -- JSON blob: sample_size, duration, assessment_method
    study_quality TEXT, -- JSON blob: institution, journal, publication_year
    clinical_significance TEXT,
    key_findings TEXT, -- JSON array
    citation TEXT,
    confidence_score REAL DEFAULT 0.5, -- 0.0-1.0 based on study quality
    created_at TEXT DEFAULT (datetime('now'))
);

-- Therapeutic Conditions (mapped from studies)
CREATE TABLE IF NOT EXISTS conditions (
    condition_id TEXT PRIMARY KEY,
    condition_name TEXT NOT NULL UNIQUE,
    category TEXT, -- anxiety, pain, sleep, epilepsy, ptsd, etc.
    description TEXT,
    recommended_cannabinoids TEXT, -- JSON array: ["CBD", "THC", "CBG"]
    recommended_terpenes TEXT, -- JSON array: ["limonene", "myrcene"]
    recommended_ratios TEXT, -- JSON: {"CBD:THC": "2:1", "terpenes": "5%"}
    evidence_count INTEGER DEFAULT 0, -- number of studies supporting this condition
    created_at TEXT DEFAULT (datetime('now'))
);

-- Cannabinoid Profiles (for strain recommendations)
CREATE TABLE IF NOT EXISTS cannabinoid_profiles (
    profile_id TEXT PRIMARY KEY,
    profile_name TEXT NOT NULL,
    thc_percent REAL,
    cbd_percent REAL,
    cbg_percent REAL,
    cbn_percent REAL,
    thcv_percent REAL,
    dominant_terpenes TEXT, -- JSON array
    effects TEXT, -- JSON array: ["relaxing", "uplifting", "pain-relief"]
    conditions_treated TEXT, -- JSON array of condition_ids
    created_at TEXT DEFAULT (datetime('now'))
);

-- Recommendation Log (track budtender queries)
CREATE TABLE IF NOT EXISTS recommendations (
    rec_id TEXT PRIMARY KEY,
    user_condition TEXT NOT NULL,
    user_severity TEXT, -- mild, moderate, severe
    user_preferences TEXT, -- JSON: delivery_method, flavor, experience_level
    recommended_profile_id TEXT,
    recommended_cannabinoids TEXT, -- JSON
    recommended_terpenes TEXT, -- JSON
    confidence_score REAL,
    studies_cited TEXT, -- JSON array of study_ids
    timestamp TEXT DEFAULT (datetime('now')),
    FOREIGN KEY (recommended_profile_id) REFERENCES cannabinoid_profiles(profile_id)
);

-- Indexes for fast queries
CREATE INDEX IF NOT EXISTS idx_studies_condition ON clinical_studies(condition);
CREATE INDEX IF NOT EXISTS idx_studies_type ON clinical_studies(study_type);
CREATE INDEX IF NOT EXISTS idx_studies_confidence ON clinical_studies(confidence_score);
CREATE INDEX IF NOT EXISTS idx_conditions_category ON conditions(category);
CREATE INDEX IF NOT EXISTS idx_recommendations_timestamp ON recommendations(timestamp);
