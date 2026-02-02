-- Dispensary tables for NeuroBotanica Budtender App

CREATE TABLE IF NOT EXISTS dispensary_profiles (
    profile_id TEXT PRIMARY KEY,
    profile_code TEXT UNIQUE,
    created_at TEXT DEFAULT (datetime('now')),
    updated_at TEXT DEFAULT (datetime('now')),
    completeness_score REAL DEFAULT 0.0,
    primary_condition TEXT,
    data TEXT -- JSON blob for profile data
);

CREATE TABLE IF NOT EXISTS dispensary_recommendations (
    recommendation_id TEXT PRIMARY KEY,
    profile_id TEXT,
    recommendations TEXT, -- JSON array of recommendations
    clinical_studies_referenced INTEGER DEFAULT 0,
    generated_at TEXT DEFAULT (datetime('now'))
);

CREATE TABLE IF NOT EXISTS dispensary_feedback (
    feedback_id TEXT PRIMARY KEY,
    recommendation_id TEXT,
    feedback TEXT, -- JSON feedback data
    submitted_at TEXT DEFAULT (datetime('now')),
    FOREIGN KEY (recommendation_id) REFERENCES dispensary_recommendations(recommendation_id)
);

CREATE TABLE IF NOT EXISTS dispensary_transactions (
    transaction_id TEXT PRIMARY KEY,
    customer_id TEXT,
    created_at TEXT DEFAULT (datetime('now')),
    total_amount REAL,
    products_json TEXT,
    notes TEXT,
    status TEXT DEFAULT 'completed'
);

CREATE TABLE IF NOT EXISTS inflammatory_profiles (
    profile_id TEXT PRIMARY KEY,
    profile_code TEXT UNIQUE,
    profile_type TEXT DEFAULT 'inflammatory',
    created_at TEXT DEFAULT (datetime('now')),
    completeness_score REAL DEFAULT 0.0,
    primary_condition TEXT,
    biomarkers TEXT, -- JSON biomarkers
    data TEXT -- JSON additional data
);