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

-- ============================================
-- NeuroBotanica Engine Tables
-- ============================================

-- Drug interaction data
CREATE TABLE IF NOT EXISTS neurobotanica_drug_interactions (
    interaction_id TEXT PRIMARY KEY,
    compound_id TEXT NOT NULL,
    drug_name TEXT NOT NULL,
    drug_class TEXT,
    cyp450_enzyme TEXT,
    interaction_type TEXT,
    severity TEXT,
    evidence_level TEXT,
    clinical_notes TEXT,
    requires_consent INTEGER,
    monitoring_parameters TEXT,
    dose_adjustment_required INTEGER,
    dose_adjustment_recommendation TEXT,
    evidence_quality TEXT,
    evidence_source TEXT,
    created_at TEXT DEFAULT (datetime('now')),
    requires_tk_consent INTEGER,
    tk_consent_id TEXT
);

-- Demographic correction factors
CREATE TABLE IF NOT EXISTS neurobotanica_demographic_factors (
    factor_id TEXT PRIMARY KEY,
    compound_id TEXT NOT NULL,
    demographic_group TEXT NOT NULL,
    adjustment_factor REAL NOT NULL,
    confidence_interval TEXT,
    evidence_basis TEXT,
    created_at TEXT DEFAULT (datetime('now'))
);

-- Synergy predictions
CREATE TABLE IF NOT EXISTS neurobotanica_synergy_predictions (
    synergy_id TEXT PRIMARY KEY,
    compound_a_id TEXT NOT NULL,
    compound_b_id TEXT NOT NULL,
    synergy_score REAL,
    synergy_mechanism TEXT,
    molecular_interaction TEXT,
    receptor_affinity_change REAL,
    pathway_modulation TEXT,
    clinical_evidence TEXT,
    computational_model TEXT,
    confidence_score REAL,
    uncertainty_range TEXT,
    requires_consent INTEGER,
    consent_verification_status TEXT,
    computation_date TEXT,
    tk_enhanced INTEGER DEFAULT 0,
    created_at TEXT DEFAULT (datetime('now'))
);

-- Plant formulations
CREATE TABLE IF NOT EXISTS neurobotanica_formulations (
    formulation_id TEXT PRIMARY KEY,
    formulation_name TEXT,
    formulation_type TEXT,
    primary_compound TEXT,
    target_condition TEXT,
    cannabinoid_ratio TEXT,
    terpene_profile TEXT,
    efficacy_score REAL,
    safety_profile TEXT,
    dosing_recommendations TEXT,
    administration_routes TEXT,
    contraindications TEXT,
    drug_interactions TEXT,
    clinical_studies TEXT,
    patient_satisfaction REAL,
    requires_consent INTEGER,
    consent_ids TEXT,
    tk_labels TEXT,
    bc_benefit_sharing_pct REAL,
    indigenous_sourcing_verified INTEGER,
    community_id TEXT,
    created_at TEXT DEFAULT (datetime('now')),
    updated_at TEXT DEFAULT (datetime('now')),
    research_status TEXT
);

-- Compound data
CREATE TABLE IF NOT EXISTS neurobotanica_compounds (
    compound_id TEXT PRIMARY KEY,
    compound_name TEXT NOT NULL,
    compound_class TEXT,
    source TEXT,
    smiles TEXT,
    molecular_weight REAL,
    structure_2d TEXT,
    structure_3d TEXT,
    therapeutic_targets TEXT,
    indications TEXT,
    bioavailability REAL,
    half_life_hours REAL,
    metabolism_pathways TEXT,
    excretion_routes TEXT,
    toxicity_profile TEXT,
    therapeutic_effects TEXT,
    requires_consent INTEGER,
    tk_consent_id TEXT,
    bc_benefit_sharing_pct REAL,
    indigenous_sourcing_verified INTEGER,
    community_id TEXT,
    patent_status TEXT,
    regulatory_status TEXT,
    created_at TEXT DEFAULT (datetime('now')),
    updated_at TEXT DEFAULT (datetime('now'))
);

-- Polysaccharides data
CREATE TABLE IF NOT EXISTS neurobotanica_polysaccharides (
    polysaccharide_id TEXT PRIMARY KEY,
    compound_id TEXT,
    polysaccharide_name TEXT NOT NULL,
    source_organism TEXT,
    chemical_structure TEXT,
    molecular_weight_range TEXT,
    solubility_profile TEXT,
    fermentation_products TEXT,
    microbiome_modulation TEXT,
    gut_brain_axis_effects TEXT,
    immune_modulation TEXT,
    anti_inflammatory_effects TEXT,
    fermentation_rate TEXT,
    scfa_production_profile TEXT,
    prebiotic_index REAL,
    clinical_evidence TEXT,
    effect_type TEXT,
    modulation_type TEXT,
    confidence_score REAL DEFAULT 0.75,
    tier_access INTEGER DEFAULT 1,
    requires_consent INTEGER,
    tk_consent_id TEXT,
    bc_benefit_sharing_pct REAL,
    indigenous_sourcing_verified INTEGER,
    community_id TEXT,
    created_at TEXT DEFAULT (datetime('now'))
);

-- Audit log for TK governance
CREATE TABLE IF NOT EXISTS omnipath_audit_log (
    audit_id INTEGER PRIMARY KEY AUTOINCREMENT,
    event_type TEXT NOT NULL,
    application TEXT NOT NULL,
    event_data TEXT,
    compound_ids TEXT,
    drug_names TEXT,
    result TEXT,
    processing_time_ms REAL,
    user_agent TEXT,
    ip_address TEXT,
    hipaa_deidentified INTEGER DEFAULT 1,
    created_at TEXT DEFAULT (datetime('now'))
);

-- Consent artifacts for TK
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
    created_at TEXT DEFAULT (datetime('now'))
);

-- Additional indexes for NeuroBotanica tables
CREATE INDEX IF NOT EXISTS idx_interactions_compound ON neurobotanica_drug_interactions(compound_id);
CREATE INDEX IF NOT EXISTS idx_interactions_drug ON neurobotanica_drug_interactions(drug_name);
CREATE INDEX IF NOT EXISTS idx_demographic_compound ON neurobotanica_demographic_factors(compound_id);
CREATE INDEX IF NOT EXISTS idx_demographic_group ON neurobotanica_demographic_factors(demographic_group);
CREATE INDEX IF NOT EXISTS idx_synergy_compounds ON neurobotanica_synergy_predictions(compound_a_id, compound_b_id);
CREATE INDEX IF NOT EXISTS idx_audit_event ON omnipath_audit_log(event_type);
CREATE INDEX IF NOT EXISTS idx_consent_status ON omnipath_consent_artifacts(consent_status);
CREATE INDEX IF NOT EXISTS idx_polysaccharides_compound ON neurobotanica_polysaccharides(compound_id);
