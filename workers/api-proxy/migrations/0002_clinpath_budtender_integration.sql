-- Migration: ClinPath Integration for Budtender Application
-- Date: February 5, 2026
-- Dependencies: Requires neurobotanica_compounds table (should exist from earlier migration)

-- ============================================================================
-- CONSENT ARTIFACTS STUB (Prevents FK failure)
-- ============================================================================

CREATE TABLE IF NOT EXISTS omnipath_consent_artifacts (
    consent_id TEXT PRIMARY KEY,
    community_id TEXT,
    consent_type TEXT,
    status TEXT DEFAULT 'pending',
    created_at DATETIME DEFAULT CURRENT_TIMESTAMP
);

-- ============================================================================
-- FORMULATIONS TABLE (Core dependency)
-- ============================================================================

CREATE TABLE IF NOT EXISTS neurobotanica_formulations (
    formulation_id TEXT PRIMARY KEY,

    -- Customer/context
    customer_id TEXT,
    dispensary_id TEXT,

    -- Compound composition
    compound_ids TEXT NOT NULL,
    cannabinoid_profile TEXT NOT NULL,
    terpene_profile TEXT NOT NULL,
    polysaccharide_profile TEXT,

    -- Predicted outcomes
    predicted_efficacy REAL,
    predicted_safety REAL,
    synergy_score REAL,
    confidence_score REAL,

    -- Target indication
    primary_indication TEXT,
    secondary_indications TEXT,

    -- Traditional knowledge integration
    tk_flag BOOLEAN DEFAULT FALSE,
    consent_id TEXT,

    -- Metadata
    created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
    updated_at DATETIME DEFAULT CURRENT_TIMESTAMP,

    FOREIGN KEY (consent_id) REFERENCES omnipath_consent_artifacts(consent_id)
);

CREATE INDEX IF NOT EXISTS idx_formulation_customer ON neurobotanica_formulations(customer_id);
CREATE INDEX IF NOT EXISTS idx_formulation_dispensary ON neurobotanica_formulations(dispensary_id);
CREATE INDEX IF NOT EXISTS idx_formulation_indication ON neurobotanica_formulations(primary_indication);

-- ============================================================================
-- CLINPATH PREDICTIONS (Therapeutic Profiling)
-- ============================================================================

CREATE TABLE IF NOT EXISTS clinpath_predictions (
    prediction_id TEXT PRIMARY KEY,
    formulation_id TEXT NOT NULL,
    patient_profile_hash TEXT NOT NULL,

    -- Therapeutic predictions (0-1 scores)
    anxiolytic_score REAL NOT NULL CHECK(anxiolytic_score BETWEEN 0 AND 1),
    antidepressant_score REAL NOT NULL CHECK(antidepressant_score BETWEEN 0 AND 1),
    sedative_score REAL NOT NULL CHECK(sedative_score BETWEEN 0 AND 1),
    analgesic_score REAL NOT NULL CHECK(analgesic_score BETWEEN 0 AND 1),

    -- Cognitive effects (-1 to +1)
    memory_impact REAL CHECK(memory_impact BETWEEN -1 AND 1),
    focus_impact REAL CHECK(focus_impact BETWEEN -1 AND 1),
    neuroprotection_score REAL CHECK(neuroprotection_score BETWEEN 0 AND 1),

    -- Side effects (0-1 probability)
    psychoactivity_risk REAL NOT NULL CHECK(psychoactivity_risk BETWEEN 0 AND 1),
    dependence_risk REAL NOT NULL CHECK(dependence_risk BETWEEN 0 AND 1),
    sedation_risk REAL NOT NULL CHECK(sedation_risk BETWEEN 0 AND 1),
    anxiety_risk REAL NOT NULL CHECK(anxiety_risk BETWEEN 0 AND 1),

    -- Confidence metrics
    overall_confidence REAL NOT NULL CHECK(overall_confidence BETWEEN 0 AND 1),
    evidence_quality TEXT NOT NULL CHECK(evidence_quality IN ('high', 'medium', 'low')),

    -- Metadata
    created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
    model_version TEXT DEFAULT 'v1.0-mvp',

    FOREIGN KEY (formulation_id) REFERENCES neurobotanica_formulations(formulation_id) ON DELETE CASCADE
);

CREATE INDEX IF NOT EXISTS idx_clinpath_formulation ON clinpath_predictions(formulation_id);
CREATE INDEX IF NOT EXISTS idx_clinpath_profile ON clinpath_predictions(patient_profile_hash);
CREATE INDEX IF NOT EXISTS idx_clinpath_confidence ON clinpath_predictions(overall_confidence);

-- ============================================================================
-- DEMOGRAPHIC FACTORS (BioPath 96% Bias Correction Integration)
-- ============================================================================

CREATE TABLE IF NOT EXISTS clinpath_demographic_factors (
    factor_id TEXT PRIMARY KEY,
    prediction_id TEXT NOT NULL,

    -- Cannabinoid metabolism polymorphisms
    cyp2c9_variant TEXT CHECK(cyp2c9_variant IN ('*1/*1', '*1/*2', '*1/*3', '*2/*2', '*2/*3', '*3/*3', NULL)),
    faah_variant TEXT CHECK(faah_variant IN ('C/C', 'C/A', 'A/A', NULL)),
    cnr1_variant TEXT,

    -- Polysaccharide metabolism
    fucosidase_variant TEXT,
    dectin1_variant TEXT CHECK(dectin1_variant IN ('WT/WT', 'WT/Y238X', 'Y238X/Y238X', 'other', NULL)),

    -- Genetic ancestry
    genetic_ancestry TEXT CHECK(genetic_ancestry IN ('african', 'asian', 'european', 'indigenous', 'admixed', 'unknown')),

    -- Bias correction metadata
    bias_correction_applied BOOLEAN DEFAULT TRUE,
    demographic_adjustment_factor REAL DEFAULT 1.0,
    equalized_odds_score REAL,

    created_at DATETIME DEFAULT CURRENT_TIMESTAMP,

    FOREIGN KEY (prediction_id) REFERENCES clinpath_predictions(prediction_id) ON DELETE CASCADE
);

CREATE INDEX IF NOT EXISTS idx_demographic_prediction ON clinpath_demographic_factors(prediction_id);
CREATE INDEX IF NOT EXISTS idx_demographic_ancestry ON clinpath_demographic_factors(genetic_ancestry);

-- ============================================================================
-- BUDTENDER RECOMMENDATIONS (Dispensary UI Integration)
-- ============================================================================

CREATE TABLE IF NOT EXISTS budtender_recommendations (
    recommendation_id TEXT PRIMARY KEY,

    -- Patient input (de-identified)
    patient_age_range TEXT CHECK(patient_age_range IN ('18-25', '26-35', '36-45', '46-55', '56-65', '65+')),
    patient_sex TEXT CHECK(patient_sex IN ('M', 'F', 'other', 'prefer_not_to_say')),
    primary_symptom TEXT NOT NULL,
    secondary_symptoms TEXT,
    contraindications TEXT,

    -- Available cultivars
    available_cultivars TEXT NOT NULL,

    -- Top recommendation
    recommended_cultivar_id TEXT NOT NULL,
    recommended_cultivar_name TEXT NOT NULL,
    predicted_efficacy REAL NOT NULL CHECK(predicted_efficacy BETWEEN 0 AND 1),
    confidence_score REAL NOT NULL CHECK(confidence_score BETWEEN 0 AND 1),

    -- Ranking
    rank INTEGER DEFAULT 1,

    -- ClinPath integration
    formulation_id TEXT,
    clinpath_prediction_id TEXT,

    -- Tracking
    budtender_id TEXT,
    dispensary_id TEXT,
    accepted BOOLEAN,

    created_at DATETIME DEFAULT CURRENT_TIMESTAMP,

    FOREIGN KEY (formulation_id) REFERENCES neurobotanica_formulations(formulation_id),
    FOREIGN KEY (clinpath_prediction_id) REFERENCES clinpath_predictions(prediction_id)
);

CREATE INDEX IF NOT EXISTS idx_recommendation_dispensary ON budtender_recommendations(dispensary_id);
CREATE INDEX IF NOT EXISTS idx_recommendation_symptom ON budtender_recommendations(primary_symptom);
CREATE INDEX IF NOT EXISTS idx_recommendation_budtender ON budtender_recommendations(budtender_id);
CREATE INDEX IF NOT EXISTS idx_recommendation_timestamp ON budtender_recommendations(created_at);

-- ============================================================================
-- TRIGGERS
-- ============================================================================

CREATE TRIGGER IF NOT EXISTS update_formulation_timestamp
AFTER UPDATE ON neurobotanica_formulations
BEGIN
    UPDATE neurobotanica_formulations
    SET updated_at = CURRENT_TIMESTAMP
    WHERE formulation_id = NEW.formulation_id;
END;

-- ============================================================================
-- VIEWS
-- ============================================================================

CREATE VIEW IF NOT EXISTS budtender_recommendation_details AS
SELECT
    br.recommendation_id,
    br.patient_age_range,
    br.patient_sex,
    br.primary_symptom,
    br.secondary_symptoms,
    br.recommended_cultivar_name,
    br.predicted_efficacy,
    br.confidence_score,
    br.created_at,
    cp.anxiolytic_score,
    cp.antidepressant_score,
    cp.sedative_score,
    cp.analgesic_score,
    cp.psychoactivity_risk,
    cp.dependence_risk,
    cp.evidence_quality,
    nf.cannabinoid_profile,
    nf.terpene_profile,
    nf.synergy_score
FROM budtender_recommendations br
LEFT JOIN clinpath_predictions cp ON br.clinpath_prediction_id = cp.prediction_id
LEFT JOIN neurobotanica_formulations nf ON br.formulation_id = nf.formulation_id;

CREATE VIEW IF NOT EXISTS high_confidence_recommendations AS
SELECT
    recommendation_id,
    primary_symptom,
    recommended_cultivar_name,
    predicted_efficacy,
    confidence_score,
    created_at
FROM budtender_recommendations
WHERE confidence_score >= 0.80
ORDER BY created_at DESC;
