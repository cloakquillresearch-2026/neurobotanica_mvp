#!/usr/bin/env python3
"""
Setup test database for NeuroBotanica MVP
Creates neurobotanica.db with hybrid schema and test data
"""
import sqlite3
import os


def create_test_database(db_path="neurobotanica.db"):
    """Create test database with hybrid schema supporting both mvp and dynamic schemas."""

    # Remove existing database
    if os.path.exists(db_path):
        os.remove(db_path)

    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Hybrid schema with columns from both schemas
    schema_sql = '''
    -- Drug interactions with BOTH schema columns
    CREATE TABLE neurobotanica_drug_interactions (
        interaction_id TEXT PRIMARY KEY,
        compound_id TEXT NOT NULL,
        drug_name TEXT NOT NULL,
        drug_class TEXT,
        cyp450_enzyme TEXT,
        interaction_type TEXT,
        severity TEXT,
        severity_level TEXT,
        clinical_effect TEXT,
        evidence_level TEXT,
        clinical_notes TEXT,
        requires_consent INTEGER,
        monitoring_parameters TEXT,
        dose_adjustment_required INTEGER,
        dose_adjustment_recommendation TEXT,
        adjustment_recommendation TEXT,
        evidence_quality TEXT,
        evidence_source TEXT,
        citations TEXT,
        created_at TEXT DEFAULT (datetime('now')),
        requires_tk_consent INTEGER,
        requires_consent_check INTEGER,
        tk_consent_id TEXT,
        consent_id TEXT
    );

    -- Demographic factors with BOTH schema columns
    CREATE TABLE neurobotanica_demographic_factors (
        factor_id TEXT PRIMARY KEY,
        compound_id TEXT,
        demographic_group TEXT,
        demographic_category TEXT,
        demographic_value TEXT,
        adjustment_factor REAL,
        cyp450_adjustment REAL,
        dosing_adjustment REAL,
        confidence_interval TEXT,
        evidence_basis TEXT,
        applicable_compounds TEXT,
        created_at TEXT DEFAULT (datetime('now'))
    );

    -- Compounds
    CREATE TABLE neurobotanica_compounds (
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

    -- Synergy predictions with consent_id
    CREATE TABLE neurobotanica_synergy_predictions (
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
        consent_id TEXT,
        created_at TEXT DEFAULT (datetime('now'))
    );

    -- Clinical studies
    CREATE TABLE neurobotanica_clinical_studies (
        study_id TEXT PRIMARY KEY,
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
        created_at TEXT DEFAULT (datetime('now'))
    );

    -- Formulations with both schemas
    CREATE TABLE neurobotanica_formulations (
        formulation_id TEXT PRIMARY KEY,
        formulation_name TEXT,
        formulation_type TEXT,
        primary_compound TEXT,
        ingredients TEXT,
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
        synergy_score REAL,
        confidence_level TEXT,
        created_at TEXT DEFAULT (datetime('now')),
        updated_at TEXT DEFAULT (datetime('now')),
        research_status TEXT
    );

    -- Polysaccharides
    CREATE TABLE neurobotanica_polysaccharides (
        polysaccharide_id TEXT PRIMARY KEY,
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
        requires_consent INTEGER,
        tk_consent_id TEXT,
        bc_benefit_sharing_pct REAL,
        indigenous_sourcing_verified INTEGER,
        community_id TEXT,
        created_at TEXT DEFAULT (datetime('now'))
    );

    -- Consent artifacts
    CREATE TABLE omnipath_consent_artifacts (
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

    -- Audit log
    CREATE TABLE omnipath_audit_log (
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

    -- Indexes
    CREATE INDEX idx_interactions_compound ON neurobotanica_drug_interactions(compound_id);
    CREATE INDEX idx_interactions_drug ON neurobotanica_drug_interactions(drug_name);
    CREATE INDEX idx_demographic_compound ON neurobotanica_demographic_factors(compound_id);
    CREATE INDEX idx_synergy_compounds ON neurobotanica_synergy_predictions(compound_a_id, compound_b_id);
    CREATE INDEX idx_audit_event ON omnipath_audit_log(event_type);
    CREATE INDEX idx_consent_status ON omnipath_consent_artifacts(consent_status);
    CREATE INDEX idx_clinical_studies_condition ON neurobotanica_clinical_studies(condition);

    -- View for SQLAlchemy compatibility
    CREATE VIEW clinical_studies AS SELECT * FROM neurobotanica_clinical_studies;
    '''

    for statement in schema_sql.split(';'):
        statement = statement.strip()
        if statement:
            cursor.execute(statement)

    # Add test data
    test_data = [
        # Compounds
        ('''INSERT INTO neurobotanica_compounds
            (compound_id, compound_name, compound_class, source, smiles, molecular_weight, requires_consent)
            VALUES
            ('cbd', 'Cannabidiol', 'cannabinoid', 'Cannabis', 'C21H30O2', 314.46, 0),
            ('thc', 'Tetrahydrocannabinol', 'cannabinoid', 'Cannabis', 'C21H30O2', 314.46, 0),
            ('curcumin', 'Curcumin', 'polyphenol', 'Turmeric', 'C21H20O6', 368.38, 0),
            ('tk_compound', 'Traditional Knowledge Compound', 'cannabinoid', 'Traditional', 'C20H30O', 290.0, 1)
        ''', ()),

        # Drug interactions
        ('''INSERT INTO neurobotanica_drug_interactions
            (interaction_id, compound_id, drug_name, drug_class, severity, severity_level,
             clinical_effect, evidence_level, citations, requires_consent_check, consent_id)
            VALUES ('cbd_warfarin', 'cbd', 'warfarin', 'anticoagulant', 'moderate', 'moderate',
                    'Increased bleeding risk', 'moderate', 'PMID:12345', 0, NULL)
        ''', ()),

        # Demographic factors
        ('''INSERT INTO neurobotanica_demographic_factors
            (factor_id, compound_id, demographic_category, demographic_value,
             adjustment_factor, cyp450_adjustment, dosing_adjustment)
            VALUES ('age_65_cbd', 'cbd', 'age', '65', 0.8, 0.9, 0.85)
        ''', ()),

        # Synergy predictions
        ('''INSERT INTO neurobotanica_synergy_predictions
            (synergy_id, compound_a_id, compound_b_id, synergy_score, confidence_score, consent_id)
            VALUES ('cbd_thc', 'cbd', 'thc', 0.75, 0.85, NULL)
        ''', ()),

        # Clinical studies
        ('''INSERT INTO neurobotanica_clinical_studies
            (study_id, condition, intervention, outcomes, confidence_score, sample_size)
            VALUES ('study_001', 'anxiety', 'cbd', 'positive', 0.85, 100)
        ''', ()),

        # Polysaccharides
        ('''INSERT INTO neurobotanica_polysaccharides
            (polysaccharide_id, polysaccharide_name, source_organism, requires_consent)
            VALUES ('beta_glucan_001', 'Beta-glucan', 'Fungi', 0)
        ''', ()),

        # Formulations
        ('''INSERT INTO neurobotanica_formulations
            (formulation_id, formulation_name, primary_compound, ingredients,
             target_condition, consent_ids, synergy_score, confidence_level)
            VALUES ('form_001', 'CBD Anxiety Relief', 'cbd', 'cbd,terpenes',
                    'anxiety', NULL, 0.75, 'moderate')
        ''', ()),
    ]

    for sql, params in test_data:
        cursor.execute(sql, params)

    conn.commit()
    conn.close()

    print(f"✓ Database created: {db_path}")
    print(f"✓ Tables created: 9")
    print(f"✓ Test data inserted: 7 records")


if __name__ == "__main__":
    import sys
    db_path = sys.argv[1] if len(sys.argv) > 1 else "neurobotanica.db"
    create_test_database(db_path)
