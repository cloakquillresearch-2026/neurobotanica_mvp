"""
Database Setup Script
Creates NeuroBotanica D1 schema tables.
"""

import sqlite3

DB_PATH = "neurobotanica.db"

def create_tables():
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()
    
    # OmniPath tables
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS omnipath_consent_artifacts (
        consent_id TEXT PRIMARY KEY,
        community_id TEXT NOT NULL,
        knowledge_type TEXT NOT NULL,
        consent_status TEXT NOT NULL,
        consent_scope TEXT NOT NULL,
        tk_labels TEXT,
        bc_labels TEXT,
        signature_rsa4096 TEXT,
        signature_ecdsa TEXT,
        manifest_hash TEXT,
        community_vote_pct REAL,
        revocation_authority TEXT,
        granted_at DATETIME NOT NULL,
        expires_at DATETIME,
        revoked_at DATETIME,
        last_verified_at DATETIME,
        created_at DATETIME DEFAULT CURRENT_TIMESTAMP
    );
    """)
    
    # Add other tables
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS neurobotanica_compounds (
        compound_id TEXT PRIMARY KEY,
        compound_name TEXT NOT NULL,
        compound_class TEXT NOT NULL,
        kingdom TEXT NOT NULL,
        molecular_formula TEXT,
        molecular_weight REAL,
        smiles TEXT,
        cas_number TEXT,
        primary_mechanisms TEXT,
        therapeutic_targets TEXT,
        bioavailability_oral REAL,
        half_life_hours REAL,
        ld50_mg_kg REAL,
        cyp450_interactions TEXT,
        contraindications TEXT,
        botanical_source TEXT,
        traditional_use TEXT,
        is_novel BOOLEAN DEFAULT 0,
        patent_protected BOOLEAN DEFAULT 0,
        consent_id TEXT,
        traditional_knowledge_flag BOOLEAN DEFAULT 0,
        community_attribution TEXT,
        requires_consent_check BOOLEAN DEFAULT 0,
        benefit_sharing_required BOOLEAN DEFAULT 0,
        available_to_computational_tier BOOLEAN DEFAULT 1,
        available_to_tk_tier BOOLEAN DEFAULT 1,
        data_source TEXT,
        created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
        updated_at DATETIME DEFAULT CURRENT_TIMESTAMP,
        FOREIGN KEY (consent_id) REFERENCES omnipath_consent_artifacts(consent_id)
    );
    """)
    
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS neurobotanica_drug_interactions (
        interaction_id TEXT PRIMARY KEY,
        compound_id TEXT NOT NULL,
        drug_name TEXT NOT NULL,
        drug_class TEXT,
        cyp450_enzyme TEXT,
        interaction_type TEXT,
        inhibition_potency TEXT,
        severity_level TEXT NOT NULL,
        clinical_effect TEXT,
        risk_description TEXT,
        monitoring_required BOOLEAN,
        monitoring_parameters TEXT,
        dose_adjustment_needed BOOLEAN,
        adjustment_recommendation TEXT,
        evidence_level TEXT,
        citations TEXT,
        requires_consent_check BOOLEAN DEFAULT 0,
        consent_id TEXT,
        created_at DATETIME DEFAULT CURRENT_TIMESTAMP
    );
    """)
    
    # Add omnipath_audit_log
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS omnipath_audit_log (
        log_id TEXT PRIMARY KEY,
        event_type TEXT NOT NULL,
        application TEXT NOT NULL,
        event_data TEXT,
        consent_id TEXT,
        policy_id TEXT,
        user_id TEXT,
        result TEXT,
        reason TEXT,
        processing_time_ms INTEGER,
        hipaa_deidentified BOOLEAN DEFAULT 0,
        timestamp DATETIME DEFAULT CURRENT_TIMESTAMP
    );
    """)
    
    # Add neurobotanica_synergy_predictions
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS neurobotanica_synergy_predictions (
        synergy_id TEXT PRIMARY KEY,
        compound_a_id TEXT NOT NULL,
        compound_b_id TEXT NOT NULL,
        synergy_score REAL NOT NULL,
        synergy_type TEXT,
        mechanism TEXT,
        is_cross_kingdom BOOLEAN,
        kingdom_a TEXT,
        kingdom_b TEXT,
        therapeutic_context TEXT,
        pathway_modulation TEXT,
        confidence_score REAL,
        evidence_basis TEXT,
        optimal_ratio TEXT,
        requires_consent BOOLEAN DEFAULT 0,
        consent_id TEXT,
        consent_verification_status TEXT,
        computation_date DATETIME DEFAULT CURRENT_TIMESTAMP
    );
    """)
    
    # Add neurobotanica_clinical_studies
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS neurobotanica_clinical_studies (
        study_id TEXT PRIMARY KEY,
        study_type TEXT NOT NULL,
        condition TEXT NOT NULL,
        intervention TEXT NOT NULL,
        outcomes TEXT,
        key_findings TEXT,
        citation TEXT,
        pubmed_id TEXT,
        doi TEXT,
        confidence_score REAL,
        sample_size INTEGER,
        study_duration_days INTEGER,
        publication_year INTEGER,
        consent_id TEXT,
        traditional_knowledge_source BOOLEAN DEFAULT 0,
        community_attribution TEXT,
        created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
        FOREIGN KEY (consent_id) REFERENCES omnipath_consent_artifacts(consent_id)
    );
    """)
    
    # Add neurobotanica_formulations
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS neurobotanica_formulations (
        formulation_id TEXT PRIMARY KEY,
        formulation_name TEXT,
        formulation_type TEXT,
        ingredients TEXT NOT NULL,
        target_condition TEXT,
        severity_level TEXT,
        synergy_score REAL,
        confidence_level REAL,
        predicted_mechanisms TEXT,
        total_dose_mg REAL,
        doses_per_day INTEGER,
        duration_days INTEGER,
        interaction_warnings TEXT,
        contraindications TEXT,
        clinical_evidence_level TEXT,
        validation_studies TEXT,
        requires_consent BOOLEAN DEFAULT 0,
        consent_ids TEXT,
        benefit_sharing_calculation TEXT,
        community_attributions TEXT,
        tk_enhanced BOOLEAN DEFAULT 0,
        non_tk_alternative_available BOOLEAN DEFAULT 0,
        pricing_tier_required TEXT,
        created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
        created_by TEXT,
        commercialization_status TEXT DEFAULT 'research'
    );
    """)
    
    # Add neurobotanica_polysaccharides
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS neurobotanica_polysaccharides (
        polysaccharide_id TEXT PRIMARY KEY,
        polysaccharide_name TEXT NOT NULL,
        kingdom TEXT NOT NULL,
        polysaccharide_type TEXT NOT NULL,
        glycosidic_linkage TEXT,
        branching_pattern TEXT,
        molecular_weight_kda REAL,
        degree_of_sulfation REAL,
        immune_modulation_score REAL,
        anti_inflammatory_potency REAL,
        gut_microbiome_effect TEXT,
        neuroprotective_mechanisms TEXT,
        fermentation_rate TEXT,
        scfa_profile TEXT,
        botanical_source TEXT,
        extraction_method TEXT,
        typical_dosage_mg TEXT,
        fucosidase_dependent BOOLEAN DEFAULT 0,
        dectin1_responsive BOOLEAN DEFAULT 0,
        consent_id TEXT,
        traditional_knowledge_flag BOOLEAN DEFAULT 0,
        community_attribution TEXT,
        indigenous_sourcing BOOLEAN DEFAULT 0,
        created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
        FOREIGN KEY (consent_id) REFERENCES omnipath_consent_artifacts(consent_id)
    );
    """)
    
    # Alter table to add missing columns if not exist
    try:
        cursor.execute("ALTER TABLE neurobotanica_drug_interactions ADD COLUMN requires_consent_check BOOLEAN DEFAULT 0")
    except sqlite3.OperationalError:
        pass  # Column already exists
    try:
        cursor.execute("ALTER TABLE neurobotanica_drug_interactions ADD COLUMN consent_id TEXT")
    except sqlite3.OperationalError:
        pass
    try:
        cursor.execute("ALTER TABLE neurobotanica_synergy_predictions ADD COLUMN consent_id TEXT")
    except sqlite3.OperationalError:
        pass
    
    # Insert test data for synergy predictions
    cursor.execute("""
    INSERT OR IGNORE INTO neurobotanica_synergy_predictions 
    (synergy_id, compound_a_id, compound_b_id, synergy_score, synergy_type, mechanism, confidence_score, requires_consent, consent_id, consent_verification_status)
    VALUES 
    ('cbd_thc_001', 'cbd', 'thc', 0.85, 'entourage', 'receptor_affinity', 0.9, 1, 'consent_001', 'verified'),
    ('cbd_curcumin_001', 'cbd', 'curcumin', 0.75, 'anti_inflammatory', 'pathway_modulation', 0.8, 0, NULL, NULL)
    """)
    
    # Insert test data for clinical studies
    cursor.execute("""
    INSERT OR IGNORE INTO neurobotanica_clinical_studies 
    (study_id, study_type, condition, intervention, outcomes, key_findings, citation, confidence_score, sample_size, publication_year)
    VALUES 
    ('cbd_thc_study_001', 'RCT', 'chronic_pain', 'cbd thc combination', 'pain_reduction', 'Synergy observed', 'J Pain 2023', 0.8, 100, 2023),
    ('cbd_study_001', 'observational', 'anxiety', 'cbd', 'anxiety_reduction', 'Effective monotherapy', 'J Anxiety 2022', 0.7, 50, 2022)
    """)
    
    # Insert test data for polysaccharides
    cursor.execute("""
    INSERT OR IGNORE INTO neurobotanica_polysaccharides 
    (polysaccharide_id, polysaccharide_name, kingdom, polysaccharide_type, fermentation_rate, scfa_profile, indigenous_sourcing, traditional_knowledge_flag, consent_id)
    VALUES 
    ('beta_glucan_001', 'Beta-Glucan', 'fungi', 'beta-glucan', 'fast', 'acetate,propionate,butyrate', 1, 1, 'consent_001'),
    ('inulin_001', 'Inulin', 'plant', 'fructan', 'moderate', 'acetate,butyrate', 0, 0, NULL)
    """)
    
    # Insert test consent data
    cursor.execute("""
    INSERT OR IGNORE INTO omnipath_consent_artifacts 
    (consent_id, community_id, knowledge_type, consent_status, consent_scope, tk_labels, bc_labels, signature_rsa4096, signature_ecdsa, manifest_hash, community_vote_pct, revocation_authority, granted_at, expires_at, revoked_at, last_verified_at)
    VALUES 
    ('consent_001', 'indigenous_community_1', 'traditional_medicine', 'active', 'research', 'cbd_thc_synergy', 'benefit_sharing', 'signature...', 'ecdsa...', 'hash...', 85.0, 'community_council', '2024-01-01', '2029-01-01', NULL, '2024-01-19')
    """)
    
    # Insert test data for demographic factors
    cursor.execute("""
    INSERT OR IGNORE INTO neurobotanica_demographic_factors 
    (factor_id, demographic_category, demographic_value, cyp450_adjustment, dosing_adjustment, applicable_compounds, evidence_basis)
    VALUES 
    ('age_65_female', 'age', '65', 0.8, 0.9, '%cbd%', 'Clinical studies show reduced metabolism'),
    ('gender_female', 'gender', 'female', 0.95, 0.95, '%cbd%', 'Hormonal differences in metabolism'),
    ('age_999_unknown', 'age', '999', 1.0, 1.0, '%', 'No adjustment for unknown demographics')
    """)
    
    # Insert test data for drug interactions
    cursor.execute("""
    INSERT OR IGNORE INTO neurobotanica_drug_interactions 
    (interaction_id, compound_id, drug_name, drug_class, cyp450_enzyme, interaction_type, inhibition_potency, severity_level, clinical_effect, risk_description, monitoring_required, monitoring_parameters, dose_adjustment_needed, adjustment_recommendation, evidence_level, citations, requires_consent_check, consent_id)
    VALUES 
    ('cbd_warfarin_001', 'cbd', 'warfarin', 'anticoagulant', 'CYP2C9', 'inhibition', 'strong', 'major', 'Increased INR', 'Risk of bleeding due to warfarin potentiation', 1, 'INR, PT', 1, 'Reduce warfarin dose by 20-30%', 'A', 'Clinical studies 2020-2023', 0, NULL),
    ('cbd_ssri_001', 'cbd', 'sertraline', 'SSRI', 'CYP2C19', 'inhibition', 'moderate', 'moderate', 'Increased SSRI levels', 'Potential serotonin syndrome', 1, 'SSRI levels, symptoms', 0, 'Monitor closely', 'B', 'Case reports 2021', 0, NULL)
    """)
    
    # Alter table to add missing columns if not exist
    try:
        cursor.execute("ALTER TABLE neurobotanica_drug_interactions ADD COLUMN requires_consent_check BOOLEAN DEFAULT 0")
    except sqlite3.OperationalError:
        pass  # Column already exists
    try:
        cursor.execute("ALTER TABLE neurobotanica_drug_interactions ADD COLUMN consent_id TEXT")
    except sqlite3.OperationalError:
        pass
    try:
        cursor.execute("ALTER TABLE neurobotanica_synergy_predictions ADD COLUMN consent_id TEXT")
    except sqlite3.OperationalError:
        pass

    conn.commit()
    conn.close()

if __name__ == "__main__":
    create_tables()