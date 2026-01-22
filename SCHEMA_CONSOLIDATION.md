# Database Schema Consolidation

## Summary

Consolidated `mvp_schema.sql` and `dynamic_schema.sql` into a single unified schema: `unified_schema.sql`

## Changes Made

### Base Schema
- Used `mvp_schema.sql` as the primary base (more production-ready)
- Retained all PRIMARY KEY constraints
- Retained all NOT NULL constraints
- Retained all DEFAULT values for timestamps
- Retained AUTOINCREMENT for audit_log
- Retained all performance indexes

### Additions from dynamic_schema.sql
- Added `neurobotanica_clinical_studies` table (was missing in mvp_schema.sql)
- Added indexes for clinical studies (condition, pubmed_id)

### Schema Differences Resolved

#### neurobotanica_drug_interactions
**mvp_schema.sql** (CHOSEN):
- Has proper PRIMARY KEY and NOT NULL constraints
- Cleaner column names: `severity`, `evidence_level`
- Better structured with `evidence_quality` and `evidence_source`

**dynamic_schema.sql** (NOT USED):
- Missing constraints
- Different column names: `severity_level`, `inhibition_potency`
- Less structured

#### neurobotanica_demographic_factors
**mvp_schema.sql** (CHOSEN):
- Compound-specific factors with `compound_id`
- Single `adjustment_factor` with confidence interval
- Cleaner model

**dynamic_schema.sql** (NOT USED):
- Category-based factors (age, gender, etc.)
- Separate `cyp450_adjustment` and `dosing_adjustment`
- More complex but less flexible

#### Tables Present in Both
All other tables follow the mvp_schema.sql structure with proper constraints and defaults.

## Migration Path

### For New Deployments
Use `unified_schema.sql` directly.

### For Existing Databases

#### From mvp_schema.sql
```sql
-- Add clinical_studies table
-- Run the neurobotanica_clinical_studies CREATE statement from unified_schema.sql
-- Add new indexes
CREATE INDEX IF NOT EXISTS idx_clinical_studies_condition ON neurobotanica_clinical_studies(condition);
CREATE INDEX IF NOT EXISTS idx_clinical_studies_pubmed ON neurobotanica_clinical_studies(pubmed_id);
```

#### From dynamic_schema.sql
```sql
-- WARNING: This requires data migration
-- Recommend: Deploy new database with unified_schema.sql and migrate data
-- Contact DBA for migration scripts
```

## Recommendations

1. **Use unified_schema.sql** for all new environments
2. **Archive old schemas**: Rename mvp_schema.sql → `archive/mvp_schema_old.sql` and dynamic_schema.sql → `archive/dynamic_schema_old.sql`
3. **Update documentation**: Update deployment docs to reference unified_schema.sql
4. **Update wrangler.toml**: Point database initialization to unified_schema.sql
5. **Test thoroughly**: Run full test suite to ensure compatibility

## Files Status

- ✅ `unified_schema.sql` - **USE THIS** (production-ready, complete)
- ⚠️ `mvp_schema.sql` - Archive (good structure, but missing clinical_studies)
- ⚠️ `dynamic_schema.sql` - Archive (incomplete, missing constraints)

## Notes

The unified schema prioritizes:
- Data integrity (PRIMARY KEY, NOT NULL constraints)
- Consistency (standardized column names)
- Performance (comprehensive indexes)
- Traditional Knowledge preservation (consent mechanisms)
- HIPAA compliance (audit logging with deidentification)
- Completeness (all necessary tables from both sources)
