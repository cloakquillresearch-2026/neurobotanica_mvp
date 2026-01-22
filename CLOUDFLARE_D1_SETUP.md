# Cloudflare D1 Database Setup Guide

This guide explains how to deploy the NeuroBotanica database schema to Cloudflare D1.

## Prerequisites

1. Cloudflare account with Workers enabled
2. Wrangler CLI installed: `npm install -g wrangler`
3. Authenticated: `wrangler login`

---

## Option 1: Using Wrangler CLI (Recommended)

### Step 1: Create D1 Database

```bash
# Create the database (if not already created)
wrangler d1 create neurobotanica

# This will output:
# database_id = "a2818297-ae0e-4514-bc1c-5dc2bb8599fc"
# (already in your wrangler.toml)
```

### Step 2: Deploy Schema to D1

The database setup script outputs pure SQL that's compatible with both SQLite and D1:

```bash
# Generate SQL from Python script
python scripts/setup_test_database.py

# Or extract just the schema from the script
python -c "
from scripts.setup_test_database import create_test_database
import sqlite3
import sys

# Print SQL commands instead of executing
# (You can modify the script to output SQL)
"
```

### Step 3: Create Schema SQL File

Create `d1_schema.sql`:

```sql
-- Copy the CREATE TABLE statements from scripts/setup_test_database.py
-- Or use the unified_schema.sql file

-- Drug interactions
CREATE TABLE IF NOT EXISTS neurobotanica_drug_interactions (
    interaction_id TEXT PRIMARY KEY,
    compound_id TEXT NOT NULL,
    drug_name TEXT NOT NULL,
    -- ... (rest of columns from script)
);

-- (Copy all other CREATE TABLE statements)
-- (Copy all CREATE INDEX statements)
-- (Copy all CREATE VIEW statements)
```

### Step 4: Execute on D1

```bash
# Execute schema
wrangler d1 execute neurobotanica --file=d1_schema.sql

# Verify tables were created
wrangler d1 execute neurobotanica --command="SELECT name FROM sqlite_master WHERE type='table'"
```

---

## Option 2: Using SQL Directly

### Create Individual Tables

```bash
# Create compounds table
wrangler d1 execute neurobotanica --command="
CREATE TABLE IF NOT EXISTS neurobotanica_compounds (
    compound_id TEXT PRIMARY KEY,
    compound_name TEXT NOT NULL,
    compound_class TEXT,
    source TEXT,
    smiles TEXT,
    molecular_weight REAL,
    -- ... rest of columns
);"

# Repeat for all tables...
```

### Insert Test Data

```bash
# Insert compounds
wrangler d1 execute neurobotanica --command="
INSERT INTO neurobotanica_compounds
(compound_id, compound_name, compound_class, source, molecular_weight, requires_consent)
VALUES
('cbd', 'Cannabidiol', 'cannabinoid', 'Cannabis', 314.46, 0),
('thc', 'Tetrahydrocannabinol', 'cannabinoid', 'Cannabis', 314.46, 0);"

# Repeat for all test data...
```

---

## Option 3: Using Migration Scripts (Production)

For production, use proper migration management:

### Create Migration File

`migrations/001_initial_schema.sql`:
```sql
-- Migration 001: Initial schema
-- Created: 2026-01-22

-- Tables
CREATE TABLE IF NOT EXISTS neurobotanica_compounds (...);
CREATE TABLE IF NOT EXISTS neurobotanica_drug_interactions (...);
-- ... etc

-- Indexes
CREATE INDEX IF NOT EXISTS idx_interactions_compound ON neurobotanica_drug_interactions(compound_id);
-- ... etc

-- Views
CREATE VIEW IF NOT EXISTS clinical_studies AS SELECT * FROM neurobotanica_clinical_studies;
```

### Execute Migration

```bash
wrangler d1 migrations create neurobotanica initial_schema
wrangler d1 migrations apply neurobotanica
```

---

## Key Schema Features for D1

### 1. Hybrid Column Names (Compatibility)

The schema includes duplicate columns for compatibility:

```sql
CREATE TABLE neurobotanica_drug_interactions (
    severity TEXT,              -- For mvp_schema.sql
    severity_level TEXT,        -- For dynamic_schema.sql (both work!)
    -- ...
);
```

### 2. Missing Features in D1

D1 (SQLite) doesn't support:
- ❌ `AUTOINCREMENT` on non-INTEGER PRIMARY KEY (use TEXT instead)
- ❌ Full-text search (use external service)
- ❌ Complex triggers (keep logic in Workers)

But **does support**:
- ✅ All column types (TEXT, INTEGER, REAL, BLOB)
- ✅ Indexes for performance
- ✅ Views for compatibility
- ✅ Foreign keys (if enabled)

### 3. Required Tables

Minimum tables for core functionality:
1. `neurobotanica_compounds` - Cannabinoids and compounds
2. `neurobotanica_drug_interactions` - Drug-drug interactions
3. `neurobotanica_demographic_factors` - Bias correction data
4. `neurobotanica_synergy_predictions` - Synergy calculations
5. `neurobotanica_formulations` - Whole plant formulations
6. `omnipath_consent_artifacts` - TK consent tracking
7. `omnipath_audit_log` - HIPAA compliance

### 4. Critical Columns Added

Recent fixes include:
- `neurobotanica_formulations.tk_enhanced` - Required by whole_plant engine
- `neurobotanica_demographic_factors.applicable_compounds` - Required by bias correction queries

---

## Production Data Loading

### Load Clinical Studies (368+ studies)

```bash
# From your data files
python scripts/load_clinical_studies.py --target=d1

# Or manually
wrangler d1 execute neurobotanica --file=data/clinical_studies_insert.sql
```

### Load Cannabinoid Compounds (63 compounds)

```bash
python scripts/load_compounds.py --target=d1

# Or from JSON
python -c "
import json
with open('data/training/neurobotanica_complete_dataset_63compounds.json') as f:
    data = json.load(f)
    # Generate INSERT statements
"
```

---

## Verify Deployment

```bash
# Check table count
wrangler d1 execute neurobotanica --command="
SELECT COUNT(*) as table_count
FROM sqlite_master
WHERE type='table';"

# Check row counts
wrangler d1 execute neurobotanica --command="
SELECT
  (SELECT COUNT(*) FROM neurobotanica_compounds) as compounds,
  (SELECT COUNT(*) FROM neurobotanica_drug_interactions) as interactions,
  (SELECT COUNT(*) FROM neurobotanica_synergy_predictions) as synergies,
  (SELECT COUNT(*) FROM neurobotanica_clinical_studies) as studies;"

# Test query
wrangler d1 execute neurobotanica --command="
SELECT compound_id, compound_name
FROM neurobotanica_compounds
LIMIT 5;"
```

---

## Test with Worker

### Simple Test Worker

`test-d1-worker.ts`:
```typescript
export default {
  async fetch(request: Request, env: any): Promise<Response> {
    // Test query
    const stmt = env.NEUROBOTANICA_DB.prepare(
      'SELECT COUNT(*) as count FROM neurobotanica_compounds'
    );
    const result = await stmt.first();

    return new Response(JSON.stringify(result), {
      headers: { 'Content-Type': 'application/json' }
    });
  }
};
```

Deploy and test:
```bash
wrangler publish
curl https://your-worker.workers.dev
```

---

## Troubleshooting

### Issue: "no such table"

```bash
# List all tables
wrangler d1 execute neurobotanica --command="
SELECT name FROM sqlite_master WHERE type='table';"

# If tables missing, re-run schema creation
```

### Issue: "no such column"

```bash
# Check table schema
wrangler d1 execute neurobotanica --command="
PRAGMA table_info(neurobotanica_formulations);"

# If column missing, add it
wrangler d1 execute neurobotanica --command="
ALTER TABLE neurobotanica_formulations
ADD COLUMN tk_enhanced INTEGER DEFAULT 0;"
```

### Issue: Query too slow

```bash
# Check indexes
wrangler d1 execute neurobotanica --command="
SELECT name FROM sqlite_master WHERE type='index';"

# Add missing index
wrangler d1 execute neurobotanica --command="
CREATE INDEX idx_formulations_condition
ON neurobotanica_formulations(target_condition);"
```

---

## Environment Configuration

Update `wrangler.toml`:

```toml
name = "neurobotanica-api"
main = "src/api/neurobotanica.py"
compatibility_date = "2024-01-19"

[[d1_databases]]
binding = "NEUROBOTANICA_DB"
database_name = "neurobotanica"
database_id = "a2818297-ae0e-4514-bc1c-5dc2bb8599fc"  # Your actual ID

[[kv_namespaces]]
binding = "NEUROBOTANICA_CACHE"
id = "586ffb1e74594be1bfa58f34b80158dd"
```

---

## Next Steps

1. ✅ Create D1 database with `wrangler d1 create`
2. ✅ Deploy schema using one of the options above
3. ✅ Insert test data or production data
4. ✅ Verify with `wrangler d1 execute` queries
5. ✅ Deploy Worker that uses the database
6. ✅ Test API endpoints

---

## Notes

- **Schema is identical** to local SQLite (neurobotanica.db)
- **Migration path**: Same SQL works on both local and D1
- **Test locally first**: Use `neurobotanica.db` for development
- **Deploy to D1**: Use `wrangler d1 execute` for production
- **Data sync**: Use scripts to keep local and D1 in sync

The database setup script (`scripts/setup_test_database.py`) creates the **exact same schema** that works in both environments!

---

## Reference Links

- [Cloudflare D1 Docs](https://developers.cloudflare.com/d1/)
- [Wrangler D1 Commands](https://developers.cloudflare.com/workers/wrangler/commands/#d1)
- [D1 Migrations](https://developers.cloudflare.com/d1/learning/migrations/)
