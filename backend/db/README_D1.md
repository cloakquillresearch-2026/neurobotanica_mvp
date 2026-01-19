Phase 1 — Cloudflare D1 Schema (Nevada Pilot)
=============================================

This folder contains the Cloudflare D1-compatible SQL schema for the Phase 1
pilot (Nevada). It is intentionally lightweight and SQLite-compatible so you
can test locally using the existing `neurobotanica_dev.db` development DB.

Files:
- `d1_schema.sql` — CREATE TABLE statements for `clinical_studies`, `conditions`, and `recommendations`.

Local apply (development):
1) Activate your virtualenv and run:

```powershell
python scripts/apply_d1_schema.py --sql backend/db/d1_schema.sql --db neurobotanica_dev.db
```

This will execute the SQL against the SQLite DB at `neurobotanica_dev.db`.

Cloudflare D1 deployment (outline):
1) Install and login with Wrangler: `npm install -g wrangler && wrangler login`
2) Create a migration from the SQL file:

   wrangler d1 migrations create init_schema --sql-file backend/db/d1_schema.sql

3) Deploy migrations to your D1 binding (replace binding name):

   wrangler d1 migrations deploy YOUR_D1_BINDING

Alternative (ad-hoc small scripts):

   wrangler d1 execute --binding YOUR_D1_BINDING --sql "$(cat backend/db/d1_schema.sql)"

Notes:
- Cloudflare D1 is SQLite-compatible but has operational differences (limits, performance). This schema avoids complex types and uses TEXT/REAL columns to remain compatible.
- Once D1 is provisioned, set the binding name in your Worker and update your code to use the D1 binding or the appropriate SDK.
