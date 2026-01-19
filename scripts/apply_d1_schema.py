"""Apply D1-compatible schema locally and print Cloudflare `wrangler d1` commands.

This script helps with Phase 1 (Nevada Pilot) where Cloudflare D1 is the
primary database. It can apply the provided SQL to a local SQLite file for
development and will print the recommended `wrangler d1` commands to run
against a Cloudflare D1 instance.

Usage (local apply):
  python scripts/apply_d1_schema.py --sql backend/db/d1_schema.sql --db neurobotanica_dev.db

To apply to Cloudflare D1 (example):
  wrangler d1 migrations create init_schema --sql-file backend/db/d1_schema.sql
  wrangler d1 migrations deploy <YOUR_D1_BINDING_NAME>

Notes:
 - You must have `wrangler` installed and authenticated to run D1 commands.
 - This script does not call `wrangler` itself; it prints the commands so
   you can review and run them in your environment.
"""
import argparse
import sqlite3
import os
import sys


def apply_to_sqlite(sql_file: str, db_path: str):
    if not os.path.exists(sql_file):
        print(f"SQL file not found: {sql_file}")
        return 1

    with open(sql_file, "r", encoding="utf-8") as fh:
        sql = fh.read()

    # Ensure directory for DB exists
    db_dir = os.path.dirname(db_path)
    if db_dir and not os.path.exists(db_dir):
        os.makedirs(db_dir, exist_ok=True)

    conn = sqlite3.connect(db_path)
    try:
        cur = conn.cursor()

        # Backup existing tables if present by renaming them with a timestamp suffix.
        tables_to_check = [
            "clinical_studies",
            "conditions",
            "cannabinoid_profiles",
            "recommendations",
        ]
        import time
        ts = int(time.time())
        for t in tables_to_check:
            cur.execute("SELECT name FROM sqlite_master WHERE type='table' AND name=?;", (t,))
            if cur.fetchone():
                backup_name = f"{t}_backup_{ts}"
                print(f"Renaming existing table {t} -> {backup_name}")
                cur.execute(f"ALTER TABLE {t} RENAME TO {backup_name};")

        cur.executescript(sql)
        conn.commit()
        print(f"Applied schema to SQLite DB: {db_path}")
        return 0
    except Exception as e:
        print(f"Failed to apply schema to SQLite: {e}")
        return 2
    finally:
        conn.close()


def print_wrangler_instructions(sql_file: str):
    print("\nCloudflare D1 migration instructions:\n")
    print("1) Create a migration from the SQL file:")
    print(f"   wrangler d1 migrations create init_schema --sql-file {sql_file}")
    print("\n2) Deploy migrations to your D1 binding (replace YOUR_D1_BINDING):")
    print("   wrangler d1 migrations deploy YOUR_D1_BINDING")
    print("\nAlternative: execute SQL directly (small scripts):")
    print(f"   wrangler d1 execute --binding YOUR_D1_BINDING --sql \"$(cat {sql_file})\"")
    print("\nMake sure `wrangler` is installed and you are authenticated (wrangler login).")


def main():
    parser = argparse.ArgumentParser(description="Apply D1 schema locally and show Wrangler commands")
    parser.add_argument("--sql", default="backend/db/d1_schema.sql", help="Path to SQL schema file")
    parser.add_argument("--db", default="neurobotanica_dev.db", help="Path to local SQLite DB to apply schema to")
    parser.add_argument("--no-local", action="store_true", help="Don't apply locally, only print wrangler instructions")
    args = parser.parse_args()

    if not args.no_local:
        code = apply_to_sqlite(args.sql, args.db)
        if code != 0:
            sys.exit(code)

    print_wrangler_instructions(args.sql)


if __name__ == "__main__":
    main()
