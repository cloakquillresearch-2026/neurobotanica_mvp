import json
import sqlite3

# Load the data from the dump
with open('neurobotanica_dump.json', 'r') as f:
    data = json.load(f)

# Connect to local SQLite to verify
conn = sqlite3.connect('neurobotanica.db')
cursor = conn.cursor()

print("Tables in database:")
for table_name, table_data in data.items():
    print(f"\n{table_name}:")
    print(f"  Columns: {len(table_data['columns'])}")
    print(f"  Rows: {len(table_data['rows'])}")

    # Show first row as sample
    if table_data['rows']:
        print(f"  Sample row: {table_data['rows'][0]}")

conn.close()