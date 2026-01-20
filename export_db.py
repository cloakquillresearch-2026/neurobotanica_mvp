import sqlite3
import json

conn = sqlite3.connect('neurobotanica.db')
cursor = conn.cursor()

# Get all tables
cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
tables = cursor.fetchall()

data = {}
for table in tables:
    table_name = table[0]
    cursor.execute(f'SELECT * FROM {table_name}')
    columns = [desc[0] for desc in cursor.description]
    rows = cursor.fetchall()
    data[table_name] = {'columns': columns, 'rows': rows}

with open('neurobotanica_dump.json', 'w') as f:
    json.dump(data, f, indent=2)

conn.close()
print('Database exported to neurobotanica_dump.json')