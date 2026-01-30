import os
import shutil
from pathlib import Path
home = Path.home()
logs_dir = home / '.wrangler' / 'logs'
if not logs_dir.exists():
    print('No logs dir:', logs_dir)
    raise SystemExit(1)
logs = sorted(logs_dir.iterdir(), key=lambda p: p.stat().st_mtime, reverse=True)
if not logs:
    print('No log files')
    raise SystemExit(1)
latest = logs[0]
dest = Path('.') / '.wrangler_log_latest.log'
shutil.copy2(latest, dest)
print('Copied', latest, 'to', dest)
