import time
import requests
import os
import sys

print('Waiting 3s for server warmup...')
time.sleep(3)

for i in range(6):
    try:
        r = requests.get('http://127.0.0.1:8002/health', timeout=5)
        print(f'--- health call {i} --- status {r.status_code}')
        print(r.text[:2000])
    except Exception as e:
        print(f'health call {i} failed: {e}')
    time.sleep(6)


def tail(path, lines=200):
    try:
        with open(path, 'rb') as f:
            f.seek(0, os.SEEK_END)
            end = f.tell()
            size = 1024
            data = b''
            while True:
                if end - size < 0:
                    f.seek(0)
                    data = f.read()
                    break
                f.seek(end - size)
                data = f.read(size) + data
                if data.count(b'\n') > lines:
                    break
                size *= 2
            return b'\n'.join(data.splitlines()[-lines:]).decode('utf-8', errors='replace')
    except Exception as e:
        return f'Error reading {path}: {e}'

print('\n==== TAIL jwks_refresher_stdout.log ====')
print(tail('jwks_refresher_stdout.log', 400))
print('\n==== TAIL jwks_refresher_stderr.log ====')
print(tail('jwks_refresher_stderr.log', 400))

# Stop uvicorn
try:
    pid = int(open('uvicorn_jwks.pid').read().strip())
    print(f'Killing PID {pid}')
    if os.name == 'nt':
        import subprocess
        subprocess.run(['taskkill','/PID',str(pid),'/F'])
    else:
        import signal
        os.kill(pid, signal.SIGTERM)
    print('Stopped')
except Exception as e:
    print('Failed to stop process:', e)
    sys.exit(1)
