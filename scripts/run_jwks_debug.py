import requests
import time

print('Waiting 3s for server warmup...')
time.sleep(3)

for i in range(8):
    try:
        r = requests.get('http://127.0.0.1:8004/internal/debug/jwks_cache', timeout=5)
        print(f'-- jwks snapshot {i} status {r.status_code} --')
        print(r.json())
    except Exception as e:
        print('jwks snapshot failed:', e)
    try:
        r2 = requests.get('http://127.0.0.1:8004/health', timeout=5)
        print(f'-- health {i} status {r2.status_code} --')
    except Exception as e:
        print('health failed:', e)
    time.sleep(10)

print('done')
