import json
from pathlib import Path
p=Path('norml_complete_200plus_studies.json')
with p.open('r',encoding='utf-8') as f:
    j=json.load(f)
print('total_studies', j.get('total_studies'))
print('total_conditions', j.get('total_conditions'))
# compounds from descriptors
p2=Path('data/processed/neurobotanica_descriptors_validated.json')
with p2.open('r',encoding='utf-8') as f:
    desc=json.load(f)
compounds=set()
for d in desc:
    name=d.get('name') or d.get('compound')
    if name:
        compounds.add(name)
print('unique_compounds', len(compounds))
# count validated dimers (load as array may be large, so iterate)
p3=Path('data/processed/validated_dimers_nextgen.json')
count=0
with p3.open('r',encoding='utf-8') as f:
    for line in f:
        if line.strip().startswith('{'):
            count+=1
print('approx_dimers_entries', count)
