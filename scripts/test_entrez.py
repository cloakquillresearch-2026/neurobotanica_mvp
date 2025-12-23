#!/usr/bin/env python3
from Bio import Entrez
import os

Entrez.email = os.environ.get('ENTREZ_EMAIL') or 'none'
print('ENTREZ_EMAIL=', Entrez.email)
try:
    handle = Entrez.esearch(db='pubmed', term='THC[Title/Abstract] AND chronic pain[Title/Abstract] AND (systematic review[Publication Type] OR meta-analysis[Publication Type])', retmax=3)
    rec = Entrez.read(handle)
    handle.close()
    print('esearch returned:', rec.keys())
    print('IdList length:', len(rec.get('IdList', [])))
    print('Ids:', rec.get('IdList', []))
except Exception as e:
    print('Entrez error:', e)
