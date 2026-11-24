#!/usr/bin/env python3
import hashlib,json
from pathlib import Path
r=Path(__file__).resolve().parents[1]
o=json.loads((r/'SOURCE_MANIFEST.json').read_text());assert o['root']==r.name
for fn,correct in o['sha256'].items():
 p=r/fn
 assert p.is_file(),str(p)
 with p.open('rb') as f:actual=hashlib.file_digest(f,'sha256').hexdigest()
 if actual!=correct:raise RuntimeError('SHA drift '+fn)
print('SOURCE MANIFEST VERIFIED',len(o['sha256']))
