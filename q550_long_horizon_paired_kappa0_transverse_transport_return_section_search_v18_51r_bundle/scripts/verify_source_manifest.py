#!/usr/bin/env python3
import hashlib,json
from pathlib import Path
r=Path(__file__).resolve().parents[1]
x=json.loads((r/'SOURCE_MANIFEST.json').read_text())
assert x['root']==r.name and x['version']=='v18.51r'
for f,h in x['files'].items():
 assert hashlib.sha256((r/f).read_bytes()).hexdigest()==h,f
print('SOURCE IDENTITY VERIFIED',len(x['files']))
