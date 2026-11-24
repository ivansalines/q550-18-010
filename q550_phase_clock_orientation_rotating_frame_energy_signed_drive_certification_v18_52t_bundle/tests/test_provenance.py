#!/usr/bin/env python3
from pathlib import Path
import hashlib,json,sys
r=Path(__file__).resolve().parents[1]
source=json.loads((r/'SOURCE_MANIFEST.json').read_text())
for n,h in source['files'].items():
 assert hashlib.sha256((r/n).read_bytes()).hexdigest()==h,n
assert source['root']==r.name
print('PROVENANCE_PASS',len(source['files']))
