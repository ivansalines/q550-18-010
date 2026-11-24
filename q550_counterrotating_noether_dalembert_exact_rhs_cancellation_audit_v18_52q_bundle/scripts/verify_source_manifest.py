#!/usr/bin/env python3
import hashlib,json
from pathlib import Path
r=Path(__file__).resolve().parents[1]
a=json.loads((r/'SOURCE_MANIFEST.json').read_text())
for name,h in a['sha256'].items():
 p=r/name
 if not p.is_file():raise SystemExit('MISSING: '+name)
 with open(p,'rb') as f:actual=hashlib.file_digest(f,'sha256').hexdigest()
 if h!=actual:raise SystemExit('SHA MISMATCH: '+name)
print('SOURCE MANIFEST PASS:',len(a['sha256']),'files',flush=True)
