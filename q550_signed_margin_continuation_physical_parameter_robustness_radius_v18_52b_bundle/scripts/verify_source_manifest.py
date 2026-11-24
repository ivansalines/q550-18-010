#!/usr/bin/env python3
"""Verify source integrity, excluding separately generated results."""
import hashlib,json,sys
from pathlib import Path
root=Path(__file__).resolve().parents[1]
m=json.loads((root/'SOURCE_MANIFEST.json').read_text())
bad=[]
for name,want in m['files'].items():
 p=root/name
 if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=want:bad.append(name)
print('Q550 v18.52b sources',len(m['files'])-len(bad),'/',len(m['files']),'verified',flush=True)
if bad:
 print('Source drift:',bad,file=sys.stderr);sys.exit(2)
