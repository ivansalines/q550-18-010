#!/usr/bin/env python3
from pathlib import Path
import hashlib,json
r=Path(__file__).resolve().parents[1]
d=json.loads((r/'SOURCE_MANIFEST.json').read_text());files=d['files']
for f,expect in files.items():
 p=r/f
 if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=expect:raise SystemExit('SOURCE_SHA_FAIL '+f)
print('SOURCE_SHA256_PASS',len(files))
