#!/usr/bin/env python3
import hashlib,json,sys
from pathlib import Path
r=Path(__file__).resolve().parents[1];d=json.loads((r/'SOURCE_MANIFEST.json').read_text())
fail=[]
for name,h in d['files'].items():
 p=r/name
 if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=h:fail.append(name)
print('SOURCE_MANIFEST_'+('PASS' if not fail else 'FAIL'),len(d['files']),fail)
sys.exit(bool(fail))
