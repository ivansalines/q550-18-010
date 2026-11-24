#!/usr/bin/env python3
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[1]
out={}
for f in sorted(ROOT.rglob('*')):
 if f.is_file() and f.parts[-2] not in ('results',) and f.name not in ('SOURCE_MANIFEST.json','OUTPUT_MANIFEST.json','PREFLIGHT.json') and '__pycache__' not in f.parts:
  with f.open('rb') as stream:out[str(f.relative_to(ROOT))]=hashlib.file_digest(stream,'sha256').hexdigest()
(ROOT/'SOURCE_MANIFEST.json').write_text(json.dumps({'version':'v18.52p','sha256':out},indent=2,sort_keys=True)+'\n')
print('SOURCE MANIFEST PASS',len(out))
