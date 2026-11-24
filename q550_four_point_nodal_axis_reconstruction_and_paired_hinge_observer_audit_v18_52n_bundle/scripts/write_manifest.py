#!/usr/bin/env python3
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[1]
out={}
for p in sorted(ROOT.rglob('*')):
 if (p.is_file() and p.name != 'OUTPUT_MANIFEST.json'
     and p.suffix not in ('.pyc',) and '__pycache__' not in p.parts):
  out[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
(ROOT/'OUTPUT_MANIFEST.json').write_text(json.dumps({'version':'v18.52n','files_sha256':out},indent=2)+'\n')
print('OUTPUT MANIFEST PASS',len(out))
