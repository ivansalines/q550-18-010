#!/usr/bin/env python3
import hashlib,json
from pathlib import Path
r=Path(__file__).resolve().parents[1]
files={p.relative_to(r).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((r/'results').glob('*')) if p.is_file()}
(r/'SHA256SUMS_OUTPUT.json').write_text(json.dumps({'root':r.name,'version':'v18.51r','file_sha256':files},indent=2)+'\n')
print('OUTPUT MANIFEST WRITTEN',len(files))
