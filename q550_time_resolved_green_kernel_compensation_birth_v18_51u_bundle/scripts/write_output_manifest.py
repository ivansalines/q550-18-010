#!/usr/bin/env python3
"""Hash the complete generated outputs after the full scientific run."""
import hashlib,json
from pathlib import Path
r=Path(__file__).resolve().parents[1]
paths=sorted(p for p in (r/'results').rglob('*') if p.is_file())
d={'version':'v18.51u','root':r.name,'file_sha256':{str(p.relative_to(r)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
(r/'SHA256SUMS_OUTPUT.json').write_text(json.dumps(d,indent=2)+'\n')
print('Q550 v18.51u outputs hashed',len(paths))
