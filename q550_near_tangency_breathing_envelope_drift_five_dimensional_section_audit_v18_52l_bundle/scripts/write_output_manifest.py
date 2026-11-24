#!/usr/bin/env python3
from pathlib import Path
import hashlib,json
root=Path(__file__).resolve().parents[1]
files=sorted(p for p in (root/'results').rglob('*') if p.is_file() and p.name!='PREFLIGHT.json')
result={'version':'v18.52l','file_sha256':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}}
(root/'OUTPUT_MANIFEST.json').write_text(json.dumps(result,indent=2)+'\n')
print('OUTPUT MANIFEST WRITTEN',len(files),'files')
