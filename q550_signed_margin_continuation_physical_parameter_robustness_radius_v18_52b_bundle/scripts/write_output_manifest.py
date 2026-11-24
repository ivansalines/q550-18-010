#!/usr/bin/env python3
import hashlib,json
from pathlib import Path
root=Path(__file__).resolve().parents[1]
files={}
for p in sorted((root/'results').glob('*')):
 if p.is_file():files[p.relative_to(root).as_posix()]=hashlib.sha256(p.read_bytes()).hexdigest()
(root/'SHA256SUMS_OUTPUT.json').write_text(json.dumps({'root':root.name,'version':'v18.52b','file_sha256':files},indent=2)+'\n')
print('Q550 v18.52b output hashes',len(files),flush=True)
