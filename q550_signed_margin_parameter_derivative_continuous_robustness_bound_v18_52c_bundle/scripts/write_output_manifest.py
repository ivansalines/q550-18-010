#!/usr/bin/env python3
"""Hash all generated result files, avoiding the B historical manifest omission."""
from pathlib import Path
import json,hashlib
root=Path(__file__).resolve().parents[1]
files={}
for p in sorted((root/'results').rglob('*')):
    if p.is_file():
        files[p.relative_to(root).as_posix()]=hashlib.sha256(p.read_bytes()).hexdigest()
(root/'SHA256SUMS_OUTPUT.json').write_text(json.dumps({'root':root.name,'version':'v18.52c','file_sha256':files},indent=2)+'\n')
print('Q550 v18.52c output hashes:',len(files),flush=True)
