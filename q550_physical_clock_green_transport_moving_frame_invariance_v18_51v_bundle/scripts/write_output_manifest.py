#!/usr/bin/env python3
"""Hash actual server outputs, never mark unexecuted checks as passed."""
import hashlib,json
from pathlib import Path
root=Path(__file__).resolve().parents[1]
paths=sorted(p for p in (root/'results').iterdir() if p.is_file() and p.name not in ('LOCAL_PREFLIGHT.json',))
files={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
obj={'version':'v18.51v','root':root.name,'file_sha256':files}
(root/'SHA256SUMS_OUTPUT.json').write_text(json.dumps(obj,indent=2)+'\n')
print('Q550 v18.51v hashed',len(files),'server outputs')
