#!/usr/bin/env python3
"""Create output SHA-256 manifest for the executed audit."""
import hashlib
import json
from pathlib import Path
r=Path(__file__).resolve().parents[1]
results=r/'results'
files={str(p.relative_to(r)):hashlib.sha256(p.read_bytes()).hexdigest()
       for p in sorted(results.iterdir()) if p.is_file()}
(r/'SHA256SUMS_OUTPUT.json').write_text(json.dumps({'version':'v18.51w','file_sha256':files},indent=2)+'\n')
print('Q550 v18.51w output SHA-256 entries:',len(files))
