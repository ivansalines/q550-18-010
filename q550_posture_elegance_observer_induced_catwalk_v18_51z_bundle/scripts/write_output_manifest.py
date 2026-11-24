#!/usr/bin/env python3
"""Hash outputs after successful tests, independent of the source manifest."""
import hashlib,json
from pathlib import Path
root=Path(__file__).resolve().parents[1]
p={str(f.relative_to(root)):hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted((root/'results').iterdir()) if f.is_file()}
(root/'SHA256SUMS_OUTPUT.json').write_text(json.dumps({'version':'v18.51z','sha256':p},indent=2)+'\n')
print('Q550 v18.51z output hashes',len(p))
