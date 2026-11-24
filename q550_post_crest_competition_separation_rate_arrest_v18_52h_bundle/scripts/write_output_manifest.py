#!/usr/bin/env python3
import hashlib,json
from pathlib import Path
root=Path(__file__).resolve().parents[1]
items={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest()
       for p in sorted((root/'results').rglob('*')) if p.is_file()}
(root/'OUTPUT_MANIFEST.json').write_text(json.dumps({'version':'v18.52h',
                             'file_sha256':items},indent=2)+'\n')
print('OUTPUT MANIFEST',len(items),'files')
