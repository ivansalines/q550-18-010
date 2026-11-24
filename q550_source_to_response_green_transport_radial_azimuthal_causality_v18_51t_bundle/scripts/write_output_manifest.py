#!/usr/bin/env python3
"""Write full result output SHA256 manifest after completed run."""
from pathlib import Path
import hashlib,json
r=Path(__file__).resolve().parents[1]
files=sorted(x for x in (r/'results').rglob('*') if x.is_file())
a={'root':r.name,'version':'v18.51t','file_sha256':{str(x.relative_to(r)):hashlib.sha256(x.read_bytes()).hexdigest() for x in files}}
(r/'SHA256SUMS_OUTPUT.json').write_text(json.dumps(a,indent=2)+'\n')
print('Q550 v18.51t output SHA256 generated',len(files))
