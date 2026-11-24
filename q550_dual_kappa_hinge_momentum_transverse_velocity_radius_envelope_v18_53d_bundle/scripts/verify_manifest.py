import json,hashlib
from pathlib import Path
r=Path(__file__).resolve().parents[1];m=json.loads((r/'SOURCE_MANIFEST.json').read_text())
for n,h in m['files'].items():
 p=r/n
 if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=h:raise RuntimeError('SOURCE_MANIFEST_FAIL '+n)
print('SOURCE_MANIFEST_PASS',len(m['files']))
