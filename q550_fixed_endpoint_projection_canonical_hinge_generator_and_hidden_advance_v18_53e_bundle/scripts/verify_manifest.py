import json,hashlib
from pathlib import Path
root=Path(__file__).resolve().parents[1]
m=json.loads((root/'SOURCE_MANIFEST.json').read_text())
for n,v in m['files'].items():
 p=root/n
 if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=v:raise SystemExit('SOURCE_MANIFEST_FAIL '+n)
print('SOURCE_MANIFEST_PASS',len(m['files']))
