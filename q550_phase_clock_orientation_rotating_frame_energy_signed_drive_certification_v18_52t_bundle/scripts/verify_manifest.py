import json,hashlib
from pathlib import Path
root=Path(__file__).resolve().parents[1]
data=json.loads((root/'SOURCE_MANIFEST.json').read_text())
for rel,expect in data['files'].items():
 f=root/rel
 if not f.is_file() or hashlib.sha256(f.read_bytes()).hexdigest()!=expect:
  raise SystemExit('SOURCE_HASH_FAIL '+rel)
print('SOURCE_MANIFEST_PASS',len(data['files']))
