from pathlib import Path
import hashlib,json
r=Path(__file__).resolve().parents[1];m=json.loads((r/'SOURCE_MANIFEST.json').read_text())
for name,expected in m['sha256'].items():
 p=r/name
 if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=expected:raise RuntimeError('SOURCE_MANIFEST_MISMATCH '+name)
print('SOURCE_MANIFEST_PASS',len(m['sha256']))
