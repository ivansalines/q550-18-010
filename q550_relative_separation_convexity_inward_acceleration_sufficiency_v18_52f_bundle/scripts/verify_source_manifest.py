import hashlib,json
from pathlib import Path
root=Path(__file__).resolve().parents[1]
manifest=json.loads((root/'SOURCE_MANIFEST.json').read_text())
for name,want in manifest['files'].items():
    p=root/name
    if not p.is_file():raise SystemExit('Missing source '+name)
    got=hashlib.sha256(p.read_bytes()).hexdigest()
    if got!=want:raise SystemExit('Source SHA-256 mismatch '+name)
print('Q550 v18.52f authenticated sources',len(manifest['files']))
