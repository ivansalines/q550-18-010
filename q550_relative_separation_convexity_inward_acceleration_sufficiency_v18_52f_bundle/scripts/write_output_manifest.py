import hashlib,json
from pathlib import Path
root=Path(__file__).resolve().parents[1]
files={p.relative_to(root).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
       for p in sorted((root/'results').rglob('*')) if p.is_file() and p.name!='OUTPUT_MANIFEST.json'}
(root/'OUTPUT_MANIFEST.json').write_text(json.dumps({'version':'v18.52f','files':files},indent=2)+'\n')
print('Q550 v18.52f output hashes',len(files))
