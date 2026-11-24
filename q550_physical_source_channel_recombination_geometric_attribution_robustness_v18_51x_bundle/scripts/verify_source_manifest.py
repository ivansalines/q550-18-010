from pathlib import Path
import hashlib,json
r=Path(__file__).resolve().parents[1]
m=json.loads((r/'SOURCE_MANIFEST.json').read_text())
assert m['root']==r.name and m['version']=='v18.51x'
for n,d in m['files'].items():
 p=r/n
 assert p.is_file() and hashlib.sha256(p.read_bytes()).hexdigest()==d,n
print('Q550 v18.51x source hashes verified:',len(m['files']))
