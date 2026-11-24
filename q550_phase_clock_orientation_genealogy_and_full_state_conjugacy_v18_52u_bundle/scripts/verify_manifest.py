import json,hashlib
from pathlib import Path
R=Path(__file__).resolve().parents[1]; d=json.loads((R/'SOURCE_MANIFEST.json').read_text())
assert d['root']==R.name and d['version']=='v18.52u'
for fn,dig in d['files'].items():
 p=R/fn
 assert p.exists(),f'MISSING SOURCE {p}'
 assert hashlib.sha256(p.read_bytes()).hexdigest()==dig,f'SOURCE DRIFT {fn}'
print('SOURCE_MANIFEST_PASS',len(d['files']))
