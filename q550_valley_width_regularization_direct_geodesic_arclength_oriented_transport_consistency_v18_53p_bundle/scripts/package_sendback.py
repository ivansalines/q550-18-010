#!/usr/bin/env python3
"""Safe, no-reintegration sendback packager (avoids prior Path/str bug)."""
import hashlib,json,zipfile
from pathlib import Path
root=Path(__file__).resolve().parents[1]
out=root.parent/(root.name+'_server_sendback.zip')

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for chunk in iter(lambda:f.read(1<<20),b''):h.update(chunk)
 return h.hexdigest()

files={str(p.relative_to(root)):sha(p) for p in sorted((root/'results').rglob('*')) if p.is_file() and p.name!='.keep'}
if not files or 'results/report.json' not in files:raise RuntimeError('INCOMPLETE_RESULTS')
(root/'OUTPUT_MANIFEST.json').write_text(json.dumps({'files':files},sort_keys=True,indent=2)+'\n')
allfiles=[p for p in sorted(root.rglob('*')) if p.is_file() and '.pycache' not in p.parts and '__pycache__' not in p.parts and p.suffix!='.pyc']
with zipfile.ZipFile(out,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=5,allowZip64=True) as z:
 for p in allfiles:
  z.write(p,arcname=str(Path(root.name)/p.relative_to(root)))
with zipfile.ZipFile(out) as z:
 if z.testzip() is not None:raise RuntimeError('SENDBACK_CRC_MISMATCH')
 pref=root.name+'/'
 for manifest in ('SOURCE_MANIFEST.json','OUTPUT_MANIFEST.json'):
  m=json.loads(z.read(pref+manifest))
  for rel,dig in m['files'].items():
   if hashlib.sha256(z.read(pref+rel)).hexdigest()!=dig:raise RuntimeError('SENDBACK_HASH_MISMATCH '+rel)
print('V1853P_SENDBACK_PASS',out,'OUTPUT',len(files),flush=True)
