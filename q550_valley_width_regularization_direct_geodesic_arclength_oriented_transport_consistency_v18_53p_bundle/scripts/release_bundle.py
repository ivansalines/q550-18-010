#!/usr/bin/env python3
"""Recreate the source manifest and distribution ZIP, preserving a single root."""
from pathlib import Path
import hashlib,json,zipfile
root=Path(__file__).resolve().parents[1]
archive=root.parent/(root.name+'.zip')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for buf in iter(lambda:f.read(1<<20),b''):h.update(buf)
 return h.hexdigest()
files={str(p.relative_to(root)):sha(p) for p in sorted(root.rglob('*')) if p.is_file()
   and not p.is_relative_to(root/'results') and p.name not in ('SOURCE_MANIFEST.json','OUTPUT_MANIFEST.json')
   and '__pycache__' not in p.parts and p.suffix!='.pyc'}
(root/'SOURCE_MANIFEST.json').write_text(json.dumps({'files':files},sort_keys=True,indent=2)+'\n')
with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=5) as z:
 for p in sorted(root.rglob('*')):
  if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.pyc' and p.name!='OUTPUT_MANIFEST.json':
   z.write(p,arcname=str(Path(root.name)/p.relative_to(root)))
with zipfile.ZipFile(archive) as z:
 assert z.testzip() is None
 prefix=root.name+'/'
 for rel,dig in files.items():
  assert hashlib.sha256(z.read(prefix+rel)).hexdigest()==dig,rel
print('V1853P_RELEASE_PASS',archive,'source_files',len(files))
