#!/usr/bin/env python3
"""Pack an authenticated scientific sendback, without changing the parent or results."""
from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
import hashlib,json,os
root=Path(__file__).resolve().parents[1]
if not (root/'results/report.json').is_file():raise RuntimeError('NO_COMPLETED_SCIENTIFIC_REPORT')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
files={str(p.relative_to(root)):sha(p) for p in sorted((root/'results').rglob('*')) if p.is_file()}
if not files:raise RuntimeError('EMPTY_RESULTS')
(root/'OUTPUT_MANIFEST.json').write_text(json.dumps({'files':files},indent=2)+'\n')
out=root.parent/(root.name+'_server_sendback.zip')
tmp=out.with_suffix('.zip.tmp')
with ZipFile(tmp,'w',ZIP_DEFLATED,compresslevel=6,allowZip64=True) as z:
 for p in sorted(root.rglob('*')):
  if p.is_file():z.write(p,arcname=str(Path(root.name)/p.relative_to(root)))
with ZipFile(tmp) as z:
 bad=z.testzip()
 if bad:raise RuntimeError('SENDBACK_CRC_FAILED '+bad)
 prefix=root.name+'/'
 for manifest in ('SOURCE_MANIFEST.json','OUTPUT_MANIFEST.json'):
  f=json.loads(z.read(prefix+manifest))['files']
  for rel,h in f.items():
   if hashlib.sha256(z.read(prefix+rel)).hexdigest()!=h:raise RuntimeError('SENDBACK_MANIFEST_ERROR '+rel)
os.replace(tmp,out)
print('SENDBACK_COMPLETE',out,'sha256',sha(out),'bytes',out.stat().st_size)
