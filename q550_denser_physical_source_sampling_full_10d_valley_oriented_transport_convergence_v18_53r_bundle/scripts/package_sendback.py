#!/usr/bin/env python3
"""Package a completed Q550 audit. No string/PosixPath concatenation."""
import zipfile,json,hashlib
from pathlib import Path
root=Path(__file__).resolve().parents[1]
out=root.parent/(root.name+'_server_sendback.zip')
sha=lambda b:hashlib.sha256(b).hexdigest()
rp=root/'results/report.json'
if not rp.is_file():raise RuntimeError('REPORT_MISSING')
report=json.loads(rp.read_text())
if report.get('scientific_decision','').startswith(('FAILED','NOT_')):
 raise RuntimeError('INCOMPLETE_PHYSICAL_AUDIT_REFUSE_SUCCESS_SENDBACK')
results={str(p.relative_to(root)):sha(p.read_bytes()) for p in sorted((root/'results').rglob('*')) if p.is_file() and p.name!='.keep'}
(root/'OUTPUT_MANIFEST.json').write_text(json.dumps({'files':results},indent=2,sort_keys=True)+'\n')
with zipfile.ZipFile(out,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=4,allowZip64=True) as z:
 for p in sorted(root.rglob('*')):
  if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.pyc':
   z.write(p,arcname=str(Path(root.name)/p.relative_to(root)))
with zipfile.ZipFile(out) as z:
 if z.testzip() is not None:raise RuntimeError('ZIP_CRC_FAIL')
 for title in ('SOURCE_MANIFEST.json','OUTPUT_MANIFEST.json'):
  for rel,h in json.loads(z.read(root.name+'/'+title))['files'].items():
   if sha(z.read(root.name+'/'+rel))!=h:raise RuntimeError('ZIP_SHA_FAIL '+rel)
print('V1853R_SENDBACK_PASS',out,len(results),flush=True)
