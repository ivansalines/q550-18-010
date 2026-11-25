#!/usr/bin/env python3
import hashlib,json,zipfile
from pathlib import Path
root=Path(__file__).resolve().parents[1];out=root.parent/(root.name+'_server_sendback.zip')
sh=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
results={str(p.relative_to(root)):sh(p) for p in sorted((root/'results').rglob('*')) if p.is_file() and p.name!='.keep'}
if not results or 'results/report.json' not in results:raise RuntimeError('EMPTY_RESULTS')
(root/'OUTPUT_MANIFEST.json').write_text(json.dumps({'files':results},indent=2,sort_keys=True)+'\n')
all_files=[p for p in sorted(root.rglob('*')) if p.is_file() and not ('__pycache__' in p.parts or p.suffix=='.pyc')]
with zipfile.ZipFile(out,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=5,allowZip64=True) as z:
 for p in all_files:z.write(p,arcname=str(Path(root.name)/p.relative_to(root)))
with zipfile.ZipFile(out) as z:
 if z.testzip() is not None:raise RuntimeError('BAD_OUTPUT_ZIP_CRC')
 for title in ('SOURCE_MANIFEST.json','OUTPUT_MANIFEST.json'):
  files=json.loads(z.read(root.name+'/'+title))['files']
  for rel,d in files.items():
   if hashlib.sha256(z.read(root.name+'/'+rel)).hexdigest()!=d:raise RuntimeError('BAD_OUTPUT_ZIP_SHA '+rel)
print('V1853Q_SENDBACK_PASS',out,len(results),'results')
