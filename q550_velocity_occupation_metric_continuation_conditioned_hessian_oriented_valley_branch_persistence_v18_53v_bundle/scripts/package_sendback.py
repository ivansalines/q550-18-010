#!/usr/bin/env python3
"""Safe, self-verified ZIP sendback; no raw Path+str concatenation."""
from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
import hashlib,json,os
root=Path(__file__).resolve().parents[1]
report=root/'results/report.json'
if not report.is_file():raise RuntimeError('NO_COMPLETED_V18_53V_RESULT')
if json.loads(report.read_text())['n_cases']!=204:raise RuntimeError('INCOMPLETE_204_CASES')
def sha(data):return hashlib.sha256(data).hexdigest()
results={str(p.relative_to(root)):sha(p.read_bytes()) for p in sorted((root/'results').rglob('*')) if p.is_file()}
if len(results)<7:raise RuntimeError('INCOMPLETE_V18_53V_OUTPUTS')
(root/'OUTPUT_MANIFEST.json').write_text(json.dumps({'files':results},indent=2)+'\n')
out=root.parent/(root.name+'_server_sendback.zip');tmp=out.with_suffix('.zip.tmp')
try:
    with ZipFile(tmp,'w',ZIP_DEFLATED,compresslevel=6,allowZip64=True) as z:
        for p in sorted(root.rglob('*')):
            if p.is_file() and '__pycache__' not in p.parts:
                z.write(p,arcname=str(Path(root.name)/p.relative_to(root)))
    with ZipFile(tmp) as z:
        if z.testzip() is not None:raise RuntimeError('SENDBACK_CRC_FAILED')
        prefix=root.name+'/'
        for manifest in ('SOURCE_MANIFEST.json','OUTPUT_MANIFEST.json'):
            for rel,digest in json.loads(z.read(prefix+manifest))['files'].items():
                if sha(z.read(prefix+rel))!=digest:raise RuntimeError('SENDBACK_MANIFEST_MISMATCH '+rel)
    os.replace(tmp,out)
except BaseException:
    if tmp.exists():tmp.unlink()
    raise
print('SENDBACK_COMPLETE',out,'sha256',sha(out.read_bytes()),'bytes',out.stat().st_size)
