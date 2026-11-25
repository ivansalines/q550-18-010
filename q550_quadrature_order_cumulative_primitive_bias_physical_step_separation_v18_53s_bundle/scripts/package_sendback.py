#!/usr/bin/env python3
"""Package completed Q550 v18.53s results (avoids Path/string arithmetic bugs)."""
import hashlib,json,zipfile
from pathlib import Path
root=Path(__file__).resolve().parents[1]
out=root.parent/(root.name+'_server_sendback.zip')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
report_path=root/'results/report.json'
if not report_path.is_file():raise RuntimeError('MISSING_REPORT')
rep=json.loads(report_path.read_text())
if rep['scientific_decision'].startswith(('FAILED','NOT_')):raise RuntimeError('FAILED_AUDIT_CANNOT_PACKAGE')
if rep['new_physical_ivps'] not in (0,6):raise RuntimeError('INCOMPLETE_PHYSICAL_REFERENCE')
if rep['new_physical_ivps']==0 and rep['scientific_decision']!='AUTHENTIC_QUADRATURE_ONLY__NO_NEW_PHYSICAL_STEP_CERTIFICATE':
    raise RuntimeError('INVALID_POSTPROCESS_ONLY_CLASSIFICATION')
files={str(p.relative_to(root)):sha(p) for p in sorted((root/'results').rglob('*')) if p.is_file() and p.name!='.keep'}
(root/'OUTPUT_MANIFEST.json').write_text(json.dumps({'files':files},sort_keys=True,indent=2)+'\n')
with zipfile.ZipFile(out,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=4,allowZip64=True) as z:
    for p in sorted(root.rglob('*')):
        if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.pyc':
            z.write(p,arcname=str(Path(root.name)/p.relative_to(root)))
with zipfile.ZipFile(out) as z:
    if z.testzip():raise RuntimeError('SENDBACK_CRC_FAIL')
    for mf in ('SOURCE_MANIFEST.json','OUTPUT_MANIFEST.json'):
        entries=json.loads(z.read(root.name+'/'+mf))['files']
        for rel,digest in entries.items():
            if hashlib.sha256(z.read(root.name+'/'+rel)).hexdigest()!=digest:
                raise RuntimeError('SENDBACK_HASH_FAIL '+rel)
print('V1853S_SENDBACK_PASS',out,'OUTPUTS',len(files))
