#!/usr/bin/env python3
"""Package entire executable genealogy + outputs with self-checking output digest."""
import hashlib,json,zipfile
from pathlib import Path
base=Path(__file__).resolve().parents[1]
files=sorted(p for p in (base/'results').rglob('*') if p.is_file() and '__pycache__' not in p.parts)
digests={str(p.relative_to(base)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
(base/'OUTPUT_MANIFEST.json').write_text(json.dumps({'files':digests},indent=2)+'\n')
out=base.parent/(base.name+'_server_sendback.zip')
with zipfile.ZipFile(out,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=7) as z:
    for p in sorted(base.rglob('*')):
        if p.is_file() and '__pycache__' not in p.parts and not p.suffix=='.pyc':
            z.write(p,arcname=base.name+'/'+str(p.relative_to(base)))
with zipfile.ZipFile(out) as z:
    bad=z.testzip()
    if bad:raise RuntimeError('SENDBACK_CRC_FAIL: '+bad)
print('SENDBACK_CREATED',out,'bytes',out.stat().st_size,'results',len(digests))
