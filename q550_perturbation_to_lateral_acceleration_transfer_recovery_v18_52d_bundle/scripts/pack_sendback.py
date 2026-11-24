#!/usr/bin/env python3
"""Pack runnable scripts, scientific report, and numerical results; exclude 88 MB parent."""
from pathlib import Path
import zipfile
ROOT=Path(__file__).resolve().parents[1]
out=ROOT.parent/(ROOT.name+'_server_sendback.zip')
with zipfile.ZipFile(out,'w',zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as z:
    for p in sorted(ROOT.rglob('*')):
        if not p.is_file() or '__pycache__' in p.parts:continue
        rel=p.relative_to(ROOT)
        if str(rel).startswith('inputs/') and p.name.endswith('.zip') and p.stat().st_size>5_000_000:continue
        if str(rel).startswith('reference/') and p.name.endswith('.zip'):continue
        z.write(p,ROOT.name+'/'+rel.as_posix())
with zipfile.ZipFile(out) as z:
    bad=z.testzip()
if bad:raise RuntimeError('ZIP CRC failure: '+bad)
print('Q550 v18.52d server sendback:',out)
