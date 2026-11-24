#!/usr/bin/env python3
"""Create an exact-name reproducible result package and validate CRC."""
import zipfile,json,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
out=ROOT.parent/(ROOT.name+'_server_sendback.zip')
# Parent zip 1.6MB is preserved, but the 84MB v18.51s lineage is referenced
# by byte-identical SHA-256 and excluded from result sendback.
with zipfile.ZipFile(out,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as z:
    for p in sorted(ROOT.rglob('*')):
        if not p.is_file() or '__pycache__' in p.parts:continue
        if p.name== 'q550_parameter_sensitivity_transport_transverse_posture_attribution_v18_51s_bundle_server_sendback.zip':continue
        z.write(p,ROOT.name+'/'+p.relative_to(ROOT).as_posix())
with zipfile.ZipFile(out) as z:
    if z.testzip() is not None:raise RuntimeError('CRC failure')
print('Q550 v18.52e sendback:',out)
