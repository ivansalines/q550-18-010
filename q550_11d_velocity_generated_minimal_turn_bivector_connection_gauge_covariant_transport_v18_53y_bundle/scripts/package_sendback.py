#!/usr/bin/env python3
"""Validated ZIP sendback; safe Path joins and CRC/SHA256 checks."""
from pathlib import Path
import zipfile,hashlib,json,sys

def pack(root,out=None):
    root=Path(root).resolve(); out=Path(out) if out else root.parent/(root.name+'_server_sendback.zip')
    source=json.loads((root/'SOURCE_MANIFEST.json').read_text())['files']
    outputs={}
    for p in sorted((root/'results').rglob('*')):
        if p.is_file():outputs[str(p.relative_to(root))]=hashlib.sha256(p.read_bytes()).hexdigest()
    if not (root/'results/report.json').is_file() or not outputs:raise RuntimeError('MISSING_REAL_RESULTS')
    manifest={'files':outputs};(root/'OUTPUT_MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
    files=list(source)+['SOURCE_MANIFEST.json','OUTPUT_MANIFEST.json']+list(outputs)
    for rel,digest in {**source,**outputs}.items():
        if hashlib.sha256((root/rel).read_bytes()).hexdigest()!=digest:raise RuntimeError('PREPACKAGE_HASH_FAILURE '+rel)
    out.parent.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(out,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=5,allowZip64=True) as z:
        for rel in files:
            z.write(root/rel,arcname=str(Path(root.name)/rel))
    with zipfile.ZipFile(out) as z:
        if (bad:=z.testzip()) is not None:raise RuntimeError('SENDBACK_CRC_FAILED '+bad)
        for rel,digest in {**source,**outputs}.items():
            if hashlib.sha256(z.read(str(Path(root.name)/rel))).hexdigest()!=digest:
                raise RuntimeError('SENDBACK_HASH_FAILURE '+rel)
    print('SENDBACK_PACKAGE_PASS',len(source),'sources',len(outputs),'outputs',out)
    return out

if __name__=='__main__':
    root=Path(__file__).resolve().parents[1]
    pack(root)
