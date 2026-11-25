#!/usr/bin/env python3
"""Package source and outputs without changing the source genealogy."""
import argparse,hashlib,json,zipfile
from pathlib import Path

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def pack(root):
    root=Path(root).resolve()
    listing=json.loads((root/'SOURCE_MANIFEST.json').read_text())['files']
    for name,digest in listing.items():
        if not (root/name).is_file() or sha(root/name)!=digest:
            raise RuntimeError('BAD_SOURCE_MANIFEST '+name)
    rdir=root/'results'
    required=('report.json','SCIENTIFIC_REPORT_SERVER.md',
        'metric_path_factorial_45d.csv','metric_path_factorial_vectors.npz',
        'metric_path_family_summary.csv','factorial_metric_path_original.png',
        'factorial_metric_path_autonomous.png')
    for name in required:
        if not (rdir/name).is_file():raise RuntimeError('MISSING_RESULT '+name)
    report=json.loads((rdir/'report.json').read_text())
    if report.get('n_cases')!=204:raise RuntimeError('INCOMPLETE_204_POINT_GRID')
    output_files={str(p.relative_to(root)):sha(p) for p in sorted(rdir.rglob('*')) if p.is_file()}
    (root/'OUTPUT_MANIFEST.json').write_text(json.dumps({'files':output_files},indent=2)+'\n')
    out=root.parent/(root.name+'_server_sendback.zip')
    if out.exists():out.unlink()
    with zipfile.ZipFile(out,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as z:
        for p in sorted(root.rglob('*')):
            if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.pyc':
                z.write(p,str(Path(root.name)/p.relative_to(root)))
    with zipfile.ZipFile(out) as z:
        bad=z.testzip()
        if bad:raise RuntimeError('SEND_BACK_CRC_FAILED '+bad)
        for m in ('SOURCE_MANIFEST.json','OUTPUT_MANIFEST.json'):
            files=json.loads(z.read(root.name+'/'+m))['files']
            for p,h in files.items():
                if hashlib.sha256(z.read(root.name+'/'+p)).hexdigest()!=h:
                    raise RuntimeError('SENDBACK_HASH_FAILED '+p)
    print('SENDBACK_PACKAGED',out,'sha256',sha(out))
    return out
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',default=str(Path(__file__).resolve().parents[1]));a=p.parse_args()
    pack(a.root)
