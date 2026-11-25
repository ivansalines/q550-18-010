#!/usr/bin/env python3
"""Validate and package an authentic v18.53x server sendback."""
import hashlib,json,zipfile
from pathlib import Path

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def pack(root):
    root=Path(root).resolve()
    sources=json.loads((root/'SOURCE_MANIFEST.json').read_text())['files']
    for rel,h in sources.items():
        if not (root/rel).is_file() or sha(root/rel)!=h:
            raise RuntimeError('SOURCE_MANIFEST_FAILURE '+rel)
    result_dir=root/'results'
    mandatory=['report.json','SCIENTIFIC_REPORT_SERVER.md',
         'metric_identification_gauge_scan.csv','fixed_physical_frame_path.csv',
         'moving_center_fixed_physical_frame.csv','moving_gauge_connection_checks.csv',
         'selected_fixed_frame_vectors.npz',
         'gauge_fixed_frame_comparison_original.png','gauge_fixed_frame_comparison_autonomous.png']
    for rel in mandatory:
        if not (result_dir/rel).is_file():raise RuntimeError('MISSING_RESULT '+rel)
    report=json.loads((result_dir/'report.json').read_text())
    if (report.get('n_parent_nodes'),report.get('n_gauge_rows'),
        report.get('n_fixed_physical_frame_rows'),
        report.get('n_reconstructed_moving_center_rows'),
        report.get('n_moving_gauge_rows'))!=(204,816,204,60,24):
        raise RuntimeError('INCOMPLETE_GAUGE_AUDIT')
    files={str(p.relative_to(root)):sha(p) for p in sorted(result_dir.rglob('*')) if p.is_file()}
    (root/'OUTPUT_MANIFEST.json').write_text(json.dumps({'files':files},indent=2)+'\n')
    destination=root.parent/(root.name+'_server_sendback.zip')
    if destination.exists():destination.unlink()
    with zipfile.ZipFile(destination,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as z:
        for p in sorted(root.rglob('*')):
            if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.pyc':
                z.write(p,str(Path(root.name)/p.relative_to(root)),compress_type=(zipfile.ZIP_STORED if p.suffix=='.zip' else zipfile.ZIP_DEFLATED))
    with zipfile.ZipFile(destination) as z:
        if (bad:=z.testzip()) is not None:raise RuntimeError('SENDBACK_ZIP_CRC_FAILED '+bad)
        for filename in ['SOURCE_MANIFEST.json','OUTPUT_MANIFEST.json']:
            manifest=json.loads(z.read(root.name+'/'+filename))['files']
            for rel,h in manifest.items():
                if hashlib.sha256(z.read(root.name+'/'+rel)).hexdigest()!=h:
                    raise RuntimeError('SENDBACK_HASH_FAILED '+rel)
    print('SENDBACK_PACKAGED',destination,'sha256',sha(destination))
    return destination
if __name__=='__main__':pack(Path(__file__).resolve().parents[1])
