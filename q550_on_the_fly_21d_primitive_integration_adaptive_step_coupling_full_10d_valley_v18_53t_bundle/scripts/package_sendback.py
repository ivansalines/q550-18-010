#!/usr/bin/env python3
"""Package only fully completed physical results; safe pathlib archive names."""
import argparse,hashlib,json,tempfile,zipfile
from pathlib import Path

def package(root,out=None):
    root=Path(root)
    report_path=root/'results/report.json'
    if not report_path.is_file():raise RuntimeError('NO_SCIENTIFIC_REPORT')
    report=json.loads(report_path.read_text())
    if report['scientific_decision'] not in (
      'AUTHENTIC_21D_ON_THE_FLY_PRIMITIVE__CONTROLLED_11D_STEP_EFFECT__NO_MATERIAL_HOLONOMY',
      'AUTHENTIC_21D_ON_THE_FLY_PRIMITIVE__SENSITIVITY_FLAGGED__NO_MATERIAL_HOLONOMY'):
        raise RuntimeError('NONPHYSICAL_OR_INCOMPLETE_CANNOT_PACKAGE')
    for key in ('completed_physical11_ivps','completed_coupled21_ivps','completed_locked10_ivps'):
        if report.get(key)!=6:raise RuntimeError('INCOMPLETE_IVPS '+key)
    required=['physical_21d_vs_11d.csv','valley_geometry_21d.csv',
              'augmented21_primitives_and_oriented_valleys.npz','SCIENTIFIC_REPORT_SERVER.md']
    if not all((root/'results'/name).is_file() for name in required):
        raise RuntimeError('REQUIRED_RESULT_MISSING')
    calc=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    files={str(p.relative_to(root)):calc(p) for p in sorted((root/'results').rglob('*')) if p.is_file() and p.name!='.keep'}
    (root/'OUTPUT_MANIFEST.json').write_text(json.dumps({'files':files},sort_keys=True,indent=2)+'\n')
    if out is None:out=root.parent/(root.name+'_server_sendback.zip')
    out=Path(out)
    with zipfile.ZipFile(out,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=4,allowZip64=True) as z:
        for p in sorted(root.rglob('*')):
            if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.pyc':
                z.write(p,arcname=str(Path(root.name)/p.relative_to(root)))
    with zipfile.ZipFile(out) as z:
        if z.testzip() is not None:raise RuntimeError('SENDBACK_CRC_FAILED')
        for manifest in ('SOURCE_MANIFEST.json','OUTPUT_MANIFEST.json'):
            entries=json.loads(z.read(root.name+'/'+manifest))['files']
            for rel,digest in entries.items():
                if hashlib.sha256(z.read(root.name+'/'+rel)).hexdigest()!=digest:
                    raise RuntimeError('SENDBACK_HASH_FAILED '+rel)
    print('V1853T_SENDBACK_PASS',out,'OUTPUTS',len(files),flush=True)
    return out

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--selftest',action='store_true')
    args=parser.parse_args()
    if args.selftest:
        with tempfile.TemporaryDirectory(prefix='q550_t_packaging_test_') as temp:
            root=Path(temp)/'q550_test_root';(root/'results').mkdir(parents=True)
            (root/'results'/'report.json').write_text(json.dumps(dict(
                 scientific_decision='AUTHENTIC_21D_ON_THE_FLY_PRIMITIVE__CONTROLLED_11D_STEP_EFFECT__NO_MATERIAL_HOLONOMY',
                 completed_physical11_ivps=6,completed_coupled21_ivps=6,completed_locked10_ivps=6)))
            for name in ('physical_21d_vs_11d.csv','valley_geometry_21d.csv',
                         'augmented21_primitives_and_oriented_valleys.npz','SCIENTIFIC_REPORT_SERVER.md'):
                (root/'results'/name).write_text('synthetic fixture\n')
            (root/'SOURCE_MANIFEST.json').write_text('{"files":{}}\n')
            package(root)
        print('V1853T_PACKAGING_SELFTEST_PASS')
    else: package(Path(__file__).resolve().parents[1])
