#!/usr/bin/env python3
"""Test the ZIP packager on a synthetic mini source tree; no physics."""
import importlib.util,hashlib,json,tempfile,zipfile
from pathlib import Path
spec=importlib.util.spec_from_file_location('q550_w_pack',Path(__file__).resolve().parents[1]/'scripts/package_sendback.py')
p=importlib.util.module_from_spec(spec);spec.loader.exec_module(p)
with tempfile.TemporaryDirectory() as tmp:
    root=Path(tmp)/'synthetic_q550_w_test'
    (root/'scripts').mkdir(parents=True);(root/'results').mkdir()
    (root/'scripts/example.py').write_text('# synthetic data only\n')
    digest=hashlib.sha256((root/'scripts/example.py').read_bytes()).hexdigest()
    (root/'SOURCE_MANIFEST.json').write_text(json.dumps({'files':{'scripts/example.py':digest}}))
    (root/'results/report.json').write_text(json.dumps({'n_cases':204}))
    required=['SCIENTIFIC_REPORT_SERVER.md','metric_path_factorial_45d.csv','metric_path_factorial_vectors.npz',
              'metric_path_family_summary.csv','factorial_metric_path_original.png','factorial_metric_path_autonomous.png']
    for name in required:(root/'results'/name).write_bytes(b'synthetic unit test fixture')
    out=p.pack(root)
    with zipfile.ZipFile(out) as z:
        assert z.testzip() is None
        assert root.name+'/results/report.json' in z.namelist()
    print('SYNTHETIC_SENDBACK_PACKAGING_PASS')
