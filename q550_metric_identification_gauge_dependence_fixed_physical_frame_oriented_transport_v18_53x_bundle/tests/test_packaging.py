#!/usr/bin/env python3
"""Regression for source validation, pathlib archive names and ZIP recovery."""
import tempfile,json,hashlib,importlib.util,zipfile
from pathlib import Path
base=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('pack_x',base/'scripts/package_sendback.py')
p=importlib.util.module_from_spec(spec);spec.loader.exec_module(p)
with tempfile.TemporaryDirectory() as tmp:
    root=Path(tmp)/'synthetic_bundle';(root/'results').mkdir(parents=True)
    (root/'s.txt').write_text('known script')
    h=hashlib.sha256((root/'s.txt').read_bytes()).hexdigest()
    (root/'SOURCE_MANIFEST.json').write_text(json.dumps({'files':{'s.txt':h}}))
    report=dict(n_parent_nodes=204,n_gauge_rows=816,n_fixed_physical_frame_rows=204,
                n_reconstructed_moving_center_rows=60,n_moving_gauge_rows=24)
    (root/'results/report.json').write_text(json.dumps(report))
    mandatory=['SCIENTIFIC_REPORT_SERVER.md','metric_identification_gauge_scan.csv',
              'fixed_physical_frame_path.csv','moving_center_fixed_physical_frame.csv',
              'moving_gauge_connection_checks.csv','selected_fixed_frame_vectors.npz',
              'gauge_fixed_frame_comparison_original.png','gauge_fixed_frame_comparison_autonomous.png']
    for name in mandatory:(root/'results'/name).write_bytes(b'synthetic-test')
    dest=p.pack(root)
    with zipfile.ZipFile(dest) as z:
        assert z.testzip() is None
        assert root.name+'/results/report.json' in z.namelist()
    print('SYNTHETIC_SENDBACK_PACKAGING_PASS')
