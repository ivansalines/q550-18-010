#!/usr/bin/env python3
"""Regression test for output directory creation and Path-safe ZIP packaging."""
import sys,tempfile,hashlib,json,zipfile
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from package_sendback import pack
with tempfile.TemporaryDirectory(prefix='v1853y_pkg_') as tmp:
    root=Path(tmp)/'synthetic_bundle';(root/'results').mkdir(parents=True)
    (root/'input.txt').write_text('reproducible synthetic datum\n')
    (root/'SOURCE_MANIFEST.json').write_text(json.dumps({'files':{'input.txt':hashlib.sha256((root/'input.txt').read_bytes()).hexdigest()}}))
    (root/'results/report.json').write_text('{"synthetic":true}\n')
    out=pack(root)
    with zipfile.ZipFile(out) as z:
        assert z.testzip() is None and z.read(root.name+'/results/report.json')
print('SYNTHETIC_SENDBACK_PACKAGING_PASS')
