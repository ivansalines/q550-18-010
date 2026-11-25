#!/usr/bin/env python3
"""Regression: complete sendback ZIP with safe Path join and no missing results dir."""
import sys, json, hashlib, tempfile, zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from package_sendback import pack
with tempfile.TemporaryDirectory(prefix='v1853z_pack_') as d:
    p=Path(d)/'dummy_full_descriptive_bundle';(p/'results').mkdir(parents=True)
    (p/'src.py').write_text('print(1)\n');(p/'results/report.json').write_text('{"test":true}\n')
    h=hashlib.sha256((p/'src.py').read_bytes()).hexdigest()
    (p/'SOURCE_MANIFEST.json').write_text(json.dumps({'files':{'src.py':h}}))
    out=pack(p)
    with zipfile.ZipFile(out) as z:
        assert z.testzip() is None
        assert z.read(p.name+'/src.py')==b'print(1)\n'
        assert z.read(p.name+'/results/report.json').startswith(b'{')
print('SYNTHETIC_SENDBACK_PACKAGING_PASS')
