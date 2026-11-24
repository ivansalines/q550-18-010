#!/usr/bin/env python3
"""Validate generated JSON/NPZ/CSV outputs; scientific failure is not a test failure."""
import csv,json
from pathlib import Path
import numpy as np
R=Path(__file__).resolve().parents[1]/'results'
report=json.loads((R/'report.json').read_text())
assert report['version']=='v18.52a'
assert report['all_passed'] and all(report['numerical_checks'].values())
assert len(report['branches'])==6
assert len(report['scientific_robustness_tests'])==4
assert 'scientific_decision' in report
with (R/'physical_parameter_robustness_sections.csv').open() as f:rows=list(csv.DictReader(f))
assert len(rows)==2*3*487,len(rows)
assert set(r['seed'] for r in rows)=={'original','autonomous'}
assert set(r['branch'] for r in rows)=={'autonomous','midpoint','original'}
with np.load(R/'physical_branch_turning_arrays.npz') as a:
    assert a['xi'].shape==(513,)
    for seed in ('original','autonomous'):
        for branch in ('autonomous','midpoint','original'):
            assert a[f'{seed}_{branch}_plane2_channels'].shape==(513,2,4)
assert (R/'physical_kappa_branch_cross_channel_turning_robustness.png').is_file()
assert (R/'SCIENTIFIC_REPORT_SERVER.md').is_file()
print('Q550 v18.52a result tests PASS',len(rows),'qualified sectional rows')
