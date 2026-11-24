#!/usr/bin/env python3
"""Verify real Q550 data endpoints, cross-seed stability and channel identities."""
from pathlib import Path
import csv,json
import numpy as np
root=Path(__file__).resolve().parents[1]
r=json.loads((root/'results/report.json').read_text())
assert r['version']=='v18.51x' and r['all_passed']
assert all(r['scientific_checks'].values())
rows=list(csv.DictReader((root/'results/channel_recombination_endpoint_ledger.csv').open()))
assert len(rows)==18
for seed in ('original','autonomous'):
 base=next(x for x in rows if x['seed']==seed and float(x['eta'])==0)
 mix=next(x for x in rows if x['seed']==seed and float(x['eta'])==.10)
 assert 88.1<float(base['end_angle_degrees'])<88.3
 assert 67.7<float(mix['end_angle_degrees'])<68.1
 assert 12<float(base['end_radial_cancellation_percent'])<12.1
 assert abs(float(mix['end_radial_cancellation_percent']))<1e-10
 assert abs(float(mix['end_total_norm'])-float(base['end_total_norm']))<1e-14
 assert float(mix['max_10d_source_additivity'])<1e-12
 assert float(mix['max_green_10d_additivity'])<1e-8
 assert float(mix['max_area_determinant_covariance_error'])<1e-12
 critical=r['seeds'][seed]['critical_symmetric_mix_for_zero_omega2_radial_at_end']
 assert .0601<critical<.0604
with np.load(root/'results/channel_recombination_robustness_arrays.npz') as a:
 assert len(a['xi'])==513
 assert np.all(np.isfinite(a['original_scan_angle']))
 assert np.all(np.isfinite(a['autonomous_scan_angle']))
print('Q550 v18.51x REAL RESULTS',len(r['scientific_checks']),'/',len(r['scientific_checks']),'PASS')
