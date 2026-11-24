#!/usr/bin/env python3
"""Check Q550 v18.51w completed on-server output and its scientific decision gates."""
import csv
import json
from pathlib import Path
import numpy as np
p=Path(__file__).resolve().parents[1]/'results'
r=json.loads((p/'report.json').read_text())
assert r['version']=='v18.51w' and r['all_passed']
assert all(r['scientific_checks'].values())
assert len(r['scientific_checks'])>=15
for seed in ('original','autonomous'):
 d=r['seeds'][seed]
 assert 88<d['angle_end_deg']<89
 assert 10<d['critical_observer_omega2_end_deg']<12
 assert d['angle_all_positive_orientation'] and d['dot_all_positive']
 assert d['total_end_norm']>.04 and 11<d['natural_radial_cancellation_end_percent']<13
 assert d['max_lab_vs_independent_t_cartesian']<1e-9
with (p/'channel_vector_section_ledger.csv').open() as f:
 rows=list(csv.DictReader(f))
assert len(rows)==1026
with np.load(p/'channel_vector_geometric_invariants.npz',allow_pickle=False) as z:
 assert z['xi'].shape==(513,)
 for seed in ('original','autonomous'):
  x=z[seed+'_mutual_angle_degrees'];a=z[seed+'_oriented_area']
  assert np.isnan(x[0]) and np.all(np.isfinite(x[1:]))
  assert np.all(np.diff(x[1:])<0) and np.all(np.diff(a)>0)
print('Q550 v18.51w result tests PASS')
