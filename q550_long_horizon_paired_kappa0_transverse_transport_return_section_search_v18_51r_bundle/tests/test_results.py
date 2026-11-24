#!/usr/bin/env python3
"""Independent saved-state and ledger checks, run after the nonlinear replay."""
import csv,json
from pathlib import Path
import numpy as np
P=Path(__file__).resolve().parents[1]/'results'
report=json.loads((P/'report.json').read_text())
assert report['version']=='v18.51r' and report['physical_ivps']==18
assert len(report['checks'])>=16 and all(report['checks'].values())
with (P/'extended_shadowing_phase_error_ledger.csv').open() as f:rows=list(csv.DictReader(f))
assert len(rows)==2*513
with np.load(P/'extended_paired_kappa_nonlinear_states.npz') as z:
 assert abs(float(z['xi'][-1])-.2)<1e-14
 for seed in ('original','autonomous'):
  a=z[f'{seed}_original_Q'];b=z[f'{seed}_autonomous_Q'];
  assert a.shape==(513,5) and b.shape==(513,5)
  assert np.linalg.norm(a[0]-b[0])<1e-12
  norms=np.linalg.norm(a-b,axis=1)
  vals=[float(r['physical_Q_gap']) for r in rows if r['seed_family']==seed]
  assert np.allclose(norms,vals,rtol=1e-7,atol=1e-11)
  assert np.all(np.isfinite(norms))
with (P/'solver_three_budget_ledger.csv').open() as f: runs=list(csv.DictReader(f))
assert len(runs)==6
with (P/'constant_lag_fit.csv').open() as f: fits=list(csv.DictReader(f))
assert len(fits)==2
print('Q550 v18.51r independent nonlinear output tests PASS')
