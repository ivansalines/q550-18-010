#!/usr/bin/env python3
"""Independent checks on computed v18.51s records and parent identity."""
import csv,json,hashlib,zipfile
from pathlib import Path
import numpy as np
R=Path(__file__).resolve().parents[1]
p=json.loads((R/'results/report.json').read_text())
assert p['version']=='v18.51s' and p['xi_end']==.2 and p['sections_per_seed']==513
assert p['all_passed'] and all(p['scientific_checks'].values())
assert p['parent_r_sha256']==hashlib.sha256(next((R/'inputs').glob('*v18_51r*.zip')).read_bytes()).hexdigest()
with (R/'results/cached_parameter_sensitivity_ledger.csv').open() as f:rows=list(csv.DictReader(f))
assert len(rows)==1026
for s in ('original','autonomous'):
 r=[x for x in rows if x['seed_family']==s]
 assert len(r)==513 and float(r[0]['xi'])==0 and abs(float(r[-1]['xi'])-.2)<1e-14
 assert float(r[-1]['plane2_squared_share'])>.9
 assert float(r[-1]['analytic_immersion_secant_absolute_error'])<1e-8
 assert float(r[-1]['secant_Q_norm'])>0
if p['mode']=='full':
 with (R/'results/variational_parameter_sensitivity_ledger.csv').open() as f:var=list(csv.DictReader(f))
 assert len(var)==1026
 with (R/'results/model_parameter_source_ledger.csv').open() as f:src=list(csv.DictReader(f))
 assert len(src)==2*(129+257)
 for s in ('original','autonomous'):
  rows=[r for r in var if r['seed_family']==s]
  assert float(rows[-1]['variational_Q_norm'])>0
  assert float(rows[-1]['variational_to_secant_Q_error'])<.02*max(1e-4,float(rows[-1]['central_secant_Q_norm']))
else:
 assert p['variational_summary'] is None
print('INDEPENDENT SAVED-RESULTS TESTS PASS; MODE:',p['mode'])
