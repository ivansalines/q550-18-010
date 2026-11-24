#!/usr/bin/env python3
"""Check completeness; physical recovery is NEVER an automatic PASS requirement."""
from pathlib import Path
import json,csv,sys
ROOT=Path(__file__).resolve().parents[1]
res=ROOT/'results'
p=json.loads((res/'report.json').read_text())
assert p['version']=='v18.52d'
assert p['parent_v18_52c_midpoints_completed']==24
assert p['all_numerical_checks_passed'],p['checks']
with (res/'fixed_tau_lateral_acceleration_ledger.csv').open() as f:rows=list(csv.DictReader(f))
assert len(rows)==1026,len(rows)
assert all(float(r['delta_acceleration'])>=0 for r in rows)
assert not p['physical_recovery_certified']
assert not p['permanent_parameter_shift_is_recovery_test']
print('Q550 v18.52d output checks:',len(p['checks']),'numerical checks; 1026 physical-time paired rows PASS')
