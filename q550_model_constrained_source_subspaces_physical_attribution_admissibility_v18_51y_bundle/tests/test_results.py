#!/usr/bin/env python3
"""Check archived real-Q550 v18.51y data and the decision boundary."""
import json,csv
from pathlib import Path
import numpy as np
r=Path(__file__).resolve().parents[1]
v=json.loads((r/'results/report.json').read_text())
checks={}
checks['report_success']=v['all_passed'] and v['version']=='v18.51y'
checks['scientific_checks_all']=all(v['scientific_checks'].values()) and len(v['scientific_checks'])>=15
rows=list(csv.DictReader((r/'results/source_support_admissibility_ledger.csv').open()))
checks['ten_mix_rows']=len(rows)==10
for s in ('original','autonomous'):
 rr=v['seeds'][s]; l=[e for e in rows if e['seed']==s]
 checks[s+'_fixed_chart_identity_only']=all((e['strict_anchored_support_admissible']=='True')==(float(e['eta'])==0) for e in l)
 checks[s+'_coupling_large_but_not_energy']=.85<rr['mean_jacobian_column_outside_v2']<1 and .80<rr['mean_jacobian_column_outside_omega2']<1
 checks[s+'_x_recombination_leaks_source']=.03<rr['critical_x_mixing_support_leak_v2_end']<.05 and .09<rr['critical_x_mixing_support_leak_omega2_end']<.12
 checks[s+'_chart_conjugacy']=rr['max_covariant_generator_defect']<1e-12 and rr['max_covariant_response_defect']<1e-12
 checks[s+'_frozen_index_control']=abs(rr['max_frozen_index_leak']-np.sin(np.deg2rad(15)))<1e-12
 checks[s+'_connection_control']=rr['max_connection_omission_relative_defect']>0.1
 checks[s+'_physical_angle_preserved']=88<rr['unmixed_response_angle_end_degrees']<89
 checks[s+'_no_fake_actuation']=rr['source_channel_independent_actuation_demonstrated'] is False
with np.load(r/'results/model_constrained_source_covariance_arrays.npz') as a:
 checks['grids_correct']=len(a['xi'])==513 and all(len(a[s+'_source_xi'])==257 for s in ('original','autonomous'))
 checks['no_nonfinite']=all(np.all(np.isfinite(a[s+'_frozen_index_leak'])) for s in ('original','autonomous'))
assert all(checks.values()),{k:v for k,v in checks.items() if not v}
print('Q550 v18.51y REAL results',len(checks),'/',len(checks),'PASS')
