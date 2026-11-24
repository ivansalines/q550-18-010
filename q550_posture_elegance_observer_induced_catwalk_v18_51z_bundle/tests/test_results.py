#!/usr/bin/env python3
"""Independent file-level consistency tests of produced v18.51z outputs."""
import json
from pathlib import Path
import numpy as np
R=Path(__file__).resolve().parents[1]
report=json.loads((R/'results/report.json').read_text())
assert report['version']=='v18.51z' and report['all_passed']
with np.load(R/'results/posture_observer_turning_arrays.npz') as data:
    checks={}
    xi=data['xi'];checks['513_sample_sections']=len(xi)==513
    for seed in ('original','autonomous'):
        p=seed+'_';lab=data[p+'lab_channel_vectors'];obs=data[p+'observer_channel_vectors']
        phi=data[p+'frame_phi_rad'];rr=data[p+'lab_rate_rad_per_xi'];cc=data[p+'corrected_lab_rate_rad_per_xi']
        mu=np.sum(lab[:,:,0]*lab[:,:,1],axis=1)
        e=np.sum(lab[:,:,3]**2,axis=1)
        checks[p+'cartesian_closure']=np.max(abs(lab[:,:,3]-lab[:,:,:3].sum(axis=2)))<2e-13
        checks[p+'observer_preserves_norms']=np.max(abs(np.linalg.norm(lab,axis=1)-np.linalg.norm(obs,axis=1)))<1e-13
        checks[p+'observer_connection']=np.nanmax(abs(rr[xi>.01]-cc[xi>.01]))<1e-9
        checks[p+'all_channels_turn_consistently']=bool(np.all(rr[xi>=.01,0]<0) and np.all(rr[xi>=.01,1]<0) and np.all(rr[xi>=.01,3]>0))
        checks[p+'total_mixed_turning_budget']=np.nanmax(abs(data[p+'total_turning'][xi>=.01]-data[p+'reconstructed_total_turning'][xi>=.01]))<2e-10
        checks[p+'positive_mutual_area']=bool(np.all(data[p+'oriented_area'][xi>=.01]>0))
        checks[p+'clock_monotone']=bool(np.all(np.diff(data[p+'tau'])>0))
        checks[p+'finite_angle_on_mask']=bool(np.all(np.isfinite(data[p+'relative_angle_rad'][xi>=.01])))
for n,v in checks.items():print('PASS' if v else 'FAIL',n)
assert all(checks.values())
print('RESULTS',len(checks),'/',len(checks),'PASS')
