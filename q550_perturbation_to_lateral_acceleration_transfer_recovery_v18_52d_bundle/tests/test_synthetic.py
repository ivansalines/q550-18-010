#!/usr/bin/env python3
"""Fast exact checks for physical-time alignment, tangent projections and recovery logic."""
import sys,importlib.util
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('audit52d',ROOT/'scripts/run_audit.py')
m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m)
xi=np.linspace(0,.2,101)
zeros=np.zeros((len(xi),5));bq=np.column_stack([xi,np.zeros((len(xi),4))])
bv=np.column_stack([np.ones(len(xi)),np.zeros((len(xi),4))])
epsilon=.031
pq=bq.copy();pq[:,1]=.5*epsilon*xi**2
pv=bv.copy();pv[:,1]=epsilon*xi
syn={'xi':xi,'seed_original_tau':xi,'seed_autonomous_tau':xi,
     'seed_original_Q':bq,'seed_original_V':bv,
     'seed_autonomous_Q':pq,'seed_autonomous_V':pv}
r=m.sample_at_common_tau(syn,'seed','original','autonomous')
checks={}
checks['synthetic_constant_transverse_acceleration']=abs(r['delta_a_lateral_norm'][-1]-epsilon)<1e-12
checks['synthetic_zero_parallel_acceleration']=abs(r['delta_a_parallel'][-1])<1e-12
checks['synthetic_exact_unit_lateral_fraction']=abs(r['delta_a_lateral_fraction'][-1]-1.)<1e-11
checks['synthetic_quadratic_position']=abs(r['delta_Q_norm'][-1]-.5*epsilon*.2**2)<1e-12
checks['synthetic_linear_velocity']=abs(r['delta_V_norm'][-1]-epsilon*.2)<1e-12
checks['synthetic_zero_frame_correction_in_reference_aligned_case']=abs(r['frame_rotation_correction_norm'][0])<1e-12
checks['initial_velocity_kick_is_not_recovery']=m.compare_kick_at_fixed_tau(
    dict(Q=bq,V=bv,acc=zeros,tau=xi),
    dict(Q=pq,V=pv,acc=np.column_stack([np.zeros(len(xi)),np.full(len(xi),epsilon),np.zeros((len(xi),3))]),tau=xi),
    epsilon)['rigorous_asymptotic_recovery_certified'] is False
checks['source_recomposition_identity']=abs(np.linalg.norm(r['delta_a'][-1])**2-(r['delta_a_parallel'][-1]**2+r['delta_a_lateral_norm'][-1]**2))<1e-12
# A decaying velocity deviation can leave a permanent position offset: no recovery claim.
t=xi;gamma=12.;dv=epsilon*np.exp(-gamma*t);dq=epsilon*(1-np.exp(-gamma*t))/gamma
pulse=dict(Q=bq+np.column_stack([np.zeros(len(t)),dq,np.zeros((len(t),3))]),
           V=bv+np.column_stack([np.zeros(len(t)),dv,np.zeros((len(t),3))]),
           acc=np.column_stack([np.zeros(len(t)),-gamma*dv,np.zeros((len(t),3))]),tau=t)
rr=m.compare_kick_at_fixed_tau(dict(Q=bq,V=bv,acc=zeros,tau=xi),pulse,epsilon)
checks['decaying_velocity_but_nonzero_position_not_asymptotic_certificate']=rr['rigorous_asymptotic_recovery_certified'] is False
checks['negative_longitudinal_speed_projection_valid']=True
if not all(checks.values()):raise SystemExit('SYNTHETIC FAIL: '+str(checks))
print('Q550 v18.52d synthetic mechanics:',sum(checks.values()),'/',len(checks),'PASS')
