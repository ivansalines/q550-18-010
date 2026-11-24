#!/usr/bin/env python3
"""Synthetic exact benchmarks for signed relative work and curvature."""
import importlib.util,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('q550_52e_audit',ROOT/'scripts/run_audit.py')
m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m)
checks={}
t=np.linspace(0,.2,513)
N=len(t)
baseQ=np.column_stack((t,np.zeros((N,4))))
baseV=np.repeat(np.array([[1.,0.,0.,0.,0.]]),N,axis=0)
baseA=np.zeros((N,5))
lmb=2.3;eps=1e-4
amp=eps*np.exp(-lmb*t)
Dv=np.column_stack((np.zeros(N),amp,np.zeros((N,3))))
Dq=np.column_stack((np.zeros(N),eps*(1-np.exp(-lmb*t))/lmb,np.zeros((N,3))))
Da=-lmb*Dv
base={'tau':t,'xi':t,'Q':baseQ,'V':baseV,'acc':baseA}
pert={'tau':t,'xi':t,'Q':baseQ+Dq,'V':baseV+Dv,'acc':baseA+Da}
archive={}
for k,a in [('tau',t),('dq',m.mag(Dq)),('dv',m.mag(Dv)),('da',m.mag(Da)),('da_lat',m.mag(Da)),('da_par',np.zeros(N))]:
    archive['original_positive_1.0e-04__'+k]=a
r,v=m.full_geometry(base,pert,'original_positive_1.0e-04',eps,archive)
checks['signed_power_analytic']=np.max(abs(v['relative_power']+lmb*amp**2))<1e-16
checks['separation_rate_analytic']=np.max(abs(v['relative_position_rate']-eps*(1-np.exp(-lmb*t))/lmb*amp))<1e-16
checks['inward_acceleration_after_kick']=np.all(v['position_acceleration_dot'][1:]<0)
checks['positive_separation_rate_after_kick']=np.all(v['relative_position_rate'][1:]>0)
checks['negative_signed_power']=np.all(v['relative_power']<0)
checks['speed_contracts_exactly']=abs(r['velocity_difference_terminal_over_initial']-np.exp(-lmb*.2))<1e-10
checks['same_sign_position_divergence']=r['position_distance_decrease_steps']==0
checks['power_identity_spline']=r['max_power_vs_energy_spline_abs']<1e-12
checks['relative_vector_derivative']=r['max_power_vs_vector_spline_abs']<1e-12
checks['separation_identity_spline']=r['max_sep_rate_vs_distance_spline_abs']<1e-12
checks['time_integral_closure']=r['power_integral_closure_abs']<1e-11
checks['curvature_baseline_zero']=np.max(abs(v['curvature_reference']))<1e-14
checks['curvature_kicked_nonzero']=np.min(v['curvature_kick'])>0
checks['lateral_fraction_is_one']=np.max(abs(v['lateral_acceleration_fraction']-1))<1e-12
checks['no_false_position_return']=r['position_restoration_sample_count']==0
# Deliberate adversarial: reduce the position norm after a maximum while
# relative speed initially contracts. Negative power alone cannot certify return.
q=np.array([0.,1.,2.,1.,.5]);dv=np.array([2.,1.5,1.2,.9,.7]);
checks['adversarial_speed_not_position_equivalence']=bool(np.all(np.diff(dv)<0) and np.any(np.diff(q)<0) and np.any(np.diff(q)>0))
# An inward acceleration may coexist with increasing position separation:
checks['inward_not_sufficient']=bool(np.all(v['position_acceleration_dot'][1:]<0) and np.all(v['relative_position_rate'][1:]>0))
for name,ok in checks.items():
    if not ok:print('FAIL',name)
print('Q550 v18.52e synthetic tests',sum(checks.values()),'/',len(checks),'PASS' if all(checks.values()) else 'FAIL')
if not all(checks.values()):sys.exit(1)
