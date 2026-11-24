#!/usr/bin/env python3
"""Light end-to-end forced tangent smoke; exact kinematic response, no model replay."""
from pathlib import Path
import sys
import numpy as np
from scipy.integrate import solve_ivp
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from run_audit import integrate_sources,angular_budget,load_import,analyze_branch
ref=load_import(ROOT/'reference'/'parent_v18_51t_run_audit.py','q550_v1852a_test_parent_t')
xi=np.linspace(0.,.2,513)
ids=np.arange(0,len(xi),2)
J=np.zeros((len(ids),10,10))
J[:,4,5]=1.;J[:,6,7]=1.;J[:,0,1]=1.
b=np.zeros((len(ids),10));b[:,5]=3.;b[:,7]=2.;b[:,1]=.1
sol,nfev=integrate_sources(xi,ids,J,b,ref)
expected_v2=1.5*xi**2;expected_omega2=1.*xi**2
checks={}
checks['source_v2_kinematic_double_integral']=np.max(np.abs(sol[:,4,0]-expected_v2))<1e-12
checks['source_omega2_kinematic_double_integral']=np.max(np.abs(sol[:,6,1]-expected_omega2))<1e-12
checks['other_channel_is_retained']=np.max(np.abs(sol[:,0,2]-.05*xi**2))<1e-12
checks['total_recomposes']=np.max(np.abs(sol[:,:,:3].sum(axis=2)-sol[:,:,3]))<1e-12
Y=np.zeros((len(xi),10))
Q=np.tile(np.array([1.,0.,1.,0.,.5]),(len(xi),1))
summary, budget, plane, _ = analyze_branch('synthetic','midpoint',xi,Y,Q,(sol,nfev),b,ref)
checks['correct_cartesian_projection']=np.max(np.abs(plane[:,0,0]-expected_v2))<1e-12
checks['correct_azimuthal_projection']=np.max(np.abs(plane[:,1,1]-expected_omega2))<1e-12
checks['exact_orthogonality']=np.abs(summary['endpoint']['mutual_angle_deg']-90)<1e-8
checks['budget_closes']=budget['wedge_reconstruction_max_abs']<1e-10
checks['no_forced_false_robustness']=summary['n_sections_mixed_margin_nonpositive']>0
print('Smoke tests:',sum(checks.values()),'/',len(checks))
for k,v in checks.items():print('PASS' if v else 'FAIL',k)
if not all(checks.values()):raise SystemExit(1)
