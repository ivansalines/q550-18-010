#!/usr/bin/env python3
"""Independent exact-law tests for source-channel splitting, including singular edge."""
from pathlib import Path
import importlib.util
import sys
import numpy as np
root=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('q550_51x',root/'scripts/run_audit.py')
m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m)
checks={}
a=np.array([3.2,1.4]);b=np.array([-.8,2.1]);c=np.array([.02,.05]);total=a+b+c
for eta in (-.5,-.2,0,.05,.10,.20,.499,.501,.6):
 M=m.transfer_matrix(eta)
 xy=np.column_stack((a,b))@M
 aa,bb=xy.T
 checks[f'sum preserved {eta}']=np.max(np.abs(aa+bb+c-total))<1e-13
 checks[f'area det law {eta}']=abs(m.signed_area(aa,bb)-np.linalg.det(M)*m.signed_area(a,b))<1e-13
 checks[f'det formula {eta}']=abs(np.linalg.det(M)-(1-2*eta))<1e-13
checks['zero-norm cancellation undefined']=np.isnan(m.cancellation_three(np.zeros((1,3)))[0])
checks['eta_half_singular']=np.linalg.matrix_rank(m.transfer_matrix(.5))==1
checks['positive_mix_orientation_unchanged']=m.signed_area(*(np.column_stack((a,b))@m.transfer_matrix(.2)).T)>0
checks['past_half_mix_orientation_reversed']=m.signed_area(*(np.column_stack((a,b))@m.transfer_matrix(.6)).T)<0
assert all(checks.values()),[k for k,v in checks.items() if not v]
print('Q550 v18.51x SYNTHETIC',len(checks),'/',len(checks),'PASS')
