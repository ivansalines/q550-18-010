#!/usr/bin/env python3
"""Small exactly-solvable Green time-kernel benchmark (not Q550 physical data)."""
import importlib.util,sys
from pathlib import Path
import numpy as np
path=Path(__file__).resolve().parents[1]/'scripts'/'run_audit.py'
spec=importlib.util.spec_from_file_location('q550_v1851u_green_synthetic',path)
m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m)
xi=np.linspace(0,.2,33);source_xi=xi[::2];n=len(xi)
J=np.zeros((10,10));J[4,5]=1.3;J[6,5]=.4;J[4,7]=-.5;J[6,7]=1.0
b=np.zeros(10);b[5]=10;b[7]=2;b[1]=.04
assert np.linalg.norm(J@J)==0
Y=np.zeros((n,10));Q=np.zeros((n,5));Q[:,2]=1
L=m.observation_covectors(Y,Q)
assert np.allclose(L[:,0,4],1.) and np.allclose(L[:,1,6],1.)
P,Jx,_,_=m.fundamental(xi,source_xi,np.repeat(J[None],len(source_xi),axis=0))
K,H,B=m.kernel_from_phi(xi,source_xi,np.repeat(b[None],len(source_xi),axis=0),P,L)
Z=m.integrated_kernel(K,xi)
t=xi[-1]
expected=np.array([.5*t*t*13,.5*t*t*(-1),0.,.5*t*t*12])
expected_az=np.array([.5*t*t*4,.5*t*t*2,0.,.5*t*t*6])
assert np.max(np.abs(Z[-1,0]-expected))<1e-10,(Z[-1,0],expected)
assert np.max(np.abs(Z[-1,1]-expected_az))<1e-10,(Z[-1,1],expected_az)
assert np.max(np.abs(K[np.triu_indices(n,k=1)]))==0
assert np.max(np.abs(K[np.arange(n),np.arange(n)]))<1e-11
# Exact two-time kernel: K_r,v2(t,u)=13*(t-u); K_r,omega2(t,u)=-(t-u).
assert abs(K[-1,0,0,0]-13*t)<1e-11
assert abs(K[-1,0,0,1]+t)<1e-11
assert np.max(np.abs(Z[:,:,:3].sum(axis=2)-Z[:,:,3]))<1e-12
frac,onset,ages=m.signed_birth(xi,Z,K[-1])
assert abs(frac[-1]-(1-12/14))<1e-11
assert ages['radial_omega2']['mean_source_lag'] is not None
assert onset['0.01'] is not None
# Regression: check orientation (projection, channel) against independent
# impulse-restart tangent IVPs. This caught the original (2,3) vs (3,2) error.
err,nrestarts=m.restart_impulse_checks(xi,Jx,source_xi,
    np.repeat(b[None],len(source_xi),axis=0),L,K,
    restart_indices=(0,8,16,24))
assert nrestarts==4 and err<1e-10, (nrestarts,err)
print('Q550 v18.51u synthetic causal Green 11/11 PASS (restart axes fixed)')
