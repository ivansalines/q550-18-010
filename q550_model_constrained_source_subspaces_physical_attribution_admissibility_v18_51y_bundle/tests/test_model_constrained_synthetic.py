#!/usr/bin/env python3
"""Synthetic model-constrained projection/chart tests, without Q550 IVPs."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import numpy as np
from run_audit import projector,split_b,rotate_chart,derivative_rotation,support_leaks_mixed,I,J,N

def error(a,b):return np.max(np.abs(np.asarray(a)-np.asarray(b)))
c={}
b=np.zeros((3,N));b[:,I]=[2,3,4];b[:,J]=[5,6,7];b[:,1]=[.1,.2,.3]
base=split_b(b)
c['source_recombination']=error(base.sum(-1),b)<1e-15
c['source_disjoint_support']=np.all(base[:,:,0][:,[a for a in range(N) if a!=I]]==0) and np.all(base[:,:,1][:,[a for a in range(N) if a!=J]]==0)
c['no_cross_channel_support']=error(base[:,:,0]@projector(J),0)==0 and error(base[:,:,1]@projector(I),0)==0
for eta in [0.,.06,.25]:
    mixed,l1,l2=support_leaks_mixed(b,eta)
    c[f'mix_source_sum_eta{eta}']=error(mixed.sum(-1),base[:,:,:2].sum(-1))<1e-14
    c[f'mix_support_zero_only_eta0_eta{eta}']=bool((max(l1.max(),l2.max())<1e-15)==(eta==0.))
    if eta:
        c[f'mix_support_analytic_eta{eta}']=abs(l1[0]-abs(eta*5)/np.hypot((1-eta)*2,eta*5))<1e-14
P1=projector(I);P2=projector(J)
J0=np.arange(N*N,dtype=float).reshape(N,N)/101
S=np.arange(N*3,dtype=float).reshape(N,3)/20
bmat=base[1]
for n,(theta,thetadot) in enumerate([(0,0),(.23,1.7),(-.11,-.8)]):
    T=rotate_chart(theta);Td=derivative_rotation(theta,thetadot);Ti=T.T
    Pt=T@P1@Ti
    bt=T@b[1];b1t=T@base[1,:,0]
    Jt=T@J0@Ti+Td@Ti
    St=T@S
    c[f'orthogonal_chart_{n}']=error(T@Ti,np.eye(N))<1e-15
    c[f'covariant_projector_{n}']=error(Pt@bt,b1t)<1e-14
    c[f'covariant_generator_{n}']=error(Jt@St+T@bmat,T@(J0@S+bmat)+Td@S)<1e-13
    c[f'covariant_recovery_{n}']=error(Ti@St,S)<1e-14
    c[f'frozen_projection_wrong_{n}']=bool((np.linalg.norm((np.eye(N)-P1)@b1t)>1e-8)==(theta!=0))
    c[f'no_connection_term_wrong_{n}']=bool((np.linalg.norm(Td@S)>1e-8)==(thetadot!=0))
assert all(c.values()),{k:v for k,v in c.items() if not v}
print('Q550 v18.51y synthetic',len(c),'/',len(c),'PASS')
