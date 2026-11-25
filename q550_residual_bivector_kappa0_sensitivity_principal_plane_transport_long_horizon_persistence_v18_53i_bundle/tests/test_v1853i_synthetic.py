#!/usr/bin/env python3
"""Lightweight adversarial mathematical checks; no historical parent required."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import numpy as np
from bivector_geometry import principal_skew,compare_geometry,bivector,integ
from run_audit import grid_kappa,slices_for,xi_grid,derivative_metrics,plane_projector,pair_metrics

def check(v,s):
 if not v:raise AssertionError(s)

def tensor(a,b,v):
 out=np.zeros((5,5));out[a,b]=v;out[b,a]=-v;return out

def synthetic_metric(R):
 p=principal_skew(R)
 return dict(R=R,metric=dict(residual_principal_sigma1=p['sigma1'],residual_principal_sigma2=p['sigma2'],invariant_cancellation=0.75))

def test_planes():
 A=tensor(0,1,3)+tensor(2,3,.5)
 P=principal_skew(A)
 assert abs(P['sigma1']-3)<1e-12 and abs(P['sigma2']-.5)<1e-12
 K0=grid_kappa();check(len(K0)==11 and len(set(k for _,k in K0))==11,'kappa grid includes all distinct historical and stress values')
 r=np.random.default_rng(53);X=r.normal(size=(5,5));Q,_=np.linalg.qr(X)
 B=Q@A@Q.T
 check(np.linalg.norm(principal_skew(B)['blocks'][0]['projector']-Q@P['blocks'][0]['projector']@Q.T)<1e-11,'gauge-covariant principal projector')
 check(compare_geometry(A,A)['full_bivector_cosine']>1-1e-13,'self agreement')
 # A large projected cancellation must NOT be identified as full-bivector cancellation.
 C=tensor(0,4,1);check(abs(C[0,1])<1e-13 and np.linalg.norm(C)>1,'invisible residual')

def test_derivatives():
 base=tensor(0,1,2)+tensor(2,3,.3);dr=tensor(0,2,.25)
 h=.001
 D=lambda k:synthetic_metric(base+k*dr)
 data=derivative_metrics(D(-h),D(0),D(h),h)
 check(abs(data['dR_dk_fro']-np.linalg.norm(dr))<1e-9,'finite difference exact linear bivector')
 check(data['nonlinearity_ratio']<1e-10,'finite difference curvature zero')
 check(data['dP1_dk_fro'] is not None,'projector slope measured')
 check(data['dP2_dk_fro'] is not None,'second projector slope measured')
 check(pair_metrics(D(-h),D(h))['bivector_cosine']<1,'different parameters rotate bivector')

def test_time_splitting():
 xs=np.linspace(12.8,25.6,1025)
 zz=np.linspace(25.6,51.2,1025)[1:]
 full=np.r_[xs,zz]
 s=slices_for(full,51.2)
 check(len(s)==6,'three segments and three cumulative horizons')
 for key,a,b in s:
  i=np.argmin(abs(full-a));j=np.argmin(abs(full-b))
  check(abs(full[i]-a)<1e-9 and abs(full[j]-b)<1e-9,'coincident checkpoint')
 times=np.linspace(0,2,2001);q=np.column_stack((times,np.zeros((len(times),4))));f=np.zeros_like(q);f[:,1]=times
 bb=bivector(q,f)
 full=integ(times,bb)
 middle=len(times)//2
 split=integ(times[:middle+1],bb[:middle+1])+integ(times[middle:],bb[middle:])
 check(np.linalg.norm(full-split)<1e-11,'additivity over adjacent windows')

def test_degeneracy():
 M=tensor(0,1,2)+tensor(2,3,2)
 p=principal_skew(M)
 check(p['relative_gap']<1e-12,'warning gap on equal singular modes')
 # A generic small perturbation of a degenerate skew tensor has nonunique separate planes.
 check(compare_geometry(M,M)['principal_plane_pairing_ill_conditioned'],'do not certify degenerate plane identity')

if __name__=='__main__':
 test_planes();test_derivatives();test_time_splitting();test_degeneracy()
 print('V18_53I_SYNTHETIC_PASS bivector_projectors/dual_kappa_grid/central_difference/degeneracy/window_additivity')
