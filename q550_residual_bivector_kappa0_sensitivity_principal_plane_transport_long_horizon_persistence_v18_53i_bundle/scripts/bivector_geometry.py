#!/usr/bin/env python3
"""Frame-independent 5D bivector geometry and numerical diagnostics.

A skew 5x5 tensor has two nonzero singular values (each repeated twice)
and a null direction.  Real Schur blocks give orthogonal oriented principal
planes; the signs of the Schur column vectors are a gauge convention, while
projectors and unoriented principal angles are invariant.
"""
from __future__ import annotations
import numpy as np
from scipy.linalg import schur

def bivector(q,f):
 q=np.asarray(q,float);f=np.asarray(f,float)
 return q[..., :,None]*f[...,None,:]-f[..., :,None]*q[...,None,:]

def integ(t,B):return np.trapezoid(np.asarray(B,float),np.asarray(t,float),axis=0)

def cancel_metric(U,C):
 nu=float(np.linalg.norm(U));nc=float(np.linalg.norm(C));net=float(np.linalg.norm(U+C))
 return dict(radial_norm=nu,phase_norm=nc,residual_norm=net,
  invariant_cancellation=(1-net/(nu+nc) if nu+nc>1e-14 else None),
  invariant_cosine=(float(np.sum(U*C)/(nu*nc)) if nu*nc>1e-15 else None))

def principal_skew(B):
 B=np.asarray(B,float)
 if B.shape!=(5,5):raise ValueError('EXPECTED_SKEW_5x5')
 if np.linalg.norm(B+B.T)>2e-10*max(1,np.linalg.norm(B)):raise ValueError('NOT_SKEW')
 B=.5*(B-B.T)
 T,Z=schur(B,output='real')
 blocks=[];i=0
 scale=max(1.0,float(np.max(abs(T))))
 while i<5:
  if i<4 and abs(T[i+1,i])>1e-10*scale:
   amp=float(T[i,i+1]);u=Z[:,i];v=Z[:,i+1];sig=abs(amp)
   plane=np.column_stack([u,v]);projector=plane@plane.T
   blocks.append(dict(sigma=sig,sign=1 if amp>0 else -1,plane=plane,projector=projector,
      oriented_unit_bivector=np.sign(amp)*(np.outer(u,v)-np.outer(v,u))))
   i+=2
  else:i+=1
 blocks.sort(key=lambda z:z['sigma'],reverse=True)
 while len(blocks)<2:
  blocks.append(dict(sigma=0.,sign=0,plane=None,projector=None,oriented_unit_bivector=None))
 reconstruct=sum((z['sigma']*z['oriented_unit_bivector'] for z in blocks if z['plane'] is not None),np.zeros((5,5)))
 error=float(np.linalg.norm(B-reconstruct)/max(1.,np.linalg.norm(B)))
 if error>2e-10:raise RuntimeError('PRINCIPAL_SCHUR_RECONSTRUCTION_FAIL: '+str(error))
 kernel=np.linalg.svd(B)[2][-1]
 s1,s2=blocks[0]['sigma'],blocks[1]['sigma']
 return dict(sigma1=s1,sigma2=s2,relative_gap=(s1-s2)/max(s1,1e-15),blocks=blocks,
  kernel=kernel,reconstruction_error=error,skew_frobenius=float(np.linalg.norm(B)))

def principal_angle_degrees(A,B):
 if A is None or B is None:return None
 sv=np.linalg.svd(np.asarray(A).T@np.asarray(B),compute_uv=False)
 return [float(np.rad2deg(np.arccos(np.clip(s,0,1)))) for s in sorted(sv,reverse=True)]

def compare_geometry(A,B,warning_gap=.003):
 a=principal_skew(A);b=principal_skew(B)
 cos=float(np.sum(A*B)/max(np.linalg.norm(A)*np.linalg.norm(B),1e-20))
 # Projection angle is invariant but individual principal planes are not
 # robust near singular-value degeneracy.  Signal that explicitly.
 unstable=(a['relative_gap']<warning_gap or b['relative_gap']<warning_gap)
 planes=[]
 for i in (0,1):
  planes.append(dict(mode=i+1,sigma_a=a['blocks'][i]['sigma'],sigma_b=b['blocks'][i]['sigma'],
     principal_angles_deg=principal_angle_degrees(a['blocks'][i]['plane'],b['blocks'][i]['plane'])))
 return dict(full_bivector_cosine=cos,relative_residual_change=float(np.linalg.norm(A-B)/max(np.linalg.norm(A),np.linalg.norm(B),1e-15)),
  principal_plane_pairing_ill_conditioned=unstable,planes=planes,
  sigma_a=[a['sigma1'],a['sigma2']],sigma_b=[b['sigma1'],b['sigma2']])
