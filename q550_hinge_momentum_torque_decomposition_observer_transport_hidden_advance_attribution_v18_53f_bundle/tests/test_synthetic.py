#!/usr/bin/env python3
"""Synthetic falsifiers — never classified as authenticated Q550 evidence."""
import math,sys,json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from torque_math import (axes,embed,force_components,material_rhs,acceleration_from_rhs,wedge,
 observer_transport,inplane_basis,tilted_plane,integral)
P=json.loads((ROOT/'inputs_config.json').read_text())['params'];K=json.loads((ROOT/'inputs_config.json').read_text())['fixed_kappa']
rng=np.random.default_rng(415053)
Y=np.column_stack([rng.uniform(-1.7,.5,80),rng.uniform(-.3,.3,80),rng.uniform(-2,2,80),rng.uniform(-1,1,80),
 rng.uniform(-2,.4,80),rng.uniform(-.5,.5,80),rng.uniform(-2,2,80),rng.uniform(-1,1,80),
 rng.uniform(-.6,.6,80),rng.uniform(-.5,.5,80)])
phase=rng.uniform(-2,2,80)
Q,V,J,g,D,force=force_components(Y,phase,K,P)
F=material_rhs(Y,phase,K,P)
A=acceleration_from_rhs(Y,F,P['m_x'])
assert np.max(abs(A-sum(force.values())))<1e-12
assert np.max(abs(F[:,0]-Y[:,1]))<1e-15
assert np.max(abs(F[:,2]-Y[:,3]))<1e-15
r1=np.exp(Y[:,0]);r2=np.exp(Y[:,4]);j1dot=r1*r1*(2*Y[:,1]*Y[:,3]+F[:,3]);j2dot=r2*r2*(2*Y[:,5]*Y[:,7]+F[:,7])
assert np.max(abs(j1dot-j2dot))<1e-13
# Independent historical support RHS agrees component-by-component with
# the material-time model on 64 off-trajectory states (clock source is tested
# separately on the authentic parent; support RHS uses positive phase clock).
import importlib.util
from types import SimpleNamespace
spec=importlib.util.spec_from_file_location('q550_old_support',ROOT/'reference/test_support_q550_oriented_coupled_bulge_cone_audit_v1.py')
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
audit=SimpleNamespace(sigmoid=lambda z:float(1/(1+np.exp(-z))),
    dU=lambda r,e:float(2*r-2.4*r*r+.8*r*r*r-e*e/r**3),
    dV_barrier=lambda x,B:float(B*(2*x-1)/(x*(1-x))))
par=SimpleNamespace(**P)
err=0.0
for j in range(64):
 old_F=old.support_rhs_float(np.r_[Y[j],phase[j]],kappa0=K,params=par,audit=audit)
 err=max(err,float(np.max(np.abs(old_F[:10] *g['x'][j]/(2*np.pi)-F[j]))))
assert err<1e-12,('HISTORICAL_OFF_ORBIT_RHS_MISMATCH',err)
# A basis rotation within the same plane gives invariant oriented moment and zero frame torque.
e,n,h=axes(-.308,[.26,-.41,.874])
for ang in [-.62,0,.298,.84]:
 a,b,ad,bd=inplane_basis(n,h,ang,.375)
 assert abs(float(wedge(Q,V,a,b)[0]-wedge(Q,V,n,h)[0]))<1e-13
 assert abs(float(observer_transport(Q,V,a,b,ad,bd)[0]))<1e-13
# A truly tilted plane in R5 has an observer-specific torque despite a fixed physical state.
y0=Y[:1];q,v=Q[:1],V[:1]
b=.17;bd=.29;nt,ht,nd,hd=tilted_plane(e,n,h,b,bd)
exact=observer_transport(q,v,nt,ht,nd,hd)[0]
step=1e-6
np_,hp_,_,_=tilted_plane(e,n,h,b+bd*step,bd)
nm_,hm_,_,_=tilted_plane(e,n,h,b-bd*step,bd)
fd=(wedge(q,v,np_,hp_)-wedge(q,v,nm_,hm_))/(2*step)
assert abs(float(exact-fd[0]))<3e-10
# A rotating observer can create apparent moment variations even at A=0.
t=np.linspace(0,1,2501);q=np.tile([1.,0.,0.,0.,1.],(len(t),1));v=np.tile([0.,1.,0.,0.,0.],(len(t),1))
beta=.4*np.sin(2*np.pi*t);betad=.4*2*np.pi*np.cos(2*np.pi*t)
nm=np.cos(beta)[:,None]*n+np.sin(beta)[:,None]*e;hm=np.broadcast_to(h,nm.shape)
nd=betad[:,None]*(-np.sin(beta)[:,None]*n+np.cos(beta)[:,None]*e)
j=wedge(q,v,nm,hm);source=observer_transport(q,v,nm,hm,nd,np.zeros_like(hm))
# Here V is an assigned constant vector but Q is held static; valid pure observer probe.
assert np.max(abs(j-j[0]-integral(t,source)))<3e-6
# Perspective falsifier: an advancing helix returns to its 2D shadow only.
a=np.linspace(0,2*np.pi,401);helix=np.column_stack([np.cos(a),np.sin(a),a])
assert np.linalg.norm(helix[-1,:2]-helix[0,:2])<1e-12
assert np.linalg.norm(helix[-1]-helix[0])>6.
print('SYNTHETIC_PASS: 80 physical force states, 64 historical source RHS probes, charge, observer transport, in-plane invariance, advancing helix')
