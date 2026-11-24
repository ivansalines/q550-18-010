#!/usr/bin/env python3
"""Synthetic algebra/model falsification controls. No genuine model integrations."""
import importlib.util,sys,math
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('aud52p',ROOT/'scripts/run_audit.py')
a=importlib.util.module_from_spec(spec);sys.modules[spec.name]=a;spec.loader.exec_module(a)
R=a.passive_nodal_basis(-.308);ang=math.degrees(math.atan(-.308))
assert ang<0 and np.allclose(R[0,1]/R[0,0],-.308)
np.random.seed(42)
q=np.random.normal(size=(13,2));p=np.random.normal(size=(13,2))
old=a.canonical_identity(q,p,R)
assert np.allclose(old['J_before'],old['J_after'],atol=1e-14)
assert np.allclose(old['norm_before'],old['norm_after'],atol=1e-14)
S=a.canonical_block(R)
E=np.eye(2);O=np.block([[np.zeros((2,2)),E],[-E,np.zeros((2,2))]])
assert np.linalg.norm(S.T@O@S-O)<1e-13
# A central oscillator: constant angular momentum, radius = 1, exact radial cancellation.
t=np.linspace(0,3*math.pi,513)
q5=np.zeros((len(t),5));v5=np.zeros_like(q5);a5=np.zeros_like(q5)
q5[:,0]=np.cos(t);q5[:,1]=np.sin(t)
v5[:,0]=-np.sin(t);v5[:,1]=np.cos(t)
a5[:,:2]=-q5[:,:2]
G,s=a.geometry(q5,v5,a5,t)
assert np.max(abs(s['j']-1))<1e-14
assert np.max(abs(s['rddot']))<1e-14
assert np.max(abs(s['T']-s['Trot']))<1e-14
assert G['relative_J_sign_flips']==0
assert G['relative_J_relative_range_to_median']<1e-13
assert not G['certified_noether_charge']
# Noncentral acceleration produces torque: exact product rule d/dt(q wedge v)=q wedge a.
q5=np.zeros((len(t),5));v5=np.zeros_like(q5);a5=np.zeros_like(q5)
q5[:,0]=1.+.3*t; q5[:,1]=.2*t*t
v5[:,0]=.3;v5[:,1]=.4*t
# q''=(0,0.4); J=(1+0.3t)*0.4t-(0.2t²)*0.3 increasing
q5[:,2]=.5; a5[:,1]=.4
G,s=a.geometry(q5,v5,a5,t)
assert G['relative_J_relative_range_to_median']>0.1
assert G['dJ_dt_vs_torque_relative_max']<2e-10
assert G['five_dimensional_bivector_squared_identity_relative_max']<1e-11
assert not G['certified_future_nocollapse']
# Native polar pair Jacobian: J=r²*w and canonical coordinate rotation leave it fixed.
y=np.zeros((len(t),10));y[:,0]=np.log(2.);y[:,2]=t;y[:,3]=1.;y[:,4]=np.log(3.);y[:,6]=t;y[:,7]=1.
qq,vv=a.extract_qv_polar(y)
assert np.max(abs(a.cross(qq[:,:2],vv[:,:2])-4.))<1e-12
assert np.max(abs(a.cross(qq[:,2:],vv[:,2:])-9.))<1e-12
# A null angular momentum toy can hit its axis: excludes unconditional collapse theorem.
assert abs(a.cross(np.array([1.,0.]),np.array([-1.,0.])))<1e-14
print('SYNTHETIC PASS: symplectic rotation, oriented nodal sign, central J, noncentral torque, 5D wedge, radial balance, kinetic chart, zero-J control')
