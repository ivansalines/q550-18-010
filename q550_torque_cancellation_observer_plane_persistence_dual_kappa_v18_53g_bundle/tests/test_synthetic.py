#!/usr/bin/env python3
"""Source-independent falsification for observer-projection and κ rescaling."""
import json,math,sys
from pathlib import Path
import numpy as np
r=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(r/'scripts'))
from torque_math import force_components,material_rhs,acceleration_from_rhs,embed,wedge,axes
from run_audit import wedge_bivector,binvariants,cscore,preview_from_previous
C=json.loads((r/'inputs_config.json').read_text());P=C['params'];K=C['kappas']
rng=np.random.default_rng(531953)
y=np.column_stack([rng.uniform(-1,.5,37),rng.normal(0,.1,37),rng.uniform(-2,2,37),rng.uniform(-.8,.8,37),rng.uniform(-1,.5,37),rng.normal(0,.1,37),rng.uniform(-2,2,37),rng.uniform(-.8,.8,37),rng.normal(0,.3,37),rng.normal(0,.2,37)])
ph=rng.uniform(-2,2,37);t=np.linspace(0,1,37)
q,v,j,g=embed(y,P['m_x'])
for kn,k in K.items():
 Q,V,J,extra,delta,f=force_components(y,ph,k,P)
 fmaterial=material_rhs(y,ph,k,P)
 assert np.max(abs(acceleration_from_rhs(y,fmaterial,P['m_x'])-sum(f.values())))<1e-11
 assert np.max(abs(q-Q))<1e-15
 # Physical counter-rotation Jdot is zero at each state, including counterfactual k.
 r1=extra['r1'];r2=extra['r2'];jd=r1*r1*(2*y[:,1]*y[:,3]+fmaterial[:,3])-r2*r2*(2*y[:,5]*y[:,7]+fmaterial[:,7])
 assert np.max(abs(jd))<1e-11
 m1=wedge_bivector(Q,f['radial_potential'],t);m2=wedge_bivector(Q,f['phase_coupling'],t)
 assert np.max(abs(m1+m1.T))<1e-13
 inv=binvariants(m1,m2)
 assert inv['bivector_cancel_score'] is not None
 e,n,h=axes(-.308,[.26,-.41,.874])
 for angle in C['observer_angles_degrees']:
  aa=np.deg2rad(angle);na=np.cos(aa)*n+np.sin(aa)*e
  direct=float(np.trapezoid(wedge(Q,f['phase_coupling'],na,h),t))
  assert abs(direct-float(na@m2@h))<1e-12
  ortho=np.eye(5);ortho[0,0]=np.cos(.4);ortho[0,1]=-np.sin(.4);ortho[1,0]=np.sin(.4);ortho[1,1]=np.cos(.4)
  m1r=ortho@m1@ortho.T;m2r=ortho@m2@ortho.T
  assert abs(inv['bivector_cancel_score']-binvariants(m1r,m2r)['bivector_cancel_score'])<1e-14
# Highly cancelling in one 2D projection, yet not in full 5D.
M=np.zeros((5,5));N=np.zeros((5,5))
M[0,1]=1;M[1,0]=-1
N[0,1]=-1;N[1,0]=1;N[2,3]=10;N[3,2]=-10
assert cscore(M[0,1],N[0,1])==1.
assert binvariants(M,N)['bivector_cancel_score']<.2
# In-plane axes rotation cannot change projected bivector.
a=.42;e,n,h=axes(-.308,[.26,-.41,.874]);na=np.cos(a)*n+np.sin(a)*h;ha=-np.sin(a)*n+np.cos(a)*h
assert abs(float(n@M@h)-float(na@M@ha))<1e-14
print('SYNTHETIC_PASS: 37 source states x 3 kappa, 7 angles, bivector frame invariance, false scalar-cancellation falsifier')

# Synthetic end-to-end analysis path, including all 3 kappa values, every
# predefined observer angle, source spotchecks and 53f baseline reconciliation.
# The fake generator is explicitly test-only; it must NOT provide scientific evidence.
from types import SimpleNamespace
import tempfile
import run_audit as runner
class SyntheticPhysical:
 def rhs_from_original(self, yy, kk, par, slope, clock, support, audit):
  yy=np.asarray(yy,float)
  f=material_rhs(yy[:10],yy[10],kk,P)[0]
  x=float(embed(yy[:10],P['m_x'])[3]['x'][0])
  return np.r_[f*.844357915384206/x,.844357915384206*clock]
class SyntheticClock:
 @staticmethod
 def cubic_clock_sign(*a):return 1.
class SyntheticDesign:
 @staticmethod
 def Design(*a):return None
fake=SimpleNamespace(omega=.844357915384206,params=SimpleNamespace(**P),slope=.1343837359721647,
 physical=SyntheticPhysical(),clock=SyntheticClock(),original_ramp=SyntheticDesign(),support=None,audit=None)
with tempfile.TemporaryDirectory(prefix='v1853g_synthetic_') as td:
 tt=np.linspace(0,2,81)
 z=np.zeros((81,11));z[:,0]=.4+.02*np.sin(tt);z[:,1]=.02*np.cos(tt)
 z[:,2]=.2+.1*tt;z[:,3]=.1;z[:,4]=-.7+.015*np.sin(tt*.8)
 z[:,5]=.012*np.cos(tt*.8);z[:,6]=.1+.04*tt;z[:,7]=.04
 z[:,8]=.02*np.sin(tt);z[:,9]=.02*np.cos(tt);z[:,10]=tt
 # Construct a synthetic path on an exactly fixed-counterrotation-charge leaf.
 J0=np.exp(2*z[0,0])*z[0,3]-np.exp(2*z[0,4])*z[0,7]
 z[:,7]=(np.exp(2*z[:,0])*z[:,3]-J0)/np.exp(2*z[:,4])
 labels=['original__reference__Y_tau','autonomous__reference__Y_tau']
 fn=Path(td)/'synthetic_native.npz'
 np.savez_compressed(fn,common_xi=tt,**{key:z for key in labels})
 ctx={key.split('_')[0]:{'theta':lambda x:0.0,'event':0.} for key in labels}
 prior={}
 obs=json.loads((r/'reference/v18_52n_authenticated_report.json').read_text())
 _,n,h=axes(C['observer_slope'],obs['pca_hidden'])
 for key in labels:
  forces=force_components(z[:,:10],np.zeros(len(tt)),K['shared'],P)[-1]
  prior[key]={'force_torque_components_integrated':{
   name:float(np.trapezoid(wedge(qfield,f,n,h),tt))
   for name,f in forces.items() for qfield in [embed(z[:,:10],P['m_x'])[0]]}}
 keep=(runner.OUT,runner.C['native_10d_sha256'],runner.C['sampled_series_expected'],runner.C['native_stride'],runner.C['source_probe_count_per_series'])
 try:
  runner.OUT=Path(td);runner.C['native_10d_sha256']=runner.sha(fn)
  runner.C['sampled_series_expected']=2;runner.C['native_stride']=1;runner.C['source_probe_count_per_series']=3
  rows,bi,metrics=runner.analyse(fn,None,fake,ctx,prior)
  assert len(rows)==2*3*7 and len(bi)==2*3
  assert metrics['max_source_field_relative_error']<1e-12
  assert metrics['max_previous_baseline_integral_abs_error']<1e-12
 finally:
  runner.OUT,runner.C['native_10d_sha256'],runner.C['sampled_series_expected'],runner.C['native_stride'],runner.C['source_probe_count_per_series']=keep
print('SYNTHETIC_PIPELINE_PASS: 42 full observer/kappa rows, 6 coordinate-free bivectors, source and genealogy checks')
