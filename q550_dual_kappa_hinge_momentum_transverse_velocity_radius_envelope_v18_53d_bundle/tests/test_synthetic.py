import importlib.util,sys,math
from pathlib import Path
import numpy as np
r=Path(__file__).resolve().parents[1];spec=importlib.util.spec_from_file_location('audit',r/'scripts/run_audit.py');a=importlib.util.module_from_spec(spec);sys.modules[spec.name]=a;spec.loader.exec_module(a)
# Different notions of return: a 3D helix can close in the 2D shadow yet advance in depth.
t=np.linspace(0,2*math.pi,33);q=np.column_stack([np.cos(t),np.sin(t),t/(2*math.pi),np.zeros_like(t),np.zeros_like(t)])
v=np.column_stack([-np.sin(t),np.cos(t),np.full_like(t,1/(2*math.pi)),np.zeros_like(t),np.zeros_like(t)])
e=np.array([1.,0,0,0,0]);w=np.array([0.,1,0,0,0]);d=np.array([0.,0,1,0,0]);dq=q[-1]-q[0]
assert np.linalg.norm([dq@e,dq@w])<1e-14 and abs(dq@d-1.)<1e-14 and np.linalg.norm(dq)>0.9
pair=a.closest_separated_pair(t,q,v,e,w,d,1.)
assert pair['screen_gap_rel']<1e-12 and pair['unseen_gap_rel']>0.9
# Kinematic wedge basis rotation in a fixed physical 2-plane is invariant.
q=np.array([[.3,.2,.1,0,0],[-.4,.25,.5,0,0]]);v=np.array([[.1,-.2,.7,0,0],[.03,-.5,.4,0,0]])
n=np.array([0.,1,0,0,0]);h=np.array([0.,0,1.,0,0]);c=math.cos(.4);s=math.sin(.4)
assert np.allclose(a.wedge(q,v,n,h),a.wedge(q,v,c*n+s*h,-s*n+c*h),rtol=0,atol=2e-16)
for ang in [-35,-17.04,0,17.04,35]:a.frame(-.308,[-.4105034129634594,.9094294538464419,-.06652079691279893],ang)
assert abs(a.rotor(1.7,.3,-.46,.93)['J']-(1.7**2*(-.46)-.3**2*.93))<1e-15
p,b=a.lineage();assert p['continuation_counts']['dual']==102
print('SYNTHETIC_PASS helix_shadow_not_full_closure wedge_covariance authentic_lineage')
