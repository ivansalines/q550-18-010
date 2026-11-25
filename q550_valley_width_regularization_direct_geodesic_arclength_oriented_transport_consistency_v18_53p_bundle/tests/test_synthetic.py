#!/usr/bin/env python3
"""All synthetic tests are mathematical toy fixtures, never physical Q550 data."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import numpy as np
from valley_convergence_math import Primitive,optimum_width,eval_curve,arcs,wedge,area_from_polygon

# Pure great-circle sampled with an arbitrary orthonormal pair in R10.
xi=np.linspace(0,1.75,251)
u=np.zeros((len(xi),10));u[:,0]=np.cos(xi);u[:,1]=np.sin(xi)
theta=arcs(u)
assert abs(np.sum(theta)-1.75)<1e-12
form=area_from_polygon(u,theta)
assert abs(form[0]-1.75)<1e-12
assert np.linalg.norm(form[1:])<1e-12
assert abs(np.linalg.norm(wedge(u[0],(u[1]-u[0])/xi[1]))-1)<1e-4
# Orthogonal changes of the ambient coordinate basis preserve scalar arclength.
rng=np.random.default_rng(4053)
rot,_=np.linalg.qr(rng.normal(size=(10,10)))
assert abs(np.sum(arcs(u@rot))-1.75)<1e-12
assert abs(np.linalg.norm(area_from_polygon(u@rot,arcs(u@rot)))-1.75)<1e-12
# Exact segmentwise global minimum of a deliberately oscillatory interpolant.
t=np.linspace(0,8,1201)
c=np.zeros((len(t),10));c[:,0]=np.cos(t);c[:,1]=np.sin(t)
p=Primitive(t,c,'linear');pp=Primitive(t,c,'pchip')
a=0.0;lo=2.03;hi=2.17
w,r=optimum_width(p,a,lo,hi)
ws=np.linspace(lo,hi,2001)
rr=np.array([np.linalg.norm(p.at(a+3*x)-p.at(a)) for x in ws])
assert lo<=w<=hi and r<=min(rr)+1e-9
wp,rp=optimum_width(pp,a,lo,hi)
assert lo<=wp<=hi and rp<1e-4
# Near-return but nonzero: analytic chain and polygonal arc agree at refinement.
t=np.linspace(0,70,2801)
c=np.zeros((len(t),10));c[:,0]=np.sin(0.1*t);c[:,1]=np.cos(0.12*t);c[:,2]=0.04*t
p=Primitive(t,c,'pchip')
for k in (141,281):
    origins=np.linspace(10,10.7,k);ww=np.full(k,3.2)
    d,s=eval_curve(p,origins,ww)
    assert np.isfinite(d['arc_chain_gradient']) and d['min_norm']>0.02
    assert abs(d['arc_geodesic']-d['arc_chain_gradient'])<3e-5
assert abs(float(np.linalg.norm(s['u'][0]))-1)<1e-12
print('V1853P_SYNTHETIC_GEODESIC_WIDTH_WEDGE_TESTS_PASS')
