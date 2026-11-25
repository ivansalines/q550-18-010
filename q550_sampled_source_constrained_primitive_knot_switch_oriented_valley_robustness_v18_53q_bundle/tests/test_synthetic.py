#!/usr/bin/env python3
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import numpy as np
from source_valley_math import Primitive,optimized_width,curve,angular_steps
x=np.linspace(0,10,201)
c=np.stack([.07*x+np.sin(x/3),np.cos(x/3),x*.1,np.zeros_like(x),np.ones_like(x)*.25,
     x*.07,x*.01,np.sin(x/4),np.cos(x/5),np.ones_like(x)*.3],axis=-1)
g=np.stack([.07+np.cos(x/3)/3,-np.sin(x/3)/3,np.full_like(x,.1),np.zeros_like(x),np.zeros_like(x),
    np.full_like(x,.07),np.full_like(x,.01),np.cos(x/4)/4,-np.sin(x/5)/5,np.zeros_like(x)],axis=-1)
for mode in ['pchip','hermite','source_integral']:
 p=Primitive(x,c,g,mode)
 assert p.eval(2.5).shape==(10,)
 assert p.eval(np.array([1.0,2.0])).shape==(2,10)
 if mode=='hermite':assert np.max(abs(p.eval(x,1)-g))<1e-12
 if mode=='source_integral':assert np.max(p.node_mismatch())<.002
 # bounded minimizer never worse than sampled candidates
 o=1.1;lo=.22;hi=.75
 m=optimized_width(p,o,lo,hi,nsearch=25)
 vals=np.array([np.linalg.norm(p.eval(o+3*w)-p.eval(o)) for w in np.linspace(lo,hi,25)])
 assert m['norm']<=vals.min()+1e-9
 t=np.linspace(1.05,1.3,31);w=np.full_like(t,.5)
 v=curve(p,t,w)
 assert abs(v['direct']-v['chain'])<0.005
 assert v['min_norm']>1e-5 and v['width_total_variation']==0.0
 assert abs(np.linalg.norm(v['u'][0])-1)<1e-12
u=np.stack([np.cos(np.linspace(0,1,100)),np.sin(np.linspace(0,1,100))],axis=1)
assert abs(float(sum(angular_steps(u)))-1)<1e-12
print('V1853Q_SYNTHETIC_SOURCE_AND_GEOMETRY_TEST_PASS')

# Helical near-return: positive residual but interior width optimum.
xx=np.linspace(.1,12,1001)
cc=np.zeros((len(xx),10));gg=np.zeros_like(cc)
cc[:,0]=np.cos(xx);cc[:,1]=np.sin(xx);cc[:,2]=.01*xx
gg[:,0]=-np.sin(xx);gg[:,1]=np.cos(xx);gg[:,2]=.01
p=Primitive(xx,cc,gg,'hermite')
opt=optimized_width(p,1.0,1.9,2.3,nsearch=51)
assert opt['label'].startswith('stationary') and 2.0<opt['w']<2.2 and opt['norm']>.01
assert opt['curvature']>0 and np.isfinite(opt['implicit_wprime'])
assert abs(opt['implicit_wprime'])<.01
print('V1853Q_HELICAL_INTERIOR_STATIONARITY_TEST_PASS')
