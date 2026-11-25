#!/usr/bin/env python3
"""Falsification controls for v18.53e; these are synthetic, not Q550 data."""
import importlib.util,math
from pathlib import Path
import numpy as np
p=Path(__file__).resolve().parents[1]/'scripts/run_audit.py'
s=importlib.util.spec_from_file_location('a',p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
e,n,h,w,d=m.axes(-.308,[.0,1.,0.],0)
assert abs(e@n)<1e-14 and abs(n@h)<1e-14
# A helix can complete a screen turn with advance exclusively in depth.
dq=3*d
assert abs(dq@e)<1e-13 and abs(dq@w)<1e-13 and np.linalg.norm(dq)>2.9
# Fixed pair: the complete displacement is independent of passive observer tilt.
D=np.linalg.norm(dq)
for a in (-35,-17.04,0,17.04,35):
 e1,n1,h1,w1,d1=m.axes(-.308,[0,1,0],a)
 screen=np.linalg.norm([dq@e1,dq@w1]);hidden=math.sqrt(max(0,D*D-screen*screen))
 assert abs(D*D-screen*screen-hidden*hidden)<1e-12
# A planar kinetic area form is invariant under simultaneous passive basis change.
q=np.array([0.5,1.,1.,0.,1.]);v=np.array([-.7,.2,1.,-1.,.4])
J=m.signed_wedge(q,v,n,h)
a=.7;nr=math.cos(a)*n+math.sin(a)*h;hr=-math.sin(a)*n+math.cos(a)*h
assert abs(J-m.signed_wedge(q,v,nr,hr))<1e-13
# A varying hinge moment does NOT mean varying J1-J2; they use different planes.
y=np.array([[math.log(2),0.1,0.2,0.4,math.log(.4),-.1,0.3,.2,0.,.03],[math.log(2),0.1,1.2,0.4,math.log(.4),-.1,0.3,.2,0.,.03]])
Q,V,Jrot=m.embed10(y,1.)
assert abs(Jrot[1]-Jrot[0])<1e-14
assert abs(m.signed_wedge(Q[1],V[1],n,h)-m.signed_wedge(Q[0],V[0],n,h))>1e-4
# The analytic optimal screen cannot exceed any predeclared screen on fixed dq.
qgap=.2*n+.5*h+.03*e
opt=math.degrees(math.atan2(-.2,.5));eo,no,ho,wo,do=m.axes(-.308,[0,1,0],opt)
assert abs(qgap@wo)<1e-13
assert abs(math.hypot(qgap@eo,qgap@wo)-abs(qgap@eo))<1e-12
print('SYNTHETIC_PASS: fixed endpoints, helix hidden advance, wedge rotation, kinetic counterexample, optimal observer')
