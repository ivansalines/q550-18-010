#!/usr/bin/env python3
"""Synthetic mathematical controls only; no physical Q550 trajectory computed."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from frame_projector_geometry import (to10,from10,svd_rank2,projector_angle,unit_rows,
 minimal_rotation,transport_curve,frame_at,exterior_lift,dynamic_projector,
 finite_difference_weights,derivative_sample,secant_product_decomposition,cosine)
rng=np.random.default_rng(1853)
# Isometric 5D bivector embedding
for _ in range(10):
 A=rng.normal(size=(5,5));A=A-A.T
 assert np.linalg.norm(from10(to10(A))-A)<1e-12
 assert abs(np.linalg.norm(A)-np.linalg.norm(to10(A)))<1e-12
# Orthogonal exterior-square action and state-space rotor
u=rng.normal(size=5);u=u/np.linalg.norm(u)
v=u+.1*rng.normal(size=5);v=v/np.linalg.norm(v)
R=minimal_rotation(u,v)
assert np.linalg.norm(R@u-v)<1e-12
assert np.linalg.norm(R.T@R-np.eye(5))<1e-12
L=exterior_lift(R)
assert np.linalg.norm(L.T@L-np.eye(10))<1e-12
# Projector secant exact identity, for different planes and nontrivial net vectors
for _ in range(5):
 M0=rng.normal(size=(3,10));M1=rng.normal(size=(3,10))
 P0,_=svd_rank2(M0);P1,_=svd_rank2(M1)
 N0=rng.normal(size=10);N1=rng.normal(size=10)
 lhs,source,drift,err=secant_product_decomposition(P0,P1,N0,N1)
 assert err<1e-13
 assert np.linalg.norm(lhs-source-drift)<1e-13
 angle,_=projector_angle(P0,P1)
 assert 0<=angle<=90.00000001
# Model-transport of minimal radial direction on an analytic circle
xi=np.linspace(0.,1.,101)
Q=np.column_stack([np.cos(.4*xi),np.sin(.4*xi),np.zeros((len(xi),3))])
C=transport_curve(xi,Q)
P0,_=svd_rank2(rng.normal(size=(3,10)))
Pd,L,R,maperr=dynamic_projector(P0,C,.35,.78)
assert maperr<1e-12
assert np.linalg.norm(Pd@Pd-Pd)<1e-12
assert np.linalg.norm(L.T@L-np.eye(10))<1e-12
# Nonuniform 3-pt stencil exact for a quadratic, 5-pt for quartic.
x=np.array([11.8,12.,12.2,12.4,12.6,12.65,12.7,12.725,12.75,12.775,12.8,12.825,12.85,12.875,12.9,12.95,13.0,13.1,13.2,13.4,13.6])
y2=x**2; y4=(x-12.8)**4
for i in range(1,len(x)-1):
 assert abs(derivative_sample(x,y2,i,3)-2*x[i])<1e-9
for i in range(len(x)):
 assert abs(derivative_sample(x,y4,i,5)-4*(x[i]-12.8)**3)<1e-7
assert cosine(np.array([1.,0.]),np.array([-1.,0.]))==-1.
print('V1853L_SYNTHETIC_GEOMETRY_TESTS_PASS')
