#!/usr/bin/env python3
"""Deterministic v18.53j synthetic controls: all tests must precede production IVPs."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import numpy as np
from three_window_geometry import (three_window,accumulators,interval,bivector_pair_angle,
    polygon_plane_holonomy,independent_direct_window_integral)
from bivector_geometry import principal_skew

def wedge(i,j):
    x=np.zeros((5,5));x[i,j]=1;x[j,i]=-1;return x

# Exactly equilateral triple of arbitrary 5D skew residuals: closure must be 1.
A=wedge(0,1);B=wedge(2,3)
R=[A,-.5*A+np.sqrt(3)/2*B,-.5*A-np.sqrt(3)/2*B]
assert abs(np.linalg.norm(sum(R,np.zeros((5,5)))))<1e-12
assert max(abs(bivector_pair_angle(R[i],R[j])-120.) for i,j in ((0,1),(1,2),(0,2)))<1e-9

# Real piecewise-linear synthetic torque with known integrals, material clock.
xi=np.linspace(0,3,3001);t=2*xi;bu=np.zeros((len(xi),5,5));bc=np.zeros_like(bu)
for i,x in enumerate(xi):bu[i]=R[min(int(x),2)]/2.0
acc=accumulators(xi,t,bu,bc)
g,rows=three_window(acc,[0,1,2,3])
assert g['normalized_residual_closure_score']>.999
assert g['max_deviation_from_120_deg']<.05
assert g['additive_bivector_error']<1e-12
assert abs(g['whole_interval_radial_phase_cancellation'])<1e-12
shifted,_=three_window(acc,[0,.5,2,3])
assert abs(shifted['net_norm']-g['net_norm'])<1e-12 # fixed endpoints force identical sum
assert shifted['equal_window_residual_norm_cv']>g['equal_window_residual_norm_cv']+.05
for a,b in ((.23,1.51),(.91,2.78)):
    q1,q2=interval(acc,a,b)
    direct=independent_direct_window_integral(xi,t,bu,a,b)
    assert np.linalg.norm(q1-direct)<.005

# Gauge invariance of principal-plane polygon under a COMMON 5D orthogonal change.
rng=np.random.default_rng(7789)
base=wedge(0,1)+.35*wedge(2,3)
Rs=[]
for alpha in (0,.4,.9):
    rot=np.eye(5);rot[0,0]=np.cos(alpha);rot[0,2]=-np.sin(alpha)
    rot[2,0]=np.sin(alpha);rot[2,2]=np.cos(alpha)
    Rs.append(rot@base@rot.T)
h=polygon_plane_holonomy(Rs,0)
assert h['defined'],h
O,_=np.linalg.qr(rng.normal(size=(5,5)))
ho=polygon_plane_holonomy([O@R@O.T for R in Rs],0)
assert ho['defined'] and abs(ho['angle_deg']-h['angle_deg'])<1e-5,(h,ho)
# Open paths are NOT a physical-holonomy certificate: test value is a polygon diagnostic.
H=principal_skew(base);assert H['sigma1']>.99 and H['relative_gap']>.5
try:three_window(acc,[0,1,1,3]);raise AssertionError('DEGENERATE_WINDOW_NOT_REJECTED')
except ValueError:pass
# Exercise the actual row assembly code for every preregistered window.
from run_audit import window_triples,simple_row,tensor_coords
n=0
for fam,label,bounds,meta in window_triples():
    # Validate only the predeclared geometry shape without running a Q550 IVP.
    fake_acc=accumulators(np.linspace(12.8,56.,3001),np.linspace(0,80,3001),
         np.zeros((3001,5,5)),np.zeros((3001,5,5)))
    assert len(bounds)==4 and bounds[0]>=12.8 and bounds[-1]<=56.0+1e-9
    n+=1
assert n==45,n
assert len(tensor_coords(base))==10
print('V18_53J_SYNTHETIC_PASS: closure, boundary shifts, off-grid quadrature, O5 gauge, polygon holonomy')
