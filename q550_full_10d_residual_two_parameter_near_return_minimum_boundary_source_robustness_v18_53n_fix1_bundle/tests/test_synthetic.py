#!/usr/bin/env python3
"""Lightweight tests (NO Q550 physical integration)."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from full_net_geometry import (check_box,make_eval,local_values,solve_minimum,
    projected_kkt,knot_clearance,unique_minima,endpoint_source_diff)

def run():
    x=np.linspace(0,12,1001)
    # Smooth synthetic 10D history in native interpolation (not a Q550 source).
    C=np.zeros((len(x),10))
    C[:,0]=np.sin(x);C[:,1]=.5*np.cos(x);C[:,2]=.02*x;C[:,3]=.007*x*x
    C[:,4]=np.cos(x*.3); C[:,5]=.02*x**1.5
    at,der=make_eval(x,C)
    p=np.array([1.,2.2]);f,g,J,N=local_values(p,at,der)
    for j in range(2):
        e=np.zeros(2);e[j]=1e-5
        fp=local_values(p+e,at,der)[0];fm=local_values(p-e,at,der)[0]
        assert abs((fp-fm)/(2e-5)-g[j])<1e-5
    assert np.allclose(g,J.T@N) and f>0
    assert np.allclose(C[0],at(x[0]))
    atp,dp=make_eval(x,C,'pchip')
    assert np.all(np.isfinite(local_values(p,atp,dp)[1]))
    orange=(.4,2.5);wrange=(1.8,2.6)
    check_box(x,orange,wrange)
    result=solve_minimum(at,der,[1.,2.1],orange,wrange)
    assert result['optimizer_success'] and result['net_norm']>0
    assert np.isfinite(result['scaled_jac_condition'])
    # Crucial counterexample: exact zero projection with a nonzero full residual.
    C2=np.zeros_like(C);C2[:,0]=np.sin(2*np.pi*x/6)
    C2[:,1]=np.cos(2*np.pi*x/6);C2[:,2]=.04*x
    a2,d2=make_eval(x,C2)
    n2,_=__import__('full_net_geometry').net_jac((1.,2.),a2,d2)
    assert np.linalg.norm(n2[:2])<1e-12 and np.linalg.norm(n2)>0.2
    try:check_box(x,(-.1,1.),wrange)
    except ValueError:pass
    else:raise AssertionError('OUTSIDE_TRACE_START_NOT_CAUGHT')
    try:check_box(x,(4.,7.),wrange)
    except ValueError:pass
    else:raise AssertionError('OUTSIDE_TRACE_END_NOT_CAUGHT')
    # Stationarity at lower bound requires positive gradient and KKT=0.
    assert projected_kkt(np.array([1.,0.]),np.array([0.,1.]),[0.,0.],[2.,2.])==0
    assert projected_kkt(np.array([-1.,0.]),np.array([0.,1.]),[0.,0.],[2.,2.])>0
    assert knot_clearance(x,1.)<=.01
    # v18.53n-fix1 regression: this exact helper raised NameError(net_j)
    # after the first genuine physical integration in the failed sendback.
    x3=np.linspace(0.,12.,241)
    rate3=np.ones(len(x3))
    B3=np.zeros((len(x3),10));B3[:,0]=2.0;B3[:,1]=-3.0
    C3=np.column_stack([2*x3,-3*x3,np.zeros((len(x3),8))])
    err,ref=endpoint_source_diff(x3,C3,B3,rate3,(1.,2.0))
    assert np.isfinite(err) and np.isfinite(ref)
    assert err<1e-10, ('ENDPOINT_SOURCE_REGRESSION',err,ref)
    # Nonconstant source: verify exact boundary identity for a smooth synthetic
    # primitive away from knots using its own piecewise-linear slopes.
    left=1.125; width=2.225; right=left+3*width
    at3,sl3=make_eval(x3,C3)
    n3,J3=__import__('full_net_geometry').net_jac((left,width),at3,sl3)
    expected=np.column_stack([sl3(right)-sl3(left),3*sl3(right)])
    assert np.allclose(J3,expected,atol=1e-12)
    assert len(unique_minima([dict(origin_xi=1.,width_xi=2.,net_norm=.1),
       dict(origin_xi=1.001,width_xi=2.0001,net_norm=.11)]))==1
    print('SYNTHETIC_FULL_10D_MINIMUM_AND_BOUNDARY_TEST_FIX1 PASS')
if __name__=='__main__':run()
