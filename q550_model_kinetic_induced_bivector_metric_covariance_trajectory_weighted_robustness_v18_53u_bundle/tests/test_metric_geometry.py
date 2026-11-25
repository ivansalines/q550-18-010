#!/usr/bin/env python3
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from metric_geometry import *

def main():
    assert np.allclose(induced(np.eye(5)),np.eye(10),atol=1e-15)
    A=np.diag([2.,3.,4.,5.,6.]); H=induced(A)
    assert np.allclose(np.diag(H),[A[i,i]*A[j,j] for i,j in PAIRS])
    rng=np.random.default_rng(53053)
    for j in range(30):
        K=np.eye(5)+.2*rng.normal(size=(5,5));L=wedge_map(K)
        G=np.eye(5)+.04*rng.normal(size=(5,5));G=G.T@G
        x=rng.normal(size=10);y=np.linalg.solve(L,x)
        lhs=y@induced(K.T@G@K)@y;rhs=x@induced(G)@x
        assert np.isclose(lhs,rhs,rtol=4e-13,atol=5e-13)
    z=np.zeros(11);z[0]=.1;z[4]=-.4;z[2]=1.1;z[6]=-1.5;z[8]=.7
    x=rng.normal(size=10)
    p=verify_pullback(z,np.eye(5),x)
    assert p['rel_error']<2e-13
    n=rng.normal(size=(19,10))
    l,o,u=oriented_polygon(n,induced(np.eye(5)))
    assert 0<l<19*np.pi and o.shape==(45,)
    S=np.diag([.1,.2,.3,.15,.25]);gg=occupied_metric(S,3)
    assert np.isclose(np.linalg.det(gg),1)
    assert np.linalg.eigvalsh(induced(gg))[0]>0
    # Canonical kinetic norm in log-polar coordinates is a pullback of Euclidean.
    J=kinetic_jacobian(z)
    assert np.allclose(J.T@J,np.diag([np.exp(.2),np.exp(.2),np.exp(-.8),np.exp(-.8),(expit(.7)*(1-expit(.7)))**2]),rtol=1e-12)
    print('SYNTHETIC_MODEL_METRIC_COVARIANCE_PASS')
if __name__=='__main__':main()
