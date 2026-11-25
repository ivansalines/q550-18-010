#!/usr/bin/env python3
"""Synthetic tests of the 11D-trajectory-induced SO(5) transport algebra."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import numpy as np
from scipy.linalg import expm
from connection_math import minrot,KinematicFrame,exterior_generator,wedge5,rot5,gauge,geodesic_polygon

def main():
    rng=np.random.default_rng(4753)
    for _ in range(16):
        a=rng.normal(size=5);a/=np.linalg.norm(a)
        b=a+.1*rng.normal(size=5);b/=np.linalg.norm(b)
        M=minrot(a,b)
        assert np.linalg.norm(M@a-b)<2e-14
        assert np.linalg.norm(M.T@M-np.eye(5))<2e-14
        A=rng.normal(size=(5,5));A=A-A.T
        assert np.linalg.norm(exterior_generator(A).T+exterior_generator(A))<5e-14
        assert np.linalg.norm((wedge5(expm(1e-7*A))-np.eye(10))/1e-7-exterior_generator(A))<5e-5
    t=np.linspace(0,2,801)
    vel=np.stack([np.cos(.31*t),np.sin(.31*t),.4*np.cos(.17*t),.2*np.sin(.17*t),np.full_like(t,.6)],axis=1)
    frame=KinematicFrame(t,vel);coarse=KinematicFrame(t,vel,4)
    assert frame.orthogonality<1e-11 and frame.tangent_error<1e-11
    assert np.linalg.norm(frame.at(np.array([.1,1.5]))-coarse.at(np.array([.1,1.5])))<1e-4
    at=.75;Q,Qp=gauge(at,t[0],t[-1]);A=frame.omega(at)
    Ag=Q@A@Q.T+Qp@Q.T
    QQ=wedge5(Q);QQp=exterior_generator(Qp@Q.T)@QQ
    assert np.linalg.norm(exterior_generator(Ag)-(QQ@exterior_generator(A)@QQ.T+QQp@QQ.T))<2e-13
    N=rng.normal(size=(80,10));a,p,u=geodesic_polygon(N)
    assert a>0 and np.linalg.norm(p)>0 and np.allclose(np.linalg.norm(u,axis=1),1)
    print('SYNTHETIC_MINIMAL_TURN_EXTERIOR_GENERATOR_GAUGE_PASS')

if __name__=='__main__':main()
