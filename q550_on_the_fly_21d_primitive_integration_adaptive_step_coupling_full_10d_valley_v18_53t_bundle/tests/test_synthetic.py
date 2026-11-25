#!/usr/bin/env python3
"""Synthetic, no native historical data or Q550 heavy IVPs."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from physics_21d import source_10d,source_array,rhs_21d,angular_polygon

class FakeM:
    PAIRS=[(i,j) for i in range(5) for j in range(i+1,5)]
    @staticmethod
    def to10(a):return np.sqrt(2)*np.array([a[i,j] for i,j in FakeM.PAIRS])

class FakeH:
    P={'m_x':1.}
    @staticmethod
    def rhs(engine,ctx,kp,t,z):
        dz=np.zeros(11)
        dz[0]=2
        dz[10]=1
        return dz
    @staticmethod
    def force_components(y,phase,kp,p):
        q=np.zeros((len(y),5));q[:,0]=1.
        pot=np.zeros_like(q);pot[:,1]=1
        fc=np.zeros_like(q)
        return q,None,None,None,None,dict(radial_potential=pot,phase_coupling=fc)
    @staticmethod
    def bivector(q,f):return q[..., :,None]*f[...,None,:]-f[..., :,None]*q[...,None,:]

CTX={'theta':lambda t:0.,'event':0.}


def run():
    y=np.zeros(21)
    v=rhs_21d(FakeH,FakeM,None,CTX,.3,0,y)
    assert v.shape==(21,)
    assert np.max(np.abs(v[:11]-FakeH.rhs(None,CTX,.3,0,y[:11])))==0
    assert np.isclose(v[11],np.sqrt(2))
    assert not np.any(v[12:])
    for k in (0,0.7,1.7):
        y[11:]=np.sin(k)*np.linspace(-100,100,10)
        # 11D RHS is exactly independent of the accumulator.
        assert np.array_equal(rhs_21d(FakeH,FakeM,None,CTX,.3,0,y)[:11],v[:11])
    sol=solve_ivp(lambda t,a:rhs_21d(FakeH,FakeM,None,CTX,.3,t,a),(0,2),np.zeros(21),
                  t_eval=[0,1,2],method='DOP853',rtol=1e-11,atol=1e-13)
    assert sol.success
    assert np.max(np.abs(sol.y[0]-[0,2,4]))<1e-12
    assert np.max(np.abs(sol.y[11]-np.sqrt(2)*np.array([0,1,2])))<1e-12
    xx=np.linspace(0,2,5)
    zz=np.zeros((len(xx),11))
    gg=source_array(FakeH,FakeM,CTX,.3,xx,zz,omega=.5)
    # expit(0)=1/2 => clock rate=1 in this synthetic model.
    assert np.allclose(gg[:,0],np.sqrt(2),atol=1e-14)
    u=np.zeros((501,10));t=np.linspace(0,1,501)
    u[:,0]=np.cos(1.3*t);u[:,1]=np.sin(1.3*t)
    pairs=[(i,j) for i in range(10) for j in range(i+1,10)]
    wedge=lambda a,b:np.stack([a[...,i]*b[...,j]-a[...,j]*b[...,i] for i,j in pairs],axis=-1)
    length,area=angular_polygon(u,wedge)
    assert np.abs(length-1.3)<1e-11, length
    assert abs(area[0]-1.3)<1e-10, area[0]
    print('SYNTHETIC_PASS 21D_PASSIVE_RHS LOCKED_SOURCE CLOCK_CONVERSION ORIENTED_ARCLENGTH')

if __name__=='__main__':run()
