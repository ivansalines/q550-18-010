#!/usr/bin/env python3
"""11D Jacobian and finite-lag configuration tangent map on authentic Q550 11D traces.

Historical model terms are sourced from authenticated v18.53f torque_math.py.
The external forcing phase is recovered by solving two linear equations in
cos(Delta), sin(Delta) from the *archived* authentic 10D source, then checked.
This is a trace-conditioned reconstruction, NOT an independent forcing replay.
The five-dimensional lifts are explicitly non-unique.
"""
from __future__ import annotations
import numpy as np
from scipy.special import expit
from scipy.interpolate import CubicSpline
from scipy.linalg import expm

PAIRS5=tuple((i,j) for i in range(5) for j in range(i+1,5))
CONFIG=(0,2,4,6,8)
VELOCITY=(1,3,5,7,9)

def biv(q,f):
    return np.sqrt(2.)*np.stack([q[...,i]*f[...,j]-q[...,j]*f[...,i] for i,j in PAIRS5],axis=-1)

def parametric_forces(z,kappa,p,torque):
    # The old source is q wedge (radial potential + phase coupling).
    q,vel,_,g=torque.embed(z[:,:10],p['m_x'])
    r1,r2,x=g['r1'],g['r2'],g['x']
    er1,et1,er2,et2=g['er1'],g['et1'],g['er2'],g['et2']
    pot=np.zeros_like(q);fc=np.zeros_like(q);fs=np.zeros_like(q)
    pot[:,:2]=-torque.dU(r1,p['e_reserve'])[:,None]*er1
    pot[:,2:4]=-torque.dU(r2,p['e_reserve'])[:,None]*er2
    fc[:,:2]=(-kappa*x*r2)[:,None]*er1
    fc[:,2:4]=(-kappa*x*r1)[:,None]*er2
    fs[:,:2]=(kappa*x*r2)[:,None]*et1
    fs[:,2:4]=(kappa*x*r1)[:,None]*et2
    return q, pot,fc,fs

def reconstruct_phase(xi,z,g10,omega,kappa,p,torque):
    xi=np.asarray(xi,float);z=np.asarray(z,float);g10=np.asarray(g10,float)
    x=expit(z[:,8]);rate=omega/x
    q,pot,fc,fs=parametric_forces(z,kappa,p,torque)
    rhs=g10/rate[:,None]-biv(q,pot)
    a=biv(q,fc);b=biv(q,fs)
    aa=np.einsum('ni,ni->n',a,a);bb=np.einsum('ni,ni->n',b,b);ab=np.einsum('ni,ni->n',a,b)
    det=aa*bb-ab*ab
    cond=det/np.maximum(aa*bb,1e-30)
    if min(cond)<=1e-11:raise RuntimeError('SOURCE_PHASE_INVERSION_NEAR_SINGULAR')
    ar=np.einsum('ni,ni->n',a,rhs);br=np.einsum('ni,ni->n',b,rhs)
    c=(bb*ar-ab*br)/det;s=(aa*br-ab*ar)/det
    strength=np.hypot(c,s)
    D=np.unwrap(np.arctan2(s,c))
    # phase is the authenticated *external* forcing, held fixed when varying z
    phase=z[:,2]+z[:,6]-p['phi0']-D
    reconstruct=biv(q,pot)+c[:,None]*a+s[:,None]*b
    error=np.linalg.norm(reconstruct-g10/rate[:,None],axis=1)
    diag=dict(max_source_10d_fit_error=float(error.max()),
              max_unit_circle_defect=float(np.max(np.abs(strength-1))),
              min_phase_basis_condition=float(cond.min()),
              omega_estimated=float(omega),
              max_external_phase_increment=float(np.max(np.abs(np.diff(phase)))))
    return CubicSpline(xi,phase,extrapolate=False),diag

def infer_omega(xi,z):
    x=expit(z[:,8]);rate=np.gradient(z[:,10],xi,edge_order=2)*x
    omega=float(np.median(rate))
    diag=dict(min_omega_probe=float(np.min(rate)),max_omega_probe=float(np.max(rate)),
              max_relative_omega_probe_spread=float(np.max(abs(rate-omega))/omega))
    return omega,diag

def config_jacobian(z,mx=1.):
    from math import sqrt
    y=np.asarray(z,float);r1,r2=np.exp(y[0]),np.exp(y[4]);t1,t2=y[2],y[6];x=expit(y[8]);g=x*(1-x)
    er1=np.array([np.cos(t1),np.sin(t1)]);et1=np.array([-np.sin(t1),np.cos(t1)])
    er2=np.array([np.cos(t2),np.sin(t2)]);et2=np.array([-np.sin(t2),np.cos(t2)])
    J=np.zeros((5,11));J[:2,0]=r1*er1;J[:2,2]=r1*et1
    J[2:4,4]=r2*er2;J[2:4,6]=r2*et2;J[4,8]=sqrt(mx)*g
    return J

def physical_velocity_jacobian(z,mx=1.):
    from math import sqrt
    y=np.asarray(z,float);r1,r2=np.exp(y[0]),np.exp(y[4]);t1,t2=y[2],y[6];x=expit(y[8]);g=x*(1-x)
    er1=np.array([np.cos(t1),np.sin(t1)]);et1=np.array([-np.sin(t1),np.cos(t1)])
    er2=np.array([np.cos(t2),np.sin(t2)]);et2=np.array([-np.sin(t2),np.cos(t2)])
    u1=r1*(y[1]*er1+y[3]*et1);u2=r2*(y[5]*er2+y[7]*et2)
    A=np.zeros((5,11));A[:2,0]=u1;A[:2,2]=np.array([-u1[1],u1[0]])
    A[:2,1]=r1*er1;A[:2,3]=r1*et1
    A[2:4,4]=u2;A[2:4,6]=np.array([-u2[1],u2[0]])
    A[2:4,5]=r2*er2;A[2:4,7]=r2*et2
    A[4,8]=sqrt(mx)*g*(1-2*x)*y[9];A[4,9]=sqrt(mx)*g
    return A

def tangent_lift(z,mx=1.,mode='fixed_cartesian_velocity'):
    J=config_jacobian(z,mx);C=J[:,CONFIG]
    S=np.linalg.solve(C,np.eye(5))
    L=np.zeros((11,5));L[CONFIG,:]=S
    if mode=='fixed_cartesian_velocity':
        D=physical_velocity_jacobian(z,mx)
        L[VELOCITY,:]=-np.linalg.solve(D[:,VELOCITY],D[:,CONFIG]@S)
        if np.linalg.norm(D@L)>2e-11:raise RuntimeError('VELOCITY_FIXED_LIFT_INVALID')
    elif mode!='fixed_generalized_velocity':raise ValueError(mode)
    if np.linalg.norm(J@L-np.eye(5))>2e-11:raise RuntimeError('LIFT_NOT_RIGHT_INVERSE')
    return L

def rhs_material_xi(z,phase,omega,kappa,p,torque):
    z=np.asarray(z,float);scale=omega/expit(z[8])
    f=torque.material_rhs(z[:10],float(phase),kappa,p)[0]
    return np.r_[scale*f,scale]

def jacobian_11d(z,phase,omega,kappa,p,torque,step=3e-5):
    J=np.zeros((11,11))
    for j in range(10):
        d=np.zeros(11);d[j]=step
        J[:,j]=(rhs_material_xi(z+d,phase,omega,kappa,p,torque)-
                rhs_material_xi(z-d,phase,omega,kappa,p,torque))/(2*step)
    # Material tau is passive and absent from the physical source.
    return J

def polar_orthogonal(F):
    U,s,Vt=np.linalg.svd(F)
    if s.min()<=1e-10 or np.linalg.det(U@Vt)<=0:raise RuntimeError('TANGENT_MAP_POLAR_SINGULAR')
    return U@Vt,s

def finite_lag_spin(xi,z,phase,omega,kappa,p,torque,times,h,mode):
    """Frozen-trajectory 11D variational midpoint map; no nonlinear 11D IVP."""
    zz=CubicSpline(xi,z,axis=0,extrapolate=False)
    times=np.asarray(times,float);out=np.empty((len(times),5,5));conditioning=[]
    for i,t in enumerate(times):
        if t+h>xi[-1]+1e-12:raise RuntimeError('TANGENT_HORIZON_OUT_OF_TRACE')
        a=zz(t);b=zz(t+h);mid=zz(t+h/2)
        ph=float(phase(t+h/2));A=jacobian_11d(mid,ph,omega,kappa,p,torque)
        M=expm(h*A);F=config_jacobian(b,p['m_x'])@M@tangent_lift(a,p['m_x'],mode)
        U,s=polar_orthogonal(F);conditioning.append(float(s.max()/s.min()))
        # Finite-h spin is a declared discrete diagnostic, NOT an exact generator.
        K=(U-U.T)/(2*h)
        # Complementary rotation does not move the velocity tangent by construction.
        V=torque.embed(a[:10],p['m_x'])[1][0];e=V/np.linalg.norm(V)
        P=np.eye(5)-np.outer(e,e)
        B=P@K@P
        out[i]=(B-B.T)/2
    return out,dict(horizon=h,lift=mode,max_projected_tangent_map_condition=float(max(conditioning)),
                    max_complement_spin_fro=float(np.max(np.linalg.norm(out,axis=(1,2)))),
                    mean_complement_spin_fro=float(np.mean(np.linalg.norm(out,axis=(1,2)))))
