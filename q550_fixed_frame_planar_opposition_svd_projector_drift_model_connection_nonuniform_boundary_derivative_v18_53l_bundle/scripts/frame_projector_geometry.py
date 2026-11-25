#!/usr/bin/env python3
"""v18.53l geometry. All 10D coordinates are Frobenius-isometric 5x5 skews.

The model connection is the MINIMAL RADIAL 5D kinetic-frame connection
associated with the actually integrated configuration q(xi):
  u = q/|q|,
  Omega = u' u.T - u u'.T,
  dR/dxi = Omega R.
This particular SO(5) lift is declared; it is NOT a unique physical law.
Its exterior-square lift acts on 5x5 bivectors by B -> R B R.T.
"""
from __future__ import annotations
import numpy as np
import math

PAIRS=[(i,j) for i in range(5) for j in range(i+1,5)]

def to10(A):
    A=np.asarray(A,float)
    if A.shape!=(5,5):raise ValueError('REQUIRE_5X5')
    if np.linalg.norm(A+A.T)>1e-9*max(1.,np.linalg.norm(A)):
        raise ValueError('NOT_SKEW')
    return np.sqrt(2.)*np.array([A[i,j] for i,j in PAIRS])

def from10(v):
    v=np.asarray(v,float)
    if v.shape!=(10,):raise ValueError('REQUIRE_10')
    A=np.zeros((5,5))
    for x,(i,j) in zip(v/np.sqrt(2.),PAIRS):A[i,j]=x;A[j,i]=-x
    return A

def svd_rank2(M):
    M=np.asarray(M,float)
    if M.shape!=(3,10):raise ValueError('THREE_BIVECTORS_REQUIRED')
    U,s,Vt=np.linalg.svd(M,full_matrices=False)
    P=Vt[:2].T@Vt[:2]
    return P,s

def projector_angle(P,Q):
    # equal-rank projector principal angles, report max and rms
    w=np.linalg.eigvalsh((P@Q@P+P@Q.T@P)/2)
    c=np.sqrt(np.clip(np.sort(w)[-2:],0,1))
    return float(np.degrees(np.arccos(np.clip(c.min(),0,1)))),float(np.linalg.norm(P-Q)/np.sqrt(2.))

def unit_rows(q):
    q=np.asarray(q,float)
    n=np.linalg.norm(q,axis=1)
    if np.any(n<1e-9):raise ValueError('DEGENERATE_PHYSICAL_Q')
    return q/n[:,None]

def minimal_rotation(u,v):
    """Minimal SO(5) rotor u -> v, well-conditioned away from antipodes."""
    u=np.asarray(u,float);v=np.asarray(v,float)
    u=u/np.linalg.norm(u);v=v/np.linalg.norm(v)
    c=float(np.dot(u,v))
    if c<-0.95:raise ValueError('RADIAL_TRANSPORT_NEAR_ANTIPODAL')
    A=np.outer(v,u)-np.outer(u,v)
    R=np.eye(len(u))+A+A@A/(1.+c)
    return R

def transport_curve(xi,Q):
    """Piecewise-minimal radial frame, a discrete covariant model lift."""
    xi=np.asarray(xi,float)
    if np.any(np.diff(xi)<=0):raise ValueError('XI_NOT_INCREASING')
    U=unit_rows(Q)
    R=np.empty((len(xi),5,5))
    R[0]=np.eye(5)
    for j in range(1,len(xi)):
        R[j]=minimal_rotation(U[j-1],U[j])@R[j-1]
    return dict(xi=xi,U=U,R=R)

def frame_at(C,x):
    xi=C['xi'];U=C['U'];R=C['R']
    if x<xi[0]-1e-9 or x>xi[-1]+1e-9:
        raise ValueError('MODEL_FRAME_OUTSIDE_TRACE')
    i=int(np.searchsorted(xi,x,side='right'))-1
    if i<0:return R[0],U[0]
    if i>=len(xi)-1:return R[-1],U[-1]
    a=float((x-xi[i])/(xi[i+1]-xi[i]))
    u=(1-a)*U[i]+a*U[i+1];u=u/np.linalg.norm(u)
    return minimal_rotation(U[i],u)@R[i],u

def exterior_lift(R):
    """10x10 exterior-square action induced by 5D orthogonal R."""
    R=np.asarray(R,float)
    if R.shape!=(5,5):raise ValueError('R_MUST_5X5')
    L=np.zeros((10,10))
    for j in range(10):
        e=np.zeros(10);e[j]=1.
        L[:,j]=to10(R@from10(e)@R.T)
    return L

def dynamic_projector(P0,C,anchor,target):
    A,u0=frame_at(C,anchor); B,u1=frame_at(C,target)
    T=B@A.T
    L=exterior_lift(T)
    P=L@P0@L.T
    return P,L,T,float(np.linalg.norm(T@u0-u1))

def cosine(a,b,floor=1e-8):
    na=np.linalg.norm(a);nb=np.linalg.norm(b)
    if min(na,nb)<floor:return None
    return float(np.clip(np.dot(a,b)/(na*nb),-1.,1.))

def finite_difference_weights(x,at):
    """Nonuniform polynomial-exact 1st derivative weights at `at`."""
    x=np.asarray(x,float);d=x-at
    if len(np.unique(x))!=len(x):raise ValueError('DUPLICATE_GRID_NODES')
    # Column j = d^k/k! in row j, derivative condition for polynomial powers.
    A=np.stack([d**k/float(math.factorial(k)) for k in range(len(x))],axis=0)
    b=np.zeros(len(x));b[1]=1.
    return np.linalg.solve(A,b)

def derivative_sample(x,y,i,count):
    x=np.asarray(x,float);y=np.asarray(y,float)
    n=len(x);half=count//2
    lo=max(0,min(n-count,i-half));idx=np.arange(lo,lo+count)
    w=finite_difference_weights(x[idx],x[i])
    return np.tensordot(w,y[idx],axes=(0,0))

def secant_product_decomposition(P0,P1,N0,N1):
    """Exact midpoint SECANT product rule, not a differentiation approximation.
       P1 N1 - P0 N0 = ((P1+P0)/2)(N1-N0) + (P1-P0)((N1+N0)/2)
    """
    Pm=.5*(P0+P1);Nm=.5*(N0+N1)
    source=Pm@(N1-N0)
    drift=(P1-P0)@Nm
    lhs=P1@N1-P0@N0
    return lhs,source,drift,float(np.linalg.norm(lhs-source-drift))
