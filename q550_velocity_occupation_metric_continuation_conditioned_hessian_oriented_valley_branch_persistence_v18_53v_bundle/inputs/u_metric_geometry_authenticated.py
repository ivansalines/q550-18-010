#!/usr/bin/env python3
"""Physical kinetic 5D exterior-square metric and explicitly noncanonical velocity-occupation families."""
import numpy as np
from scipy.special import expit

PAIRS=tuple((i,j) for i in range(5) for j in range(i+1,5))
PAIRS45=tuple((i,j) for i in range(10) for j in range(i+1,10))

def induced(G):
    G=np.asarray(G,float)
    if G.shape!=(5,5):raise ValueError('EXPECTED_5x5_METRIC')
    H=np.array([[G[i,k]*G[j,l]-G[i,l]*G[j,k] for k,l in PAIRS] for i,j in PAIRS])
    return (H+H.T)/2

def wedge_map(J):
    J=np.asarray(J,float)
    if J.shape!=(5,5):raise ValueError('EXPECTED_5x5_JACOBIAN')
    return np.array([[J[i,k]*J[j,l]-J[i,l]*J[j,k] for k,l in PAIRS] for i,j in PAIRS])

def kinetic_jacobian(z,mx=1.):
    z=np.asarray(z,float)
    if z.shape!=(11,):raise ValueError('EXPECTED_11D_STATE')
    if not mx>0:raise ValueError('NONPOSITIVE_MODEL_MASS')
    r1,r2=np.exp(z[0]),np.exp(z[4]);t1,t2=z[2],z[6]
    x=expit(z[8]);g=x*(1-x)
    J=np.zeros((5,5))
    J[:2,0]=r1*np.array([np.cos(t1),np.sin(t1)])
    J[:2,1]=r1*np.array([-np.sin(t1),np.cos(t1)])
    J[2:4,2]=r2*np.array([np.cos(t2),np.sin(t2)])
    J[2:4,3]=r2*np.array([-np.sin(t2),np.cos(t2)])
    J[4,4]=np.sqrt(mx)*g
    return J

def physical_velocity(z,mx=1.):
    z=np.asarray(z,float)
    if z.ndim!=2 or z.shape[1]!=11:raise ValueError('EXPECTED_Nx11_TRACE')
    r1=np.exp(z[:,0]);r2=np.exp(z[:,4]);t1=z[:,2];t2=z[:,6]
    v1=z[:,1];v2=z[:,5];w1=z[:,3];w2=z[:,7]
    x=expit(z[:,8]);xd=x*(1-x)*z[:,9]
    V=np.column_stack([r1*(v1*np.cos(t1)-w1*np.sin(t1)),
                       r1*(v1*np.sin(t1)+w1*np.cos(t1)),
                       r2*(v2*np.cos(t2)-w2*np.sin(t2)),
                       r2*(v2*np.sin(t2)+w2*np.cos(t2)),
                       np.sqrt(mx)*xd])
    if not np.all(np.isfinite(V)):raise RuntimeError('BAD_KINETIC_VELOCITY')
    return V

def occupation_covariance(xi,z,mx,epoch='preformal'):
    xi=np.asarray(xi,float);z=np.asarray(z,float)
    if len(xi)!=len(z) or np.any(np.diff(xi)<=0):raise ValueError('BAD_TRACE')
    select=(xi<=25.6+1e-10) if epoch=='preformal' else np.ones(len(xi),bool)
    t=xi[select];v=physical_velocity(z[select],mx)
    # Physical time dtau/dxi proportional to 1/x in the authenticated clock.
    # The overall omega factor cancels when normalizing the covariance trace.
    x=expit(z[select,8]); wt=1/x
    if len(t)<10:raise RuntimeError('INSUFFICIENT_EPOCH')
    dt=np.diff(t)
    A=np.sum(dt[:,None,None]*.5*(wt[1:,None,None]*v[1:,:,None]*v[1:,None,:]
                +wt[:-1,None,None]*v[:-1,:,None]*v[:-1,None,:]),axis=0)
    tr=np.trace(A)
    if tr<=1e-15:raise RuntimeError('ZERO_KINETIC_OCCUPATION')
    return (A+A.T)/(2*tr)

def occupied_metric(S,strength):
    if strength<0:raise ValueError('NEGATIVE_STRENGTH')
    G=np.eye(5)+float(strength)*S
    sign,logdet=np.linalg.slogdet(G)
    if sign<=0:raise RuntimeError('NON_SPD_METRIC')
    G=G*np.exp(-logdet/5) # volume normalization, explicitly not model-unique
    if np.linalg.eigvalsh(G)[0]<=0:raise RuntimeError('NON_SPD_METRIC')
    return G

def metric_norm(x,H):
    x=np.asarray(x,float)
    p=np.einsum('...i,ij,...j->...',x,H,x)
    if np.any(p<=0):raise RuntimeError('ZERO_OR_NEGATIVE_METRIC_NORM')
    return np.sqrt(p)

def canonical_sqrt(H):
    e,V=np.linalg.eigh(H)
    if e.min()<=0:raise RuntimeError('NON_SPD_INDUCED_METRIC')
    return (V*np.sqrt(e))@V.T

def oriented_polygon(N,H):
    Q=canonical_sqrt(H)
    q=np.asarray(N)@Q.T
    u=q/np.linalg.norm(q,axis=1)[:,None]
    a=u[:-1];b=u[1:]
    angle=2*np.arctan2(np.linalg.norm(a-b,axis=1),np.linalg.norm(a+b,axis=1))
    f=np.where(angle<1e-5,1+angle**2/6+7*angle**4/360,
                 angle/np.maximum(np.sin(angle),1e-300))
    wedge=np.array([np.sum((a[:,i]*b[:,j]-a[:,j]*b[:,i])*f) for i,j in PAIRS45])
    return float(angle.sum()),wedge,u

def verify_pullback(z,Gcart,x,mx=1.):
    J=kinetic_jacobian(z,mx);L=wedge_map(J)
    Gmat=J.T@Gcart@J
    Hmat=induced(Gmat)
    Hq=induced(Gcart)
    xmat=np.linalg.solve(L,np.asarray(x,float))
    a=xmat@Hmat@xmat;b=np.asarray(x,float)@Hq@np.asarray(x,float)
    return dict(rel_error=abs(a-b)/max(abs(b),1e-20), jac_cond=np.linalg.cond(J),
                metric_min_eig=float(np.linalg.eigvalsh(Hmat)[0]))
