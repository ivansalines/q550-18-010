"""Exact kinetic embedding of historical Q550 material state into T(R^5).

Coordinates (solver material chart):
  y=(l1,v1,t1,w1,l2,v2,t2,w2,z,vz),  r=exp(l), x=sigmoid(z).
  T=0.5*(r1^2*(v1^2+w1^2)+r2^2*(v2^2+w2^2)+m_x*(x*(1-x)*vz)^2).

Global Cartesian configuration and physical tau velocities:
  Q=(r1 cos t1,r1 sin t1,r2 cos t2,r2 sin t2,sqrt(m_x)*x)
  V=(r1*(v1 e_r+w1 e_t),r2*(v2 e_r+w2 e_t),sqrt(m_x)*g*vz).

Sasaki metric of flat kinetic base represented as pullback of
   H(y)=(Q, tau_ref * V),   G(y)=DH(y)^T DH(y).
This is Riemannian on material 10D chart for r>0 and 0<x<1; tau_ref
is a DECLARED physical scaling convention, not a discovered invariant.
"""
import numpy as np

def sigmoid(z):
    return 1./(1.+np.exp(-z)) if z>=0 else np.exp(z)/(1.+np.exp(z))

def embedding_and_jac(y,mx,tau_ref=1.):
    y=np.asarray(y,dtype=float)[:10]
    if not(mx>0 and tau_ref>0):raise ValueError('Need positive mx and tau_ref')
    H=np.zeros(10);J=np.zeros((10,10))
    for off,k in ((0,0),(4,2)):
        l,v,theta,w=y[off:off+4]
        r=np.exp(l);c=np.cos(theta);s=np.sin(theta)
        er=np.array([c,s]);et=np.array([-s,c])
        p=r*er;vel=r*(v*er+w*et)
        H[k:k+2]=p;H[5+k:7+k]=tau_ref*vel
        J[k:k+2,off]=p
        J[k:k+2,off+2]=r*et
        J[5+k:7+k,off]=tau_ref*vel
        J[5+k:7+k,off+2]=tau_ref*r*(v*et-w*er)
        J[5+k:7+k,off+1]=tau_ref*r*er
        J[5+k:7+k,off+3]=tau_ref*r*et
    x=sigmoid(y[8]);g=x*(1-x);sx=np.sqrt(mx)
    H[4]=sx*x;H[9]=tau_ref*sx*g*y[9]
    J[4,8]=sx*g
    J[9,8]=tau_ref*sx*g*(1-2*x)*y[9]
    J[9,9]=tau_ref*sx*g
    return H,J

def kinetic_energy(y,mx):
    H,_=embedding_and_jac(y,mx)
    return .5*float(H[5:]@H[5:])

def material_metric(y,mx,tau_ref=1.):
    _,J=embedding_and_jac(y,mx,tau_ref)
    return J.T@J

def metric_deformation(F,start,end,mx,tau_ref=1.):
    _,J0=embedding_and_jac(start,mx,tau_ref)
    _,J1=embedding_and_jac(end,mx,tau_ref)
    B=J1@F@np.linalg.inv(J0)
    C=B.T@B
    stretches=np.linalg.svd(B,compute_uv=False)
    return C,stretches,B

def nonlinear_chart_gauge(F,G0,G1,y0,y1):
    """Nonlinear local chart psi_i(y)=y_i + a_i sin(y_i), positive derivative."""
    a=np.linspace(.025,.11,10)
    D0=np.diag(1+a*np.cos(np.asarray(y0)[:10]));D1=np.diag(1+a*np.cos(np.asarray(y1)[:10]))
    Fp=D1@F@np.linalg.inv(D0)
    G0p=np.linalg.inv(D0).T@G0@np.linalg.inv(D0)
    G1p=np.linalg.inv(D1).T@G1@np.linalg.inv(D1)
    return Fp,G0p,G1p
