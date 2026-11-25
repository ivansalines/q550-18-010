#!/usr/bin/env python3
"""Q550 v18.53f — explicit material-time 5D forces and observer torque.

Flat kinetic 5D embedding is independently pinned to historical v18.49h.
The first ten material-time equations agree with the archived physical source;
the last native NPZ column is material tau, NOT the external phase.
"""
from __future__ import annotations
import math
import numpy as np

def unit(v):
 v=np.asarray(v,dtype=float);m=np.linalg.norm(v)
 if not m>1.e-14:raise ValueError('DEGENERATE_OBSERVER')
 return v/m

def axes(slope,hidden):
 e=unit([1,slope,0,0,0]);n=unit([-slope,1,0,0,0]);h=unit([0,0,*hidden])
 assert np.max(abs(np.stack([e,n,h])@np.stack([e,n,h]).T-np.eye(3)))<1.e-13
 return e,n,h

def wedge(Q,V,n,h):
 return np.sum(Q*np.asarray(n),axis=-1)*np.sum(V*np.asarray(h),axis=-1)-np.sum(Q*np.asarray(h),axis=-1)*np.sum(V*np.asarray(n),axis=-1)

def sigmoid(z):
 z=np.asarray(z,float)
 return np.where(z>=0,1/(1+np.exp(-np.clip(z,-700,700))),np.exp(np.clip(z,-700,700))/(1+np.exp(np.clip(z,-700,700))))

def dU(r,e):return 2*r-2.4*r*r+.8*r*r*r-e*e/(r*r*r)

def dV_barrier(x,B):return B*(2*x-1)/(x*(1-x))

def embed(Y,mx=1.):
 a=np.atleast_2d(np.asarray(Y,float));r1=np.exp(a[:,0]);r2=np.exp(a[:,4]);t1=a[:,2];t2=a[:,6];v1=a[:,1];v2=a[:,5];w1=a[:,3];w2=a[:,7]
 x=sigmoid(a[:,8]);g=x*(1-x);xd=g*a[:,9]
 er1=np.column_stack([np.cos(t1),np.sin(t1)]);et1=np.column_stack([-np.sin(t1),np.cos(t1)])
 er2=np.column_stack([np.cos(t2),np.sin(t2)]);et2=np.column_stack([-np.sin(t2),np.cos(t2)])
 Q=np.column_stack([r1[:,None]*er1,r2[:,None]*er2,np.sqrt(mx)*x]);V=np.column_stack([r1[:,None]*(v1[:,None]*er1+w1[:,None]*et1),r2[:,None]*(v2[:,None]*er2+w2[:,None]*et2),np.sqrt(mx)*xd]);J=r1*r1*w1-r2*r2*w2
 return Q,V,J,dict(r1=r1,r2=r2,x=x,g=g,xd=xd,er1=er1,et1=et1,er2=er2,et2=et2)

def force_components(Y,phase,kappa,p):
 """Physical 5D acceleration separated by literal terms in documented RHS.

 Phase is the external phase of the historical return protocol. Its rate is
 NOT an extra instantaneous force at frozen state; only Delta is present.
 """
 a=np.atleast_2d(np.asarray(Y,float));phase=np.broadcast_to(np.asarray(phase,float),len(a));Q,V,J,g=embed(a,p['m_x'])
 r1=g['r1'];r2=g['r2'];x=g['x'];er1=g['er1'];er2=g['er2'];et1=g['et1'];et2=g['et2'];xd=g['xd']
 D=a[:,2]+a[:,6]-p['phi0']-phase
 def blank():return np.zeros_like(Q)
 pot=blank();pot[:,:2]=-dU(r1,p['e_reserve'])[:,None]*er1;pot[:,2:4]=-dU(r2,p['e_reserve'])[:,None]*er2
 coupling=blank();c1=(-kappa*x*r2*np.cos(D))[:,None]*er1+(kappa*x*r2*np.sin(D))[:,None]*et1
 c2=(-kappa*x*r1*np.cos(D))[:,None]*er2+(kappa*x*r1*np.sin(D))[:,None]*et2
 coupling[:,:2]=c1;coupling[:,2:4]=c2
 radial_diss=blank();radial_diss[:,:2]=(-p['gamma_rho']*r1*a[:,1])[:,None]*er1;radial_diss[:,2:4]=(-p['gamma_rho']*r2*a[:,5])[:,None]*er2
 x_barrier=blank();x_barrier[:,4]=-dV_barrier(x,p['barrier_B'])/np.sqrt(p['m_x'])
 x_activity=blank();x_activity[:,4]=p['mu_x']*r1*r2*np.sin(D)/np.sqrt(p['m_x'])
 x_diss=blank();x_diss[:,4]=-p['gamma_x']*xd/np.sqrt(p['m_x'])
 out=dict(radial_potential=pot,phase_coupling=coupling,radial_dissipation=radial_diss,
   x_barrier=x_barrier,x_activity=x_activity,x_dissipation=x_diss)
 return Q,V,J,g,D,out

def material_rhs(Y,phase,kappa,p):
 a=np.atleast_2d(np.asarray(Y,float));Q,V,J,g,D,C=force_components(a,phase,kappa,p)
 r1=g['r1'];r2=g['r2'];x=g['x'];xd=g['xd'];gg=g['g'];sin=np.sin(D);cos=np.cos(D)
 out=np.zeros((len(a),10))
 out[:,0]=a[:,1];out[:,2]=a[:,3];out[:,4]=a[:,5];out[:,6]=a[:,7];out[:,8]=a[:,9]
 out[:,1]=(r1*a[:,3]**2-dU(r1,p['e_reserve'])-kappa*x*r2*cos-p['gamma_rho']*r1*a[:,1])/r1-a[:,1]**2
 out[:,3]=kappa*x*r2*sin/r1-2*a[:,1]*a[:,3]
 out[:,5]=(r2*a[:,7]**2-dU(r2,p['e_reserve'])-kappa*x*r1*cos-p['gamma_rho']*r2*a[:,5])/r2-a[:,5]**2
 out[:,7]=kappa*x*r1*sin/r2-2*a[:,5]*a[:,7]
 xdd=(-dV_barrier(x,p['barrier_B'])+p['mu_x']*r1*r2*sin-p['gamma_x']*xd)/p['m_x']
 out[:,9]=xdd/gg-(1-2*x)*a[:,9]**2
 return out

def acceleration_from_rhs(Y,F,mx=1.):
 a=np.atleast_2d(np.asarray(Y,float));f=np.atleast_2d(np.asarray(F,float));Q,V,J,g=embed(a,mx);A=np.zeros_like(Q)
 for i,(off,vel,om,er,et,r) in enumerate(((0,1,3,g['er1'],g['et1'],g['r1']),(1,5,7,g['er2'],g['et2'],g['r2']))):
  co=i*2;rd=r*(a[:,vel]**2+f[:,vel]-a[:,om]**2);td=r*(2*a[:,vel]*a[:,om]+f[:,om])
  A[:,co:co+2]=rd[:,None]*er+td[:,None]*et
 x=g['x'];A[:,4]=np.sqrt(mx)*g['g']*((1-2*x)*a[:,9]**2+f[:,9])
 return A

def inplane_basis(n,h,alpha,alpha_dot):
 c=np.cos(alpha);s=np.sin(alpha)
 a=c*n+s*h;b=-s*n+c*h
 ad=alpha_dot*b;bd=-alpha_dot*a
 return a,b,ad,bd

def tilted_plane(e,n,h,beta,beta_dot):
 c=np.cos(beta);s=np.sin(beta)
 return c*n+s*e,h,beta_dot*(-s*n+c*e),np.zeros(5)

def observer_transport(Q,V,n,h,ndot,hdot):
 return wedge(Q,V,ndot,h)+wedge(Q,V,n,hdot)

def integral(t,x):
 t=np.asarray(t,float);x=np.asarray(x,float);dx=(x[1:]+x[:-1])*np.diff(t)/2
 return np.r_[0,np.cumsum(dx)]
