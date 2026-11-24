#!/usr/bin/env python3
"""Exact formulas for archived Q550 physical-time model; no ODE integration.

For Y=(log r1,v1,theta1,omega1,log r2,v2,theta2,omega2,z,vz,Theta):
 r_i'=r_i v_i, x=sigmoid(z), x'=x(1-x)vz, Theta'=x.
Input rows can be shape (...,11). Values are physical-time point functions.
"""
import math
from dataclasses import dataclass
import numpy as np

@dataclass(frozen=True)
class Params:
 phi0:float=.28
 e_reserve:float=.26
 gamma_rho:float=.3
 barrier_B:float=.48
 mu_x:float=.5
 gamma_x:float=.02
 m_x:float=1.2

def sigmoid(z):
 z=np.asarray(z,dtype=float)
 return np.where(z>=0,1/(1+np.exp(-np.clip(z,-700,700))),np.exp(np.clip(z,-700,700))/(1+np.exp(np.clip(z,-700,700))))

def geometry(Y,p,phase=None):
 y=np.asarray(Y,float)
 r1=np.exp(y[...,0]);r2=np.exp(y[...,4]);x=sigmoid(y[...,8]);g=x*(1-x)
 return (r1,r2,x,g,r1*y[...,1],r2*y[...,5],g*y[...,9],y[...,2]+y[...,6]-p.phi0-(y[...,10] if phase is None else np.asarray(phase,float)))

def U(r,e):
 return .2*r*r*((r-2)**2+1)+e*e/(2*r*r)

def dU(r,e):
 return 2*r-2.4*r*r+.8*r*r*r-e*e/(r*r*r)

def energy_terms(Y,k,p,phase=None):
 y=np.asarray(Y,float)
 r1,r2,x,g,dr1,dr2,dx,D=geometry(y,p,phase=phase)
 T_rad=.5*(dr1**2+dr2**2)
 T_ang=.5*(r1*r1*y[...,3]**2+r2*r2*y[...,7]**2)
 T_x=.5*p.m_x*dx*dx
 U_rad=U(r1,p.e_reserve)+U(r2,p.e_reserve)
 V_bar=-p.barrier_B*np.log(x*(1-x))
 V_c=k*x*r1*r2*np.cos(D)
 E=T_rad+T_ang+T_x+U_rad+V_bar+V_c
 J1=r1*r1*y[...,3];J2=r2*r2*y[...,7];Jdiff=J1-J2;Jsum=J1+J2
 # Ordinary energy identity: E' = -D_rho-D_x + W_x + W_Theta.
 D_rho=p.gamma_rho*(dr1**2+dr2**2)
 D_x=p.gamma_x*dx**2
 W_x=r1*r2*(p.mu_x*np.sin(D)+k*np.cos(D))*dx
 W_Theta=k*x*x*r1*r2*np.sin(D)
 E_dot=-D_rho-D_x+W_x+W_Theta
 # Legendre transformed energy for frame theta_i-Theta/2: H=E - x*(J1+J2)/2.
 H=E-.5*x*Jsum
 H_dot=-D_rho-D_x+W_x-.5*dx*Jsum
 rhoR=np.hypot(r1,r2)
 return locals()

def exact_energy_derivative(Y,F,k,p):
 """Independently evaluate dE/dt = grad(E) . F, any physical-time F."""
 y=np.asarray(Y,float); f=np.asarray(F,float)
 r1,r2,x,g,dr1,dr2,dx,D=geometry(y,p)
 # Differentiated coordinates (NOT the original-model RHS inserted into the chain rule).
 rdr1=r1*f[...,0]; rdr2=r2*f[...,4]
 ddr1=r1*(f[...,0]*y[...,1]+f[...,1]); ddr2=r2*(f[...,4]*y[...,5]+f[...,5])
 dxd=g*f[...,8]
 ddxd=g*(f[...,9]+(1-2*x)*y[...,9]*f[...,8])
 dD=f[...,2]+f[...,6]-f[...,10]
 val=dr1*ddr1+dr2*ddr2
 val+=r1*rdr1*y[...,3]**2+(r1*r1)*y[...,3]*f[...,3]
 val+=r2*rdr2*y[...,7]**2+(r2*r2)*y[...,7]*f[...,7]
 val+=p.m_x*dx*ddxd
 val+=dU(r1,p.e_reserve)*rdr1+dU(r2,p.e_reserve)*rdr2
 val+=(-p.barrier_B*(1/x-1/(1-x)))*dxd
 val+=k*(dxd*r1*r2*np.cos(D)+x*(rdr1*r2+r1*rdr2)*np.cos(D)-x*r1*r2*np.sin(D)*dD)
 return val

def exact_H_derivative(Y,F,k,p):
 y=np.asarray(Y,float);f=np.asarray(F,float)
 r1,r2,x,g,*_=geometry(y,p)
 dx=g*y[...,9];dxd=g*f[...,8]
 Jsum=r1*r1*y[...,3]+r2*r2*y[...,7]
 dJsum=r1*r1*(2*y[...,3]*f[...,0]+f[...,3])+r2*r2*(2*y[...,7]*f[...,4]+f[...,7])
 return exact_energy_derivative(y,f,k,p)-.5*(dxd*Jsum+x*dJsum)

def a_priori_bound(k,p):
 """Elementary conservative comparison: E+86.4 >= all listed positive terms.
 r^2-.8r^3+.2r^4 = .1r^4+r^2+[.1r^4-.8r^3];
 square bracket >= -43.2 at r=6. Assumes x in (0,1), B>0, m_x>0,
 gamma_rho,gamma_x>=0, |k|<2 (hence quadratic coefficient positive).
 """
 c=1-abs(k)/2
 if c<=0 or p.e_reserve<=0 or p.barrier_B<=0 or p.m_x<=0 or p.gamma_rho<0 or p.gamma_x<0:
  return dict(status='NOT_APPLICABLE',reason='parameter conditions failed')
 A=math.hypot(p.mu_x,k)
 # r1^2 r2^2 <= (r1^4+r2^4)/2 <= 5 G; dx^2 <= 2G/m_x
 # |dx*r1*r2| <= .5*(dx^2+r1^2r2^2) <= (2.5+1/m_x)*G
 K=A*(2.5+1/p.m_x)+abs(k)/(2*c)
 return dict(status='FINITE_TIME_GRONWALL_BOUND',energy_shift=86.4,quadratic_coefficient=c,
             K=K,source_norm=A,proof='G(t)<=G(0)*exp(K*(t-t0)), G=E+86.4',
             radial_bound='r_i(t)>=e_reserve/sqrt(2*G(0)*exp(K*t))',
             charge_bound='r1^2+r2^2>=Jdiff^2/(2*G(0)*exp(K*t))',
             uniform_all_time_bound_proved=False)
