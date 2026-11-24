#!/usr/bin/env python3
"""Exact Q550 signed-phase energetics; no IVP solved.

All formulas refer to physical material time tau, not support arclength.
The global phase Theta(t) is external to the recorded 10D state.
"""
import numpy as np
from energy_math import energy_terms, exact_energy_derivative, Params, sigmoid


def signed_power(y,k,p,phase,phase_dot):
    """True E balance when external global phase has ANY differentiable rate."""
    d=energy_terms(y,k,p,phase=phase)
    w_phase=k*d['x']*d['r1']*d['r2']*np.sin(d['D'])*np.asarray(phase_dot,float)
    return -d['D_rho']-d['D_x']+d['W_x']+w_phase


def frame_quantity(y,k,p,phase,sigma):
    d=energy_terms(y,k,p,phase=phase)
    return d['E']-sigma*.5*d['x']*d['Jsum']


def frame_ideal_power(y,k,p,phase,sigma):
    d=energy_terms(y,k,p,phase=phase)
    return -d['D_rho']-d['D_x']+d['W_x']-sigma*.5*d['dx']*d['Jsum']


def phase_work_per_unit(y,k,p,phase):
    d=energy_terms(y,k,p,phase=phase)
    return k*d['x']*d['r1']*d['r2']*np.sin(d['D'])


def physical_rhs_with_signed_external_phase(y,k,p,phase,sigma=None,rate=None):
    """Independent archived mechanical 11D RHS for a prescribed external Theta.

    sigma selects Theta_dot = sigma*x; alternatively rate sets arbitrary Theta_dot.
    Not used to generate any physical orbit.
    """
    y=np.asarray(y,float)
    r1,r2=np.exp(y[0]),np.exp(y[4]);v1,v2=y[1],y[5]
    o1,o2=y[3],y[7];x=float(sigmoid(y[8]));g=x*(1-x)
    dth=y[2]+y[6]-p.phi0-float(phase)
    sd,cd=np.sin(dth),np.cos(dth);ke=k*x
    dr1=r1*v1;dr2=r2*v2;dx=g*y[9]
    a1=r1*o1**2-(2*r1-2.4*r1**2+.8*r1**3-p.e_reserve**2/r1**3)-ke*r2*cd-p.gamma_rho*dr1
    a2=r2*o2**2-(2*r2-2.4*r2**2+.8*r2**3-p.e_reserve**2/r2**3)-ke*r1*cd-p.gamma_rho*dr2
    dx2=(-(-p.barrier_B*(1/x-1/(1-x)))+p.mu_x*r1*r2*sd-p.gamma_x*dx)/p.m_x
    if rate is None:
        if sigma not in (-1,1):raise ValueError('Provide sign or explicit phase rate')
        rate=sigma*x
    return np.array([v1,a1/r1-v1**2,o1,ke*(r2/r1)*sd-2*v1*o1,
                     v2,a2/r2-v2**2,o2,ke*(r1/r2)*sd-2*v2*o2,
                     y[9],dx2/g-(1-2*x)*y[9]**2,float(rate)],float)


def exact_identity_probe(y,k,p,phase,rate):
    F=physical_rhs_with_signed_external_phase(y,k,p,phase,rate=rate)
    # The last entry of state is phase, not material time.
    y11=np.asarray(y,float).copy();y11[10]=phase
    chain=float(exact_energy_derivative(y11,F,k,p))
    balance=float(signed_power(y11,k,p,phase,rate))
    return chain-balance


def integral(t,v):
    t=np.asarray(t,float);v=np.asarray(v,float)
    if len(t)<3 or not np.all(np.diff(t)>0):raise ValueError('Non-monotone or short time grid')
    return np.r_[0.,np.cumsum(.5*(v[1:]+v[:-1])*np.diff(t))]


def summary_distribution(x):
    x=np.asarray(x,float)
    return {'min':float(np.min(x)),'max':float(np.max(x)),
            'rms':float(np.sqrt(np.mean(x*x))), 'median':float(np.median(x))}
