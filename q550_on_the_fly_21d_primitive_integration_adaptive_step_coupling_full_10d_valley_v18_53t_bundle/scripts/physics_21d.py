#!/usr/bin/env python3
"""Exact functional extension of the authenticated Q550 11D RHS to 11D+10D.

The last ten components are a passive integral. They never enter the physical
11D RHS. The instantaneous 10D source is the same radial-potential plus
phase-coupling bivector used by v18.53r/s, in physical xi clock.
"""
from __future__ import annotations
import numpy as np


def source_10d(h, m, context, kappa, xi, z11, tau_rate=None):
    """Return g(xi,z) as 10D Euclidean skew coordinates; no quadrature."""
    z11=np.asarray(z11,float)
    if z11.shape!=(11,):raise ValueError('EXPECTED_11D_STATE')
    phase=float(context['theta'](context['event']+float(xi)))
    q,_,_,_,_,forces=h.force_components(z11[:10][None,:],np.array([phase]),kappa,h.P)
    f=forces['radial_potential'][0]+forces['phase_coupling'][0]
    biv=h.bivector(q[0],f)
    if tau_rate is None:
        raise ValueError('PASS_AUTHENTIC_MATERIAL_CLOCK_RATE_EXPLICITLY')
    answer=float(tau_rate)*m.to10(biv)
    if answer.shape!=(10,) or not np.all(np.isfinite(answer)):
        raise RuntimeError('BAD_10D_SOURCE')
    return answer


def source_array(h,m,context,kappa,xi,z11,omega):
    """Batch evaluation, same literal force components as source_10d."""
    xi=np.asarray(xi,float);z11=np.asarray(z11,float)
    if z11.shape!=(len(xi),11):raise ValueError('INVALID_PHYSICAL_TRACE')
    from scipy.special import expit
    phase=np.array([context['theta'](context['event']+float(a)) for a in xi])
    q,_,_,_,_,force=h.force_components(z11[:,:10],phase,kappa,h.P)
    ff=force['radial_potential']+force['phase_coupling']
    B=h.bivector(q,ff)
    # Keep the authenticated conversion rather than assume a public PAIRS table.
    ret=np.stack([m.to10(mat) for mat in B],axis=0)
    ret *= (omega/expit(z11[:,8]))[:,None]
    if ret.shape!=(len(xi),10) or not np.all(np.isfinite(ret)):
        raise RuntimeError('BAD_SOURCE_TRACE')
    return ret


def rhs_21d(h,m,engine,context,kappa,xi,y21):
    y21=np.asarray(y21,float)
    if y21.shape!=(21,):raise ValueError('EXPECTED_21D_STATE')
    dz=h.rhs(engine,context,kappa,xi,y21[:11])
    dg=source_10d(h,m,context,kappa,xi,y21[:11],tau_rate=dz[10])
    return np.concatenate((dz,dg))


def angular_polygon(u, wedge):
    """Direct, metric-defined arclength and oriented polygon integral."""
    u=np.asarray(u,float)
    norm=np.linalg.norm(u,axis=1)
    if u.ndim!=2 or u.shape[1]!=10 or np.any(norm<1e-12):
        raise ValueError('UNDEFINED_DIRECTION')
    u=u/norm[:,None]
    # Stable on small adjacent angular steps.
    angles=2*np.arctan2(np.linalg.norm(u[1:]-u[:-1],axis=1),
                        np.linalg.norm(u[1:]+u[:-1],axis=1))
    factor=np.where(angles<1e-6,1+angles*angles/6,angles/np.maximum(np.sin(angles),1e-30))
    poly=np.sum(wedge(u[:-1],u[1:])*factor[:,None],axis=0)
    return float(np.sum(angles)),poly
