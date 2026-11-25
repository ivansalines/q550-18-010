#!/usr/bin/env python3
"""Quadrature operations on authenticated, nested Q550 xi/source samples.

Independent of model integration. All input arrays are numerical data, never treated
as an oracle for the exact continuous-time source integral.
"""
from __future__ import annotations
import numpy as np
from scipy.interpolate import PchipInterpolator


def assert_history(x,g):
    x=np.asarray(x,dtype=float)
    g=np.asarray(g,dtype=float)
    if x.ndim!=1 or g.ndim!=2 or g.shape[0]!=len(x) or len(x)<5:
        raise ValueError('INVALID_HISTORY_SHAPE')
    if np.any(~np.isfinite(x)) or np.any(~np.isfinite(g)) or np.any(np.diff(x)<=0):
        raise ValueError('INVALID_HISTORY_GRID')
    return x,g


def cumulative_trapezoid(x,g,c0=None):
    x,g=assert_history(x,g)
    c0=np.zeros(g.shape[1]) if c0 is None else np.asarray(c0,dtype=float)
    if c0.shape!=(g.shape[1],):raise ValueError('INVALID_C0')
    panels=(g[:-1]+g[1:])*(np.diff(x)[:,None]*0.5)
    return np.vstack([c0,c0+np.cumsum(panels,axis=0)])


def baseline_quadrature(x,g,block,method='trapezoid',c0=None):
    """Cumulative integral at coarse endpoints, panel grouping exactly nested.

    Composite Simpson requires even block size and equal subpanel widths in each
    Simpson pair (not globally uniform coarse knots).
    """
    x,g=assert_history(x,g)
    block=int(block)
    if block<1 or (len(x)-1)%block:raise ValueError('NOT_NESTED')
    n=(len(x)-1)//block
    c0=np.zeros(g.shape[1]) if c0 is None else np.asarray(c0,dtype=float)
    if c0.shape!=(g.shape[1],):raise ValueError('INVALID_C0')
    if method=='trapezoid':
        p=(g[:-1]+g[1:])*(np.diff(x)[:,None]*.5)
        panels=p.reshape(n,block,g.shape[1]).sum(axis=1)
    elif method=='simpson':
        if block%2:raise ValueError('SIMPSON_REQUIRES_EVEN_BLOCK')
        k=np.arange(n)[:,None]*block+np.arange(0,block,2)[None,:]
        h1=x[k+1]-x[k]
        h2=x[k+2]-x[k+1]
        if not np.allclose(h1,h2,rtol=0,atol=1e-10):raise ValueError('SIMPSON_REQUIRES_EQUAL_PAIRS')
        panels=((x[k+2]-x[k])[...,None]*(g[k]+4*g[k+1]+g[k+2])/6).sum(axis=1)
    else:raise ValueError('UNKNOWN_QUADRATURE')
    return np.vstack([c0,c0+np.cumsum(panels,axis=0)])


def source_pchip_integral(x,g,c0=None):
    """Integrate a componentwise shape-preserving interpolant of sampled SOURCE.

    C1 source interpolant, not an exact primitive of the continuous model.
    """
    x,g=assert_history(x,g)
    c0=np.zeros(g.shape[1]) if c0 is None else np.asarray(c0,dtype=float)
    q=PchipInterpolator(x,g,axis=0,extrapolate=False)
    anti=q.antiderivative()
    return np.asarray(c0+anti(x)-anti(x[0]),float)


def max_norm(x):
    x=np.asarray(x,float)
    return float(np.max(np.linalg.norm(x,axis=-1)))


def decompose_nested_2x4x(x2,g2,x4,g4,c0=None):
    """Exact telescoping difference for trapezoidal cumulative quadratures.

    T4(g4)-T2(g2) = [T4(g4)-T2(g4 restricted)]
                    + [T2(g4 restricted)-T2(g2)] .
    Left term: sample-density effect on SAME 4x source. Right: changes of
    source at identical 2x knots (from different physical integrations).
    This is an accounting identity, not causal attribution or an error bound.
    """
    x2,g2=assert_history(x2,g2);x4,g4=assert_history(x4,g4)
    if len(x4)!=2*(len(x2)-1)+1 or not np.allclose(x4[::2],x2,atol=1e-12,rtol=0):
        raise ValueError('NON_NESTED_PHYSICAL_GRIDS')
    a=cumulative_trapezoid(x4,g4,c0)[::2]
    b=cumulative_trapezoid(x2,g4[::2],c0)
    c=cumulative_trapezoid(x2,g2,c0)
    residual=a-c-(a-b)-(b-c)
    return dict(total=a-c,quadrature_on_common_g4=a-b,physical_sample_change=b-c,
                identity_max_error=max_norm(residual))


def norm_rows(c):
    return float(np.max(np.linalg.norm(c,axis=1))),float(np.sqrt(np.mean(np.sum(c*c,axis=1))))
