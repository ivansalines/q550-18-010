#!/usr/bin/env python3
"""10D full residual geometry; no SVD observer projection in objective.

For C(xi) the accumulated full bivector and p=(o,w), the exactly telescoped
three-window sum is N=C(o+3w)-C(o). Inside native piecewise-linear source cells,
J=[g_end-g_start,3*g_end]. For f=1/2||N||^2, grad=J^T N and
H=J^T J (only between knots). H=J^T J is NOT a smooth physical Hessian
certificate: the trace interpolant is merely piecewise linear.
"""
from __future__ import annotations
import numpy as np
from scipy.optimize import least_squares, minimize
from scipy.interpolate import PchipInterpolator

def check_box(xi, orange, wrange):
    xi=np.asarray(xi,float)
    assert xi.ndim==1 and len(xi)>3 and np.all(np.isfinite(xi)) and np.all(np.diff(xi)>0)
    lo,hi=map(float,orange); wl,wh=map(float,wrange)
    if not (lo<hi and 0<wl<wh):raise ValueError('INVALID_PARAM_BOX')
    if lo<xi[0]-1e-10 or hi+3*wh>xi[-1]+1e-10:
        raise ValueError(f'ENDPOINT_OUTSIDE_TRACE: box=[{lo},{hi+3*wh}], trace=[{xi[0]},{xi[-1]}]')
    return dict(trace_start=float(xi[0]),trace_end=float(xi[-1]),box_start=lo,box_end=hi+3*wh)

def make_eval(xi, primitive, kind='linear'):
    x=np.asarray(xi,float); c=np.asarray(primitive,float)
    if c.ndim!=2 or c.shape!=(len(x),10):raise ValueError('PRIMITIVE_NOT_10D')
    if not np.all(np.isfinite(c)) or not np.all(np.diff(x)>0):raise ValueError('INVALID_PRIMITIVE')
    if kind=='pchip':
        itp=PchipInterpolator(x,c,axis=0,extrapolate=False)
        def at(t):
            if t<x[0]-1e-10 or t>x[-1]+1e-10:raise ValueError('ENDPOINT_OUTSIDE_TRACE')
            return np.asarray(itp(float(np.clip(t,x[0],x[-1]))))
        der=itp.derivative()
        def slope(t):
            if t<x[0]-1e-10 or t>x[-1]+1e-10:raise ValueError('ENDPOINT_OUTSIDE_TRACE')
            return np.asarray(der(float(np.clip(t,x[0],x[-1]))))
    elif kind=='linear':
        def idx(t):
            if t<x[0]-1e-10 or t>x[-1]+1e-10:raise ValueError('ENDPOINT_OUTSIDE_TRACE')
            k=int(np.searchsorted(x,float(np.clip(t,x[0],x[-1])),side='right')-1)
            return max(0,min(len(x)-2,k))
        def at(t):
            k=idx(t); a=(float(t)-x[k])/(x[k+1]-x[k]);return c[k]*(1-a)+c[k+1]*a
        def slope(t):
            k=idx(t);return (c[k+1]-c[k])/(x[k+1]-x[k])
    else: raise ValueError('UNSUPPORTED_INTERPOLANT')
    return at,slope

def net_jac(p,at,slope):
    o,w=map(float,p);e=o+3*w
    n=at(e)-at(o)
    g0=slope(o);g1=slope(e)
    j=np.column_stack((g1-g0,3*g1))
    return n,j

def local_values(p,at,slope):
    n,j=net_jac(p,at,slope)
    return .5*float(n@n),j.T@n,j,n

def knot_clearance(x,t):
    idx=int(np.searchsorted(x,t))
    return float(min(abs(x[max(0,idx-1)]-t),abs(x[min(len(x)-1,idx)]-t)))

def projected_kkt(grad,p,lb,ub,tol=1e-7):
    # min f, lower bound requires grad>=0, upper bound requires grad<=0
    v=np.asarray(grad,float).copy()
    for i in (0,1):
        if p[i]<=lb[i]+tol: v[i]=min(v[i],0.)
        elif p[i]>=ub[i]-tol:v[i]=max(v[i],0.)
    return float(np.linalg.norm(v))

def starts_from_grid(xi,primitive,origins,widths,orange,wrange,n=16,minsep=.035):
    at,slope=make_eval(xi,primitive)
    scores=[]
    for o in origins:
        for w in widths:
            f,_,_,_=local_values((o,w),at,slope)
            scores.append((f,float(o),float(w)))
    scores.sort()
    selected=[]
    for f,o,w in scores:
        if all(np.hypot((o-a)/(orange[1]-orange[0]),(w-b)/(wrange[1]-wrange[0]))>minsep for a,b in selected):
            selected.append((o,w))
        if len(selected)>=n:break
    for test in [(14.392245,12.79964),(14.4,12.8),(15.2,12.80949),
                 (13.6,12.81336),(12.8,12.8),(15.4,12.8)]:
        if all(np.linalg.norm(np.asarray(test)-np.asarray(q))>1e-5 for q in selected):selected.append(test)
    return selected,scores

def solve_minimum(at,slope,start,orange,wrange):
    lb=np.array([orange[0],wrange[0]],float); ub=np.array([orange[1],wrange[1]],float)
    init=np.clip(np.asarray(start,float),lb+1e-10,ub-1e-10)
    def fun(p):return net_jac(p,at,slope)[0]
    def jac(p):return net_jac(p,at,slope)[1]
    r=least_squares(fun,init,jac=jac,bounds=(lb,ub),method='trf',
       max_nfev=250,ftol=1e-13,xtol=1e-13,gtol=1e-13)
    p=np.asarray(r.x,float)
    f,g,j,net=local_values(p,at,slope)
    sv=np.linalg.svd(j,compute_uv=False)
    scaled=np.linalg.svd(j@np.diag(ub-lb),compute_uv=False)
    return dict(origin_xi=float(p[0]),width_xi=float(p[1]),end_xi=float(p[0]+3*p[1]),
        net_norm=float(np.linalg.norm(net)),objective=f,
        gradient_origin=float(g[0]),gradient_width=float(g[1]),gradient_norm=float(np.linalg.norm(g)),
        constrained_kkt_norm=projected_kkt(g,p,lb,ub),
        jac_sigma_max=float(sv[0]),jac_sigma_min=float(sv[-1]),
        jac_condition=float(sv[0]/max(sv[-1],1e-300)),
        scaled_jac_sigma_min=float(scaled[-1]),
        scaled_jac_condition=float(scaled[0]/max(scaled[-1],1e-300)),
        gauss_newton_hessian_min_eigen=float(scaled[-1]**2),
        on_box_boundary=bool(np.any(np.r_[p-lb,ub-p]<1e-6)),
        optimizer_success=bool(r.success),optimizer_nfev=int(r.nfev),
        initial_origin_xi=float(init[0]),initial_width_xi=float(init[1]),
        optimizer_message=str(r.message))

def unique_minima(rows,dtol=.003,wtol=.0008):
    result=[]
    for r in sorted(rows,key=lambda x:x['net_norm']):
        if not any(abs(r['origin_xi']-p['origin_xi'])<dtol and
                   abs(r['width_xi']-p['width_xi'])<wtol for p in result):result.append(r)
    return result

def endpoint_source_diff(primitive_xi,primitive,B,tau_rate,p):
    """Linear-primitive derivative vs interpolated continuous-source endpoint.
    Disagreement is a quadrature/stencil diagnostic, not a source of physics.
    """
    at,slope=make_eval(primitive_xi,primitive)
    o,w=map(float,p);e=o+3*w
    def source(t):
        k=np.clip(int(np.searchsorted(primitive_xi,t,side='right')-1),0,len(primitive_xi)-2)
        a=np.clip((t-primitive_xi[k])/(primitive_xi[k+1]-primitive_xi[k]),0.,1.)
        return np.asarray(B[k])*(1-a)+np.asarray(B[k+1])*a
    from_rate=lambda t:float(np.interp(t,primitive_xi,tau_rate))
    g0=source(o)*from_rate(o);g1=source(e)*from_rate(e)
    _,J=net_j(p,at,slope)
    diff=np.column_stack((g1-g0,3*g1))
    return float(np.linalg.norm(J-diff)),float(np.linalg.norm(diff))
