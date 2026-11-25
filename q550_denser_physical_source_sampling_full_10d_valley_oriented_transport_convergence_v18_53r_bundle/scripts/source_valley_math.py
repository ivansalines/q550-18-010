#!/usr/bin/env python3
"""Source-constrained post-processing geometry; no Q550 physical integration."""
import numpy as np
from scipy.interpolate import PchipInterpolator, CubicHermiteSpline
from scipy.optimize import brentq
I,J=np.triu_indices(10,1)

def wedge(a,b): return a[...,I]*b[...,J]-a[...,J]*b[...,I]

def angular_steps(u):
    return 2*np.arctan2(np.linalg.norm(np.diff(u,axis=0),axis=1),np.linalg.norm(u[1:]+u[:-1],axis=1))

class Primitive:
    def __init__(self,xi,c,source,kind):
        self.x=np.asarray(xi,float)
        self.c=np.asarray(c,float)
        self.source=np.asarray(source,float)
        self.kind=kind
        if self.c.shape!=(len(xi),10) or self.source.shape!=self.c.shape or not np.all(np.diff(self.x)>0):
            raise ValueError('BAD_SOURCE_PRIMITIVE_SHAPE')
        if kind=='pchip': self.f=PchipInterpolator(self.x,self.c,axis=0,extrapolate=False)
        elif kind=='hermite': self.f=CubicHermiteSpline(self.x,self.c,self.source,axis=0,extrapolate=False)
        elif kind=='source_integral':
            dx=np.diff(self.x)
            self.source_c=np.vstack([self.c[0],self.c[0]+np.cumsum((self.source[:-1]+self.source[1:])*.5*dx[:,None],axis=0)])
            self.slope=(self.source[1:]-self.source[:-1])/dx[:,None]
            self.f=None
        else: raise ValueError('BAD_PRIMITIVE_MODE')
    def indices(self,t):
        t=np.asarray(t,float)
        if not np.all(np.isfinite(t)) or np.any(t<self.x[0]-1e-10) or np.any(t>self.x[-1]+1e-10):
            raise ValueError('EVALUATION_OUTSIDE_AUTHENTIC_TRACE')
        return np.clip(np.searchsorted(self.x,t,side='right')-1,0,len(self.x)-2)
    def eval(self,t,order=0):
        t=np.asarray(t,float)
        k=self.indices(t)
        if self.f is not None: return np.asarray(self.f.derivative(order)(t),float)
        a=np.expand_dims(t-self.x[k],-1)
        if order==0:return self.source_c[k]+self.source[k]*a+.5*self.slope[k]*a*a
        if order==1:return self.source[k]+self.slope[k]*a
        if order==2:return self.slope[k]
        return np.zeros_like(self.source[k])
    def knot_clearance(self,t):
        t=np.asarray(t,float);k=self.indices(t)
        return np.minimum(np.abs(t-self.x[k]),np.abs(t-self.x[k+1]))
    def node_mismatch(self):
        if self.kind=='source_integral':return np.linalg.norm(self.source_c-self.c,axis=1)
        return np.zeros(len(self.x))
    def node_source_derivative_mismatch(self):
        return np.linalg.norm(self.eval(self.x,1)-self.source,axis=1)

def obj(p,o,w):
    e=o+3*w
    N=p.eval(e)-p.eval(o)
    g=p.eval(e,1)
    return float(.5*(N@N)),float(N@g)

def optimized_width(p,o,wlo,whi,nsearch=25):
    """Search all sampled derivative sign brackets; compare all roots and endpoints."""
    candidates=[]
    def add(w,label):
        f,s=obj(p,o,w);candidates.append((f,w,label,s))
    add(wlo,'lower_bound');add(whi,'upper_bound')
    grid=np.linspace(wlo,whi,nsearch)
    station=np.array([obj(p,o,float(w))[1] for w in grid])
    for j in range(len(grid)-1):
        a=float(grid[j]);b=float(grid[j+1]);sa=station[j];sb=station[j+1]
        if abs(sa)<1e-14: add(a,'stationary_sample')
        if sa*sb<0:
            w=brentq(lambda x:obj(p,o,x)[1],a,b,xtol=1e-12,rtol=1e-14)
            add(float(w),'stationary_root')
    if abs(station[-1])<1e-14:add(float(grid[-1]),'stationary_sample')
    # Extra gate: stationary candidates can be missed if two roots hide in a subinterval.
    # This is a stated grid-bracketing limitation, not a global-continuum certificate.
    candidates=sorted(candidates,key=lambda z:z[0])
    f,w,label,s=candidates[0]
    N=p.eval(o+3*w)-p.eval(o);ge=p.eval(o+3*w,1);he=p.eval(o+3*w,2);go=p.eval(o,1)
    curvature=float(9*(ge@ge+N@he))
    denom=3*(ge@ge+N@he)
    numer=(ge-go)@ge+N@he
    if label.startswith('stationary') and curvature>1e-10 and abs(denom)>1e-10:
        wprime=float(-numer/denom)
    else:wprime=float('nan')
    runnerup=float(candidates[1][0]-f) if len(candidates)>1 else float('nan')
    return dict(w=float(w),norm=float(np.sqrt(2*f)),label=label,num_candidates=len(candidates),
                score_margin=runnerup,stationarity=float(s),curvature=curvature,
                implicit_wprime=wprime,bound=int(label.endswith('bound')),
                clearance=float(p.knot_clearance(o+3*w)))

def curve(p,o,w,implicit=None):
    o=np.asarray(o,float);w=np.asarray(w,float);e=o+3*w
    N=p.eval(e)-p.eval(o);r=np.linalg.norm(N,axis=1)
    if np.any(r<=1e-12):raise ValueError('UNDEFINED_10D_UNIT_DIRECTION')
    u=N/r[:,None]
    dw=np.gradient(w,o,edge_order=2)
    if implicit is not None:
        ww=np.asarray(implicit,float)
        dw=np.where(np.isfinite(ww),ww,dw)
    g0=p.eval(o,1);ge=p.eval(e,1)
    dN=ge*(1+3*dw[:,None])-g0
    radial=np.sum(u*dN,axis=1)
    dU=(dN-u*radial[:,None])/r[:,None]
    speed=np.linalg.norm(dU,axis=1)
    angles=angular_steps(u)
    fac=1/np.sinc(angles/np.pi)
    poly=np.sum(wedge(u[:-1],u[1:])*fac[:,None],axis=0)
    tang=np.trapezoid(wedge(u,dU),o,axis=0)
    direct=float(np.sum(angles));chain=float(np.trapezoid(speed,o))
    if len(u)<3:raise ValueError('TOO_FEW_GRID_NODES')
    fd=np.gradient(u,o,axis=0,edge_order=2)
    return dict(direct=direct,chain=chain,arc_diff=float(direct-chain),
                fd_length=float(np.trapezoid(np.linalg.norm(fd,axis=1),o)),
                oriented_cos=float(poly@tang/max(np.linalg.norm(poly)*np.linalg.norm(tang),1e-100)),
                oriented_area_difference=float(np.linalg.norm(poly-tang)),
                oriented_polygon=poly,oriented_tangent=tang,
                width_range=float(np.ptp(w)),width_total_variation=float(np.sum(np.abs(np.diff(w)))),
                width_sign_changes=int(np.sum(np.diff(w)[1:]*np.diff(w)[:-1]<0)),
                max_angular_step=float(np.max(angles)),min_norm=float(np.min(r)),
                max_speed=float(np.max(speed)),directional_speed=speed,min_endpoint_knot_clearance=float(np.min(p.knot_clearance(e))),
                u=u,norm=r,wd=dw)
