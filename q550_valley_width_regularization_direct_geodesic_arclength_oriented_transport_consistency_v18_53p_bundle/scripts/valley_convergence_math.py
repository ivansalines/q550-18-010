#!/usr/bin/env python3
"""Numerical geometry of an authenticated 10D cumulative bivector primitive.
No physical IVPs or SVD observer transport are computed here.
"""
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import minimize_scalar
I,J=np.triu_indices(10,1)

class Primitive:
    def __init__(self,xi,c,mode):
        self.x=np.asarray(xi,float); self.c=np.asarray(c,float);self.mode=mode
        assert self.c.shape==(len(self.x),10) and np.all(np.diff(self.x)>0)
        self.g=np.diff(self.c,axis=0)/np.diff(self.x)[:,None]
        self.p=PchipInterpolator(self.x,self.c,axis=0,extrapolate=False) if mode=='pchip' else None
    def idx(self,t):
        t=np.asarray(t,float)
        if np.any(t<self.x[0]-1e-10) or np.any(t>self.x[-1]+1e-10):raise ValueError('OUTSIDE_AUTHENTIC_TRACE')
        return np.clip(np.searchsorted(self.x,t,side='right')-1,0,len(self.x)-2)
    def at(self,t):
        t=np.asarray(t,float)
        if self.p is not None:return np.asarray(self.p(t),float)
        k=self.idx(t)
        return self.c[k]+np.expand_dims(t-self.x[k],-1)*self.g[k]
    def slope(self,t):
        if self.p is not None:return np.asarray(self.p.derivative()(t),float)
        return self.g[self.idx(t)]
    def clearance(self,t):
        t=np.asarray(t,float);k=self.idx(t)
        return np.minimum(abs(t-self.x[k]),abs(t-self.x[k+1]))

def optimum_width(p,o,lo,hi):
    """PWL: exact global segment-wise optimum. PCHIP: bracketed multi-minimum search."""
    co=p.at(o); emin=o+3*lo;emax=o+3*hi
    if p.mode=='linear':
        best=(np.inf,np.nan)
        for k in range(int(p.idx(emin)),int(p.idx(emax))+1):
            a=max(emin,p.x[k]);b=min(emax,p.x[k+1]);g=p.g[k]
            if b<a:continue
            v=p.c[k]+g*(a-p.x[k])-co
            step=np.clip(-float(v@g)/max(float(g@g),1e-200),0,b-a)
            vv=v+step*g;score=float(vv@vv)
            if score<best[0]:best=(score,a+step)
        return (best[1]-o)/3.,np.sqrt(best[0])
    def obj(w):
        v=p.at(o+3*w)-co
        return float(v@v)
    grid=np.linspace(lo,hi,17);vals=[obj(w) for w in grid]
    cand=[(vals[0],lo),(vals[-1],hi)]
    for k in range(1,16):
        if vals[k]<=vals[k-1] and vals[k]<=vals[k+1]:
            z=minimize_scalar(obj,bounds=(grid[k-1],grid[k+1]),method='bounded',options={'xatol':1e-12})
            if not z.success:raise RuntimeError('PCHIP_OPTIMIZER_FAILED')
            cand.append((z.fun,z.x))
    v,w=min(cand)
    return float(w),float(np.sqrt(v))

def wedge(a,b):return a[...,I]*b[...,J]-a[...,J]*b[...,I]

def arcs(u):
    return 2*np.arctan2(np.linalg.norm(np.diff(u,axis=0),axis=1),np.linalg.norm(u[1:]+u[:-1],axis=1))

def area_from_polygon(u,theta):
    # theta/sin(theta) = 1/sinc(theta/pi), with the correct limit at zero.
    fac=1/np.sinc(theta/np.pi)
    return np.sum(wedge(u[:-1],u[1:])*fac[:,None],axis=0)

def eval_curve(p,o,w):
    o=np.asarray(o);w=np.asarray(w);end=o+3*w
    n=p.at(end)-p.at(o);r=np.linalg.norm(n,axis=1)
    if np.min(r)<1e-12:raise RuntimeError('UNDEFINED_DIRECTION')
    u=n/r[:,None];wd=np.gradient(w,o,edge_order=2)
    wd_p=PchipInterpolator(o,w).derivative()(o)
    gs=p.slope(o);ge=p.slope(end)
    def derive(v):
        dn=ge*(1+3*v[:,None])-gs
        radial=np.einsum('ij,ij->i',u,dn)
        du=(dn-u*radial[:,None])/r[:,None]
        return np.linalg.norm(du,axis=1),wedge(u,du),radial
    speed,form,radial=derive(wd)
    speed_p,form_p,radial_p=derive(wd_p)
    theta=arcs(u);polygon=area_from_polygon(u,theta)
    form_trap=np.trapezoid(form,o,axis=0)
    form_p_trap=np.trapezoid(form_p,o,axis=0)
    form_norm=np.linalg.norm(form,axis=1)
    inner_center=form@form[len(o)//2]/np.maximum(form_norm*form_norm[len(o)//2],1e-100)
    inner_adjacent=np.einsum('ij,ij->i',form[1:],form[:-1])/np.maximum(form_norm[1:]*form_norm[:-1],1e-100)
    area_cos=float(polygon@form_trap/max(float(np.linalg.norm(polygon)*np.linalg.norm(form_trap)),1e-100))
    fd_unit=np.gradient(u,o,axis=0,edge_order=2)
    endpoint_angle=arcs(np.stack([u[0],u[-1]]))[0]
    wd2=np.gradient(wd,o,edge_order=2)
    stats=dict(
      arc_geodesic=float(sum(theta)),arc_chain_gradient=float(np.trapezoid(speed,o)),
      arc_chain_pchip_width=float(np.trapezoid(speed_p,o)),
      arc_fd_unit=float(np.trapezoid(np.linalg.norm(fd_unit,axis=1),o)),
      arc_geodesic_minus_chain=float(sum(theta)-np.trapezoid(speed,o)),
      arc_geodesic_minus_pchip_width=float(sum(theta)-np.trapezoid(speed_p,o)),
      geodesic_endpoint=float(endpoint_angle),max_arc_step=float(np.max(theta)),
      center_speed=float(speed[len(o)//2]),min_norm=float(np.min(r)),
      center_norm=float(r[len(o)//2]),
      area_polygon_norm=float(np.linalg.norm(polygon)),
      area_polygon_chain_form_cos=area_cos,
      min_oriented_turn_form_cos_to_center=float(np.min(inner_center)),
      min_adjacent_turn_form_cos=float(np.min(inner_adjacent)),
      num_turn_form_direction_flips=int(np.sum(inner_adjacent<0)),
      area_polygon_minus_chain=float(np.linalg.norm(polygon-form_trap)),
      area_polygon_minus_pchip_width=float(np.linalg.norm(polygon-form_p_trap)),
      radial_fundamental_theorem_defect=float(np.trapezoid(radial,o)-(r[-1]-r[0])),
      width_range=float(np.ptp(w)),width_total_variation=float(np.sum(np.abs(np.diff(w)))),
      max_width_step=float(np.max(abs(np.diff(w)))),
      max_abs_wprime=float(np.max(abs(wd))),
      max_abs_wdoubleprime=float(np.max(abs(wd2))),
      width_derivative_sign_changes=int(np.sum(wd[1:]*wd[:-1]<0)),
      n_origin_on_knots=int(np.sum(p.clearance(o)<1e-9)),
      n_end_on_knots=int(np.sum(p.clearance(end)<1e-9)),
      closest_end_knot=float(np.min(p.clearance(end))),
      closest_origin_knot=float(np.min(p.clearance(o))),
      min_adjacent_direction_cos=float(np.min(np.einsum('ij,ij->i',u[1:],u[:-1]))))
    series=dict(origin=o,width=w,end=end,r=r,u=u,speed=speed,radial=radial,
                wdot=wd,wdot_pchip=wd_p,arc_steps=theta,polygon_area=polygon,
                chain_area=form_trap)
    return stats,series
