#!/usr/bin/env python3
"""10D, SVD-gauge-free orientation of the integrated bivector residual.

All exact derivative identities below concern the chosen C(xi) interpolant.
"""
from __future__ import annotations
import numpy as np
from scipy.optimize import minimize_scalar

I,J=np.triu_indices(10,k=1)

def wedge(a,b):
    a=np.asarray(a,float);b=np.asarray(b,float)
    return a[I]*b[J]-a[J]*b[I]

def cos(a,b):
    a=np.asarray(a,float);b=np.asarray(b,float)
    den=float(np.linalg.norm(a)*np.linalg.norm(b))
    return float(np.dot(a,b)/den) if den>1e-25 else float('nan')

def angle_deg(a,b):
    c=cos(a,b)
    return float(np.rad2deg(np.arccos(np.clip(c,-1,1)))) if np.isfinite(c) else float('nan')

def net_geometry(o,w,at,slope,wdot=0.):
    o=float(o);w=float(w);e=o+3*w
    net=np.asarray(at(e)-at(o),float)
    r=float(np.linalg.norm(net))
    g0=np.asarray(slope(o),float);g1=np.asarray(slope(e),float)
    jo=g1-g0;jw=3*g1
    dnet=jo+float(wdot)*jw
    if r<=1e-15:
        raise ValueError('UNDEFINED_DIRECTION_AT_NEAR_ZERO_RESIDUAL')
    u=net/r
    radial=float(u@dnet)
    du=(dnet-u*radial)/r
    form=wedge(u,du)
    du_norm=float(np.linalg.norm(du))
    if abs(np.linalg.norm(form)-du_norm)>2e-10*max(1.,du_norm):
        raise AssertionError('WEDGE_NORM_IDENTITY_FAILED')
    return dict(o=o,w=w,e=e,net=net,r=r,u=u,g0=g0,g1=g1,jo=jo,jw=jw,
        dnet=dnet,radial=radial,du=du,ang_speed=du_norm,form=form,
        start_radial=float(-u@g0),end_radial=float(u@g1*(1+3*wdot)),
        end_source_radial_fixed=float(u@g1),start_source_radial_fixed=float(-u@g0),
        g0_cos=cos(net,g0),g1_cos=cos(net,g1),
        end_transverse_norm=float(np.linalg.norm(g1-u*(u@g1))),
        start_transverse_norm=float(np.linalg.norm(g0-u*(u@g0))))

def optimize_local_width(o,at,base_w,lo,hi):
    lo=max(float(lo),float(base_w)-.05)
    hi=min(float(hi),float(base_w)+.05)
    if hi<=lo:raise ValueError('INVALID_LOCAL_WIDTH_BAND')
    def objective(w):
        n=at(float(o)+3*float(w))-at(float(o))
        return float(n@n)
    x=minimize_scalar(objective,bounds=(lo,hi),method='bounded',
        options={'xatol':1e-12,'maxiter':300})
    if not x.success:raise RuntimeError('WIDTH_OPTIMIZATION_FAILED')
    return float(x.x),bool(min(x.x-lo,hi-x.x)<5e-5)

def record(g,kind,index,ref,xi,source=None):
    f=g['form']
    record=dict(curve=kind,index=int(index),origin_xi=g['o'],width_xi=g['w'],end_xi=g['e'],
        norm_10d=g['r'],radial_derivative=g['radial'],directional_angular_speed_per_xi=g['ang_speed'],
        direction_angle_from_minimum_deg=angle_deg(g['u'],ref['u']),
        turning_form_cos_to_minimum=cos(f,ref['form']),
        start_source_radial=g['start_radial'],end_source_radial=g['end_radial'],
        boundary_radial_sum_abs_error=abs(g['radial']-(g['start_radial']+g['end_radial'])),
        signed_unit_dot_endpoint_source_start=g['g0_cos'],
        signed_unit_dot_endpoint_source_end=g['g1_cos'],
        start_source_transverse_norm=g['start_transverse_norm'],
        end_source_transverse_norm=g['end_transverse_norm'],
        net_10d_components=' '.join(f'{v:.15g}' for v in g['net']),
        unit_10d_components=' '.join(f'{v:.15g}' for v in g['u']),
        unit_tangent_10d_components=' '.join(f'{v:.15g}' for v in g['du']),
        oriented_turning_45d_components=' '.join(f'{v:.15g}' for v in g['form']))
    if source is not None:
        src0,src1=source(g['o']),source(g['e'])
        record['start_physical_source_interpolation_defect']=float(np.linalg.norm(g['g0']-src0))
        record['end_physical_source_interpolation_defect']=float(np.linalg.norm(g['g1']-src1))
        record['relative_two_source_defect']=float(np.sqrt(np.linalg.norm(g['g0']-src0)**2+np.linalg.norm(g['g1']-src1)**2)/max(np.sqrt(np.linalg.norm(src0)**2+np.linalg.norm(src1)**2),1e-30))
    return record

def derive_curves(xi,at,slope,optimal_origin,optimal_width,halfwidth=0.35,npoints=141,width_range=(12.7,12.91),physical_source=None):
    origins=np.linspace(float(optimal_origin)-halfwidth,float(optimal_origin)+halfwidth,int(npoints))
    if len(origins)<7 or len(origins)%2!=1:raise ValueError('REQUIRE_ODD_CURVE_POINTS')
    origin0=float(optimal_origin);w0=float(optimal_width)
    if origins[0]<xi[0] or origins[-1]+3*(w0+.05)>xi[-1]:
        raise ValueError('CURVE_ENDPOINT_OUTSIDE_TRACE')
    ref=net_geometry(origin0,w0,at,slope)
    fixed=[net_geometry(o,w0,at,slope) for o in origins]
    opt=[optimize_local_width(o,at,w0,*width_range) for o in origins]
    ww=np.array([t[0] for t in opt])
    wdot=np.gradient(ww,origins,edge_order=2)
    valley=[net_geometry(o,w,at,slope,wp) for o,w,wp in zip(origins,ww,wdot)]
    fixed_rows=[record(g,'fixed_width',i,ref,xi,physical_source) for i,g in enumerate(fixed)]
    valley_rows=[record(g,'valley_follow',i,ref,xi,physical_source) for i,g in enumerate(valley)]
    for i,(f,v) in enumerate(zip(fixed_rows,valley_rows)):
        f.update(width_slope=0.0,width_on_bound=False)
        v.update(width_slope=float(wdot[i]),width_on_bound=bool(opt[i][1]))
    return dict(fixed=fixed,valley=valley,fixed_rows=fixed_rows,valley_rows=valley_rows,
        fixed_unit=np.stack([g['u'] for g in fixed]),valley_unit=np.stack([g['u'] for g in valley]),
        fixed_forms=np.stack([g['form'] for g in fixed]),valley_forms=np.stack([g['form'] for g in valley]),
        origins=origins,widths=ww,ref=ref,optimization_bounds_count=sum(t[1] for t in opt))

def curve_stats(curve,kind):
    gs=curve['fixed'] if kind=='fixed_width' else curve['valley']
    origins=curve['origins'];center=len(gs)//2
    angular=np.array([g['ang_speed'] for g in gs]);norm=np.array([g['r'] for g in gs]);
    u=np.stack([g['u'] for g in gs]);form=np.stack([g['form'] for g in gs]);
    f0=gs[center]['form'];l=center-6;r=center+6
    integral=float(np.trapz(angular,origins)) if hasattr(np,'trapz') else float(np.trapezoid(angular,origins))
    return dict(curve=kind,n_samples=len(gs),min_norm=float(np.min(norm)),max_norm=float(np.max(norm)),
        minimum_sample_origin_xi=float(origins[np.argmin(norm)]),
        reference_norm=float(curve['ref']['r']),
        reference_angular_speed=float(angular[center]),
        orientation_path_length_rad=integral,
        endpoint_unit_angle_rad=float(np.arccos(np.clip(u[0]@u[-1],-1,1))),
        normalized_endpoint_turn_cos=cos(form[l],form[r]),
        min_turning_form_10d_norm=float(np.min(np.linalg.norm(form,axis=1))),
        max_turning_form_10d_norm=float(np.max(np.linalg.norm(form,axis=1))),
        num_width_bounds=int(curve['optimization_bounds_count'] if kind=='valley' else 0))

def fd_derivative_check(o,w,at,slope,h,xi):
    a=net_geometry(o,w,at,slope)
    knots=np.asarray(xi)
    d0=min(np.min(abs(knots-o)),np.min(abs(knots-(o+3*w))))
    if d0<h*1.5:
        return dict(evaluated=False,knot_clearance=float(d0),max_abs_defect=None,relative_defect=None)
    netp=at(o+h+3*w)-at(o+h)
    netm=at(o-h+3*w)-at(o-h)
    diff=(netp-netm)/(2*h)
    absd=float(np.linalg.norm(a['jo']-diff))
    return dict(evaluated=True,knot_clearance=float(d0),max_abs_defect=absd,
        relative_defect=absd/max(float(np.linalg.norm(a['jo'])),1e-30))
