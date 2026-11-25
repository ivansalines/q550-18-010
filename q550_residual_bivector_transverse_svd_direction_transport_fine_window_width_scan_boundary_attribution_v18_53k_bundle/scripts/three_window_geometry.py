#!/usr/bin/env python3
"""Fixed, inspectable geometry for Q550 v18.53j.

All residuals are complete 5x5 skew tensors. The additional 'holonomy' is a
GEODESICALLY CLOSED Grassmannian triangle of principal 2-planes constructed
from 3 window-integrated bivectors, not a physical return map or a Noether
charge. Its O(2)-frame gauge cancels under the cyclic overlap construction.
"""
from __future__ import annotations
import math
import numpy as np
from scipy.integrate import cumulative_trapezoid
from bivector_geometry import principal_skew,cancel_metric


def accumulators(xi,t,Bu,Bc):
    """Trapezoidal source integrals using the source material clock t=tau(xi)."""
    xi=np.asarray(xi,float);t=np.asarray(t,float)
    bu=np.asarray(Bu,float);bc=np.asarray(Bc,float)
    if xi.ndim!=1 or t.shape!=xi.shape or bu.shape!=(len(xi),5,5) or bc.shape!=bu.shape:
        raise ValueError('BAD_INTEGRAL_INPUT_SHAPE')
    if np.any(np.diff(xi)<=0) or np.any(np.diff(t)<=0):raise ValueError('NONMONOTONE_CLOCK')
    cu=cumulative_trapezoid(bu,t,axis=0,initial=0.)
    cc=cumulative_trapezoid(bc,t,axis=0,initial=0.)
    return {'xi':xi,'t':t,'u':cu,'c':cc}


def at(acc,x):
    xi=acc['xi'];x=float(x)
    if x<xi[0]-1e-8 or x>xi[-1]+1e-8:raise ValueError('REQUEST_OUTSIDE_INTEGRATED_DOMAIN')
    k=int(np.searchsorted(xi,x,side='left'))
    if k==0:return acc['u'][0],acc['c'][0]
    if k==len(xi):return acc['u'][-1],acc['c'][-1]
    alpha=np.clip((x-xi[k-1])/(xi[k]-xi[k-1]),0,1)
    return tuple((1-alpha)*acc[key][k-1]+alpha*acc[key][k] for key in ('u','c'))


def interval(acc,a,b):
    if b<=a:raise ValueError('NONPOSITIVE_WINDOW')
    u0,c0=at(acc,a);u1,c1=at(acc,b)
    u=u1-u0;c=c1-c0
    return .5*(u-u.T),.5*(c-c.T)


def equal_norm_cv(norms):
    m=float(np.mean(norms))
    return float(np.std(norms)/m) if m>1e-15 else None


def bivector_pair_angle(A,B):
    na=float(np.linalg.norm(A));nb=float(np.linalg.norm(B))
    if na*nb<=1e-18:return None
    cos=float(np.clip(np.sum(A*B)/(na*nb),-1,1))
    return float(np.degrees(np.arccos(cos)))


def _closest_frame(next_basis,previous_frame):
    # Frame-independent minimal principal-angle geodesic transport between
    # adjacent 2-planes.  It becomes nonunique as a principal angle nears 90deg.
    M=np.asarray(next_basis).T@np.asarray(previous_frame)
    U,s,Vt=np.linalg.svd(M)
    return np.asarray(next_basis)@(U@Vt),float(np.min(s))


def polygon_plane_holonomy(bivectors,mode,cosine_warning=.08,gap_warning=.003):
    """Transport frames around P1->P2->P3->P1, with final edge geodesic.

    Angle is the absolute O(2)-gauge-invariant trace angle; not the holonomy of
    a physically closed orbit. Orthogonal projectors, not arbitrary Schur axes,
    are the underlying observable.
    """
    if len(bivectors)!=3:raise ValueError('EXPECTED_THREE_BIVECTORS')
    pd=[principal_skew(R) for R in bivectors]
    if any(z['relative_gap']<gap_warning for z in pd):
        return dict(defined=False,reason='NEAR_EQUAL_SINGULAR_STRENGTHS',angle_deg=None,
                    min_adjacent_overlap_cosine=None,max_principal_angle_deg=None)
    bases=[z['blocks'][mode]['plane'] for z in pd]
    if any(b is None for b in bases):
        return dict(defined=False,reason='MISSING_PRINCIPAL_TWO_PLANE',angle_deg=None,
                    min_adjacent_overlap_cosine=None,max_principal_angle_deg=None)
    transported=bases[0].copy()
    overlap=[]
    for next_basis in [bases[1],bases[2],bases[0]]:
        transported,c=_closest_frame(next_basis,transported);overlap.append(c)
    # Both transported and initial lie in first plane; a gauge invariant
    # rotation magnitude follows from the trace of the return overlap.
    W=bases[0].T@transported
    orth_err=float(np.linalg.norm(W.T@W-np.eye(2)))
    if orth_err>1e-8:raise RuntimeError('CLOSED_TRANSPORT_NOT_ORTHOGONAL')
    det=float(np.linalg.det(W))
    angle=float(np.degrees(np.arccos(np.clip(np.trace(W)/2,-1,1)))) if det>0 else None
    m=float(min(overlap))
    return dict(defined=bool(m>=cosine_warning and det>0),
        reason=('OK_GEODESIC_POLYGON_ONLY' if m>=cosine_warning and det>0 else 'AMBIGUOUS_GEODESIC_OR_ORIENTATION'),
        angle_deg=angle if m>=cosine_warning and det>0 else None,
        min_adjacent_overlap_cosine=m,max_principal_angle_deg=float(np.degrees(np.arccos(np.clip(m,0,1)))),
        closed_frame_orthogonality_error=orth_err,closed_frame_determinant=det)


def three_window(acc,bounds,cosine_warning=.08,gap_warning=.003):
    bounds=np.asarray(bounds,float)
    if bounds.shape!=(4,) or np.any(np.diff(bounds)<=0):raise ValueError('EXPECTED_FOUR_INCREASING_BOUNDARIES')
    rows=[]
    for a,b in zip(bounds[:-1],bounds[1:]):
        u,c=interval(acc,a,b);rows.append(dict(U=u,C=c,R=u+c))
    Rs=[d['R'] for d in rows]
    norms=[float(np.linalg.norm(R)) for R in Rs]
    net=sum(Rs,np.zeros((5,5)))
    net_direct=sum(interval(acc,bounds[0],bounds[-1]),np.zeros((5,5)))
    add_error=float(np.linalg.norm(net-net_direct)/max(1,np.linalg.norm(net_direct)))
    if add_error>1e-10:raise RuntimeError('WINDOW_ADDITIVITY_FAILURE')
    pair_angles=[bivector_pair_angle(Rs[0],Rs[1]),bivector_pair_angle(Rs[1],Rs[2]),bivector_pair_angle(Rs[0],Rs[2])]
    strength=sum(norms)
    closure_score=(1-float(np.linalg.norm(net))/strength) if strength>1e-14 else None
    Csum=sum((d['C'] for d in rows),np.zeros((5,5)))
    Usum=sum((d['U'] for d in rows),np.zeros((5,5)))
    holo=[polygon_plane_holonomy(Rs,i,cosine_warning,gap_warning) for i in range(2)]
    arr=np.stack([R.reshape(-1) for R in Rs],axis=0)
    gram=arr@arr.T;evals=np.sort(np.linalg.eigvalsh(gram))
    gram_min_fraction=float(max(0,evals[0])/max(np.sum(evals),1e-24))
    data=dict(xi_start=float(bounds[0]),xi_end=float(bounds[-1]),
        left_boundary_xi=float(bounds[1]),right_boundary_xi=float(bounds[2]),
        window_lengths_xi=[float(v) for v in np.diff(bounds)],
        norms=norms,net_norm=float(np.linalg.norm(net)),
        normalized_residual_closure_score=closure_score,
        equal_window_residual_norm_cv=equal_norm_cv(norms),
        pairwise_bivector_angles_deg=pair_angles,
        max_deviation_from_120_deg=max(abs(a-120) for a in pair_angles) if all(a is not None for a in pair_angles) else None,
        whole_interval_radial_phase_cancellation=cancel_metric(Usum,Csum)['invariant_cancellation'],
        gram_smallest_eigenvalue_fraction=gram_min_fraction,
        additive_bivector_error=add_error,
        plane1_geodesic_polygon_holonomy=holo[0],plane2_geodesic_polygon_holonomy=holo[1])
    return data,rows


def independent_direct_window_integral(xi,t,field,a,b):
    """Independent check of off-grid boundary interpolation by direct trapz."""
    from scipy.integrate import trapezoid
    xi=np.asarray(xi);t=np.asarray(t);field=np.asarray(field)
    mask=(xi>a)&(xi<b)
    points=np.r_[a,xi[mask],b]
    vals=[]
    for value in points:
        k=int(np.searchsorted(xi,value))
        if k==0:tt=t[0];vv=field[0]
        elif k==len(xi):tt=t[-1];vv=field[-1]
        else:
            fac=float((value-xi[k-1])/(xi[k]-xi[k-1]))
            tt=(1-fac)*t[k-1]+fac*t[k];vv=(1-fac)*field[k-1]+fac*field[k]
        vals.append((tt,vv))
    return trapezoid(np.stack([z[1] for z in vals]),
                     np.array([z[0] for z in vals]),axis=0)
