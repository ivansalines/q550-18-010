#!/usr/bin/env python3
"""Observer-dependent 2-parameter roots of a fixed 2D projection in Lambda^2 R^5.

The three-window NET has only two material-clock endpoints. For numerically
piecewise-linear cumulative C(xi), the derivatives below are exact on open grid
cells, not at knots. Input tensors and the frame use Frobenius-isometric 10D.
"""
from __future__ import annotations
import numpy as np
from scipy.optimize import least_squares

PAIRS=[(a,b) for a in range(5) for b in range(a+1,5)]

def to10(a):
    a=np.asarray(a,dtype=float)
    if a.shape!=(5,5):raise ValueError('EXPECTED_5X5')
    if np.linalg.norm(a+a.T)>1e-8*max(1.,np.linalg.norm(a)):
        raise ValueError('NOT_SKEW')
    return np.sqrt(2.)*np.array([a[i,j] for i,j in PAIRS],dtype=float)

def primitive_from_acc(acc):
    q=np.asarray(acc['u'])+np.asarray(acc['c'])
    return np.asarray(acc['xi'],float),np.stack([to10(x) for x in q])

def primitive_at(xi,primitive,x):
    if x < xi[0]-1e-10 or x > xi[-1]+1e-10:raise ValueError('ENDPOINT_OUTSIDE_TRACE')
    x=float(np.clip(x,xi[0],xi[-1]))
    k=int(np.searchsorted(xi,x,side='right')-1)
    k=max(0,min(len(xi)-2,k))
    a=(x-xi[k])/(xi[k+1]-xi[k])
    return primitive[k]+a*(primitive[k+1]-primitive[k])

def primitive_slope(xi,primitive,x):
    # One-sided at knots; report knot clearance as a separate condition.
    if x < xi[0]-1e-10 or x > xi[-1]+1e-10:raise ValueError('SLOPE_OUTSIDE_TRACE')
    k=int(np.searchsorted(xi,x,side='right')-1)
    k=max(0,min(len(xi)-2,k))
    return (primitive[k+1]-primitive[k])/(xi[k+1]-xi[k])

def validate_trace_coverage(xi,origin_range,width_range,atol=1e-10):
    """Reject optimizer boxes that could ask for unrecorded xi values.

    The three equal windows telescope to endpoints [o,o+3*w].
    For a rectangular box, its earliest start is origin_range[0] and its
    latest endpoint is origin_range[1]+3*width_range[1].
    Never cure an out-of-domain trial by silently clipping the physics.
    """
    x=np.asarray(xi,float)
    if x.ndim!=1 or len(x)<2 or not np.all(np.diff(x)>0):
        raise ValueError('INVALID_TRACE_GRID')
    lo_o,hi_o=map(float,origin_range)
    lo_w,hi_w=map(float,width_range)
    if not (lo_o<hi_o and 0<lo_w<hi_w):
        raise ValueError('INVALID_PARAMETER_BOX')
    if lo_o < x[0]-atol:
        raise ValueError(f'OPTIMIZER_BOX_START_BELOW_TRACE: {lo_o:.12g} < {x[0]:.12g}')
    if hi_o+3*hi_w > x[-1]+atol:
        raise ValueError(f'OPTIMIZER_BOX_END_ABOVE_TRACE: {hi_o+3*hi_w:.12g} > {x[-1]:.12g}')
    return {'xi_min':float(x[0]),'xi_max':float(x[-1]),
            'box_start_min':lo_o,'box_end_max':hi_o+3*hi_w}

def near_knot_distance(xi,x):
    k=int(np.searchsorted(xi,x))
    return float(min(abs(x-xi[max(0,k-1)]),abs(x-xi[min(len(xi)-1,k)])))

def basis_from_triplet(acc,o,w):
    from three_window_geometry import interval
    rows=[]
    for j in range(3):
        a,b=o+j*w,o+(j+1)*w
        u,c=interval(acc,a,b)
        rows.append(to10(u+c))
    matrix=np.stack(rows)
    U,S,Vt=np.linalg.svd(matrix,full_matrices=False)
    base=Vt[:2].T.copy()
    for j in range(2):
        pivot=int(np.argmax(np.abs(base[:,j])))
        if base[pivot,j]<0:base[:,j]*=-1
    relative_gap=float((S[1]-S[2])/max(S[1],1e-15))
    return base,S,relative_gap

def residual(xi,c,o,w):
    s,e=float(o),float(o+3*w)
    return primitive_at(xi,c,e)-primitive_at(xi,c,s)

def projected(xi,c,basis,point):
    o,w=map(float,point)
    return basis.T@residual(xi,c,o,w)

def boundary_jacobian(xi,c,basis,point):
    o,w=map(float,point)
    g0=basis.T@primitive_slope(xi,c,o)
    g1=basis.T@primitive_slope(xi,c,o+3*w)
    return np.column_stack((g1-g0,3*g1))

def determinant_from_edges(gstart,gend):
    """det([gend-gstart,3*gend]) = -3 det([gstart,gend])."""
    return -3*float(np.linalg.det(np.column_stack((gstart,gend))))

def singular_diagnostics(jac,scaled_columns=None):
    s=np.linalg.svd(jac,compute_uv=False)
    ss=np.linalg.svd(jac@np.diag(scaled_columns),compute_uv=False) if scaled_columns is not None else s
    return {'det':float(np.linalg.det(jac)),'sigma_max':float(s[0]),
       'sigma_min':float(s[-1]),'condition_number':float(s[0]/max(s[-1],1e-300)),
       'scaled_sigma_min':float(ss[-1]),
       'scaled_condition_number':float(ss[0]/max(ss[-1],1e-300))}

def fixed_grid(xi,c,basis,origins,widths):
    # Endpoint additivity makes the 2D map cheap to scan, no new IVPs.
    grid=[]
    for o in origins:
        for w in widths:
            f=projected(xi,c,basis,(o,w))
            grid.append((float(o),float(w),float(np.linalg.norm(f)),float(f[0]),float(f[1])))
    return grid

def spaced_seeds(grid,n,orange,wrange):
    sorted_points=sorted(grid,key=lambda r:r[2]);seeds=[]
    for o,w,_,_,_ in sorted_points:
        if all(np.hypot((o-u)/(orange[1]-orange[0]),(w-v)/(wrange[1]-wrange[0]))>0.06
               for u,v in seeds):seeds.append((o,w))
        if len(seeds)>=n:break
    for p in [(14.4,12.8),(15.2,12.809),(13.6,12.813)]:
        if p not in seeds:seeds.append(p)
    return seeds

def optimize_projected_root(xi,c,basis,start,orange,wrange):
    validate_trace_coverage(xi,orange,wrange)
    lb=np.asarray((orange[0],wrange[0]),float)
    ub=np.asarray((orange[1],wrange[1]),float)
    initial=np.clip(np.asarray(start,dtype=float),lb+1e-9,ub-1e-9)
    res=least_squares(lambda p:projected(xi,c,basis,p),initial,
        jac=lambda p:boundary_jacobian(xi,c,basis,p),bounds=(lb,ub),
        method='trf',ftol=1e-13,xtol=1e-13,gtol=1e-13,max_nfev=180)
    p=np.asarray(res.x,float)
    f=projected(xi,c,basis,p)
    J=boundary_jacobian(xi,c,basis,p)
    s,e=p[0],p[0]+3*p[1]
    raw=singular_diagnostics(J,(orange[1]-orange[0],wrange[1]-wrange[0]))
    gs=basis.T@primitive_slope(xi,c,s)
    ge=basis.T@primitive_slope(xi,c,e)
    # Domain-aware diagnostic finite differences; the optimizer is allowed
    # to converge exactly on a native trace boundary. At such points use
    # a one-sided stencil instead of requesting xi outside the measured trace.
    def valid_trial(q):
        return (xi[0]-1e-12<=q[0]<=xi[-1]+1e-12 and
                xi[0]-1e-12<=q[0]+3*q[1]<=xi[-1]+1e-12)
    fdp=[]
    for axis in (0,1):
        d=np.zeros(2);d[axis]=1e-5
        plus,minus=p+d,p-d
        if valid_trial(plus) and valid_trial(minus):
            derivative=(projected(xi,c,basis,plus)-projected(xi,c,basis,minus))/(2e-5)
        elif valid_trial(plus):
            derivative=(projected(xi,c,basis,plus)-f)/1e-5
        elif valid_trial(minus):
            derivative=(f-projected(xi,c,basis,minus))/1e-5
        else:
            raise ValueError('FINITE_DIFFERENCE_OUTSIDE_TRACE')
        fdp.append(derivative)
    fd=np.column_stack(fdp)
    raw.update(origin_xi=float(s),width_xi=float(p[1]),end_xi=float(e),
      projected_f1=float(f[0]),projected_f2=float(f[1]),
      projected_norm=float(np.linalg.norm(f)),
      full_bivector_norm=float(np.linalg.norm(residual(xi,c,s,p[1]))),
      complement_bivector_norm=float(np.linalg.norm(residual(xi,c,s,p[1])-basis@f)),
      end_start_edge_source_det=determinant_from_edges(gs,ge),
      boundary_det_identity_error=float(abs(np.linalg.det(J)-determinant_from_edges(gs,ge))),
      jacobian_fd_norm_difference=float(np.linalg.norm(J-fd)),
      origin_knot_clearance_xi=near_knot_distance(xi,s),
      end_knot_clearance_xi=near_knot_distance(xi,e),
      optimizer_success=bool(res.success),optimizer_message=str(res.message),
      optimizer_nfev=int(res.nfev),optimizer_cost=float(res.cost),
      on_search_boundary=bool(np.min(np.r_[p-lb,ub-p])<1e-6),
      original_seed_origin=float(initial[0]),original_seed_width=float(initial[1]))
    return raw

def unique_solutions(results,o_tol=2e-3,w_tol=1e-3):
    unique=[]
    for r in sorted(results,key=lambda d:d['projected_norm']):
        if not any(abs(r['origin_xi']-e['origin_xi'])<o_tol and abs(r['width_xi']-e['width_xi'])<w_tol
                   for e in unique):unique.append(r)
    return unique
