#!/usr/bin/env python3
"""Q550 v18.53k: 10D bivector SVD without asserting a physical 3rd axis.

The 5x5 skew upper triangle is multiplied by sqrt(2) to give the
Frobenius-isometric embedding into R^10.  M is a 3x10 array of
interval-integrated bivectors.  Its right singular directions live in
bivector-coordinate space, NOT physical 5D state space.
"""
from __future__ import annotations
import numpy as np

PAIRS=[(i,j) for i in range(5) for j in range(i+1,5)]

def to10(B):
    B=np.asarray(B,dtype=float)
    if B.shape!=(5,5):raise ValueError('SKEW_5X5_REQUIRED')
    if np.linalg.norm(B+B.T)>1e-10*max(1.,np.linalg.norm(B)):
        raise ValueError('ANTISYMMETRY_FAIL')
    return np.sqrt(2.)*np.array([B[i,j] for i,j in PAIRS],dtype=float)

def from10(v):
    v=np.asarray(v,dtype=float)
    if v.shape!=(10,):raise ValueError('VECTOR10_REQUIRED')
    B=np.zeros((5,5))
    for val,(i,j) in zip(v/np.sqrt(2),PAIRS):B[i,j]=val;B[j,i]=-val
    return B

def metrics(M,rank_floor=1e-8,gap_floor=0.02):
    M=np.asarray(M,dtype=float)
    if M.shape!=(3,10) or not np.all(np.isfinite(M)):
        raise ValueError('EXPECTED_THREE_10D_BIVECTORS')
    U,S,Vt=np.linalg.svd(M,full_matrices=False)
    net=M.sum(axis=0)
    norm=float(np.linalg.norm(net)); norms=np.linalg.norm(M,axis=1)
    strength=float(sum(norms));s1,s2,s3=map(float,S)
    gap=(s2-s3)/max(s2,1e-15)
    well_defined=bool(s3/max(s1,1e-15)>rank_floor and gap>gap_floor)
    # A gauge-invariant VECTOR, unlike the signed singular vector v3.
    projection=np.dot(net,Vt[2])
    T=projection*Vt[2]
    P=net-T
    score=float(1.-norm/strength) if strength>1e-15 else None
    angles=[]
    for i,j in ((0,1),(1,2),(0,2)):
        cosine=float(np.dot(M[i],M[j])/max(norms[i]*norms[j],1e-30))
        angles.append(float(np.degrees(np.arccos(np.clip(cosine,-1,1)))))
    return dict(S=S,V3=Vt[2].copy(),transverse=T,planar=P,net=net,
        rank3_strength_ratio=s3/max(s1,1e-15),rank23_gap=gap,
        direction_defined=well_defined,
        net_norm=norm,norms=norms,
        transverse_norm=float(np.linalg.norm(T)),planar_norm=float(np.linalg.norm(P)),
        transverse_share_squared=float(np.dot(T,T)/norm**2) if norm>1e-15 else None,
        cancellation_score=score,max_120_deviation_deg=float(max(abs(a-120.) for a in angles)),
        angles=angles,amplitude_cv=float(np.std(norms)/np.mean(norms)) if strength>1e-15 else None)

def transported_sign_frame(rows,reference_width=12.8,minimum_adjacent_abs_cos=0.25):
    """Phase-align v3 along sampled widths. Independent for every origin/seed/kappa.

    For v3 -> -v3, T=(net.v3)v3 is identical. The signed coefficient
    net.v3 is a CONDITIONAL diagnostic, only accepted where the frame
    can be continued through all adjacent overlap thresholds.
    """
    rows=sorted(rows,key=lambda r:r['width'])
    anchor=min(range(len(rows)),key=lambda i:abs(rows[i]['width']-reference_width))
    if abs(rows[anchor]['width']-reference_width)>1e-9:raise ValueError('NO_ANCHOR_WIDTH')
    vecs=[None]*len(rows);ok=[False]*len(rows);overlaps=[None]*len(rows)
    a=rows[anchor]
    if a['m']['direction_defined']:
        v=a['m']['V3'].copy()
        if np.dot(v,a['m']['net'])<0:v=-v
        vecs[anchor]=v;ok[anchor]=True;overlaps[anchor]=1.0
    for direction in (1,-1):
        i=anchor+direction
        while 0<=i<len(rows):
            prev=i-direction
            if not ok[prev] or not rows[i]['m']['direction_defined']:
                ok[i]=False;i+=direction;continue
            v=rows[i]['m']['V3'].copy()
            dot=float(np.dot(v,vecs[prev]));overlaps[i]=abs(dot)
            if abs(dot)<minimum_adjacent_abs_cos:
                ok[i]=False;i+=direction;continue
            vecs[i]=v if dot>=0 else -v
            ok[i]=True;i+=direction
    for i,row in enumerate(rows):
        row['transport_ok']=ok[i]
        row['adjacent_abs_cos']=overlaps[i]
        row['signed_transverse_coefficient']=float(np.dot(row['m']['net'],vecs[i])) if ok[i] else None
        row['transported_v3']=vecs[i]
        row['v3_projector_cosine_to_anchor']=float(abs(np.dot(row['m']['V3'],vecs[anchor]))) if ok[anchor] and row['m']['direction_defined'] else None
        row['v3_principal_angle_to_anchor_deg']=float(np.degrees(np.arccos(np.clip(row['v3_projector_cosine_to_anchor'],0,1)))) if row['v3_projector_cosine_to_anchor'] is not None else None
        ta=row['m']['transverse'];tb=rows[anchor]['m']['transverse']
        row['transverse_oriented_cosine_to_anchor']=(float(np.dot(ta,tb)/(np.linalg.norm(ta)*np.linalg.norm(tb)))
                    if np.linalg.norm(ta)>1e-12 and np.linalg.norm(tb)>1e-12 else None)
    for i,row in enumerate(rows):
        # Differentiation of a vector T(w) is invariant to the arbitrary sign of v3.
        if i==0 or i==len(rows)-1:row['dtransverse_dwidth_norm']=None;continue
        w_left,w_right=rows[i-1]['width'],rows[i+1]['width']
        row['dtransverse_dwidth_norm']=float(np.linalg.norm(rows[i+1]['m']['transverse']-rows[i-1]['m']['transverse'])/(w_right-w_left))
    return rows

def boundary_derivative(Bu,Bc,xi,tau_rate,w_start,win_width):
    """FTC prediction d/dw of integrated sum for fixed origin: 3 B_tau(o+3w) tau'(o+3w).
    The integrand is a material-clock torque bivector; interpolated on xi.
    """
    endpoint=w_start+3.*win_width
    xi=np.asarray(xi,float)
    if endpoint<xi[0]-1e-9 or endpoint>xi[-1]+1e-9:raise ValueError('ENDPOINT_OUT_OF_RANGE')
    k=int(np.searchsorted(xi,endpoint))
    k=min(max(k,1),len(xi)-1)
    t=float(np.clip((endpoint-xi[k-1])/(xi[k]-xi[k-1]),0,1))
    val=(1-t)*(Bu[k-1]+Bc[k-1])+t*(Bu[k]+Bc[k])
    speed=(1-t)*tau_rate[k-1]+t*tau_rate[k]
    return 3.*to10(val)*speed

def derivative_diagnostic(rows,xi,Bu,Bc,tau_rate):
    """Independent sampled d(sum)/dw versus FTC endpoint source, no new physical law."""
    for i,row in enumerate(rows):
        if i==0 or i==len(rows)-1:row['boundary_derivative_error']=None;continue
        prev,nxt=rows[i-1],rows[i+1]
        fd=(nxt['m']['net']-prev['m']['net'])/(nxt['width']-prev['width'])
        pred=boundary_derivative(Bu,Bc,xi,tau_rate,row['origin'],row['width'])
        row['boundary_derivative_error']=float(np.linalg.norm(fd-pred)/max(1.,np.linalg.norm(fd),np.linalg.norm(pred)))
        row['boundary_derivative_fd_norm']=float(np.linalg.norm(fd))
        row['boundary_derivative_pred_norm']=float(np.linalg.norm(pred))
    return rows
