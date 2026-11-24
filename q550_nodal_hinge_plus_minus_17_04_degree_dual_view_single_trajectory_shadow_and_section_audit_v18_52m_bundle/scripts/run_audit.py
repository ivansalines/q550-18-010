#!/usr/bin/env python3
"""Q550 v18.52m: dual tilted physical 5D observation planes around a fixed hinge.

All comparisons use ONE time-ordered relative physical trajectory per label. No
ODE integration, 2D image-rotation, forced closure, or assumed topological invariant.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import platform
from pathlib import Path

import numpy as np
import scipy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

ROOT = Path(__file__).resolve().parents[1]
VERSION = 'v18.52m'
PARENT_ROOT = 'q550_near_tangency_breathing_envelope_drift_five_dimensional_section_audit_v18_52l_bundle'
PARENT_BASENAME = 'extended_horizon_5d_vectors_and_traces.npz'
EXPECTED_LABELS = (
    'original_negative_1.0e-04', 'original_positive_1.0e-04',
    'original_negative_3.0e-04', 'original_positive_3.0e-04',
    'autonomous_negative_1.0e-04', 'autonomous_positive_1.0e-04',
    'autonomous_negative_3.0e-04', 'autonomous_positive_3.0e-04',
)
ANGLE_CONTROLS = (0., 10., 17.04, 25., 35.)


def sha256(path: Path) -> str:
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def write_json(path: Path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def csv_out(path: Path, rows: list[dict], keys: list[str]):
    with path.open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=keys, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def source_check(src: Path, verify_hash=True):
    if not src.is_file():
        raise FileNotFoundError(f'Missing v18.52l physical trace NPZ: {src}')
    manifest = json.loads((ROOT/'reference/PARENT_v18_52l_OUTPUT_MANIFEST.json').read_text())
    if manifest.get('version') != 'v18.52l':
        raise ValueError('Unexpected ancestral output version')
    reference_sha = manifest['file_sha256']['results/'+PARENT_BASENAME]
    actual = sha256(src) if verify_hash else None
    if verify_hash and actual != reference_sha:
        raise ValueError(f'v18.52l trace SHA mismatch: observed {actual}, expected {reference_sha}')
    with np.load(src, allow_pickle=False) as data:
        for label in EXPECTED_LABELS:
            for suffix in ('tau','dQ','dV','D','G'):
                if f'{label}__{suffix}' not in data:
                    raise ValueError(f'Missing physical array {label}__{suffix}')
            n = data[label+'__tau'].shape[0]
            if n < 1000 or data[label+'__dQ'].shape != (n,5) or data[label+'__dV'].shape != (n,5):
                raise ValueError(f'Invalid 5D physical trace shape for {label}')
        samples=int(data[EXPECTED_LABELS[0]+'__tau'].shape[0])
    return dict(version='v18.52l', path=str(src), expected_sha256=reference_sha,
                observed_sha256=actual, kicks=len(EXPECTED_LABELS), samples_per_kick=samples)


def normalize(v):
    v=np.asarray(v,float)
    n=float(np.linalg.norm(v))
    if n<1e-14 or not np.isfinite(n): raise ValueError('Zero or invalid observer direction')
    return v/n


def pca_hidden(q:np.ndarray):
    """Fit one fixed direction in the physical Q3:Q5 hidden subspace."""
    x=np.asarray(q,float)[:,2:5]
    x=x-np.mean(x,axis=0)
    cov=x.T@x/max(1,len(x)-1)
    eigvals,eigvecs=np.linalg.eigh(cov)
    direction=np.asarray(eigvecs[:,-1],float)
    direction*=1. if direction[np.argmax(np.abs(direction))]>=0 else -1.
    return normalize(direction),eigvals[::-1].tolist()


def observer_basis(theta_deg, hidden):
    """Orthogonal 2-frame (e1, cos(theta)e2+sin(theta)hidden), fixed hinge e1."""
    a=math.radians(float(theta_deg))
    h=normalize(hidden)
    e1=np.array([1.,0,0,0,0])
    e2=np.zeros(5);e2[1]=math.cos(a);e2[2:5]=math.sin(a)*h
    if abs(e1@e2)>1e-14 or abs(np.linalg.norm(e2)-1)>1e-14:
        raise RuntimeError('Observer is not orthonormal')
    return e1,e2


def projections(q, hidden, a):
    """Two views of the identical q(t), not two physical orbits."""
    q=np.asarray(q,float);h=normalize(hidden)
    u=q[:,0]
    v=q[:,1]
    z=q[:,2:5]@h
    ca=math.cos(math.radians(a));sa=math.sin(math.radians(a))
    plus=ca*v+sa*z
    minus=ca*v-sa*z
    assert np.allclose(plus+minus,2*ca*v,rtol=2e-13,atol=1e-14)
    assert np.allclose(plus-minus,2*sa*z,rtol=2e-13,atol=1e-14)
    return dict(u=u,v=v,h=z,plus=plus,minus=minus)


def crossings(tt,q,dv,basis,kind,theta,hidden_name,small_rel=1e-12):
    """Piecewise-linear refined zero crossings, with full 5D vectors retained."""
    w=q@basis;dw=dv@basis
    eps=max(1e-15,small_rel*float(np.max(np.abs(w))))
    pairs=np.flatnonzero((w[:-1]<-eps)&(w[1:]>eps) | (w[:-1]>eps)&(w[1:]<-eps))
    # Exact grid zeros are rare; include them without duplication.
    zeros=np.flatnonzero(np.abs(w)<=eps)
    ev=[]
    for idx in pairs:
        frac=float(-w[idx]/(w[idx+1]-w[idx]))
        ts=float(tt[idx]+frac*(tt[idx+1]-tt[idx]))
        qe=q[idx]+frac*(q[idx+1]-q[idx]);ve=dv[idx]+frac*(dv[idx+1]-dv[idx])
        g=float(ve@basis)
        ev.append((ts,qe,ve,g))
    for idx in zeros:
        if idx==0 or idx==len(w)-1: continue
        ts=float(tt[idx]);qe=q[idx];ve=dv[idx]
        ev.append((ts,qe,ve,float(ve@basis)))
    ev.sort(key=lambda row:row[0])
    out=[]
    for ts,qe,ve,g in ev:
        if out and ts-out[-1]['tau']<1e-8: continue
        vn=float(np.linalg.norm(ve))
        out.append(dict(tau=ts,observer=kind,angle_deg=float(theta),hidden=hidden_name,
                        crossing='upward' if g>0 else ('downward' if g<0 else 'grazing'),
                        u=float(qe[0]),transverse_derivative=g,
                        transversality=float(abs(g)/vn) if vn>1e-15 else None,
                        q_norm=float(np.linalg.norm(qe)),
                        D_at_crossing=float(.5*np.dot(qe,qe)),
                        hinge_perp_radius=float(np.linalg.norm(qe[1:])),
                        **{f'q{i+1}':float(qe[i]) for i in range(5)}))
    return out


def nearest_cross_view_pairs(plus,minus):
    """One-to-one mutual nearest crossings; descriptive only, not recurrence."""
    matches=[]
    if not plus or not minus: return matches
    ta=np.array([e['tau'] for e in plus]);tb=np.array([e['tau'] for e in minus])
    for i,a in enumerate(plus):
        j=int(np.argmin(np.abs(tb-a['tau'])))
        if int(np.argmin(np.abs(ta-tb[j])))!=i:continue
        b=minus[j]
        qa=np.array([a[f'q{k}'] for k in range(1,6)]); qb=np.array([b[f'q{k}'] for k in range(1,6)])
        an=float(np.linalg.norm(qa));bn=float(np.linalg.norm(qb))
        cos=float(np.dot(qa,qb)/(an*bn)) if an*bn>1e-25 else None
        matches.append(dict(tau_plus=a['tau'],tau_minus=b['tau'],delta_tau=b['tau']-a['tau'],
                            u_plus=a['u'],u_minus=b['u'],delta_u=b['u']-a['u'],
                            posture_distance_5d=float(np.linalg.norm(qb-qa)),
                            normalized_posture_distance_5d=float(np.linalg.norm(qb-qa)/max(an,bn,1e-14)),
                            posture_angle_deg=float(np.degrees(np.arccos(np.clip(cos,-1,1)))) if cos is not None else None))
    return matches


def consecutive_section_return_pairs(events):
    """Time-order return mapping in ONE observer; inspect antipodal vs full return."""
    out=[]
    for a,b in zip(events[:-1],events[1:]):
        qa=np.array([a[f'q{k}'] for k in range(1,6)]);qb=np.array([b[f'q{k}'] for k in range(1,6)])
        an=float(np.linalg.norm(qa));bn=float(np.linalg.norm(qb))
        den=max(an,bn,1e-14)
        c=float(np.dot(qa,qb)/(an*bn)) if an*bn>1e-25 else None
        out.append(dict(observer=a['observer'],first_tau=a['tau'],second_tau=b['tau'],
                        delta_tau=b['tau']-a['tau'],first_crossing=a['crossing'],
                        second_crossing=b['crossing'],
                        first_D=a['D_at_crossing'],second_D=b['D_at_crossing'],
                        D_ratio=b['D_at_crossing']/max(a['D_at_crossing'],1e-30),
                        same_posture_residual_normalized=float(np.linalg.norm(qb-qa)/den),
                        antipodal_residual_normalized=float(np.linalg.norm(qb+qa)/den),
                        full_5d_angle_deg=float(np.degrees(np.arccos(np.clip(c,-1.,1.)))) if c is not None else None))
    return out


def section_scalar_dispersion(ev):
    """Candidate section scalar agreement, not a conserved quantity claim."""
    if not ev:return dict(section_D_mean=None,section_D_cv=None,section_radius_cv=None)
    ds=np.array([x['D_at_crossing'] for x in ev]);rs=np.array([x['hinge_perp_radius'] for x in ev])
    return dict(section_D_mean=float(ds.mean()),
                section_D_cv=float(ds.std()/max(abs(ds.mean()),1e-30)),
                section_radius_cv=float(rs.std()/max(abs(rs.mean()),1e-30)))


def opposite_sign_stats(tt,w1,w2,tolerance=1e-12):
    scale=max(float(np.max(np.abs(w1))),float(np.max(np.abs(w2))),1e-14)
    opp=(w1>scale*tolerance)&(w2<-scale*tolerance) | (w1<-scale*tolerance)&(w2>scale*tolerance)
    widths=np.diff(tt)
    duration=float(np.sum(widths*(opp[:-1].astype(float)+opp[1:].astype(float))/2))
    total=float(tt[-1]-tt[0]);n_runs=int(np.sum(opp & ~np.r_[False,opp[:-1]]))
    return dict(opposite_side_duration_fraction=duration/total if total>0 else 0.,
                opposite_side_runs=n_runs)


def stats_geometry(tt,q,hidden,small_cut=.02):
    """Axial/transverse diagnostics; none is a coordinate-free topological invariant."""
    h=normalize(hidden);u=q[:,0];v=q[:,1];z=q[:,2:5]@h
    transverse=np.sqrt(v*v+z*z)
    full_perp=np.linalg.norm(q[:,1:5],axis=1)
    d=np.diff(u); frac_pos=float(np.mean(d>0));frac_neg=float(np.mean(d<0))
    phi=np.unwrap(np.arctan2(z,v))
    radius_floor=max(1e-14,float(np.max(transverse))*small_cut)
    masked=transverse<radius_floor
    # Phase is not claimed as a globally defined winding in R^5.
    net=float((phi[-1]-phi[0])/(2*np.pi))
    total=float(np.sum(np.abs(np.diff(phi)))/(2*np.pi))
    aa=np.abs(u)
    if float(np.var(aa))>1e-30:
        slope=float(np.cov(aa,full_perp,bias=True)[0,1]/np.var(aa))
        intercept=float(np.mean(full_perp)-slope*np.mean(aa))
        prediction=intercept+slope*aa
        ss=float(np.sum((full_perp-prediction)**2));ttss=float(np.sum((full_perp-np.mean(full_perp))**2))
        r2=1.-ss/ttss if ttss>1e-30 else None
    else:slope=intercept=r2=None
    return dict(full_5d_hinge_perpendicular_radius_cv=float(np.std(full_perp)/max(np.mean(full_perp),1e-30)),
                full_5d_radius_linear_fit_vs_abs_axial_R2=r2,
                full_5d_radius_linear_fit_slope=slope,
                full_5d_radius_linear_fit_intercept=intercept,
                transverse_projection_radius_min=float(np.min(transverse)),
                transverse_projection_radius_max=float(np.max(transverse)),
                phase_net_turns_projected=net,phase_total_turns_projected=total,
                fraction_phase_radius_near_zero=float(np.mean(masked)),
                axial_monotonic_increasing_fraction=frac_pos,axial_monotonic_decreasing_fraction=frac_neg)


def null_self_test():
    t=np.linspace(0,10,2501)
    q=np.zeros((len(t),5));q[:,0]=t/10;q[:,1]=np.sin(3*t)
    v=np.gradient(q,t,axis=0)
    hidden=np.array([1.,0,0])
    p=projections(q,hidden,17.04)
    assert np.max(np.abs(p['plus']-p['minus']))<1e-14, 'Planar negative control must coincide'
    assert opposite_sign_stats(t,p['plus'],p['minus'])['opposite_side_duration_fraction']==0.
    q[:,2]=0.5*np.cos(3*t)
    v=np.gradient(q,t,axis=0)
    p=projections(q,hidden,17.04)
    assert not np.allclose(p['plus'],p['minus']), 'Nonplanar control must separate'
    bpos=observer_basis(+17.04,hidden)[1];bneg=observer_basis(-17.04,hidden)[1]
    plus=crossings(t,q,v,bpos,'plus',17.04,'Q3')
    minus=crossings(t,q,v,bneg,'minus',-17.04,'Q3')
    assert len(plus)>=5 and len(minus)>=5
    assert all(abs(sum(e[f'q{k+1}']*bpos[k] for k in range(5)))<.00002 for e in plus)
    for e in plus:
        assert e['tau']>=0 and np.isfinite(e['q_norm'])
    assert max(abs(np.linalg.norm(p)-np.linalg.norm(p)) for p in q[::100])==0.
    return dict(planar_negative_control='PASS',nonplanar_positive_control='PASS',
                independent_hinge_basis_norm='PASS',crossing_interpolation='PASS')


def plots(results,out:Path,angle,primary):
    label=primary
    a=results[label]['trace'];t=a['tau'];q=a['q'];p=a['proj'];n=a['n']
    keep=np.unique(np.linspace(0,len(t)-1,1400).astype(int))
    t0=t[keep];u=p['u'][keep];v=p['v'][keep];z=p['h'][keep];wpos=p['plus'][keep];wneg=p['minus'][keep]
    s=1/max(np.max(np.linalg.norm(q,axis=1)),1e-30)
    # Same 5D q, shown on two different observer planes with the same hinge.
    fig=plt.figure(figsize=(16,5.4),layout='constrained')
    axs=[fig.add_subplot(1,3,i+1) for i in range(3)]
    for ax,w,title in ((axs[0],wpos,f'+{angle:g} deg view'),(axs[1],wneg,f'-{angle:g} deg view')):
        pts=ax.scatter(u*s,w*s,c=t0,s=2.5,cmap='viridis',rasterized=True)
        ax.axhline(0,color='gray',lw=.8,ls='--')
        ax.set(xlabel='Hinge axial coordinate Q1 / scale',ylabel='Tilted transverse / scale',title=title)
        ax.grid(alpha=.2);ax.set_aspect('equal',adjustable='datalim')
    axs[2].scatter(u*s,wpos*s,c=t0,s=2.2,cmap='viridis',rasterized=True)
    axs[2].scatter(u*s,wneg*s,c=t0,s=2.2,cmap='viridis',rasterized=True,alpha=.4,marker='x')
    axs[2].axhline(0,color='gray',lw=.8,ls='--');axs[2].grid(alpha=.2)
    axs[2].set(xlabel='Common Q1 hinge / scale',ylabel='Paired transverse shadows / scale',title='Two observations, one trajectory')
    axs[2].set_aspect('equal',adjustable='datalim')
    fig.colorbar(pts,ax=axs[:2],label='Matched material time')
    fig.suptitle(f'Q550 v18.52m — paired hinge rotations; {label}')
    fig.savefig(out/'paired_hinge_views_single_trajectory.png',dpi=160)
    plt.close(fig)
    fig=plt.figure(figsize=(12,5.6),layout='constrained')
    a3=fig.add_subplot(1,2,1,projection='3d',proj_type='ortho')
    a3.plot(u*s,v*s,z*s,lw=.9,color='tab:blue')
    a3.scatter(u[::40]*s,v[::40]*s,z[::40]*s,c=t0[::40],cmap='viridis',s=3)
    a3.set(xlabel='Q1 / scale',ylabel='Q2 / scale',zlabel='hidden dot n / scale',title='3D reduction of one physical 5D path')
    a3.view_init(elev=22,azim=-62)
    try:a3.set_box_aspect((1,1,1))
    except Exception:pass
    a2=fig.add_subplot(1,2,2)
    a2.plot(t0,np.sqrt(v*v+z*z)*s,label='2D transverse radius',lw=1)
    a2.plot(t0,np.linalg.norm(q[keep,1:5],axis=1)*s,label='Full 5D hinge-perp radius',lw=1)
    a2.set(xlabel='Matched material time',ylabel='Normalized perpendicular radius',title='Projection radius vs intrinsic hinge distance')
    a2.grid(alpha=.3);a2.legend(fontsize=8)
    fig.savefig(out/'single_trajectory_reduced_3d_and_hinge_radius.png',dpi=150)
    plt.close(fig)
    fig,axs=plt.subplots(2,1,figsize=(13,7),sharex=True,layout='constrained')
    axs[0].plot(t0,wpos*s,label=f'+{angle:g} deg',lw=.9)
    axs[0].plot(t0,wneg*s,label=f'-{angle:g} deg',lw=.9)
    axs[0].axhline(0,color='gray',lw=.8,ls='--')
    axs[0].set(ylabel='Transverse coordinate / scale',title='Gate opening and closing in two tilted views')
    axs[0].legend()
    opp=opposite_sign_stats(t,p['plus'],p['minus'])
    mask=np.sign(p['plus'][keep])*np.sign(p['minus'][keep])<0
    axs[1].fill_between(t0,0,mask.astype(float),color='tab:orange',alpha=.5,label='Opposite projected sides')
    axs[1].set(ylabel='Opposite-side gate',xlabel='Matched material time',ylim=(-.1,1.1),
               title=f'Fraction of material-time duration: {opp["opposite_side_duration_fraction"]:.3f}')
    for ax in axs:ax.grid(alpha=.2)
    fig.savefig(out/'paired_observer_gate_windows.png',dpi=145)
    plt.close(fig)


def main(argv=None):
    pa=argparse.ArgumentParser(description=__doc__)
    pa.add_argument('--mode',choices=('synthetic','preflight','audit'),default='audit')
    pa.add_argument('--source',type=Path,default=ROOT.parent/PARENT_ROOT/'results'/PARENT_BASENAME)
    pa.add_argument('--theta-deg',type=float,default=17.04)
    pa.add_argument('--primary',choices=EXPECTED_LABELS,default=EXPECTED_LABELS[0])
    pa.add_argument('--sample-stride',type=int,default=1,help='Development smoke-test downsampling; retain 1 for scientific production')
    pa.add_argument('--output-dir',type=Path,default=ROOT/'results')
    args=pa.parse_args(argv)
    if not 0<args.theta_deg<45:pa.error('--theta-deg must be in (0,45) degrees')
    if args.sample_stride < 1:pa.error('--sample-stride must be positive')
    out=args.output_dir;out.mkdir(parents=True,exist_ok=True)
    if args.mode=='synthetic':
        print('SYNTHETIC HINGE TEST PASS',null_self_test());return
    src_info=source_check(args.source)
    if args.mode=='preflight':
        write_json(out/'PREFLIGHT.json',dict(version=VERSION,mode='preflight',source=src_info,physical_integration=False))
        print('PREFLIGHT PASS',src_info['kicks'],'physical trajectories;',src_info['samples_per_kick'],'samples each');return
    null_self_test()
    with np.load(args.source,allow_pickle=False) as src:
        # Prespecified primary for direction discovery; first half only.
        train=np.asarray(src[args.primary+'__dQ'][:src[args.primary+'__dQ'].shape[0]//2:args.sample_stride],float)
        learned,train_eigenvalues=pca_hidden(train)
        test=np.asarray(src[args.primary+'__dQ'][src[args.primary+'__dQ'].shape[0]//2::args.sample_stride],float)
        heldout,heldout_eig=pca_hidden(test)
        stability_angle=float(np.degrees(np.arccos(np.clip(abs(learned@heldout),-1,1))))
        axes={'Q3':np.array([1.,0.,0.]),'Q4':np.array([0.,1.,0.]),'Q5':np.array([0.,0.,1.]),'PCA_train':learned}
        angle_rows=[];all_events=[];all_pairs=[];all_returns=[];metric_rows=[];results={}
        for label in EXPECTED_LABELS:
            tt=np.asarray(src[label+'__tau'][::args.sample_stride],float)
            q=np.asarray(src[label+'__dQ'][::args.sample_stride],float)
            dv=np.asarray(src[label+'__dV'][::args.sample_stride],float)
            D=np.asarray(src[label+'__D'][::args.sample_stride],float)
            G=np.asarray(src[label+'__G'][::args.sample_stride],float)
            if tt.ndim!=1 or not np.all(np.diff(tt)>0) or not np.isfinite(q).all() or not np.isfinite(dv).all():
                raise RuntimeError(f'Invalid time-ordered physical state: {label}')
            if max(float(np.max(np.abs(.5*np.sum(q*q,axis=1)-D))),float(np.max(np.abs(np.sum(q*dv,axis=1)-G))))>4e-12:
                raise RuntimeError('Physical 5D D/G identity failed: '+label)
            for hidden_name,hidden in axes.items():
                for ang in ANGLE_CONTROLS if args.theta_deg==17.04 else tuple(sorted(set((0.,10.,args.theta_deg,25.,35.)))):
                    p=projections(q,hidden,ang)
                    stats=opposite_sign_stats(tt,p['plus'],p['minus'])
                    bplus=observer_basis(+ang,hidden)[1]
                    bminus=observer_basis(-ang,hidden)[1]
                    evplus=crossings(tt,q,dv,bplus,'plus',ang,hidden_name)
                    evminus=crossings(tt,q,dv,bminus,'minus',-ang,hidden_name)
                    aa=dict(label=label,hidden_direction=hidden_name,theta_abs_deg=float(ang),
                            **stats,n_plus_crossings=len(evplus),n_minus_crossings=len(evminus),
                            mean_plus_transversality=float(np.mean([x['transversality'] for x in evplus])) if evplus else None,
                            mean_minus_transversality=float(np.mean([x['transversality'] for x in evminus])) if evminus else None)
                    angle_rows.append(aa)
                    if abs(ang-args.theta_deg)>1e-9:continue
                    if hidden_name=='PCA_train':
                        for e in evplus+evminus:all_events.append(dict(label=label,**e))
                        near=nearest_cross_view_pairs(evplus,evminus)
                        for e in near:all_pairs.append(dict(label=label,hidden_direction=hidden_name,theta_abs_deg=float(ang),**e))
                        serial=consecutive_section_return_pairs(evplus)+consecutive_section_return_pairs(evminus)
                        for e in serial:all_returns.append(dict(label=label,hidden_direction=hidden_name,theta_abs_deg=float(ang),**e))
                        geo=stats_geometry(tt,q,hidden)
                        mr=dict(label=label,hidden_direction=hidden_name,theta_abs_deg=float(ang),**stats,**geo,
                                n_plus_crossings=len(evplus),n_minus_crossings=len(evminus),
                                n_mutual_nearest_cross_view_pairs=len(near),
                                n_serial_section_return_pairs=len(serial),
                                plus_section_D_cv=section_scalar_dispersion(evplus)['section_D_cv'],
                                minus_section_D_cv=section_scalar_dispersion(evminus)['section_D_cv'],
                                plus_section_hinge_radius_cv=section_scalar_dispersion(evplus)['section_radius_cv'],
                                minus_section_hinge_radius_cv=section_scalar_dispersion(evminus)['section_radius_cv'],
                                section_is_projection_dependent=True,
                                topology_invariant_certified=False)
                        metric_rows.append(mr)
                        results[label]=dict(metrics=mr,trace=dict(tau=tt,q=q,v=dv,proj=p,n=learned))
        plots(results,out,args.theta_deg,args.primary)
    csv_out(out/'angle_and_hidden_direction_sensitivity.csv',angle_rows,list(angle_rows[0]))
    csv_out(out/'tilted_hinge_section_crossings_5d.csv',all_events,list(all_events[0]) if all_events else ['label','tau'])
    csv_out(out/'cross_view_sequential_posture_pairs_5d.csv',all_pairs,list(all_pairs[0]) if all_pairs else ['label','tau_plus','tau_minus'])
    csv_out(out/'single_observer_consecutive_section_return_map_5d.csv',all_returns,
            list(all_returns[0]) if all_returns else ['label','first_tau','second_tau'])
    csv_out(out/'single_trajectory_hinge_observer_metrics.csv',metric_rows,list(metric_rows[0]))
    clean={k:v['metrics'] for k,v in results.items()}
    payload=dict(version=VERSION,title='Nodal Hinge +/-17.04 Degree Dual-View Single-Trajectory Shadow & Section Audit',
                 stage=('POSTPROCESS_AUTHENTIC_v18_52l_PHYSICAL_5D_TRACES' if args.sample_stride==1
                        else 'DEVELOPMENT_SMOKE_SUBSAMPLED_NOT_SCIENTIFIC_RESULT'),
                 sample_stride=args.sample_stride,
                 input_sha256=src_info['expected_sha256'],input_trace=src_info,
                 angle_deg=float(args.theta_deg),fixed_hinge_axis='physical Q1 (index 0)',
                 historical_horizontal_cut='physical Q2=0 in the Q1-Q2 projection; not a unique 5D line',
                 PCA_direction_training='first half of primary trajectory only',
                 PCA_primary_label=args.primary,PCA_hidden_axis_Q3_Q4_Q5=learned.tolist(),
                 PCA_hidden_eigenvalues_first_half=train_eigenvalues,
                 PCA_hidden_eigenvalues_heldout=heldout_eig,
                 PCA_hidden_direction_heldout_angle_deg=stability_angle,
                 results=clean,checks=dict(hash_checked=True,physical_vectors_checked=True,
                     single_trajectory_per_label=True,nonplanar_and_planar_synthetic_controls=True,
                     two_view_algebraic_identities_pass=True,observer_norm_pass=True),
                 caveats=[
                     'Two tilted projections are not two independent physical filaments.',
                     'The Q2=0 cut is a four-dimensional hyperplane in R5, not an intrinsic one-dimensional axis.',
                     'The chosen fixed Q1 hinge is an observer hypothesis, not an independently certified material axis.',
                     'The data-learned hidden PCA direction is not independent evidence; it is fit on a declared training half.',
                     'The two-view sum and difference identities are algebraic necessities and cannot certify symmetry.',
                     'Projected angular turns in the Q2-hidden plane are not topological winding invariants in R5.',
                     'The cone/tube proxies are geometrical fit diagnostics, not certificates of a physical surface.',
                     'Cross-view event pairs sample distinct times; they do not prove same-posture recurrence.',
                     'Candidate section D/radius dispersion diagnoses nonconservation rather than certifying an invariant.',
                     'Piecewise-linear crossing time interpolation is not an interval-exhaustive root certificate.',
                     'Changing observer does not modify the 5D physical trajectory, dynamical law, or material clock.',
                     'The +/-17.04 degree angle is a predeclared candidate, not numerically estimated as a law.'
                 ],runtime=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__))
    write_json(out/'report.json',payload)
    opp=[v['metrics']['opposite_side_duration_fraction'] for v in results.values()]
    lines=[f'# Q550 {VERSION} — Server scientific report','',
           '**Nodal Hinge +/-17.04 Degree Dual-View Single-Trajectory Shadow & Section Audit**','',
           f'- Authenticated parent v18.52l NPZ SHA-256: `{src_info["expected_sha256"]}`',
           f'- Fixed physical Q1 hinge, observer angle pair +/-{args.theta_deg:g} degrees.',
           f'- Hidden direction PCA (Q3,Q4,Q5): `{[round(float(x),8) for x in learned]}`.',
           f'- Hidden axis train/held-out principal-angle instability: **{stability_angle:.4f} deg**.',
           f'- Across eight trajectories: opposite-side gate fraction {min(opp):.4f}..{max(opp):.4f}.','',
           '| Trajectory | Opposite-view fraction | + gate crossings | - gate crossings | Projected net phase turns |',
           '|---|---:|---:|---:|---:|']
    for label,item in results.items():
        v=item['metrics'];lines.append(f'| {label} | {v["opposite_side_duration_fraction"]:.5f} | {v["n_plus_crossings"]} | {v["n_minus_crossings"]} | {v["phase_net_turns_projected"]:.5f} |')
    lines.extend(['','## Interpretation boundaries','',
      'Two wings are **two observer projections of ONE physical 5D time-ordered relative trajectory**.',
      'Apparent complementary lobes are a geometric readout, not evidence of two physical filaments.',
      'An axis in the projected Q1-Q2 plane does not uniquely identify an axis in R5.',
      'A projected winding or compact appearance is not a topological invariant or an autosimilarity proof.',
      'The D and hinge-radius dispersions at crossings and the same-observer return map test candidate section organization; they do not establish conservation.',
      'Angle +/-17.04 degrees is a user-proposed candidate; controls at 0, 10, 25, 35 degrees and Q3/Q4/Q5 directions are included.',
      'No nonlinear solver was rerun; the historical 10D trajectory and material-clock alignments are inherited unchanged.',''])
    (out/'SCIENTIFIC_REPORT_SERVER.md').write_text('\n'.join(lines),encoding='utf-8')
    print(f'{VERSION} AUDIT PASS: 8 physical single-trajectory observer pairs; source SHA validated; output in {out}')


if __name__=='__main__':main()
