#!/usr/bin/env python3
"""Q550 v18.53x — gauge identification versus fixed physical-frame geometry.

The physical 5D Cartesian frame and its exterior-square 10D coordinates are fixed.
Gauge frames are *orthonormal representations* of the same SPD metric, produced
only from SO(5) rotations; Lambda^2 SO(5) acts on the whitened 10D bivectors.
No 11D/21D equations are changed, and no new IVP is integrated.
Science-facing strings and comments are in English.
"""
from __future__ import annotations
import argparse, csv, hashlib, importlib.util, io, json, math, sys, zipfile
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import minimize_scalar

ROOT=Path(__file__).resolve().parents[1]
CFG=json.loads((ROOT/'inputs_config.json').read_text())
OUT=ROOT/'results'
OUT.mkdir(parents=True,exist_ok=True)
spec=importlib.util.spec_from_file_location('u_authenticated_geometry',ROOT/'inputs/u_metric_geometry_authenticated.py')
mg=importlib.util.module_from_spec(spec);sys.modules[spec.name]=mg;spec.loader.exec_module(mg)
PAIRS45=mg.PAIRS45

def sha(b:bytes)->str:return hashlib.sha256(b).hexdigest()

def verify_zip(z:zipfile.ZipFile):
    if (bad:=z.testzip()) is not None:raise RuntimeError('ZIP_CRC_FAILED '+bad)
    roots={s.split('/')[0] for s in z.namelist() if '/' in s}
    if len(roots)!=1:raise RuntimeError('PARENT_ROOT_AMBIGUOUS')
    root=next(iter(roots))+'/'
    checks={}
    for n in ['SOURCE_MANIFEST.json','OUTPUT_MANIFEST.json']:
        manifest=json.loads(z.read(root+n))['files']
        for rel,digest in manifest.items():
            if sha(z.read(root+rel))!=digest:raise RuntimeError('ANCESTOR_HASH_FAILURE '+n+' '+rel)
        checks[n]=len(manifest)
    return root,checks

def load_parent(light=False):
    parent=(ROOT/'inputs'/CFG['parent_filename']).read_bytes()
    if sha(parent)!=CFG['parent_sha256']:raise RuntimeError('WRONG_W_PARENT_HASH')
    helper=(ROOT/'inputs/u_metric_geometry_authenticated.py').read_bytes()
    if sha(helper)!=CFG['u_helper_sha256']:raise RuntimeError('WRONG_U_GEOMETRY_HASH')
    with zipfile.ZipFile(io.BytesIO(parent)) as zw:
        wr,wcount=verify_zip(zw)
        wreport=json.loads(zw.read(wr+'results/report.json'))
        if wreport['n_cases']!=CFG['n_parent_nodes'] or 'AUTHENTIC_THREE_FACTOR' not in wreport['scientific_decision']:
            raise RuntimeError('UNEXPECTED_W_SCIENCE')
        w_rows=list(csv.DictReader(io.StringIO(zw.read(wr+'results/metric_path_factorial_45d.csv').decode())))
        with np.load(io.BytesIO(zw.read(wr+'results/metric_path_factorial_vectors.npz')),allow_pickle=False) as f:
            polygons={name:np.asarray(f[name]).copy() for name in f.files}
        vblob=zw.read(wr+'inputs/v18_53v_authenticated_server_sendback.zip')
    with zipfile.ZipFile(io.BytesIO(vblob)) as zv:
        vr,vcount=verify_zip(zv)
        vrows=list(csv.DictReader(io.StringIO(zv.read(vr+'results/metric_lambda_continuation.csv').decode())))
        vcurves=list(csv.DictReader(io.StringIO(zv.read(vr+'results/metric_lambda_valley_curves.csv').decode())))
        ublob=zv.read(vr+'inputs/v18_53u_authenticated_server_sendback.zip')
    with zipfile.ZipFile(io.BytesIO(ublob)) as zu:
        ur,ucount=verify_zip(zu)
        if sha(zu.read(ur+'scripts/metric_geometry.py'))!=CFG['u_helper_sha256']:
            raise RuntimeError('U_SCRIPT_MISMATCH')
        tblob=zu.read(ur+'inputs/v18_53t_authenticated_server_sendback.zip')
    with zipfile.ZipFile(io.BytesIO(tblob)) as zt:
        tr,tcount=verify_zip(zt)
        if not light:
            with np.load(io.BytesIO(zt.read(tr+'results/augmented21_primitives_and_oriented_valleys.npz')),allow_pickle=False) as f:
                traces={}
                for seed in ['original','autonomous']:
                    for rule in ['shared','original','autonomous']:
                        tag=f'{seed}__{rule}'
                        traces[tag]=tuple(np.asarray(f[tag+'__'+k]).copy() for k in ('xi','z21','C21','g21'))
        else:traces={}
    if len(w_rows)!=204 or len(vrows)!=204 or len(vcurves)!=28764:
        raise RuntimeError('INCOMPLETE_W_V_PARENT')
    if len([k for k in polygons if k.endswith('__P011')])!=204:
        raise RuntimeError('MISSING_W_POLYGONS')
    paths={}
    for row in vcurves:
        key=(row['seed'],row['kappa_rule'],row['epoch'],float(row['lambda_value']))
        paths.setdefault(key,[]).append((float(row['origin']),float(row['optimized_width'])))
    if len(paths)!=204 or any(len(p)!=141 for p in paths.values()):
        raise RuntimeError('BAD_V_WIDTH_PATHS')
    return w_rows,vrows,polygons,paths,traces,{'w':wcount,'v':vcount,'u':ucount,'t':tcount}

def exterior_square(A):
    """General exterior-square matrix for square A, in ordered i<j basis."""
    n=A.shape[0]
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    return np.array([[A[i,k]*A[j,l]-A[i,l]*A[j,k] for k,l in pairs] for i,j in pairs])

def rotation5(i,j,angle):
    c,s=np.cos(angle),np.sin(angle)
    R=np.eye(5);R[i,i]=c;R[j,j]=c;R[i,j]=-s;R[j,i]=s
    return R

def frame5(family,lam):
    if family=='symmetric':return np.eye(5)
    if family=='intraplane':return rotation5(0,1,0.10*lam)
    if family=='mixed':return rotation5(0,4,0.08*lam)@rotation5(1,3,-0.06*lam)
    if family=='countermixed':return rotation5(0,4,-0.08*lam)@rotation5(1,3,0.06*lam)
    raise ValueError('UNKNOWN_GAUGE_FAMILY '+family)

def angular_deg(a,b):
    na,nb=np.linalg.norm(a),np.linalg.norm(b)
    if na<1e-14 or nb<1e-14:return float('nan')
    return float(np.degrees(np.arccos(np.clip(np.dot(a,b)/(na*nb),-1,1))))

def save_csv(path,rows):
    if not rows:raise RuntimeError('EMPTY_CSV_'+str(path))
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def trace_primitive(traces,seed,rule):
    xi,z,C,g=traces[seed+'__'+rule]
    return CubicHermiteSpline(xi,C,g,axis=0,extrapolate=False),xi,z

def get_path(prox,samples):
    origins=np.array([a for a,w in samples],float)
    widths=np.array([w for a,w in samples],float)
    N=np.array([prox(a+3*w)-prox(a) for a,w in zip(origins,widths)])
    if not np.isfinite(N).all():raise RuntimeError('NONFINITE_10D_PATH')
    return origins,widths,N

def width_minimizer(prox,H,a,bounds,grid=21):
    lo,hi=bounds
    ws=np.linspace(lo,hi,grid)
    val=lambda w:float((prox(a+3*w)-prox(a))@H@(prox(a+3*w)-prox(a)))
    vals=np.array([val(w) for w in ws])
    trials=[(float(vals[0]),float(ws[0])),(float(vals[-1]),float(ws[-1]))]
    inds=[i for i in range(1,grid-1) if vals[i]<=vals[i-1] and vals[i]<=vals[i+1]]
    if not inds:inds=[int(np.argmin(vals[1:-1]))+1]
    for i in inds:
        fit=minimize_scalar(val,bounds=(ws[i-1],ws[i+1]),method='bounded',options={'xatol':2e-11})
        trials.append((float(fit.fun),float(fit.x)))
    return min(trials,key=lambda x:x[0])[1]

def moving_path(prox,center,H,reference_width,n):
    aa=np.linspace(center-CFG['halfspan'],center+CFG['halfspan'],n)
    bounds=(reference_width-0.05,reference_width+0.05)
    ww=np.array([width_minimizer(prox,H,float(a),bounds) for a in aa])
    NN=np.array([prox(a+3*w)-prox(a) for a,w in zip(aa,ww)])
    return aa,ww,NN

def moving_frame_test(N,H,amplitude):
    """A moving SO(5) gauge creates naive extra rotation, removed by transport.

    Frame R_i maps canonical whitened u_i to v_i=R_i u_i. Pullback with R_i.T
    is precisely the discrete covariant frame identification, not a theorem
    about the physical connection of the 11D dynamics.
    """
    W=mg.canonical_sqrt(H)
    u=N@W.T
    u=u/np.linalg.norm(u,axis=1)[:,None]
    s=np.linspace(0,1,len(u))
    RR=np.array([mg.wedge_map(rotation5(0,4,amplitude*np.sin(2*np.pi*x))
                  @rotation5(1,3,amplitude*.5*np.sin(4*np.pi*x))) for x in s])
    v=np.einsum('nij,nj->ni',RR,u)
    na,poly_na,_=mg.oriented_polygon(v,np.eye(10))
    truth,poly_truth,_=mg.oriented_polygon(u,np.eye(10))
    pulled=np.einsum('nji,nj->ni',RR,v)
    co,poly_cov,_=mg.oriented_polygon(pulled,np.eye(10))
    # Local pair-wise connection: transport next v_i+1 back into the frame at i.
    pair_errors=[]
    for i in range(len(u)-1):
        next_to_i=RR[i]@RR[i+1].T@v[i+1]
        pair_errors.append(np.linalg.norm(RR[i].T@next_to_i-u[i+1]))
    diagnostics=dict(raw_arc=na,covariant_arc=co,fixed_arc=truth,
                raw_poly_norm=float(np.linalg.norm(poly_na)),
                max_connection_pair_error=float(max(pair_errors)),
                covariant_arc_error=float(abs(co-truth)),
                covariant_polygon_error=float(np.linalg.norm(poly_cov-poly_truth)))
    return diagnostics, (u,v,poly_na,poly_truth)

def plot_results(gauge_rows,physical_rows):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for seed in ['original','autonomous']:
        fig,axes=plt.subplots(1,2,figsize=(11,4.5),constrained_layout=True)
        for epoch in ['preformal','full_formal']:
            for family in CFG['gauge_families']:
                tab=sorted([r for r in gauge_rows if r['seed']==seed and r['kappa_rule']=='shared' and r['epoch']==epoch and r['gauge']==family],key=lambda r:r['lambda_value'])
                axes[0].plot([r['lambda_value'] for r in tab],[r['gauge_fixed_center_angle_deg'] for r in tab],label=f'{family} / {epoch}',lw=1.4)
            phy=sorted([r for r in physical_rows if r['seed']==seed and r['kappa_rule']=='shared' and r['epoch']==epoch],key=lambda r:r['lambda_value'])
            axes[1].plot([r['lambda_value'] for r in phy],[r['fixed_frame_physical_path_angle_deg'] for r in phy],marker='.',label=epoch)
        axes[0].set(ylabel='Cross-lambda 45D angle in chosen gauge [deg]',xlabel='Metric strength lambda',title='Same metric, different orthonormal frames')
        axes[1].set(ylabel='Fixed physical-frame 45D path angle [deg]',xlabel='Metric strength lambda',title='Physical 10D path, canonical frame')
        for ax in axes:ax.grid(alpha=.2);ax.legend(fontsize=7)
        fig.suptitle(f'Q550 v18.53x: gauge dependence and physical-frame separation ({seed}, shared)')
        fig.savefig(OUT/f'gauge_fixed_frame_comparison_{seed}.png',dpi=135)
        plt.close(fig)

def main(mode):
    wrows,vrows,polygons,paths,traces,ancestry=load_parent(light=(mode=='preflight'))
    print('AUTHENTIC_W_V_U_T_ANCESTRY_PASS',ancestry,flush=True)
    if mode=='preflight':return
    selected=wrows if mode=='full' else [next(r for r in wrows if r['seed']=='original' and r['kappa_rule']=='shared' and r['epoch']=='preformal' and float(r['lambda_value'])==lam) for lam in [0.,1.,4.]]
    gauge_rows=[];physical_rows=[];moving_rows=[];variable_rows=[]
    err_cov=0.;err_same=0.;err_fixed=0.;err_moving=0.;max_connection=0.;artifacts={}
    cache={}
    for n,row in enumerate(selected,1):
        seed=row['seed'];rule=row['kappa_rule'];epoch=row['epoch'];lam=float(row['lambda_value'])
        tag=f'{seed}__{rule}';code=f'{tag}__{epoch}__lam{int(round(lam/.25)):02d}'
        P0=polygons[code+'__P000'];Pw=polygons[code+'__P011'];Pm=polygons[code+'__P111']
        if tag not in cache:
            cache[tag]=trace_primitive(traces,seed,rule)
        prox,xi,z=cache[tag]
        if (tag,epoch) not in cache:
            S=mg.occupation_covariance(xi,z,CFG['m_x'],'preformal' if epoch=='preformal' else 'full')
            cache[tag,epoch]=S
        H=mg.induced(mg.occupied_metric(cache[tag,epoch],lam))
        physical_path=paths[(seed,rule,epoch,lam)]
        origins,widths,N=get_path(prox,physical_path)
        normal_arc,physical_poly,_=mg.oriented_polygon(N,np.eye(10))
        weighted_arc,replayed_poly,_=mg.oriented_polygon(N,H)
        err_fixed=max(err_fixed,float(np.linalg.norm(replayed_poly-Pw)))
        if mode=='full' and err_fixed>CFG['parent_polygon_tolerance']:
            raise RuntimeError('HISTORICAL_POLYGON_REPLAY_FAILURE '+code)
        fixed_angle=angular_deg(P0,physical_poly)
        physical_rows.append(dict(seed=seed,kappa_rule=rule,epoch=epoch,lambda_value=lam,
            fixed_frame_physical_path_angle_deg=fixed_angle,
            parent_metric_angle_deg=angular_deg(P0,Pw),
            fixed_frame_physical_arc=float(normal_arc),
            parent_metric_arc=float(weighted_arc),
            historical_poly_replay_error=float(np.linalg.norm(replayed_poly-Pw)),
            fixed_frame_representation='mass-normalized Cartesian physical coordinates',
            time_endpoint_formal=bool(max(a+3*w for a,w in physical_path)>CFG['formal_clock_boundary'])))
        for family in CFG['gauge_families']:
            O=frame5(family,lam)
            R10=mg.wedge_map(O)
            R45=exterior_square(R10)
            if np.linalg.norm(R10.T@R10-np.eye(10))>5e-13 or np.linalg.norm(R45.T@R45-np.eye(45))>5e-13:
                raise RuntimeError('GAUGE_NOT_ORTHOGONAL')
            Qw=R45@Pw;Qm=R45@Pm
            rec=R45.T@Qw
            cov_err=float(np.linalg.norm(rec-Pw))
            err_cov=max(err_cov,cov_err)
            invariant_err=abs(angular_deg(Pw,Pm)-angular_deg(Qw,Qm))
            err_same=max(err_same,invariant_err)
            if mode=='full' and cov_err>CFG['gauge_covariance_tolerance']:
                raise RuntimeError('GAUGE_PULLBACK_FAILURE '+code)
            gauge_rows.append(dict(seed=seed,kappa_rule=rule,epoch=epoch,lambda_value=lam,gauge=family,
                 gauge_fixed_center_angle_deg=angular_deg(P0,Qw),
                 gauge_moving_center_angle_deg=angular_deg(P0,Qm),
                 parent_fixed_center_angle_deg=angular_deg(P0,Pw),
                 parent_moving_center_angle_deg=angular_deg(P0,Pm),
                 common_lambda_moving_vs_fixed_deg=angular_deg(Qw,Qm),
                 common_lambda_invariant_error_deg=invariant_err,
                 covariant_pullback_error=cov_err,
                 norm45_invariant_error=float(abs(np.linalg.norm(Qw)-np.linalg.norm(Pw))),
                 gauge_frame_max_rotation_deg=float(np.degrees(2*np.arcsin(np.clip(np.linalg.norm(O-np.eye(5),2)/2,0,1))))))
        # The moving-center reconstruction is intentionally separately budgeted.
        if lam in CFG['moving_path_lambdas']:
            # The canonical reference baseline width and moving minimum were authenticated by w.
            width0=float(row['baseline_width'])
            mov=float(row['moving_origin'])
            ns=CFG['n_samples'] if mode=='full' else 41
            aam,wwm,Nm=moving_path(prox,mov,H,width0,ns)
            _,Qphysical,_=mg.oriented_polygon(Nm,np.eye(10))
            _,Qmetric,_=mg.oriented_polygon(Nm,H)
            metric_err=float(np.linalg.norm(Qmetric-Pm)) if mode=='full' else float('nan')
            if mode=='full' and metric_err>CFG['parent_moving_polygon_tolerance']:
                raise RuntimeError('HISTORICAL_MOVING_POLYGON_REPLAY_FAILURE '+code+' '+str(metric_err))
            if mode=='full':err_moving=max(err_moving,metric_err)
            moving_rows.append(dict(seed=seed,kappa_rule=rule,epoch=epoch,lambda_value=lam,
                  physical_frame_moving_center_angle_deg=angular_deg(P0,Qphysical),
                  metric_frame_moving_center_angle_deg=angular_deg(P0,Pm),
                  physical_frame_vs_fixed_center_angle_deg=angular_deg(physical_poly,Qphysical) if mode=='full' else float('nan'),
                  replay_metric_polygon_error=metric_err, moving_canonical_path_arc=float(mg.oriented_polygon(Nm,np.eye(10))[0])))
            if lam==CFG['gauge_moving_frame_lambda']:
                for name,Htest in [('canonical',np.eye(10)),('weighted',H)]:
                    diag,(u,v,pna,pco)=moving_frame_test(Nm,Htest,CFG['moving_gauge_amplitude_rad'])
                    max_connection=max(max_connection,diag['covariant_polygon_error'])
                    variable_rows.append(dict(seed=seed,kappa_rule=rule,epoch=epoch,metric=name,lambda_value=lam,**diag,
                        additional_naive_arc=diag['raw_arc']-diag['fixed_arc']))
                    if mode=='full' and diag['covariant_polygon_error']>CFG['gauge_covariance_tolerance']:
                        raise RuntimeError('MOVING_FRAME_CONNECTION_FAILURE '+code)
                    if mode=='full' and seed=='original' and rule=='shared':
                        artifacts[f'{epoch}__{name}__fixed_u']=u
                        artifacts[f'{epoch}__{name}__gauge_v']=v
                        artifacts[f'{epoch}__{name}__raw_polygon']=pna
                        artifacts[f'{epoch}__{name}__covariant_polygon']=pco
        if n%12==0 or n==len(selected):print('GAUGE_FRAME_NODES',n,'/',len(selected),'MAX_PULLBACK',f'{err_cov:.2e}',flush=True)
    if mode=='smoke':
        print('AUTHENTIC_V1853X_SMOKE_PASS',len(gauge_rows),'GAUGES',len(moving_rows),'MOVING',flush=True)
        return
    if len(gauge_rows)!=816 or len(physical_rows)!=204 or len(moving_rows)!=60 or len(variable_rows)!=24:
        raise RuntimeError('INCOMPLETE_GAUGE_DATA')
    save_csv(OUT/'metric_identification_gauge_scan.csv',gauge_rows)
    save_csv(OUT/'fixed_physical_frame_path.csv',physical_rows)
    save_csv(OUT/'moving_center_fixed_physical_frame.csv',moving_rows)
    save_csv(OUT/'moving_gauge_connection_checks.csv',variable_rows)
    np.savez_compressed(OUT/'selected_fixed_frame_vectors.npz',**artifacts)
    plot_results(gauge_rows,physical_rows)
    max_angle=max(r['gauge_fixed_center_angle_deg'] for r in gauge_rows)
    min_angle=min(r['gauge_fixed_center_angle_deg'] for r in gauge_rows)
    max_fixed=max(r['fixed_frame_physical_path_angle_deg'] for r in physical_rows)
    metric_vs_physical_at_lam4={}
    for epoch in ['preformal','full_formal']:
        subset=[r for r in physical_rows if r['epoch']==epoch and r['lambda_value']==4.]
        metric_vs_physical_at_lam4[epoch]=dict(physical_mean_deg=float(np.mean([r['fixed_frame_physical_path_angle_deg'] for r in subset])),
              metric_mean_deg=float(np.mean([r['parent_metric_angle_deg'] for r in subset])))
    decision='GAUGE_DEPENDENT_CROSS_METRIC_45D_ANGLES__FIXED_PHYSICAL_FRAME_DIAGNOSTIC__NO_MATERIAL_HOLONOMY'
    report=dict(version=CFG['version'],scientific_decision=decision,
         n_parent_nodes=204,n_gauge_rows=816,n_fixed_physical_frame_rows=204,
         n_reconstructed_moving_center_rows=60,n_moving_gauge_rows=24,
         n_new_physical_ivps=0,ancestry=ancestry,
         max_constant_gauge_pullback_error=err_cov,
         max_constant_gauge_common_metric_angle_error_deg=err_same,
         max_historical_fixed_center_replay_error=err_fixed,
         max_historical_moving_center_replay_error=err_moving,
         max_moving_gauge_connection_error=max_connection,
         min_cross_metric_gauge_angle_deg=min_angle,max_cross_metric_gauge_angle_deg=max_angle,
         max_fixed_frame_physical_path_angle_deg=max_fixed,
         lambda4_group_means=metric_vs_physical_at_lam4,
         limits=[
          'Any cross-lambda 45D angle requires an identification of orthonormal metric frames; it is not intrinsically gauge invariant.',
          'All demonstration gauge frames are induced by SO(5) rotations of mass-normalized Cartesian coordinates, not arbitrary SO(45) matrices.',
          'The fixed-physical-frame diagnostic uses a common 5D Cartesian kinetic frame and its induced Euclidean 10D and 45D bases; this is a declared physical frame, not a unique coordinate-free material transport law.',
          'The moving xi-dependent frame uses an explicitly imposed pure-gauge connection to test consistency; its corrected recovery is a covariance identity, not independently generated material holonomy.',
          'The velocity-occupation metrics are model-motivated sensitivity metrics and not alternate fundamental kinetic metrics.',
          'The parent physical trajectories are reused without new IVPs; interpolation and finite quadrature remain numerical approximations.',
          'The window endpoints extend past xi=25.6 in the formal clock continuation; no material 11D recurrence is certified.'
         ])
    (OUT/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    (OUT/'SCIENTIFIC_REPORT_SERVER.md').write_text('# Q550 v18.53x — Gauge identification and fixed physical-frame audit\n\n'
       +decision+'\n\n'
       +f'204 parent nodes; 816 gauge tests; 204 fixed physical frame curves; 60 moving-center curves; 24 xi-dependent-gauge controls.\n\n'
       +f'Gauge pullback maximum {err_cov:.4e}. Fixed-path replay maximum {err_fixed:.4e}. Moving-center replay maximum {err_moving:.4e}.\n\n'
       +'Cross-lambda angles under independent orthonormal gauge frames are convention-dependent. Corrected pure-gauge transport recovers the fixed physical frame geometrically. Neither is a material holonomy certificate.\n')
    print('V1853X',decision,'GAUGES',len(gauge_rows),'PHYSICAL',len(physical_rows),'MOVING',len(moving_rows),'MAX_PULLBACK',err_cov,flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--preflight',action='store_true')
    parser.add_argument('--smoke',action='store_true')
    parser.add_argument('--full',action='store_true')
    a=parser.parse_args()
    main('full' if a.full else 'smoke' if a.smoke else 'preflight')
