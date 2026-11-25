#!/usr/bin/env python3
"""Q550 v18.53y: trajectory-induced minimal-turn SO(5) frame and Lambda^2 transport.

The kinematic transport is derived from the authenticated mass-normalized 5D
velocity along a Q550 11D trajectory. Its minimal-twist condition is a stated
geometric choice, NOT a unique physical material connection or a new 11D RHS.
All scientific strings and comments are in English.
"""
from __future__ import annotations
import argparse,csv,hashlib,importlib.util,io,json,sys,tempfile,zipfile
from pathlib import Path
from contextlib import ExitStack
import numpy as np
from connection_math import KinematicFrame, wedge5, exterior_generator, gauge, geodesic_polygon, angle_deg

ROOT=Path(__file__).resolve().parents[1]
CFG=json.loads((ROOT/'inputs_config.json').read_text())
OUT=ROOT/'results';OUT.mkdir(exist_ok=True)
PARENT=ROOT/'inputs'/CFG['parent_filename']
def sha(b):return hashlib.sha256(b).hexdigest()

def read_parent_zip(z):
    if z.testzip() is not None:raise RuntimeError('PARENT_X_CRC_FAILED')
    roots={n.split('/')[0] for n in z.namelist() if '/' in n}
    if len(roots)!=1:raise RuntimeError('AMBIGUOUS_PARENT_ROOT')
    root=next(iter(roots))+'/'
    counts={}
    for man in ('SOURCE_MANIFEST.json','OUTPUT_MANIFEST.json'):
        manifest=json.loads(z.read(root+man))['files']
        for name,digest in manifest.items():
            if sha(z.read(root+name))!=digest:raise RuntimeError('X_MANIFEST_MISMATCH '+name)
        counts[man]=len(manifest)
    report=json.loads(z.read(root+'results/report.json'))
    if report.get('n_parent_nodes')!=204 or 'GAUGE_DEPENDENT' not in report['scientific_decision']:
        raise RuntimeError('UNEXPECTED_X_REPORT')
    fixed=list(csv.DictReader(io.StringIO(z.read(root+'results/fixed_physical_frame_path.csv').decode())))
    return root,counts,fixed

def import_authentic_x(stack,full):
    if not PARENT.is_file() or sha(PARENT.read_bytes())!=CFG['parent_sha256']:
        raise RuntimeError('PARENT_X_SHA256_MISMATCH')
    stage=Path(stack.enter_context(tempfile.TemporaryDirectory(prefix='q550_v1853y_x_')))
    with zipfile.ZipFile(PARENT) as z:
        root,counts,fixed=read_parent_zip(z)
        for src in ('scripts/gauge_frame_audit.py','inputs/u_metric_geometry_authenticated.py',
                    'inputs/v18_53w_authenticated_server_sendback.zip','inputs_config.json'):
            dest=stage/src;dest.parent.mkdir(parents=True,exist_ok=True)
            dest.write_bytes(z.read(root+src))
    spec=importlib.util.spec_from_file_location('authenticated_x',stage/'scripts/gauge_frame_audit.py')
    x=importlib.util.module_from_spec(spec);sys.modules[spec.name]=x;spec.loader.exec_module(x)
    wrows,vrows,polygons,paths,traces,ancestry=x.load_parent(light=not full)
    if len(fixed)!=204 or len(vrows)!=204:raise RuntimeError('PARENT_NODE_COUNT')
    return x,paths,traces,fixed,{'x':counts,**ancestry}

def record_csv(path,rows):
    if not rows:raise RuntimeError('NO_OUTPUT_ROWS')
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)

def pull_to_initial(N,time,frame):
    R=frame.at(time)
    # R(t): columns of the transported velocity-aligned physical frame.
    # A 10D bivector is pulled back by Lambda^2(R(t)^T).
    return np.array([wedge5(mat.T)@n for mat,n in zip(R,N)])

def moving_gauge_check(N,ts,frame,amplitude):
    R=frame.at(ts)
    G=np.array([gauge(t,frame.t[0],frame.t[-1],amplitude)[0] for t in ts])
    fixed=pull_to_initial(N,ts,frame)
    gauged=np.array([wedge5(q)@n for q,n in zip(G,N)])
    reconstructed=np.array([wedge5((q@r).T)@n for q,r,n in zip(G,R,gauged)])
    err=float(np.max(np.linalg.norm(fixed-reconstructed,axis=1)))
    a1,p1,_=geodesic_polygon(fixed)
    a2,p2,_=geodesic_polygon(reconstructed)
    # Exact gauge connection transformation, sampled at an authentic physical event.
    t=float(ts[len(ts)//2]);Q,Qp=gauge(t,frame.t[0],frame.t[-1],amplitude)
    A=frame.omega(t)
    if np.linalg.norm(A+A.T)>2e-12:raise RuntimeError('KINEMATIC_GENERATOR_NON_SKEW')
    Ag=Q@A@Q.T+Qp@Q.T
    Q10=wedge5(Q)
    Qp10=exterior_generator(Qp@Q.T)@Q10
    L=exterior_generator(A);Lg=exterior_generator(Ag)
    generr=float(np.linalg.norm(Lg-(Q10@L@Q10.T+Qp10@Q10.T)))
    return dict(max_gauge_pullback_error=err,geo_arc_error=abs(a1-a2),
                polygon_error=float(np.linalg.norm(p1-p2)),
                generator_covariance_error=generr,
                physical_omega_frobenius=float(np.linalg.norm(A)),
                correction_is_pure_gauge_identity=True)

def main(mode):
    with ExitStack() as stack:
        x,paths,traces,fixed,anc=import_authentic_x(stack,mode!='preflight')
        print('AUTHENTIC_X_W_V_U_T_ANCESTRY_PASS',anc,flush=True)
        if mode=='preflight':return
        parent_fixed={(r['seed'],r['kappa_rule'],r['epoch'],float(r['lambda_value'])):r for r in fixed}
        selected=sorted(paths,key=lambda k:(k[0],k[1],k[2],k[3]))
        if mode=='smoke':selected=[k for k in selected if k[:3]==('original','shared','preformal') and k[3] in (0.,1.,4.)]
        if len(selected)!=(3 if mode=='smoke' else 204):raise RuntimeError('WRONG_SELECTION_SIZE')
        frames={};prox={};stores={};rows=[];gauge_rows=[];refinement=[];replay_error=0.
        base={};max_gauge=0.;min_speed=1e9
        for n,key in enumerate(selected,1):
            seed,rule,epoch,lam=key;tag=f'{seed}__{rule}'
            if tag not in frames:
                xi,z,C,g=traces[tag]
                V=x.mg.physical_velocity(z[:,:11],CFG['m_x'])
                speed=np.linalg.norm(V,axis=1)
                min_speed=min(min_speed,float(speed.min()))
                if speed.min()<CFG['min_allowed_speed']:raise RuntimeError('VELOCITY_FRAME_SINGULAR')
                frames[tag]=KinematicFrame(xi,V,1)
                prox[tag]=x.trace_primitive(traces,seed,rule)[0]
            frame=frames[tag]
            aa,ww,NN=x.get_path(prox[tag],paths[key]);n_samples=len(aa)
            raw_arc,raw_poly,raw_u=geodesic_polygon(NN)
            if lam==0.:base[seed,rule,epoch,'fixed']=raw_poly
            canon_angle=angle_deg(base[seed,rule,epoch,'fixed'],raw_poly)
            recorded=float(parent_fixed[key]['fixed_frame_physical_path_angle_deg'])
            err=float(abs(canon_angle-recorded)) if np.isfinite(canon_angle) else float('inf')
            replay_error=max(replay_error,err)
            if err>CFG['fixed_parent_angle_tolerance_deg']:raise RuntimeError('X_CANONICAL_REPLAY_BREACH '+str(key)+' '+str(err))
            for lift in CFG['lifts']:
                times=aa if lift=='start' else aa+3*ww if lift=='end' else aa+1.5*ww
                if np.min(times)<frame.t[0] or np.max(times)>frame.t[-1]:raise RuntimeError('LIFT_OUTSIDE_TRACE')
                Nm=pull_to_initial(NN,times,frame)
                arc,poly,u=geodesic_polygon(Nm)
                if lam==0.:base[seed,rule,epoch,lift]=poly
                rel_angle=angle_deg(base[seed,rule,epoch,lift],poly)
                rows.append(dict(seed=seed,kappa_rule=rule,epoch=epoch,lambda_value=lam,lift=lift,
                    canonical_fixed_angle_deg=canon_angle,dynamic_frame_angle_deg=rel_angle,
                    canonical_arc_rad=raw_arc,dynamic_frame_arc_rad=arc,
                    dynamic_minus_canonical_arc_rad=arc-raw_arc,
                    canonical_polygon_norm=float(np.linalg.norm(raw_poly)),
                    dynamic_polygon_norm=float(np.linalg.norm(poly)),
                    min_10d_residual_norm=float(np.min(np.linalg.norm(NN,axis=1))),
                    max_10d_residual_norm=float(np.max(np.linalg.norm(NN,axis=1))),
                    min_5d_velocity=float(np.min(np.linalg.norm(frame.v,axis=1))),
                    max_frame_orthogonality=frame.orthogonality,
                    max_frame_tangent_error=frame.tangent_error,
                    parent_fixed_angle_replay_error_deg=err,
                    temporal_reference=lift,
                    endpoint_beyond_preformal=bool((aa+3*ww).max()>CFG['formal_clock_boundary']),
                    n_path_samples=n_samples))
                if lam in CFG['selected_lambdas'] and mode=='full':
                    code=f'{tag}__{epoch}__lam{int(round(lam/.25)):02d}__{lift}'
                    stores[code+'__times']=times
                    stores[code+'__normalized10']=u
                    stores[code+'__polygon45']=poly
                if lam==CFG['gauge_test_lambda'] and lift in CFG['gauge_test_lifts']:
                    chk=moving_gauge_check(NN,times,frame,CFG['gauge_amplitude_rad'])
                    max_gauge=max(max_gauge,chk['max_gauge_pullback_error'])
                    gauge_rows.append(dict(seed=seed,kappa_rule=rule,epoch=epoch,lift=lift,
                                           lambda_value=lam,**chk))
                    if max(chk['polygon_error'],chk['generator_covariance_error'],chk['max_gauge_pullback_error'])>CFG['gauge_tolerance']:
                        raise RuntimeError('GAUGE_COVARIANCE_FAILURE '+str(key)+' '+lift)
            if (mode=='full' and lam==4.) or (mode=='smoke' and lam==4.):
                coarse=KinematicFrame(frame.t,frame.v,CFG['refinement_stride'])
                for lift in CFG['lifts']:
                    times=aa if lift=='start' else aa+3*ww if lift=='end' else aa+1.5*ww
                    P1=pull_to_initial(NN,times,frame)
                    Pc=pull_to_initial(NN,times,coarse)
                    arc1,poly1,_=geodesic_polygon(P1)
                    arcc,polyc,_=geodesic_polygon(Pc)
                    refinement.append(dict(seed=seed,kappa_rule=rule,epoch=epoch,lambda_value=lam,lift=lift,
                       stride_coarse=CFG['refinement_stride'],
                       max_10d_coordinate_frame_difference=float(np.max(np.linalg.norm(P1-Pc,axis=1))),
                       covariant_arclength_refinement_error_rad=abs(arc1-arcc),
                       covariant_polygon_refinement_angle_deg=angle_deg(poly1,polyc),
                       coarse_frame_orthogonality=coarse.orthogonality))
            if n%12==0 or n==len(selected):
                print('V1853Y_TRANSPORT',n,'/',len(selected),'LIFTS',len(rows),'PARENT_ERROR',f'{replay_error:.3e}',flush=True)
        if mode=='smoke':
            print('V1853Y_AUTHENTIC_SMOKE_PASS',len(rows),'REPLAY_ERROR',replay_error,'MAX_GAUGE',max_gauge,flush=True)
            return
        if len(rows)!=612 or len(gauge_rows)!=24 or len(refinement)!=36:
            raise RuntimeError('INCOMPLETE_TRANSPORT_OR_REFINEMENT')
        record_csv(OUT/'dynamic_frame_oriented_transport.csv',rows)
        record_csv(OUT/'moving_frame_gauge_covariance.csv',gauge_rows)
        record_csv(OUT/'frame_resolution_convergence.csv',refinement)
        np.savez_compressed(OUT/'selected_dynamic_10d_45d_transport.npz',**stores)
        from plot_transport import make_plots
        make_plots(rows,OUT)
        by_lift={lift:max(r['dynamic_frame_angle_deg'] for r in rows if r['lift']==lift) for lift in CFG['lifts']}
        max_change={lift:max(abs(r['dynamic_minus_canonical_arc_rad']) for r in rows if r['lift']==lift) for lift in CFG['lifts']}
        report=dict(version=CFG['version'],scientific_decision='AUTHENTIC_11D_TRAJECTORY_INDUCED_MINIMAL_TURN_KINEMATIC_CONNECTION__GAUGE_COVARIANT_TRANSPORT__NO_UNIQUE_MATERIAL_CONNECTION',
                    n_parent_nodes=204,n_lift_rows=len(rows),n_moving_gauge_tests=len(gauge_rows),
                    n_refinement_checks=len(refinement),n_new_physical_ivps=0,
                    ancestor_manifest_counts=anc,
                    minimum_5d_speed=min_speed,
                    maximum_parent_x_fixed_frame_replay_error_deg=replay_error,
                    maximum_gauge_pullback_error=max_gauge,
                    maximum_generator_gauge_error=max(r['generator_covariance_error'] for r in gauge_rows),
                    maximum_polygon_gauge_error=max(r['polygon_error'] for r in gauge_rows),
                    maximum_frame_refinement_polygon_angle_deg=max(r['covariant_polygon_refinement_angle_deg'] for r in refinement),
                    maximum_dynamic_frame_angle_by_temporal_lift_deg=by_lift,
                    maximum_absolute_dynamic_minus_fixed_arclength_by_lift_rad=max_change,
                    limits=[
                     'The 5D tangent is reconstructed from authenticated mass-normalized 11D phase trajectories; its derivative uses cubic spline differentiation of stored velocities, not a newly integrated 11D tangent map.',
                     'Omega=e_prime wedge e is the unique minimum-Frobenius-norm skew generator satisfying Omega e=e_prime. This is a kinematic no-twist convention, not a physically unique material connection.',
                     'A 10D residual is an integral over [a,a+3w] and has no intrinsically assigned time. Start, midpoint, and end lifts are distinct, explicitly tested conventions.',
                     'Gauge recovery and generator covariance are controlled algebraic/numerical identities; they cannot alone certify material holonomy.',
                     'The canonical 10D metric is used for the dynamical-frame transport; noncanonical velocity-occupation metric sensitivity remains a separate question.',
                     'All endpoints of the near-return windows extend beyond xi=25.6 in the formal clock continuation. Neither 11D material recurrence nor unique frame transport is certified.'
                    ])
        (OUT/'report.json').write_text(json.dumps(report,indent=2)+'\n')
        (OUT/'SCIENTIFIC_REPORT_SERVER.md').write_text(
            '# Q550 v18.53y — 11D-generated kinematic frame, 10D induced transport and gauge checks\n\n'
            +report['scientific_decision']+'\n\n'
            +f"Authentic parent nodes: 204; temporal-lift observations: {len(rows)}; gauge tests: {len(gauge_rows)}; refinement checks: {len(refinement)}.\n\n"
            +f"Minimum model 5D speed: {min_speed:.10g}. Parent x replay error: {replay_error:.4e} deg. Gauge pullback maximum: {max_gauge:.4e}.\n\n"
            +f"Maximum dynamical-frame angle by temporal lift (degrees): {by_lift}.\n\n"
            +'The frame is a model-trajectory-induced minimum-turn SO(5) convention, not a unique material transport law or evidence of 11D material holonomy.\n')
        print('V1853Y',report['scientific_decision'],'ROWS',len(rows),'MAX_GAUGE',max_gauge,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--preflight',action='store_true');p.add_argument('--smoke',action='store_true');p.add_argument('--full',action='store_true')
    a=p.parse_args();main('full' if a.full else 'smoke' if a.smoke else 'preflight')
