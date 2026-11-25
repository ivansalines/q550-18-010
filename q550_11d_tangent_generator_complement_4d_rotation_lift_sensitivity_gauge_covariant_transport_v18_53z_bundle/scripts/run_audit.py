#!/usr/bin/env python3
"""Q550 v18.53z: model-conditioned 11D tangent complement vs minimal-turn frame.

Authenticated 11D state/source traces are immutable. The external phase is
recovered from the physical 10D source. The tangent map is a midpoint finite-lag
linearization of the *historical material RHS*, not an independently integrated
11D orbit. A horizontal lift and its alternative are declared. No material
connection/holonomy is asserted.
"""
from __future__ import annotations
import argparse,csv,hashlib,importlib.util,io,json,sys,tempfile,zipfile
from contextlib import ExitStack
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicHermiteSpline, CubicSpline
from tangent_source import infer_omega,reconstruct_phase,finite_lag_spin,jacobian_11d,rhs_material_xi
from complement_frame import ComplementFrame,pullback,skew
from connection_math import geodesic_polygon,angle_deg,wedge5,gauge,exterior_generator

ROOT=Path(__file__).resolve().parents[1]
CFG=json.loads((ROOT/'inputs_config.json').read_text())
PARENT=ROOT/'inputs'/CFG['parent_filename'];OUT=ROOT/'results';OUT.mkdir(exist_ok=True)
SHA=lambda data:hashlib.sha256(data).hexdigest()

def import_py(path,name):
    p=Path(path);sys.path.insert(0,str(p.parent))
    spec=importlib.util.spec_from_file_location(name,p)
    m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m)
    return m

def csvout(path,rows):
    if not rows:raise RuntimeError('NO_RESULTS_TO_SAVE')
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)

def parent_import(stack,read_traces):
    if not PARENT.is_file() or SHA(PARENT.read_bytes())!=CFG['parent_sha256']:
        raise RuntimeError('Y_PARENT_HASH_MISMATCH')
    stage=Path(stack.enter_context(tempfile.TemporaryDirectory(prefix='q550_v1853z_y_')))
    with zipfile.ZipFile(PARENT) as z:
        if (bad:=z.testzip()) is not None:raise RuntimeError('Y_PARENT_CRC_FAILED '+bad)
        names=z.namelist();roots={n.split('/')[0] for n in names if '/' in n}
        if len(roots)!=1:raise RuntimeError('AMBIGUOUS_Y_PARENT_ROOT')
        root=next(iter(roots))+'/'
        counts={}
        for man in ('SOURCE_MANIFEST.json','OUTPUT_MANIFEST.json'):
            mf=json.loads(z.read(root+man))['files']
            for rel,dig in mf.items():
                if SHA(z.read(root+rel))!=dig:raise RuntimeError('Y_PARENT_MANIFEST_FAIL '+rel)
            counts[man]=len(mf)
        report=json.loads(z.read(root+'results/report.json'))
        if report.get('n_parent_nodes')!=204 or 'AUTHENTIC_11D_TRAJECTORY' not in report['scientific_decision']:
            raise RuntimeError('NOT_AUTHENTIC_V1853Y_REPORT')
        rows=list(csv.DictReader(io.StringIO(z.read(root+'results/dynamic_frame_oriented_transport.csv').decode())))
        for rel in ('scripts/physical_frame_transport.py','scripts/connection_math.py',
                    'inputs_config.json','inputs/v18_53x_authenticated_server_sendback.zip'):
            dest=stage/rel;dest.parent.mkdir(parents=True,exist_ok=True)
            dest.write_bytes(z.read(root+rel))
    y=import_py(stage/'scripts/physical_frame_transport.py','q550_v1853z_authentic_y')
    x,paths,traces,fixed,anc=y.import_authentic_x(stack,read_traces)
    if len(rows)!=612 or len(paths)!=204:raise RuntimeError('Y_PARENT_INCOMPLETE')
    return y,x,paths,traces,rows,{'y':counts,**anc}

def pin_torque():
    torq=ROOT/'reference/torque_math_v18_53f_authenticated.py'
    if SHA(torq.read_bytes())!=CFG['torque_source_sha256']:
        raise RuntimeError('HISTORICAL_53F_TORQUE_SOURCE_CHANGED')
    return import_py(torq,'q550_v1853z_historical_torque')

def phase_and_source(tag,xi,z,g,torq):
    omega,omdiag=infer_omega(xi,z)
    phase,phdiag=reconstruct_phase(xi,z,g,omega,CFG['kappa_rules'][tag.split('__')[1]],
                                   CFG['parameters'],torq)
    if phdiag['max_source_10d_fit_error']>CFG['source_reconstruction_tolerance']:
        raise RuntimeError('SOURCE_RECONSTRUCTION_FAILED '+tag+': '+str(phdiag))
    if phdiag['max_unit_circle_defect']>CFG['phase_unit_circle_tolerance']:
        raise RuntimeError('PHASE_UNIT_CIRCLE_FAILED '+tag)
    return phase,omega,{**omdiag,**phdiag}

def derivative_probe(xi,z,phase,omega,kappa,torq):
    spline=CubicSpline(xi,z,axis=0,extrapolate=False)
    idx=np.linspace(80,len(xi)-81,15,dtype=int)
    mismatch=[];jerr=[];passive=[]
    for i in idx:
        t=float(xi[i]);src=rhs_material_xi(z[i],phase(t),omega,kappa,CFG['parameters'],torq)
        gd=spline.derivative()(t);mismatch.append(float(np.linalg.norm(src-gd)/max(1,np.linalg.norm(src))))
        J=jacobian_11d(z[i],phase(t),omega,kappa,CFG['parameters'],torq,3e-5)
        J2=jacobian_11d(z[i],phase(t),omega,kappa,CFG['parameters'],torq,1.5e-5)
        jerr.append(float(np.linalg.norm(J-J2)/max(1,np.linalg.norm(J2))))
        passive.append(float(np.linalg.norm(J[:,10])))
    return dict(max_trajectory_rhs_vs_spline_relative=float(max(mismatch)),
                max_jacobian_half_stencil_relative=float(max(jerr)),
                max_passive_tau_jacobian_column=float(max(passive)))

def gauge_check(N,ts,frame,torq,tau=0.22):
    R=frame.at(ts);fixed=pullback(N,ts,frame)
    Q=np.array([gauge(t,frame.t[0],frame.t[-1],tau)[0] for t in ts])
    rotated=np.array([wedge5(q)@n for q,n in zip(Q,N)])
    rec=np.array([wedge5((q@r).T)@v for q,r,v in zip(Q,R,rotated)])
    t=float(ts[len(ts)//2]);qq,qp=gauge(t,frame.t[0],frame.t[-1],tau)
    # Algebraic covariance of the induced generator is tested for an arbitrary
    # skew 5D physical-plus-complement generator at this authentic time.
    baseline=frame.base.omega(t)
    # For this test the intrinsic B is recovered from finite difference of R_z.
    dt=1e-4;R0=frame.at(np.array([t]))[0]
    rp=frame.at(np.array([t+dt]))[0];rm=frame.at(np.array([t-dt]))[0]
    Omega=skew(((rp-rm)/(2*dt))@R0.T)
    Omegag=qq@Omega@qq.T+qp@qq.T
    Q10=wedge5(qq);Q10prime=exterior_generator(qp@qq.T)@Q10
    gen=exterior_generator(Omega);geng=exterior_generator(Omegag)
    generr=float(np.linalg.norm(geng-(Q10@gen@Q10.T+Q10prime@Q10.T)))
    return dict(max_gauge_pullback_error=float(np.max(np.linalg.norm(fixed-rec,axis=1))),
                generator_covariance_error=generr,
                polygon_gauge_recovery_error=float(np.linalg.norm(geodesic_polygon(fixed)[1]-geodesic_polygon(rec)[1])))

def run(mode):
    with ExitStack() as stack:
        torq=pin_torque()
        y,x,paths,traces,yrows,ancestry=parent_import(stack,mode!='preflight')
        print('AUTHENTIC_Y_X_W_V_U_T_CHAIN_PASS',ancestry,flush=True)
        if mode=='preflight':return
        ymap={(r['seed'],r['kappa_rule'],r['epoch'],float(r['lambda_value']),r['lift']):r for r in yrows}
        tags=sorted(traces)
        if mode=='smoke':tags=['original__shared']
        main_rows=[];source_rows=[];gauge_rows=[];refine_rows=[];stores={};max_parent=0.;base45={}
        modes=CFG['lifts']
        for ti,tag in enumerate(tags,1):
            seed,rule=tag.split('__');xi,z,C,g=traces[tag]
            xi=np.asarray(xi,float);z=np.asarray(z,float)
            if len(xi)!=8961 or z.shape!=(len(xi),11):raise RuntimeError('WRONG_11D_SAMPLES')
            phase,omega,diag=phase_and_source(tag,xi,z,g,torq)
            jacdiag=derivative_probe(xi,z,phase,omega,CFG['kappa_rules'][rule],torq)
            if jacdiag['max_jacobian_half_stencil_relative']>CFG['jacobian_stencil_relative_tolerance']:
                raise RuntimeError('JACOBIAN_STENCIL_UNSTABLE '+tag)
            source_rows.append(dict(seed=seed,kappa_rule=rule,**diag,**jacdiag))
            V=x.mg.physical_velocity(z,CFG['parameters']['m_x'])
            frame=y.KinematicFrame(xi,V)
            prox=CubicHermiteSpline(xi,C,g,axis=0,extrapolate=False)
            h=CFG['horizon_xi'];nn=CFG['spin_nodes_smoke'] if mode=='smoke' else CFG['spin_nodes']
            grid=np.linspace(float(xi[0]),float(xi[-1])-h,nn)
            engines={}
            for lift in CFG['horizontal_lifts']:
                B,met=finite_lag_spin(xi,z,phase,omega,CFG['kappa_rules'][rule],
                             CFG['parameters'],torq,grid,h,lift)
                engines[lift]=ComplementFrame(grid,B,frame)
                source_rows.append(dict(seed=seed,kappa_rule=rule,kind='spin',**met,
                                        max_frame_orthogonality=engines[lift].max_orthogonality))
            usekeys=sorted((k for k in paths if k[:2]==(seed,rule)),key=lambda k:(k[2],k[3]))
            if mode=='smoke':usekeys=[k for k in usekeys if k[2]=='preformal' and k[3] in (0.,1.,4.)]
            for key in usekeys:
                _,_,epoch,lam=key
                aa,ww,NN=x.get_path(prox,paths[key])
                for lift in modes:
                    times=aa if lift=='start' else aa+3*ww if lift=='end' else aa+1.5*ww
                    baseline=y.pull_to_initial(NN,times,frame)
                    arc0,pol0,_=geodesic_polygon(baseline)
                    key5=(*key,lift)
                    if lam==0.:base45[seed,rule,epoch,lift,'minimal']=pol0
                    refer=angle_deg(base45[seed,rule,epoch,lift,'minimal'],pol0)
                    target=float(ymap[key5]['dynamic_frame_angle_deg'])
                    err=abs(refer-target);max_parent=max(max_parent,err)
                    if err>CFG['parent_y_replay_tolerance_deg']:
                        raise RuntimeError('Y_PARENT_ANGULAR_REPLAY_MISMATCH '+str(key5)+' '+str(err))
                    for horiz_lift,engine in engines.items():
                        moved=pullback(NN,times,engine)
                        arc,poly,u=geodesic_polygon(moved)
                        if lam==0.:base45[seed,rule,epoch,lift,horiz_lift]=poly
                        dynamic=angle_deg(base45[seed,rule,epoch,lift,horiz_lift],poly)
                        row=dict(seed=seed,kappa_rule=rule,epoch=epoch,lambda_value=lam,
                                 temporal_lift=lift,horizontal_lift=horiz_lift,
                                 minimal_turn_lambda_angle_deg=refer,
                                 model_tangent_lambda_angle_deg=dynamic,
                                 model_vs_minimal_same_lambda_polygon_angle_deg=angle_deg(pol0,poly),
                                 minimal_turn_arclength_rad=arc0,model_tangent_arclength_rad=arc,
                                 model_minus_minimal_arc_rad=arc-arc0,
                                 model_polygon_norm=float(np.linalg.norm(poly)),
                                 model_min_residual_norm=float(np.min(np.linalg.norm(NN,axis=1))),
                                 remains_formal_clock=bool(np.max(aa+3*ww)>CFG['formal_clock_boundary']))
                        main_rows.append(row)
                        if lam in CFG['selected_lambdas'] and mode=='full':
                            name=f'{tag}__{epoch}__lambda{int(round(4*lam)):02d}__{lift}__{horiz_lift}'
                            stores[name+'__u10d']=u;stores[name+'__polygon45d']=poly
                        if lam==CFG['gauge_test_lambda'] and lift in ('start','end') and (mode=='full' or mode=='smoke'):
                            check=gauge_check(NN,times,engine,torq)
                            gauge_rows.append(dict(seed=seed,kappa_rule=rule,epoch=epoch,lambda_value=lam,
                                                   temporal_lift=lift,horizontal_lift=horiz_lift,**check))
                            if max(check.values())>CFG['gauge_tolerance']:
                                raise RuntimeError('MODEL_FRAME_GAUGE_FAILURE '+str(key))
                if lam==CFG['gauge_test_lambda'] and (mode=='full' or mode=='smoke'):
                    # Horizon-sensitivity direct field comparison at selected true times.
                    probe=np.linspace(float(aa[0]),float(aa[-1]),(5 if mode=='smoke' else CFG['horizon_probe_nodes']))
                    B1,di1=finite_lag_spin(xi,z,phase,omega,CFG['kappa_rules'][rule],CFG['parameters'],
                                       torq,probe,h,'fixed_cartesian_velocity')
                    B2,di2=finite_lag_spin(xi,z,phase,omega,CFG['kappa_rules'][rule],CFG['parameters'],
                                       torq,probe,h/2,'fixed_cartesian_velocity')
                    refine_rows.append(dict(seed=seed,kappa_rule=rule,epoch=epoch,lambda_value=lam,
                       h=h,h_half=h/2,max_spin_horizon_absolute_difference=float(np.max(np.linalg.norm(B1-B2,axis=(1,2)))),
                       relative_difference_to_mean_spin=float(np.max(np.linalg.norm(B1-B2,axis=(1,2)))/max(1e-12,float(np.mean(np.linalg.norm(B2,axis=(1,2)))))),
                       mean_spin_norm=float(np.mean(np.linalg.norm(B2,axis=(1,2))))))
            print('V1853Z_TANGENT_CASE',tag,'PATHS',len(usekeys),'RESULT_ROWS',len(main_rows),'Y_REPLAY',f'{max_parent:.3e}',flush=True)
        if mode=='smoke':
            print('V1853Z_AUTHENTIC_SMOKE_PASS',len(main_rows),'GAUGE',len(gauge_rows),'HORIZON',len(refine_rows),'MAX_Y_REPLAY',max_parent,flush=True)
            return
        assert len(main_rows)==204*3*2,(len(main_rows),'incomplete case count')
        assert len(gauge_rows)==6*2*2*2,(len(gauge_rows),'incomplete gauge count')
        assert len(refine_rows)==6*2,(len(refine_rows),'incomplete horizon count')
        csvout(OUT/'model_tangent_complement_transport.csv',main_rows)
        csvout(OUT/'authentic_source_and_jacobian_checks.csv',source_rows)
        csvout(OUT/'model_frame_gauge_covariance.csv',gauge_rows)
        csvout(OUT/'finite_lag_horizon_sensitivity.csv',refine_rows)
        np.savez_compressed(OUT/'selected_tangent_transport_10d_45d.npz',**stores)
        by_lift={lift:max(r['model_tangent_lambda_angle_deg'] for r in main_rows if r['temporal_lift']==lift) for lift in modes}
        max_lift=max(r['model_vs_minimal_same_lambda_polygon_angle_deg'] for r in main_rows)
        min_source=max(r['max_source_10d_fit_error'] for r in source_rows if r.get('max_source_10d_fit_error') is not None)
        report=dict(version=CFG['version'],scientific_decision='AUTHENTIC_TRACE_CONDITIONED_11D_TANGENT_COMPLEMENT_DIAGNOSTIC__HORIZONTAL_LIFT_AND_FINITE_LAG_SENSITIVE__NO_UNIQUE_MATERIAL_CONNECTION',
                    n_historical_paths=204,n_transport_rows=len(main_rows),n_gauge_tests=len(gauge_rows),n_horizon_tests=len(refine_rows),n_new_nonlinear_11d_ivps=0,
                    ancestor_manifest_counts=ancestry,
                    max_parent_y_angle_replay_error_deg=max_parent,
                    max_source_10d_fit_error=min_source,
                    max_gauge_recovery_error=max(r['max_gauge_pullback_error'] for r in gauge_rows),
                    max_gauge_generator_error=max(r['generator_covariance_error'] for r in gauge_rows),
                    max_model_vs_minimal_polygon_angle_deg=max_lift,
                    max_tangent_metric_horizon_spin_difference=max(r['max_spin_horizon_absolute_difference'] for r in refine_rows),
                    max_model_tangent_lambda_angle_by_temporal_lift_deg=by_lift,
                    limitations=[
                    'The model 11D tangent matrix is reconstructed from the historical material-time RHS and the forcing phase inferred from authenticated source samples, not from a new independent forcing dataset.',
                    'A midpoint finite-lag tangent map approximates a short variational flow. Its horizon dependence is explicitly tested but this is not an exact continuous-time variational integrator.',
                    'Mapping a 5D configuration perturbation into 11D state space requires a horizontal lift. Fixed generalized velocity and fixed Cartesian velocity are distinct defensible lifts. Any detected rotation is conditional on that choice.',
                    'Projecting the 5D tangent polar rotation into the complement of the kinetic velocity leaves e_prime unchanged but does not establish a unique material spin.',
                    'The SO(5) gauge covariance test is an algebraic consistency test, not material holonomy.',
                    'All near-return window endpoints remain in formal clock continuation beyond xi=25.6; no 11D material closure is claimed.' ])
        (OUT/'report.json').write_text(json.dumps(report,indent=2)+'\n')
        (OUT/'SCIENTIFIC_REPORT_SERVER.md').write_text(
            '# Q550 v18.53z — 11D tangent-induced complementary rotation\n\n'
            +report['scientific_decision']+'\n\n'
            +f"Authentic paths: 204. Complement transport rows: {len(main_rows)}. Gauge tests: {len(gauge_rows)}. Horizon probes: {len(refine_rows)}.\n\n"
            +f"Maximum parent y replay difference: {max_parent:.6e} deg. Maximum model versus minimal frame polygon angle: {max_lift:.6g} deg.\n\n"
            +'The force-conditioned complement is a lift-dependent kinematic tangent diagnostic, not a certified material holonomy.\n')
        try:
            from plot_transport import plots
            plots(main_rows,OUT)
        except Exception as e:
            raise RuntimeError('FIGURE_GENERATION_FAILURE '+repr(e))
        print('V1853Z',report['scientific_decision'],'ROWS',len(main_rows),'MAX_SOURCE',min_source,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--preflight',action='store_true');p.add_argument('--smoke',action='store_true');p.add_argument('--full',action='store_true')
    a=p.parse_args();run('full' if a.full else 'smoke' if a.smoke else 'preflight')
