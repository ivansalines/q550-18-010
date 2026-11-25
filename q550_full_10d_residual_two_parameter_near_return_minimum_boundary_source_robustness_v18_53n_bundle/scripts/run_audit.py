#!/usr/bin/env python3
"""v18.53n: min of full 10D residual, without observer-plane null projection.

Heavy Q550 physical IVPs only on user's server. No synthetic substitution.
"""
from __future__ import annotations
import argparse,csv,hashlib,importlib.util,io,json,sys,tempfile,traceback,zipfile
from collections import defaultdict
from pathlib import Path
import numpy as np
from full_net_geometry import (check_box,make_eval,local_values,starts_from_grid,
   solve_minimum,unique_minima,knot_clearance,endpoint_source_diff)
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results';OUT.mkdir(parents=True,exist_ok=True)
CFG=json.loads((ROOT/'inputs_config.json').read_text())
PARENT=ROOT/'inputs/v18_53m_fix1_authenticated_server_sendback.zip'

def sha_bytes(s):return hashlib.sha256(s).hexdigest()
def sha(p):return sha_bytes(Path(p).read_bytes())
def output_json(name,x): (OUT/name).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def output_csv(name,rows):
    if not rows:return
    fields=list(dict.fromkeys(q for r in rows for q in r.keys()))
    with (OUT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)

def authenticate_parent_m(tmp):
    if sha(PARENT)!=CFG['parent_m_sha256']:raise RuntimeError('PARENT_M_SHA256_MISMATCH')
    with zipfile.ZipFile(PARENT) as z:
        if z.testzip() is not None:raise RuntimeError('PARENT_M_CRC_FAILED')
        members=z.namelist();roots={p.split('/')[0] for p in members if '/' in p}
        if len(roots)!=1:raise RuntimeError('PARENT_M_MULTIPLE_ROOTS')
        prefix=next(iter(roots))+'/'
        if any(Path(p).is_absolute() or '..' in Path(p).parts for p in members):raise RuntimeError('UNSAFE_PARENT_PATH')
        verified={}
        for mn in ('SOURCE_MANIFEST.json','OUTPUT_MANIFEST.json'):
            manifest=json.loads(z.read(prefix+mn))
            for rel,h in manifest['files'].items():
                if sha_bytes(z.read(prefix+rel))!=h:raise RuntimeError('PARENT_M_MANIFEST_MISMATCH: '+rel)
            verified[mn]=len(manifest['files'])
        report=json.loads(z.read(prefix+'results/report.json'))
        if (report['scientific_decision']!=CFG['parent_m_required_decision'] or
              report['actual_Q550_ivps']!=10):
            raise RuntimeError('PARENT_M_NOT_AUTHENTIC_PHYSICAL_RUN')
        prior=list(csv.DictReader(io.StringIO(z.read(prefix+'results/deduplicated_projected_root_candidates.csv').decode())))
        z.extractall(tmp)
    return Path(tmp)/next(iter(roots)),prior,verified,report

def import_m(root):
    sys.path.insert(0,str(root/'scripts'))
    s=importlib.util.spec_from_file_location('q550_parent_v1853mfix1',root/'scripts/run_audit.py')
    m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m)
    return m

def case_primitive_comparison(xi,primitive,rows,seed,kn):
    at,der=make_eval(xi,primitive)
    out=[]
    for r in rows:
        if r['seed']!=seed or r['kappa_rule']!=kn:continue
        if str(r['solved_projected_root']).lower()!='true':continue
        p=(float(r['origin_xi']),float(r['width_xi']))
        _,_,_,net=local_values(p,at,der)
        measured=float(np.linalg.norm(net));prior=float(r['full_bivector_norm'])
        err=abs(measured-prior)
        if err>CFG['max_parent_m_primitive_replay_abs']:
            raise RuntimeError('PARENT_M_FULL_NET_REPLAY_MISMATCH '+str(err))
        out.append(dict(seed=seed,kappa_rule=kn,anchor_origin_xi=float(r['anchor_origin_xi']),
            parent_origin_xi=p[0],parent_width_xi=p[1],parent_full_net_norm=prior,
            reconstructed_full_net_norm=measured,full_net_replay_abs_error=err))
    if len(out)!=6:raise RuntimeError(f'EXPECTED_SIX_PARENT_PROJECTED_ROOTS: {seed} {kn} {len(out)}')
    return out

def analyze_case(seed,kn,xi,acc,b10,rate,prior_rows,results):
    _,primitive = results['m'].primitive_from_acc(acc)
    base=case_primitive_comparison(xi,primitive,prior_rows,seed,kn)
    results['parent_baseline'].extend(base)
    bounds=[CFG['origin_range_xi'],CFG['width_range_xi']]
    o_range,w_range=bounds
    at,der=make_eval(xi,primitive)
    starts,grid=starts_from_grid(xi,primitive,CFG['grid_origins_xi'],CFG['grid_widths_xi'],
        o_range,w_range,CFG['start_points_per_case'],CFG['min_separation_normalized'])
    for v,o,w in grid:
        results['grid'].append(dict(seed=seed,kappa_rule=kn,origin_xi=o,width_xi=w,
            net_norm=float(np.sqrt(2*v)),objective=v))
    runs=[]
    for start in starts:
        r=solve_minimum(at,der,start,o_range,w_range)
        r['seed']=seed;r['kappa_rule']=kn
        results['multistart'].append(r)
        runs.append(r)
    uniq=unique_minima(runs)
    if not uniq:raise RuntimeError('NO_OPTIMIZATION_RESULTS')
    best=uniq[0]
    if not best['optimizer_success']:raise RuntimeError('BEST_MINIMIZER_FAILED')
    bestp=[best['origin_xi'],best['width_xi']]
    _,gradient,j,net=local_values(bestp,at,der)
    pchip_at,pchip_der=make_eval(xi,primitive,'pchip')
    smoothed=solve_minimum(pchip_at,pchip_der,bestp,o_range,w_range)
    src_diff,src_ref=endpoint_source_diff(xi,primitive,b10,rate,bestp)
    best['start_knot_clearance_xi']=knot_clearance(xi,bestp[0])
    best['end_knot_clearance_xi']=knot_clearance(xi,bestp[0]+3*bestp[1])
    best['material_source_vs_quadrature_derivative_error']=src_diff
    best['material_source_relative_error']=src_diff/max(src_ref,1e-30)
    best['smooth_interpolation_net_norm_at_linear_min']=float(np.linalg.norm(pchip_at(best['end_xi'])-pchip_at(best['origin_xi'])))
    best['pchip_reoptimized_net_norm']=smoothed['net_norm']
    best['pchip_reoptimized_origin_xi']=smoothed['origin_xi']
    best['pchip_reoptimized_width_xi']=smoothed['width_xi']
    best['pchip_minimum_shift_scaled']=float(np.linalg.norm((np.array([smoothed['origin_xi'],smoothed['width_xi']])-bestp)/
       np.array([o_range[1]-o_range[0],w_range[1]-w_range[0]])))
    best['interior_stationarity_gradient_pass']=bool(not best['on_box_boundary'] and best['gradient_norm']<1e-5 and
       min(best['start_knot_clearance_xi'],best['end_knot_clearance_xi'])>CFG['minimum_knot_clearance_xi_for_smooth_jacobian'])
    best['observer_independent_objective']='FULL_10D_FROBENIUS_NET_NORM'
    best['numerical_closure_screen_only']=bool(best['net_norm']<CFG['min_full_norm_material_closure_screen'])
    best['source_clock_after_xi_25p6']=bool(best['end_xi']>CFG['continued_clock_boundary_xi'])
    best['multistart_near_best_count']=sum(abs(r['origin_xi']-bestp[0])<.003 and
        abs(r['width_xi']-bestp[1])<.0008 for r in runs)
    best['parent_best_projected_root_full_net_norm']=min(x['parent_full_net_norm'] for x in base)
    best['gain_over_parent_projected_root']=best['parent_best_projected_root_full_net_norm']-best['net_norm']
    results['best'].append(best.copy())
    for r in uniq:
        results['minima'].append(dict(r,local_minimum_status='PIECEWISE_LINEAR_TRACE_CANDIDATE'))
    return primitive,best,base

def figure(best):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        for seed in ('original','autonomous'):
            rows=[r for r in best if r['seed']==seed]
            fig,ax=plt.subplots(figsize=(8,5))
            for r in rows:
                ax.scatter(r['origin_xi'],r['width_xi'],s=50+400*r['net_norm'],
                           label=f"{r['kappa_rule']}: {r['net_norm']:.5f}")
            ax.set(xlabel='Origin xi',ylabel='Window width w',
                title='v18.53n — full 10D net residual minima ('+seed+')')
            ax.legend(loc='best',fontsize=8);fig.tight_layout()
            fig.savefig(OUT/f'full_10d_minima_{seed}.png',dpi=145);plt.close(fig)
    except Exception as exc:
        (OUT/'PLOT_WARNING.txt').write_text(str(exc))

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--preflight',action='store_true')
    p.add_argument('--base',type=Path,default=ROOT.parent)
    p.add_argument('--parent-v1851s',type=Path)
    p.add_argument('--native-npz',type=Path)
    p.add_argument('--require-parent',action='store_true')
    p.add_argument('--require-native',action='store_true')
    p.add_argument('--require-dynamic',action='store_true')
    a=p.parse_args()
    report=dict(version='v18.53n',scientific_decision='NOT_EVALUATED',actual_Q550_ivps=0,
        refined_Q550_ivps=0,objective_is_unprojected_full_10d=True,
        native_trace_tail_is_formal_continuation_after_xi_25p6=True)
    failure=0
    try:
        with tempfile.TemporaryDirectory(prefix='q550_v1853n_m_') as md:
            mroot,prior_rows,pins,pr=authenticate_parent_m(md)
            report['authenticated_parent_m_manifests']=pins
            m=import_m(mroot)
            with tempfile.TemporaryDirectory(prefix='q550_v1853n_l_') as ld:
                lroot,old_l,lc=m.authenticate_l(ld)
                l=m.import_l(lroot)
                with tempfile.TemporaryDirectory(prefix='q550_v1853n_k_') as kd:
                    kroot,krows,krep=l.verified_parent(kd)
                    k=l.import_k(kroot)
                    j=k.load_j();k.verify_j_parent()
                    h=j.hist_module()
                    if any(not callable(getattr(h,fn,None)) for fn in
                        ('loader','solve','validated_rhs','primary_metrics','native_53g_reconciliation')):
                        raise RuntimeError('INCOMPLETE_PHYSICAL_GENEALOGY')
                    if a.preflight:
                        if a.require_parent or a.require_native or a.require_dynamic:
                            raise RuntimeError('PHYSICAL_FLAGS_FORBIDDEN_IN_PREFLIGHT')
                        report['scientific_decision']='PREFLIGHT_PASS__NO_NEW_Q550_INTEGRATION'
                        report['parent_m_decision']=pr['scientific_decision']
                        print('V1853N_GENEALOGY_PREFLIGHT_PASS',pins,flush=True)
                        return 0
                    loader=h.loader()
                    parent=loader.find_parent(a.parent_v1851s)
                    native=a.native_npz or loader.find_npz('10d',a.base)
                    if not parent or not Path(parent).is_file() or sha(parent)!=CFG['original_source_parent_sha256']:
                        raise RuntimeError('ORIGINAL_PHYSICAL_PARENT_SHA_OR_FILE_MISSING')
                    if not native or not Path(native).is_file() or sha(native)!=CFG['native_10d_sha256']:
                        raise RuntimeError('NATIVE_10D_SHA_OR_FILE_MISSING')
                    from three_window_geometry import accumulators
                    from torque_math import sigmoid
                    staged,engine,contexts,metadata=loader.stage_original(parent)
                    try:
                        data=dict(m=m,grid=[],multistart=[],minima=[],best=[],parent_baseline=[],
                            quadrature=[],refinement=[],ivps=[])
                        with np.load(native,allow_pickle=False) as f:
                            xi,nh=j.xi_grid(f['common_xi'])
                            report['trace_parameter_box']=check_box(xi,CFG['origin_range_xi'],CFG['width_range_xi'])
                            for seed in ('original','autonomous'):
                                ctx=contexts[seed]
                                perr,cerr,history=h.native_53g_reconciliation(f,xi[:nh],seed,ctx,engine.omega)
                                z0=np.asarray(history[0],float)
                                for kn in ('shared','original','autonomous'):
                                    kp=float(k.K[kn]);z,info=h.solve(engine,ctx,kp,xi,z0)
                                    report['actual_Q550_ivps']+=1
                                    rhs_err=h.validated_rhs(engine,ctx,kp,xi,z)
                                    d=h.primary_metrics(xi,z,ctx,kp)
                                    if d['metric']['charge_range']>k.CFG['three_window']['charge_abs_drift_tolerance']:
                                        raise RuntimeError('PHYSICAL_CHARGE_DRIFT')
                                    acc=accumulators(xi,d['t'],d['B_U'],d['B_C'])
                                    ten=np.stack([m.to10(v) for v in d['B_U']+d['B_C']])
                                    rate=engine.omega/np.asarray(sigmoid(z[:,8]),float)
                                    primitive,best,baseline=analyze_case(seed,kn,xi,acc,ten,rate,prior_rows,data)
                                    use_idx=np.arange(0,len(xi),2)
                                    if use_idx[-1]!=len(xi)-1:use_idx=np.r_[use_idx,len(xi)-1]
                                    coarse=accumulators(xi[use_idx],d['t'][use_idx],d['B_U'][use_idx],d['B_C'][use_idx])
                                    _,c2=m.primitive_from_acc(coarse)
                                    at2,d2=make_eval(xi[use_idx],c2)
                                    q=solve_minimum(at2,d2,(best['origin_xi'],best['width_xi']),
                                        CFG['origin_range_xi'],CFG['width_range_xi'])
                                    fv=at2(best['end_xi'])-at2(best['origin_xi'])
                                    data['quadrature'].append(dict(seed=seed,kappa_rule=kn,
                                        full_net_at_primary_min_coarse=float(np.linalg.norm(fv)),
                                        coarse_reoptimized_net_norm=q['net_norm'],
                                        coarse_minimum_origin_xi=q['origin_xi'],
                                        coarse_minimum_width_xi=q['width_xi'],
                                        coarse_minimum_parameter_shift=float(np.linalg.norm(np.array([q['origin_xi'],q['width_xi']])-np.array([best['origin_xi'],best['width_xi']])))))
                                    if kn in ('original','autonomous'):
                                        ztight,_=h.solve(engine,ctx,kp,xi,z0,tight=True)
                                        report['actual_Q550_ivps']+=1
                                        report['refined_Q550_ivps']+=1
                                        zdev=float(np.max(np.abs(z-ztight)/np.maximum(1.,np.abs(ztight))))
                                        if zdev>CFG['max_refined_scaled_trajectory_deviation']:
                                            raise RuntimeError('TIGHT_INTEGRATOR_TRAJECTORY_DIVERGENCE')
                                        dt=h.primary_metrics(xi,ztight,ctx,kp)
                                        a3=accumulators(xi,dt['t'],dt['B_U'],dt['B_C'])
                                        _,c3=m.primitive_from_acc(a3)
                                        at3,der3=make_eval(xi,c3)
                                        refined=solve_minimum(at3,der3,(best['origin_xi'],best['width_xi']),
                                            CFG['origin_range_xi'],CFG['width_range_xi'])
                                        data['refinement'].append(dict(seed=seed,kappa_rule=kn,
                                            max_relative_scaled_trajectory_error=zdev,
                                            refined_net_at_primary_min=float(np.linalg.norm(at3(best['end_xi'])-at3(best['origin_xi']))),
                                            refined_minimum_net_norm=refined['net_norm'],
                                            refined_minimum_origin_xi=refined['origin_xi'],
                                            refined_minimum_width_xi=refined['width_xi'],
                                            refined_minimum_parameter_shift=float(np.linalg.norm(np.array([refined['origin_xi'],refined['width_xi']])-np.array([best['origin_xi'],best['width_xi']])))))
                                    data['ivps'].append(dict(seed=seed,kappa_rule=kn,kappa0=kp,
                                        nfev=info['nfev'],max_rhs_error=rhs_err,charge_drift=d['metric']['charge_range'],
                                        ancestral_53g_reference_error=perr,ancestral_clock_error=cerr))
                                    print('V1853N_PHYSICAL_CASE',seed,kn,'ivps',report['actual_Q550_ivps'],flush=True)
                        if report['actual_Q550_ivps']!=10 or len(data['best'])!=6:
                            raise RuntimeError('INCOMPLETE_FULL_10D_SCAN')
                        for key,fn in [('grid','full_10d_objective_grid.csv'),('multistart','multistart_minima.csv'),
                            ('minima','deduplicated_local_minima.csv'),('best','best_full_10d_minimum_per_case.csv'),
                            ('parent_baseline','parent_m_projected_root_full_net_reconciliation.csv'),
                            ('quadrature','decimated_quadrature_full_minimum_sensitivity.csv'),
                            ('refinement','tight_integrator_full_minimum_sensitivity.csv'),
                            ('ivps','genuine_integrated_ivps.csv')]:output_csv(fn,data[key])
                        figure(data['best'])
                        overall=min(data['best'],key=lambda q:q['net_norm'])
                        closure=any(r['numerical_closure_screen_only'] for r in data['best'])
                        report['summary']=dict(n_grid_rows=len(data['grid']),n_multistart_runs=len(data['multistart']),
                            n_distinct_local_minima=len(data['minima']),min_full_net_norm=overall['net_norm'],
                            best_case=dict(seed=overall['seed'],kappa_rule=overall['kappa_rule'],
                                origin_xi=overall['origin_xi'],width_xi=overall['width_xi'],
                                end_xi=overall['end_xi']),
                            min_parent_root_replay_error=max(r['full_net_replay_abs_error'] for r in data['parent_baseline']),
                            best_full_vs_m_projected_gains=[dict(seed=r['seed'],kappa_rule=r['kappa_rule'],
                                gain=r['gain_over_parent_projected_root']) for r in data['best']],
                            max_pchip_reoptimized_origin_width_scaled_shift=max(r['pchip_minimum_shift_scaled'] for r in data['best']),
                            max_quadrature_minimum_parameter_shift=max(r['coarse_minimum_parameter_shift'] for r in data['quadrature']),
                            max_refined_minimum_parameter_shift=max(r['refined_minimum_parameter_shift'] for r in data['refinement']),
                            all_bests_in_formal_continuation=all(r['source_clock_after_xi_25p6'] for r in data['best']))
                        report['scientific_decision']=('AUTHENTIC_FULL_10D_NUMERICAL_CLOSURE_SCREEN_REQUIRES_INDEPENDENT_PHYSICAL_CERTIFICATION'
                            if closure else 'AUTHENTIC_FULL_10D_NONZERO_NEAR_RETURN_MINIMUM__NO_MATERIAL_CLOSURE')
                        report['scientific_limits']=[
                            'A full 10D residual norm near zero is neither a full 11D state return nor a physical/material holonomy certificate.',
                            'The accumulator depends only on endpoints: the three interior window boundaries are algebraically irrelevant.',
                            'The model clock xi>25.6 is a formal physical-model continuation, not independent observation.',
                            'Piecewise-linear interpolation makes objective locally quadratic but is not a smooth physical Hessian certificate.',
                            'PCHIP is an interpolation sensitivity test, not an alternative physical differential equation.',
                            'No theorem of global minimization is certified by finite grids and multistart optimizations.' ]
                    finally:staged.cleanup()
    except Exception as exc:
        failure=2;report['scientific_decision']='FAILED_OR_INCOMPLETE__NO_SCIENTIFIC_CERTIFICATE'
        report['error']=str(exc);report['traceback']=traceback.format_exc()
    finally:
        output_json('report.json',report)
        (OUT/'SCIENTIFIC_REPORT_SERVER.md').write_text('# Q550 v18.53n — Full 10D Residual Minimum Audit\n\n'
            '**Scientific decision:** `'+report['scientific_decision']+'`\n\n'
            'Actual Q550 physical IVPs: '+str(report['actual_Q550_ivps'])+'\n\n'
            'Full bivector near-return is not a full material/physical orbit return.\n\n'
            '```json\n'+json.dumps(report.get('summary',{}),indent=2)+'\n```\n\n'
            'Error: '+str(report.get('error','none'))+'\n')
        print('V1853N',report['scientific_decision'],'IVPS',report['actual_Q550_ivps'],flush=True)
    if a.require_dynamic and not report['scientific_decision'].startswith('AUTHENTIC_FULL_10D_'):
        return max(3,failure)
    return failure
if __name__=='__main__':sys.exit(main())
