#!/usr/bin/env python3
"""Q550 v18.53v fix1: aligned-parent anchor regression and sendback recovery.

No 11D/21D IVPs, no rerun of the 204 metric-continuation nodes.
Historical v18.53u/v18.53t source ZIPs and their manifests remain immutable.
All comments, scientific labels and generated records are in English.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,math,os,subprocess,sys,zipfile
from pathlib import Path
from collections import Counter
import numpy as np
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import minimize_scalar

HERE=Path(__file__).resolve().parent
ORIGINAL='q550_velocity_occupation_metric_continuation_conditioned_hessian_oriented_valley_branch_persistence_v18_53v_bundle'
ROOT=HERE.parent/ORIGINAL
SOURCE=ROOT/'scripts/metric_continuation.py'
CFG_PATH=ROOT/'inputs_config.json'

def require(cond,msg):
    if not cond:raise RuntimeError(msg)

def read_csv(p):
    require(p.is_file(),'MISSING_ORIGINAL_RESULT '+str(p))
    with p.open(newline='',encoding='utf-8') as f:return list(csv.DictReader(f))

def write_csv(p,rows):
    require(bool(rows),'EMPTY_OUTPUT '+str(p))
    with p.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,list(rows[0]));w.writeheader();w.writerows(rows)

def archived_u_best_width(prox,H,a,wlo,whi,grid=25):
    """Functionally matches authenticated v18.53u run_audit.py best_width.

    It is intentionally separate from v18.53v best_width; this allows a
    historical-algorithm check on a common origin anchor and interval.
    """
    def residual(w):return prox(a+3*w)-prox(a)
    def val(w):
        N=residual(w)
        return float(N@H@N)
    ws=np.linspace(wlo,whi,grid)
    vals=np.array([val(w) for w in ws])
    candidates=[(float(vals[0]),float(ws[0])),(float(vals[-1]),float(ws[-1]))]
    for i in range(1,len(ws)-1):
        if vals[i]<=vals[i-1] and vals[i]<=vals[i+1]:
            fit=minimize_scalar(val,bounds=(ws[i-1],ws[i+1]),method='bounded',
                                options={'xatol':1e-11,'maxiter':120})
            candidates.append((float(fit.fun),float(fit.x)))
    if len(candidates)==2:
        j=int(np.argmin(vals))
        fit=minimize_scalar(val,bounds=(ws[max(0,j-1)],ws[min(len(ws)-1,j+1)]),
                            method='bounded',options={'xatol':1e-11})
        candidates.append((float(fit.fun),float(fit.x)))
    return min(candidates,key=lambda p:p[0])[1]

def load_historical_curves(v,epochs):
    """Read manifest-authenticated U NPZ/T center from unchanged ancestor ZIP."""
    parent=ROOT/'inputs'/v.CFG['parent_filename']
    with zipfile.ZipFile(parent) as z:
        require(z.testzip() is None,'BAD_HISTORICAL_U_CRC')
        root=next(iter({n.split('/')[0] for n in z.namelist() if '/' in n}))+'/'
        counts,_=v.verify_zip(z,root)
        Udata=z.read(root+'results/metric_oriented_10d_data.npz')
        Ucfg=json.loads(z.read(root+'inputs_config.json'))
        oldscript=z.read(root+'scripts/run_audit.py').decode()
        require('def best_width' in oldscript and 'grid=25' in oldscript,
                'HISTORICAL_U_WIDTH_ALGORITHM_UNEXPECTED')
        require(Ucfg['curve_points']==141 and Ucfg['width_grid']==25,
                'HISTORICAL_U_CURVE_GRID_UNEXPECTED')
        tdata=z.read(root+'inputs/v18_53t_authenticated_server_sendback.zip')
    with zipfile.ZipFile(io.BytesIO(tdata)) as z:
        require(z.testzip() is None,'BAD_HISTORICAL_T_CRC')
        root=next(iter({n.split('/')[0] for n in z.namelist() if '/' in n}))+'/'
        v.verify_zip(z,root)
        rows=list(csv.DictReader(io.StringIO(z.read(root+'results/valley_geometry_21d.csv').decode())))
        centers={(r['seed'],r['kappa_rule']):(float(r['origin_center']),float(r['width_center']))
                 for r in rows if r['variant']=='on_the_fly_coupled21'}
    require(len(centers)==6,'HISTORICAL_T_SIX_CENTERS_MISSING')
    with np.load(io.BytesIO(Udata),allow_pickle=False) as f:
        out={k:np.array(f[k],copy=True) for k in f.files if k.endswith(('__N10','__origins'))}
    return out,centers,Ucfg,counts

def calculate(run_all=True):
    require(ROOT.is_dir(),'ORIGINAL_V18_53V_FOLDER_NOT_FOUND '+str(ROOT))
    require(SOURCE.is_file() and CFG_PATH.is_file(),'ORIGINAL_SOURCE_FILES_MISSING')
    subprocess.run([sys.executable,str(ROOT/'scripts/verify_manifest.py')],cwd=ROOT,check=True)
    sys.path.insert(0,str(ROOT/'scripts'))
    import metric_continuation as v
    data,old,ancestry,ureport=v.load_authentic_parent()
    old_map={(r['seed'],r['kappa_rule'],r['metric']):r for r in old}
    required=['metric_lambda_continuation.csv','metric_lambda_valley_curves.csv',
              'metric_lambda_parent_u_replay.csv','metric_lambda_family_summary.csv',
              'metric_lambda_oriented_45d.npz','metric_continuation_original.png',
              'metric_continuation_autonomous.png']
    require(all((ROOT/'results'/p).is_file() for p in required),
            'INCOMPLETE_EXISTING_204_NODE_RESULTS: ' + ', '.join(p for p in required if not (ROOT/'results'/p).is_file()))
    nodes=read_csv(ROOT/'results/metric_lambda_continuation.csv')
    curves=read_csv(ROOT/'results/metric_lambda_valley_curves.csv')
    rawreplay=read_csv(ROOT/'results/metric_lambda_parent_u_replay.csv')
    family=read_csv(ROOT/'results/metric_lambda_family_summary.csv')
    require(len(nodes)==204 and len(curves)==28764 and len(rawreplay)==36 and len(family)==12,
            f'UNEXPECTED_CASE_COUNTS {len(nodes)}, {len(curves)}, {len(rawreplay)}, {len(family)}')
    lams=sorted(set(float(x['lambda_value']) for x in nodes))
    require(len(lams)==17 and abs(lams[0])<1e-12 and abs(lams[-1]-4)<1e-12,
            'UNEXPECTED_LAMBDA_GRID')
    keys=[(r['seed'],r['kappa_rule'],r['epoch'],float(r['lambda_value'])) for r in nodes]
    require(len(set(keys))==204,'DUPLICATE_LAMBDA_NODE')
    with np.load(ROOT/'results/metric_lambda_oriented_45d.npz',allow_pickle=False) as z:
        require(len(z.files)==204 and all(z[k].shape==(45,) and np.isfinite(z[k]).all() for k in z.files),
                'INCOMPLETE_ORIENTED_45D_POLYGONS')
    if not run_all:
        keys_to_run={('original','shared')};lam_to_run={0.,1.,3.};epoch_to_run={'preformal'}
    else:
        keys_to_run={(x['seed'],x['kappa_rule']) for x in nodes};lam_to_run={0.,1.,3.};epoch_to_run={'preformal','full_formal'}
    chosen={(r['seed'],r['kappa_rule'],r['epoch'],float(r['lambda_value'])):r for r in nodes
            if (r['seed'],r['kappa_rule']) in keys_to_run and float(r['lambda_value']) in lam_to_run
            and r['epoch'] in epoch_to_run}
    require(len(chosen)==(36 if run_all else 3),'MISSING_HISTORICAL_REPLAY_NODES')
    archived,centers,Ucfg,ucounts=load_historical_curves(v, v.CFG['epochs'])
    result=[]
    worst_archive=0.
    bycase={}
    for seed,rule in sorted(keys_to_run):
        xi,z,C,g=data[seed+'__'+rule]
        prox=CubicHermiteSpline(xi,C,g,axis=0,extrapolate=False)
        oldbase=old_map[(seed,rule,'kinetic_exact')]
        a_center=float(oldbase['min_origin'])
        w_center=float(oldbase['min_width'])
        t_center,t_width=centers[(seed,rule)]
        for epoch in sorted(epoch_to_run):
            S=v.mg.occupation_covariance(xi,z,v.CFG['m_x'], 'preformal' if epoch=='preformal' else 'full')
            for lam in sorted(lam_to_run):
                metric='kinetic_exact' if lam==0 else f'occupation_{epoch}_lambda{lam:g}'
                oldrow=old_map[(seed,rule,metric)]
                newrow=chosen[(seed,rule,epoch,lam)]
                G,Gp,H,Hp=v.analytic_metric(S,lam)
                ukey=seed+'__'+rule+'__'+metric
                old_N=archived[ukey+'__N10']
                old_origins=archived[ukey+'__origins']
                require(old_N.shape==(141,10) and old_origins.shape==(141,), 'HISTORICAL_U_N_SHAPE')
                require(abs(old_origins[70]-t_center)<1e-11,'HISTORICAL_U_T_CENTER_DRIFT')
                old_arc=float(v.mg.oriented_polygon(old_N,H)[0])
                archive_gap=old_arc-float(oldrow['curve_arclength'])
                worst_archive=max(worst_archive,abs(archive_gap))
                require(abs(archive_gap)<3e-10,'HISTORICAL_U_ARCHIVED_ARCLENGTH_MISMATCH')
                # Replay U's 25-point width search on the SAME anchors as v's 21-point curve.
                orig=np.linspace(a_center-v.CFG['origin_halfspan'],a_center+v.CFG['origin_halfspan'],141)
                wl=w_center-v.CFG['width_halfband'];wh=w_center+v.CFG['width_halfband']
                widths=[archived_u_best_width(prox,H,float(a),wl,wh,grid=25) for a in orig]
                N=np.array([prox(a+3*w)-prox(a) for a,w in zip(orig,widths)])
                historical_algorithm_same_anchor=float(v.mg.oriented_polygon(N,H)[0])
                v_arc=float(newrow['arclength'])
                raw_archive_delta=v_arc-float(oldrow['curve_arclength'])
                comparable_delta=v_arc-historical_algorithm_same_anchor
                origin_delta=float(newrow['origin'])-float(oldrow['min_origin'])
                width_delta=float(newrow['width'])-float(oldrow['min_width'])
                norm_delta=float(newrow['minimum_norm'])-float(oldrow['min_norm'])
                result.append(dict(seed=seed,kappa_rule=rule,epoch=epoch,lambda_value=lam,
                    old_parent_t_center=t_center,new_v_curve_center=a_center,
                    center_shift=a_center-t_center,
                    archived_u_arc=float(oldrow['curve_arclength']),
                    v_arc=v_arc,raw_arc_difference_noncomparable=raw_archive_delta,
                    historical_u_algorithm_replayed_at_v_anchor=historical_algorithm_same_anchor,
                    matched_anchor_arc_delta=comparable_delta,
                    min_origin_delta=origin_delta,min_width_delta=width_delta,
                    min_norm_delta=norm_delta,archived_u_internal_arc_error=archive_gap))
                print('ALIGNED_U_REPLAY',seed,rule,epoch,lam,
                      'raw_arc_delta',f'{raw_archive_delta:.9g}',
                      'matched_arc_delta',f'{comparable_delta:.9g}',flush=True)
    max_matched=max(abs(x['matched_anchor_arc_delta']) for x in result)
    max_param=max(max(abs(x['min_origin_delta']),abs(x['min_width_delta']),abs(x['min_norm_delta'])) for x in result)
    max_raw=max(abs(x['raw_arc_difference_noncomparable']) for x in result)
    require(max_matched<1e-5,'MATCHED_ANCHOR_U_ARC_BREACH '+str(max_matched))
    require(max_param<v.CFG['relative_u_audit_tolerance'], 'PARENT_U_MINIMUM_REGRESSION_BREACH '+str(max_param))
    print('ALIGNED_PARENT_U_REPLAY_PASS',len(result),'max_matched_arc',max_matched,
          'max_minimum',max_param,'max_unaligned_arc',max_raw,flush=True)
    if not run_all:return None
    # Preserve the original failed-regression CSV; add a new strictly comparable audit.
    write_csv(ROOT/'results/metric_lambda_parent_u_matched_anchor.csv',result)
    all_pos=all(r['positive_hessian']=='True' for r in nodes)
    all_internal=all(r['is_internal']=='True' for r in nodes)
    all_stat=all(r['stationarity_ok']=='True' for r in nodes)
    multi=max([float(r['competing_distance']) for r in nodes if r['independent_search']=='True'] or [0.])
    branch=multi<=v.CFG['max_branch_vs_global_param_diff']
    category=('REGULAR_CONDITIONED_NUMERIC_METRIC_CONTINUATION__NO_MATERIAL_HOLONOMY'
            if all_pos and all_internal and all_stat and branch else
            'METRIC_CONTINUATION_NUMERIC_BRANCH_AMBIGUITY_OR_CONDITIONING_WARNING__NO_MATERIAL_HOLONOMY')
    positive_norm=all(float(r['minimum_norm'])>v.CFG['zero_screen'] for r in nodes)
    require(positive_norm, 'ZERO_OR_UNRESOLVED_MINIMUM_NORM')
    report=dict(version='v18.53v-fix1',scientific_decision=category,
        mode='AUTHENTIC_U_T_POSTPROCESS_NO_NEW_PHYSICAL_IVPS__POST_RUN_ANCHOR_ALIGNMENT_RECOVERY',
        n_cases=len(nodes),n_metric_families=2,n_physical_cases=6,n_lambda_values=17,
        n_geometry_rows=len(curves),n_ancestral_u_replays=len(result),
        parent_manifest_counts=ancestry,
        max_historical_unchanged_U_arclength_internal_error=worst_archive,
        max_original_unaligned_arclength_difference_not_a_regression=max_raw,
        max_parent_u_matched_anchor_arclength_difference=max_matched,
        max_parent_u_minimum_fields_difference=max_param,
        recovered_without_repeating_204_continuation_nodes=True,
        all_interior_minima=all_internal,all_positive_scaled_hessians=all_pos,
        all_stationary_minima=all_stat,independent_branch_screen_pass=branch,
        max_multistart_param_displacement=multi,
        min_scaled_hessian_eigenvalue=min(float(x['hessian_scaled_min_eig']) for x in nodes),
        max_scaled_hessian_condition=max(float(x['hessian_scaled_cond']) for x in nodes),
        max_origin_shift=max(abs(float(x['origin_shift_from_canonical'])) for x in nodes),
        max_oriented_polygon_angle_deg=max(float(x['oriented_polygon_angle_from_zero_deg']) for x in nodes),
        max_implicit_origin_derivative_secant_mismatch=max(abs(float(x['fd_error_d_origin_d_lambda'])) for x in nodes if x.get('fd_error_d_origin_d_lambda') not in (None,'nan')),
        max_implicit_width_derivative_secant_mismatch=max(abs(float(x['fd_error_d_width_d_lambda'])) for x in nodes if x.get('fd_error_d_width_d_lambda') not in (None,'nan')),
        parent_anchor_explanation=(
            'v18.53u archived 141-point curve is centered on the historical v18.53t curve origin; '
            'v18.53v 141-point curve is centered on the v18.53u optimized canonical minimum. '
            'The original absolute arclength difference therefore compares distinct paths. '
            'Fix1 recomputes the historical u width-search algorithm at the v anchor and compares identical observation windows.'),
        limitations=[
            'The canonical kinetic metric is Euclidean in the historical mass-whitened Cartesian embedding.',
            'Velocity-occupation-weighted metrics are model-conditioned sensitivity probes, not fundamental physical metrics.',
            'Analytic Hessians and implicit derivatives refer to the archived CubicHermiteSpline interpolant, not the continuum flow.',
            'The oriented 45D polygon is compared using a convention-dependent SPD symmetric square-root identification.',
            'Far endpoints remain in formal clock continuation beyond xi=25.6.',
            'Regularity of this metric-conditioned numerical branch does not prove unique global continuation or 11D recurrence, material memory or holonomy.'
        ])
    (ROOT/'results/report.json').write_text(json.dumps(report,indent=2)+'\n')
    (ROOT/'results/SCIENTIFIC_REPORT_SERVER.md').write_text(
        '# Q550 v18.53v-fix1 — Authenticated server results and matched-anchor regression\n\n'
        f'Classification: `{category}`\n\n'
        f'- Authenticated metric-continuation nodes reused: {len(nodes)}; new IVPs: 0.\n'
        f'- Historical U minimum replay maximum delta: {max_param:.9g}.\n'
        f'- Raw U/V arclength comparison (different anchors, non-regression): {max_raw:.9g} rad.\n'
        f'- Matched-anchor historical U algorithm replay: {max_matched:.9g} rad.\n'
        f'- Historical archived U arclength internal check: {worst_archive:.9g} rad.\n'
        f'- Minimum scaled Hessian eigenvalue: {report["min_scaled_hessian_eigenvalue"]:.9g}.\n'
        f'- Maximum scaled Hessian condition: {report["max_scaled_hessian_condition"]:.9g}.\n'
        f'- Internal minima: {all_internal}; positive Hessians: {all_pos}; stationary: {all_stat}; '
        f'independent branch screen: {branch}.\n\n'
        'The original `metric_lambda_parent_u_replay.csv` is preserved as a diagnostic of the formerly mismatched anchors. '
        'Use `metric_lambda_parent_u_matched_anchor.csv` for the regression.\n\n'
        'No material holonomy or full 11D return is claimed.\n')
    print('V18_53V_FIX1_POSTPROCESS_COMPLETE',category,'NODES',len(nodes),flush=True)
    subprocess.run([sys.executable,str(ROOT/'scripts/package_sendback.py')],cwd=ROOT,check=True)
    return report

def self_test():
    require(ROOT.is_dir() and SOURCE.is_file(),'ORIGINAL_V18_53V_FOLDER_MISSING')
    subprocess.run([sys.executable,str(ROOT/'scripts/verify_manifest.py')],cwd=ROOT,check=True)
    sys.path.insert(0,str(ROOT/'scripts'))
    import metric_continuation as v
    data,old,ancestry,_=v.load_authentic_parent()
    archived,centers,Ucfg,ucounts=load_historical_curves(v,v.CFG['epochs'])
    mp={(r['seed'],r['kappa_rule'],r['metric']):r for r in old}
    seed,rule='original','shared'
    xi,z,C,g=data[seed+'__'+rule]
    prox=CubicHermiteSpline(xi,C,g,axis=0,extrapolate=False)
    base=mp[(seed,rule,'kinetic_exact')]
    oc=float(base['min_origin']);wc=float(base['min_width'])
    S=v.mg.occupation_covariance(xi,z,v.CFG['m_x'],'preformal')
    for lam in (0.,1.,3.):
        metric='kinetic_exact' if lam==0 else f'occupation_preformal_lambda{lam:g}'
        G,Gp,H,Hp=v.analytic_metric(S,lam)
        Nold=archived[f'{seed}__{rule}__{metric}__N10']
        oldarc=float(v.mg.oriented_polygon(Nold,H)[0])
        require(abs(oldarc-float(mp[(seed,rule,metric)]['curve_arclength']))<3e-10,'U_ARCHIVE_ARC_SELF_TEST')
        vgeom,_=v.path_geometry(prox,H,None,oc,wc,[(oc-.35,oc+.35),(wc-.05,wc+.05)])
        origins=np.linspace(oc-.35,oc+.35,141)
        widths=[archived_u_best_width(prox,H,float(a),wc-.05,wc+.05,25) for a in origins]
        N=np.array([prox(a+3*w)-prox(a) for a,w in zip(origins,widths)])
        aligned=float(v.mg.oriented_polygon(N,H)[0])
        diff=float(vgeom['arclength']-aligned)
        print('SELFTEST_MATCHED_ANCHOR',lam,'U_archive_arc',oldarc,
              'V_anchor_arc',vgeom['arclength'],'matched_delta',diff,flush=True)
        require(abs(diff)<1e-5,'SELFTEST_MATCHED_ANCHOR_ARC_MISMATCH')
    print('SELFTEST_ANCHOR_ALIGNMENT_PASS',flush=True)

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--self-test',action='store_true',help='Checks three real U/V replay cases without touching the outputs')
    p.add_argument('--recover',action='store_true',help='Recover existing complete 204-node results and package the sendback')
    a=p.parse_args()
    require(a.self_test!=a.recover,'Choose exactly one of --self-test or --recover')
    if a.self_test:self_test()
    else:calculate(run_all=True)

if __name__=='__main__':main()
