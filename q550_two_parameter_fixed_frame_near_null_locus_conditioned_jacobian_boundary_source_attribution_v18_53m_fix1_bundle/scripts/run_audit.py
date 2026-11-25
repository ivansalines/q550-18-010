#!/usr/bin/env python3
"""Q550 v18.53m: two-parameter fixed 10D-bivector observer-plane zeros.

Physical integrations are staged solely from authenticated historical Q550
source code and native data on the user's server; no absent source is simulated.
"""
from __future__ import annotations
import argparse, csv, hashlib, importlib.util, io, json, sys, tempfile, traceback, zipfile
from pathlib import Path
from collections import defaultdict
import numpy as np
from near_null_geometry import (to10,primitive_from_acc,primitive_slope,projected,
    residual,boundary_jacobian,basis_from_triplet, fixed_grid,spaced_seeds,
    validate_trace_coverage,
    optimize_projected_root, unique_solutions, singular_diagnostics)
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results';OUT.mkdir(exist_ok=True)
CFG=json.loads((ROOT/'inputs_config.json').read_text())
PARENT=ROOT/'inputs/v18_53l_authenticated_server_sendback.zip'
REVIEW=ROOT/'inputs/v18_53l_independent_review.zip'

def sha_bytes(data):return hashlib.sha256(data).hexdigest()
def sha(path):return sha_bytes(Path(path).read_bytes())
def save_json(name,obj):(OUT/name).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')
def csv_rows(name,rows):
    if not rows:return
    keys=list(dict.fromkeys(k for r in rows for k in r.keys()))
    with (OUT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)

def authenticate_l(temp):
    if sha(PARENT)!=CFG['parent_l_sha256']:raise RuntimeError('PARENT_L_SHA256_MISMATCH')
    if sha(REVIEW)!=CFG['parent_l_independent_review_sha256']:raise RuntimeError('PARENT_L_REVIEW_SHA256_MISMATCH')
    with zipfile.ZipFile(PARENT) as z:
        if z.testzip() is not None:raise RuntimeError('PARENT_L_CRC_FAIL')
        roots={n.split('/')[0] for n in z.namelist() if '/' in n}
        if len(roots)!=1:raise RuntimeError('PARENT_L_INVALID_ROOTS')
        prefix=next(iter(roots))+'/'
        for n in z.namelist():
            p=Path(n)
            if p.is_absolute() or '..' in p.parts:raise RuntimeError('UNSAFE_ZIP_PATH')
        counts={}
        for which in ('SOURCE_MANIFEST.json','OUTPUT_MANIFEST.json'):
            files=json.loads(z.read(prefix+which))['files']
            for rel,digest in files.items():
                if sha_bytes(z.read(prefix+rel))!=digest:raise RuntimeError('PARENT_L_MANIFEST_FAIL '+rel)
            counts[which]=len(files)
        report=json.loads(z.read(prefix+'results/report.json'))
        if report['scientific_decision']!=CFG['parent_l_required_decision'] or report['new_actual_Q550_ivps']!=10:
            raise RuntimeError('PARENT_L_MISSING_PHYSICAL_AUDIT')
        table=list(csv.DictReader(io.StringIO(z.read(prefix+'results/three_plane_comparison.csv').decode())))
        if len(table)!=756:raise RuntimeError('PARENT_L_756_ROWS_REQUIRED')
        z.extractall(temp)
        root=Path(temp)/next(iter(roots))
    with zipfile.ZipFile(REVIEW) as zr:
        if zr.testzip() is not None:raise RuntimeError('PARENT_L_INDEPENDENT_REVIEW_CRC_FAIL')
    return root,table,counts

def import_l(path):
    sys.path.insert(0,str(path/'scripts'))
    spec=importlib.util.spec_from_file_location('q550_v1853m_historical_parent_l',path/'scripts/run_audit.py')
    m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m)
    return m

def interp_source(xi,B,tau_rate,x):
    idx=np.clip(int(np.searchsorted(xi,x,side='right')-1),0,len(xi)-2)
    a=np.clip((x-xi[idx])/(xi[idx+1]-xi[idx]),0.,1.)
    return to10((1-a)*B[idx]+a*B[idx+1])*float((1-a)*tau_rate[idx]+a*tau_rate[idx+1])

def baseline_check(xi,primitive,old,seed,kn,anchors):
    out=[]
    for o in anchors:
        for w in (12.775,12.8,12.825):
            key=(seed,kn,round(o,6),round(w,6))
            if key not in old:raise RuntimeError('MISSING_L_BASELINE_ROW '+str(key))
            net=residual(xi,primitive,o,w)
            v=np.array([float(old[key][f'net_{j:02d}']) for j in range(10)])
            err=float(np.max(np.abs(net-v)))
            if err>CFG['max_parent_net_replay_abs_error']:
                raise RuntimeError('PARENT_L_REPLAY_DISAGREEMENT '+str(key)+' '+str(err))
            out.append(dict(seed=seed,kappa_rule=kn,origin_xi=o,width_xi=w,max_component_abs_error=err))
    return out

def analyze_trace(seed,kn,xi,acc,B,tau_rate,previous,allgrids,allcand,allroots,allbounds,allanchors,quadrature):
    origins=CFG['grid_origins_xi'];widths=CFG['grid_widths_xi']
    orange=CFG['origin_range_xi'];wrange=CFG['width_range_xi']
    x,c=primitive_from_acc(acc)
    bases={}
    for ao in CFG['anchor_origin_xi']:
        basis,S,gap=basis_from_triplet(acc,ao,CFG['anchor_width_xi'])
        bases[str(ao)]=basis
        allanchors.append(dict(seed=seed,kappa_rule=kn,anchor_origin_xi=ao,
            anchor_width_xi=CFG['anchor_width_xi'],singular_1=float(S[0]),
            singular_2=float(S[1]),singular_3=float(S[2]),rank2_relative_gap=gap,
            well_conditioned_frame=bool(gap>=CFG['projector_relative_gap_min'])))
    for ao in CFG['anchor_origin_xi']:
        basis=bases[str(ao)]
        grid=fixed_grid(x,c,basis,origins,widths)
        for o,w,n,f1,f2 in grid:
            allgrids.append(dict(seed=seed,kappa_rule=kn,anchor_origin_xi=ao,
                origin_xi=o,width_xi=w,projected_norm=n,component_1=f1,component_2=f2))
        best=spaced_seeds(grid,CFG['root_starts_per_anchor'],orange,wrange)
        candidates=[]
        for origin0,width0 in best:
            r=optimize_projected_root(x,c,basis,(origin0,width0),orange,wrange)
            r.update(seed=seed,kappa_rule=kn,anchor_origin_xi=ao,
                frame_rank2_gap=next(a['rank2_relative_gap'] for a in allanchors if a['seed']==seed and a['kappa_rule']==kn and a['anchor_origin_xi']==ao))
            root=(r['projected_norm']<CFG['root_absolute_residual_tolerance'] and r['optimizer_success']
                   and not r['on_search_boundary'])
            conditioned=(r['sigma_min']>=CFG['root_minimum_sigma'] and
                   r['condition_number']<=CFG['root_condition_number_limit']
                   and r['frame_rank2_gap']>=CFG['projector_relative_gap_min'])
            r.update(solved_projected_root=bool(root),conditioned_projected_root=bool(root and conditioned),
                scientific_scope='OBSERVER_PROJECTED_ZERO_ONLY_NOT_FULL_BIVECTOR_CLOSURE')
            candidates.append(r)
        distinct=unique_solutions(candidates)
        # Avoid multiplying several starts converging to the same identical root.
        for r in distinct:
            point=(r['origin_xi'],r['width_xi'])
            r.update(**{f'at_other_anchor_{str(other).replace(".","p")}_norm':float(np.linalg.norm(projected(x,c,bb,point)))
                        for other,bb in ((f,b) for f,b in [(p,bases[str(p)]) for p in CFG['anchor_origin_xi']])})
            gs=basis.T@primitive_slope(x,c,point[0]);ge=basis.T@primitive_slope(x,c,point[0]+3*point[1])
            actual_o=interp_source(xi,B,tau_rate,point[0]);actual_e=interp_source(xi,B,tau_rate,point[0]+3*point[1])
            Jsrc=np.column_stack((basis.T@(actual_e-actual_o),3*basis.T@actual_e))
            r['material_point_source_vs_trapz_endpoint_jacobian_norm']=float(np.linalg.norm(Jsrc-boundary_jacobian(x,c,basis,point)))
            r['boundary_start_primitive_slope_norm']=float(np.linalg.norm(gs));r['boundary_end_primitive_slope_norm']=float(np.linalg.norm(ge))
            # SO(2) frame rotations must leave a projected root and its condition number intact.
            angle=.637;rot=np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])
            rotated=basis@rot
            r['global_fixed_SO2_gauge_norm_error']=float(abs(np.linalg.norm(projected(x,c,rotated,point))-r['projected_norm']))
            r['global_fixed_SO2_gauge_jacobian_singular_error']=float(np.max(np.abs(np.linalg.svd(boundary_jacobian(x,c,rotated,point),compute_uv=False)-np.linalg.svd(boundary_jacobian(x,c,basis,point),compute_uv=False))))
            allroots.append(r)
        for r in candidates:
            allcand.append(dict(**{k:v for k,v in r.items() if not k.startswith('at_other_anchor_')}))
    return bases,distinct

def derivative_source_checks(seed,kn,xi,acc,B,tau_rate,bases,roots,boundary_rows):
    x,c=primitive_from_acc(acc)
    for ao,b in bases.items():
        eligible=[r for r in roots if r['seed']==seed and r['kappa_rule']==kn and str(r['anchor_origin_xi'])==ao]
        for r in eligible:
            o,w=float(r['origin_xi']),float(r['width_xi']);e=o+3*w
            p=np.array([o,w])
            gs=b.T@primitive_slope(x,c,o);ge=b.T@primitive_slope(x,c,e)
            jac=boundary_jacobian(x,c,b,p)
            jac_alt=np.column_stack((ge-gs,3*ge))
            det_alt=-3*float(np.linalg.det(np.column_stack((gs,ge))))
            boundary_rows.append(dict(seed=seed,kappa_rule=kn,anchor_origin_xi=ao,
                origin_xi=o,width_xi=w,projected_norm=r['projected_norm'],
                exact_endpoint_jacobian_identity_error=float(np.linalg.norm(jac-jac_alt)),
                exact_endpoint_determinant_identity_error=float(abs(np.linalg.det(jac)-det_alt)),
                endpoint_material_source_jacobian_gap=float(np.linalg.norm(jac-np.column_stack((
                   b.T@(interp_source(xi,B,tau_rate,e)-interp_source(xi,B,tau_rate,o)),
                   3*b.T@interp_source(xi,B,tau_rate,e))))),
                start_endpoint_knot_clearance_xi=r['origin_knot_clearance_xi'],
                final_endpoint_knot_clearance_xi=r['end_knot_clearance_xi']))

def save_diagnostics_plot(grids,roots):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for seed in ('original','autonomous'):
        pts=[r for r in grids if r['seed']==seed and r['kappa_rule']=='shared' and r['anchor_origin_xi']==14.4]
        if not pts:continue
        xx=sorted({r['origin_xi'] for r in pts});ww=sorted({r['width_xi'] for r in pts})
        data={(r['origin_xi'],r['width_xi']):r['projected_norm'] for r in pts}
        z=np.array([[data[(o,w)] for w in ww] for o in xx]);fig,ax=plt.subplots(figsize=(8,5.5),layout='constrained')
        mesh=ax.pcolormesh(ww,xx,np.log10(np.maximum(z,1e-12)),shading='nearest')
        fig.colorbar(mesh,ax=ax,label='log10 ||fixed projected net||')
        rr=[r for r in roots if r['seed']==seed and r['kappa_rule']=='shared' and r['anchor_origin_xi']==14.4 and r['solved_projected_root']]
        if rr:ax.scatter([r['width_xi'] for r in rr],[r['origin_xi'] for r in rr],s=85,marker='x',c='white',label='projected-zero candidate');ax.legend()
        ax.set(xlabel='Equal window width w (xi)',ylabel='Origin o (xi)',title='Q550 v18.53m | '+seed+' | fixed anchor o=14.4')
        fig.savefig(OUT/f'fixed_frame_2parameter_map_{seed}.png',dpi=140);plt.close(fig)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--preflight',action='store_true')
    parser.add_argument('--base',type=Path,default=ROOT.parent)
    parser.add_argument('--parent-v1851s',type=Path)
    parser.add_argument('--native-npz',type=Path)
    parser.add_argument('--require-parent',action='store_true')
    parser.add_argument('--require-native',action='store_true')
    parser.add_argument('--require-dynamic',action='store_true')
    a=parser.parse_args()
    report={'version':'v18.53m-fix1','scientific_decision':'NOT_EVALUATED',
            'actual_Q550_ivps':0,'refined_Q550_ivps':0,
            'observer_projected_zero_is_not_full_closure':True,
            'continued_clock_after_xi_25p6':True}
    code=0
    try:
        with tempfile.TemporaryDirectory(prefix='q550_v1853m_stage_') as td:
            lroot,oldrows,counts=authenticate_l(td)
            report['parent_l_verified_manifests']=counts
            report['parent_l_sha256']=CFG['parent_l_sha256']
            l=import_l(lroot)
            with tempfile.TemporaryDirectory(prefix='q550_v1853m_k_') as kd:
                kroot,krows,kr=l.verified_parent(kd)
                k=l.import_k(kroot)
                j=k.load_j();k.verify_j_parent()
                h=j.hist_module()
                required=('loader','solve','validated_rhs','primary_metrics','native_53g_reconciliation')
                if any(not callable(getattr(h,method,None)) for method in required):
                    raise RuntimeError('PARENT_L_DYNAMIC_GENEALOGY_UNAVAILABLE')
                if a.preflight:
                    if a.require_dynamic or a.require_parent or a.require_native:
                        raise RuntimeError('PREFLIGHT_WITH_PHYSICAL_REQUIRED_FLAGS_NOT_ALLOWED')
                    report['scientific_decision']='PREFLIGHT_PASS__NO_NEW_Q550_INTEGRATIONS'
                    print('V1853M_PREFLIGHT_SOURCE_GENEALOGY_PASS',counts,flush=True)
                    return 0
                loader=h.loader()
                par=loader.find_parent(a.parent_v1851s)
                native=a.native_npz or loader.find_npz('10d',a.base)
                report['original_parent_authenticated']=bool(par and Path(par).is_file() and sha(par)==CFG['original_source_parent_sha256'])
                report['native_10d_authenticated']=bool(native and Path(native).is_file() and sha(native)==CFG['native_10d_sha256'])
                if not report['original_parent_authenticated'] or not report['native_10d_authenticated']:
                    raise RuntimeError('ORIGINAL_PHYSICAL_PARENT_OR_NATIVE_DATA_ABSENT_OR_SHA_MISMATCH')
                from three_window_geometry import accumulators
                from torque_math import sigmoid
                l_index={(r['seed'],r['kappa_rule'],round(float(r['origin_xi']),6),round(float(r['width_xi']),6)):r for r in oldrows}
                stage,engine,contexts,metadata=loader.stage_original(par)
                try:
                    grids=[];candidates=[];roots=[];bounds=[];anchors=[];quadrature=[];ivps=[];refinements=[];baseline=[]
                    with np.load(native,allow_pickle=False) as f:
                        xi,hcount=j.xi_grid(f['common_xi'])
                        trace_box=validate_trace_coverage(xi,CFG['origin_range_xi'],CFG['width_range_xi'])
                        report['validated_trace_parameter_box']=trace_box
                        print('V1853M_FIX1_VALIDATED_TRACE_DOMAIN',trace_box,flush=True)
                        for seed in ('original','autonomous'):
                            ctx=contexts[seed]
                            referr,clockerr,archived=h.native_53g_reconciliation(f,xi[:hcount],seed,ctx,engine.omega)
                            z0=np.asarray(archived[0],float)
                            for kn in ('shared','original','autonomous'):
                                kp=float(k.K[kn]);z,info=h.solve(engine,ctx,kp,xi,z0)
                                report['actual_Q550_ivps']+=1
                                rhs_err=h.validated_rhs(engine,ctx,kp,xi,z)
                                d=h.primary_metrics(xi,z,ctx,kp)
                                if d['metric']['charge_range']>k.CFG['three_window']['charge_abs_drift_tolerance']:
                                    raise RuntimeError('MODEL_CHARGE_DRIFT')
                                acc=accumulators(xi,d['t'],d['B_U'],d['B_C'])
                                B=d['B_U']+d['B_C']
                                tau_rate=engine.omega/np.asarray(sigmoid(z[:,8]),float)
                                xx,primitive=primitive_from_acc(acc)
                                baseline.extend(baseline_check(xx,primitive,l_index,seed,kn,CFG['anchor_origin_xi']))
                                bases,newroots=analyze_trace(seed,kn,xi,acc,B,tau_rate,l_index,grids,candidates,roots,bounds,anchors,quadrature)
                                derivative_source_checks(seed,kn,xi,acc,B,tau_rate,bases,roots,bounds)
                                step=np.arange(0,len(xi),2)
                                if step[-1]!=len(xi)-1:step=np.r_[step,len(xi)-1]
                                coarse=accumulators(xi[step],d['t'][step],d['B_U'][step],d['B_C'][step])
                                xc,cc=primitive_from_acc(coarse)
                                case_roots=[r for r in roots if r['seed']==seed and r['kappa_rule']==kn]
                                for r in case_roots:
                                    bb=bases[str(r['anchor_origin_xi'])]
                                    p=np.array([r['origin_xi'],r['width_xi']])
                                    fcoarse=projected(xc,cc,bb,p)
                                    entry=dict(seed=seed,kappa_rule=kn,anchor_origin_xi=r['anchor_origin_xi'],
                                        primary_origin_xi=p[0],primary_width_xi=p[1],
                                        coarse_projected_norm_at_primary_root=float(np.linalg.norm(fcoarse)))
                                    if r['solved_projected_root'] and r['conditioned_projected_root']:
                                        alt=optimize_projected_root(xc,cc,bb,p,CFG['origin_range_xi'],CFG['width_range_xi'])
                                        entry.update(coarse_reoptimized_origin_xi=alt['origin_xi'],
                                           coarse_reoptimized_width_xi=alt['width_xi'],
                                           coarse_root_norm=alt['projected_norm'],
                                           coarse_root_shift_xi=float(np.linalg.norm(np.array([alt['origin_xi'],alt['width_xi']])-p)),
                                           coarse_root_remains_resolved=bool(alt['projected_norm']<CFG['root_absolute_residual_tolerance']))
                                    quadrature.append(entry)
                                if kn in ('original','autonomous'):
                                    zref,_=h.solve(engine,ctx,kp,xi,z0,tight=True)
                                    report['actual_Q550_ivps']+=1;report['refined_Q550_ivps']+=1
                                    z_error=float(np.max(np.abs(z-zref)/np.maximum(1.,np.abs(zref))))
                                    if z_error>CFG['max_refined_scaled_trajectory_deviation']:
                                        raise RuntimeError('TIGHT_INTEGRATOR_REFINEMENT_FAILED')
                                    ref=h.primary_metrics(xi,zref,ctx,kp)
                                    aa=accumulators(xi,ref['t'],ref['B_U'],ref['B_C']);xr,cr=primitive_from_acc(aa)
                                    for r in case_roots:
                                        bb=bases[str(r['anchor_origin_xi'])];pt=np.array([r['origin_xi'],r['width_xi']])
                                        metric=dict(seed=seed,kappa_rule=kn,anchor_origin_xi=r['anchor_origin_xi'],
                                            primary_origin_xi=pt[0],primary_width_xi=pt[1],
                                            max_scaled_integrator_trajectory_error=z_error,
                                            refined_projected_norm_at_primary_root=float(np.linalg.norm(projected(xr,cr,bb,pt))))
                                        if r['solved_projected_root'] and r['conditioned_projected_root']:
                                            found=optimize_projected_root(xr,cr,bb,pt,CFG['origin_range_xi'],CFG['width_range_xi'])
                                            metric.update(refined_root_norm=found['projected_norm'],
                                                refined_root_shift_xi=float(np.linalg.norm(np.array([found['origin_xi'],found['width_xi']])-pt)),
                                                refined_root_remains_resolved=bool(found['projected_norm']<CFG['root_absolute_residual_tolerance']))
                                        refinements.append(metric)
                                ivps.append(dict(seed=seed,kappa_rule=kn,kappa0=kp,nfev=info['nfev'],
                                    maximum_rhs_error=rhs_err,charge_drift=d['metric']['charge_range'],
                                    parent_53g_reconciliation_error=referr,material_clock_error=clockerr))
                                print('V1853M_PHYSICAL_TRACE',seed,kn,'IVPS',report['actual_Q550_ivps'],flush=True)
                    if report['actual_Q550_ivps']!=10 or len(grids)!=2*3*3*len(CFG['grid_origins_xi'])*len(CFG['grid_widths_xi']):
                        raise RuntimeError('UNEXPECTED_PHYSICAL_OR_GRID_COUNT')
                    csv_rows('two_parameter_fixed_frame_grid.csv',grids)
                    csv_rows('multistart_optimizer_candidates.csv',candidates)
                    csv_rows('deduplicated_projected_root_candidates.csv',roots)
                    csv_rows('boundary_source_jacobian_attribution.csv',bounds)
                    csv_rows('anchor_basis_conditioning.csv',anchors)
                    csv_rows('parent_l_net_reconciliation.csv',baseline)
                    csv_rows('decimated_quadrature_root_sensitivity.csv',quadrature)
                    csv_rows('tight_integrator_root_sensitivity.csv',refinements)
                    csv_rows('genuine_integrated_ivps.csv',ivps)
                    valid=[r for r in roots if r['solved_projected_root']]
                    conditional=[r for r in valid if r['conditioned_projected_root']]
                    # Any converged observer-plane zero is only a projected candidate.
                    report['summary']=dict(n_2d_map_samples=len(grids),n_unfiltered_optimizations=len(candidates),
                        n_distinct_optimizer_candidates=len(roots),n_projected_root_candidates=len(valid),
                        n_conditioned_projected_root_candidates=len(conditional),
                        minimum_projected_norm=min(r['projected_norm'] for r in roots),
                        minimum_full_net_norm_at_projected_roots=min((r['full_bivector_norm'] for r in valid),default=None),
                        max_boundary_jacobian_identity_error=max(r['exact_endpoint_jacobian_identity_error'] for r in bounds),
                        max_parent_l_net_reconciliation_error=max(r['max_component_abs_error'] for r in baseline),
                        max_SO2_gauge_norm_error=max(r['global_fixed_SO2_gauge_norm_error'] for r in roots),
                        candidate_condition_numbers=[dict(seed=r['seed'],kappa_rule=r['kappa_rule'],anchor=r['anchor_origin_xi'],
                            origin=r['origin_xi'],width=r['width_xi'],condition=r['condition_number'],
                            full_net_norm=r['full_bivector_norm']) for r in conditional],
                        anchor_selection_is_observer_dependent=True)
                    save_diagnostics_plot(grids,roots)
                    report['scientific_decision']=('AUTHENTIC_TWO_PARAMETER_PROJECTED_ROOT_CANDIDATES_ONLY__NO_MATERIAL_CLOSURE'
                        if valid else 'AUTHENTIC_TWO_PARAMETER_NO_RESOLVED_PROJECTED_ZERO__NEAR_NULL_ONLY')
                    report['scientific_limits']=[
                      'A zero of a chosen fixed 2D bivector-space projection does NOT imply the 10D net bivector is zero.',
                      'The SVD anchor is a predeclared observer constructed from three windows and is tested against alternate fixed anchors.',
                      'The three-window sum is additive; in a fixed plane both Jacobian columns depend entirely on the two endpoint source values.',
                      'A zero need not be robust under choice of global fixed plane or quadrature.',
                      'No full material return, hidden mechanical flip, exact temporal holonomy, or universal near-null law is proved.',
                      'xi>25.6 belongs to the continued Q550 model clock, not an independently measured history.' ]
                finally:
                    stage.cleanup()
    except Exception as exc:
        code=2;report['scientific_decision']='FAILED_OR_INCOMPLETE__NO_SCIENTIFIC_CERTIFICATE'
        report['failure']=str(exc);report['traceback']=traceback.format_exc()
    finally:
        save_json('report.json',report)
        (OUT/'SCIENTIFIC_REPORT_SERVER.md').write_text('# Q550 v18.53m-fix1 — Two-Parameter Fixed-Frame Near-Null Audit\n\n'
           'Decision: `'+report['scientific_decision']+'`\n\n'
           'Actual physical Q550 integrations: '+str(report['actual_Q550_ivps'])+'\n\n'
           'Two-dimensional projected zeros are not full ten-dimensional bivector closure.\n\n'
           '```json\n'+json.dumps(report.get('summary',{}),indent=2)+'\n```\n\n'
           'Error: '+str(report.get('failure','none'))+'\n')
        print('V1853M',report['scientific_decision'],'ACTUAL_IVPS',report['actual_Q550_ivps'],flush=True)
    if a.require_dynamic and not report['scientific_decision'].startswith('AUTHENTIC_TWO_PARAMETER_'):
        return max(3,code)
    return code
if __name__=='__main__':sys.exit(main())
