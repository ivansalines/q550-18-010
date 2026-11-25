#!/usr/bin/env python3
"""Q550 v18.53k — Fine-width SVD transverse-vector geometry.

Genuinely reintegrates the authenticated Q550 model on two historical seeds
and three kappa rules, then varies integration windows without changing the
underlying trajectory within a seed/kappa. Detects SVD direction transport,
possible oriented reversals, and the exact endpoint dependence predicted by
FTC; never labels a bivector-coordinate SVD vector a material axis.
"""
from __future__ import annotations
import argparse,csv,hashlib,importlib.util,io,json,sys,zipfile,traceback
from collections import defaultdict
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];REF=ROOT/'reference';OUT=ROOT/'results';OUT.mkdir(exist_ok=True)
sys.path.insert(0,str(ROOT/'scripts'))
from fine_svd_geometry import to10,metrics,transported_sign_frame,derivative_diagnostic
from three_window_geometry import accumulators,interval,independent_direct_window_integral
from torque_math import sigmoid
CFG=json.loads((ROOT/'inputs_config.json').read_text());K=CFG['kappas'];F=CFG['fine_svd_scan'];JCFG=CFG['three_window']

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(name,obj): (OUT/name).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')
def dump_rows(name,rows):
    if not rows:return
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with (OUT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)

def load_j():
    p=REF/'v18_53j_original_run_audit.py'
    s=importlib.util.spec_from_file_location('q550_v1853k_historical_parent_j',p)
    mod=importlib.util.module_from_spec(s);sys.modules[s.name]=mod;s.loader.exec_module(mod)
    return mod

def verify_j_parent():
    """Authenticate entire j sendback and check review numeric provenance."""
    for name,digest in CFG['v1853k_parent_pins'].items():
        if sha(ROOT/name)!=digest:raise RuntimeError('V1853J_PARENT_RESOURCE_SHA_MISMATCH '+name)
    zpath=ROOT/'inputs/v18_53j_authenticated_server_sendback.zip'
    with zipfile.ZipFile(zpath) as z:
        if z.testzip() is not None:raise RuntimeError('V1853J_ZIP_CRC_FAILURE')
        roots={n.split('/')[0] for n in z.namelist() if '/' in n}
        if len(roots)!=1:raise RuntimeError('V1853J_MULTIPLE_ROOTS')
        prefix=next(iter(roots))+'/'
        ledgers={}
        for mname in ('SOURCE_MANIFEST.json','OUTPUT_MANIFEST.json'):
            m=json.loads(z.read(prefix+mname));files=m['files'];fail=[]
            for rel,expected in files.items():
                if hashlib.sha256(z.read(prefix+rel)).hexdigest()!=expected:fail.append(rel)
            if fail:raise RuntimeError('V1853J_PARENT_MANIFEST_FAIL '+mname+':'+','.join(fail))
            ledgers[mname]=len(files)
        report=json.loads(z.read(prefix+'results/report.json'))
        if report['scientific_decision']!='AUTHENTIC_THREE_WINDOW_SHIFT_WIDTH_SCAN__WINDOW_SENSITIVE_OR_NONCLOSING' or report['actual_completed_Q550_ivps']!=10:
            raise RuntimeError('V1853J_REQUIRED_PHYSICAL_OUTPUT_MISSING')
        coarse=list(csv.DictReader(io.StringIO(z.read(prefix+'results/window_shift_width_triplets.csv').decode())))
        components=list(csv.DictReader(io.StringIO(z.read(prefix+'results/all_residual_bivector_components.csv').decode())))
        if len(coarse)!=120 or len(components)!=810:raise RuntimeError('V1853J_EXPECTED_NUMBER_OF_ROWS')
    review=json.loads((REF/'v18_53j_independent_svd_review_report.json').read_text())
    if review['original_server_sendback_sha256']!=sha(zpath) or review['triplets_reconstructed']!=270:
        raise RuntimeError('INDEPENDENT_SVD_REVIEW_WRONG_PARENT')
    return {'source_manifest_entries':ledgers['SOURCE_MANIFEST.json'],
      'output_manifest_entries':ledgers['OUTPUT_MANIFEST.json'],
      'server_sendback_sha256':sha(zpath),'independent_review_matches_parent':True},coarse,review

def validate_sources(j,coarse,review):
    h=j.hist_module()
    old,old_windows=j.validate_genealogy(h)
    target=review['baseline_original_seed_shared_kappa']
    baseline=next(r for r in coarse if r['seed']=='original' and r['kappa_rule']=='shared' and
        abs(float(r['xi_start'])-12.8)<1e-9 and abs(float(r['width_xi'])-12.8)<1e-9)
    if abs(target['near_closure_score']-float(baseline['cancellation_three_window']))>1e-12:
        raise RuntimeError('PARENT_SVD_REVIEW_RECONCILIATION_FAIL')
    return h,old_windows

def make_row(seed,kn,origin,w,acc,rank_floor,gap_floor):
    bounds=[origin,origin+w,origin+2*w,origin+3*w]
    parts=[]
    for a,b in zip(bounds[:-1],bounds[1:]):
        U,C=interval(acc,a,b);parts.append(to10(U+C))
    M=np.stack(parts)
    m=metrics(M,rank_floor,gap_floor)
    whole=to10(sum(interval(acc,bounds[0],bounds[-1])))
    error=float(np.linalg.norm(m['net']-whole))
    if error>1e-9:raise RuntimeError('NONADDITIVE_THREE_WINDOW_BIVECTORS')
    return {'seed':seed,'kn':kn,'origin':float(origin),'width':float(w),'m':m,
        'xi_start':float(origin),'xi_end':float(origin+3*w),'add_error':error}

def compact_row(row):
    m=row['m']
    d=dict(seed=row['seed'],kappa_rule=row['kn'],origin_xi=row['origin'],width_xi=row['width'],
        xi_start=row['xi_start'],xi_end=row['xi_end'],
        cancellation_score=m['cancellation_score'],R_total_norm=m['net_norm'],
        rank2_net_norm=m['planar_norm'],rank3_net_norm=m['transverse_norm'],
        rank3_share_of_net_squared=m['transverse_share_squared'],
        singular_value1=m['S'][0],singular_value2=m['S'][1],singular_value3=m['S'][2],
        rank3_strength_ratio=m['rank3_strength_ratio'],rank23_gap=m['rank23_gap'],
        direction_defined=m['direction_defined'],
        max_120_deviation_deg=m['max_120_deviation_deg'],equal_norm_cv=m['amplitude_cv'],
        signed_rank3_coefficient=row.get('signed_transverse_coefficient'),
        transport_defined=row.get('transport_ok'),adjacent_abs_v3_cosine=row.get('adjacent_abs_cos'),
        v3_plane_angle_to_anchor_deg=row.get('v3_principal_angle_to_anchor_deg'),
        transverse_oriented_cosine_to_anchor=row.get('transverse_oriented_cosine_to_anchor'),
        dtransverse_dwidth_norm=row.get('dtransverse_dwidth_norm'),
        endpoint_derivative_scaled_error=row.get('boundary_derivative_error'),
        endpoint_fd_derivative_norm=row.get('boundary_derivative_fd_norm'),
        endpoint_predicted_derivative_norm=row.get('boundary_derivative_pred_norm'),
        integral_additivity_error=row['add_error'])
    for label,v in [('net',m['net']),('transverse',m['transverse']),('v3_projector',m['V3'])]:
        for i,component in enumerate(v):d[f'{label}_{i:02d}']=float(component)
    return d

def sign_candidates(rows,zero_floor=1e-8):
    out=[]
    for a,b in zip(rows[:-1],rows[1:]):
        x,y=a['signed_transverse_coefficient'],b['signed_transverse_coefficient']
        if x is None or y is None or min(abs(x),abs(y))<zero_floor:continue
        if x*y>=0:continue
        out.append(dict(seed=a['seed'],kappa_rule=a['kn'],origin_xi=a['origin'],
            width_left=a['width'],width_right=b['width'],left_signed_coefficient=x,right_signed_coefficient=y,
            left_rank3_ratio=a['m']['rank3_strength_ratio'],right_rank3_ratio=b['m']['rank3_strength_ratio'],
            left_v3_neighbor_overlap=a['adjacent_abs_cos'],right_v3_neighbor_overlap=b['adjacent_abs_cos'],
            classification='SAMPLED_SIGN_CHANGE_CANDIDATE_NOT_ROOT_CERTIFICATE'))
    return out

def reconcile_coarse(rows,coarse):
    index={(r['seed'],r['kappa_rule'],round(float(r['xi_start']),8),round(float(r['width_xi']),8)):r for r in coarse}
    compare=[]
    for x in rows:
        key=(x['seed'],x['kn'],round(x['origin'],8),round(x['width'],8))
        if key not in index:continue
        old=index[key];m=x['m']
        absdiff=max(abs(m['net_norm']-float(old['R_total_norm'])),
            abs(m['cancellation_score']-float(old['cancellation_three_window'])),
            abs(m['max_120_deviation_deg']-float(old['max_120_deviation_deg'])))
        compare.append(dict(seed=x['seed'],kappa_rule=x['kn'],origin_xi=x['origin'],width_xi=x['width'],max_abs_difference=absdiff))
    if len(compare)!=len(CFG['kappas'])*2*4*3:raise RuntimeError('INCOMPLETE_J_COARSE_RECONCILIATION')
    if max(r['max_abs_difference'] for r in compare)>F['reference_reconciliation_abs_tolerance']:
        raise RuntimeError('V1853J_COARSE_RECONCILIATION_DRIFT')
    return compare

def coarse_grid_check(acc,rows):
    """Every-other-grid quadrature sensitivity (not solver-independent)."""
    s=np.arange(0,len(acc['xi']),2)
    if s[-1]!=len(acc['xi'])-1:s=np.r_[s,len(acc['xi'])-1]
    # Recover accumulated values at odd grid using even-cumulative: not possible
    # from cumulative values alone without source array. Parent supplies sources
    # separately; checked in run context below.
    return s

def plot(all_rows):
    try:
        import matplotlib;matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:return
    for seed in ('original','autonomous'):
        chosen=sorted([r for r in all_rows if r['seed']==seed and r['kn']=='shared' and abs(r['origin']-14.4)<1e-8],key=lambda a:a['width'])
        w=np.array([r['width'] for r in chosen]);score=[r['m']['cancellation_score'] for r in chosen]
        p=[r['m']['planar_norm'] for r in chosen];tr=[r['m']['transverse_norm'] for r in chosen]
        fig,ax=plt.subplots(figsize=(9,5.5),layout='constrained')
        ax.plot(w,score,marker='.',label='near-closure score')
        ax.axvline(12.8,color='grey',ls='--',label='v18.53j sampled width')
        ax.axhline(.9,color='grey',ls=':')
        ax.set_xlabel('three equal windows, each width Δξ');ax.set_ylabel('cancellation score');ax.grid(alpha=.2)
        ax.set_title(f'Q550 v18.53k fine-width near-closure ({seed}, shared κ₀)');ax.legend()
        fig.savefig(OUT/f'fine_width_near_closure_{seed}.png',dpi=150);plt.close(fig)
        fig,ax=plt.subplots(figsize=(9,5.5),layout='constrained')
        ax.plot(w,tr,marker='.',label='third SVD direction contribution |T|')
        ax.plot(w,p,marker='.',label='best rank-2 contribution |P|')
        ax.axvline(12.8,color='grey',ls='--');ax.set_xlabel('three-window width Δξ')
        ax.set_ylabel('Frobenius norm');ax.set_title(f'Q550 v18.53k transverse residual decomposition ({seed})')
        ax.grid(alpha=.2);ax.legend();fig.savefig(OUT/f'transverse_residual_{seed}.png',dpi=150);plt.close(fig)
        fig,ax=plt.subplots(figsize=(9,5.5),layout='constrained')
        ax.plot(w,[r['v3_principal_angle_to_anchor_deg'] for r in chosen],marker='.',label='unoriented angle of v3 line to 12.8')
        ax.set_xlabel('three-window width Δξ');ax.set_ylabel('principal line angle (deg)')
        ax.set_title(f'Q550 v18.53k v3 line geometry ({seed}); NOT physical-axis rotation')
        ax.grid(alpha=.2);ax.legend();fig.savefig(OUT/f'v3_line_angle_{seed}.png',dpi=150);plt.close(fig)

def main():
    p=argparse.ArgumentParser();p.add_argument('--base',type=Path,default=ROOT.parent)
    p.add_argument('--native-npz',type=Path);p.add_argument('--parent-v1851s',type=Path)
    p.add_argument('--preflight',action='store_true')
    p.add_argument('--require-parent',action='store_true');p.add_argument('--require-native',action='store_true')
    p.add_argument('--require-dynamic',action='store_true')
    a=p.parse_args();exit_code=0;temp=None
    report=dict(version='v18.53k',scientific_decision='UNCLASSIFIED',actual_Q550_ivps=0,
        actual_refined_ivps=0,parent_status='MISSING',native_status='MISSING',
        output_not_physical_axis_or_closed_orbit=True,
        temporal_range_from_real_history=[12.8,25.6],continued_clock_beyond_xi25_6=True)
    try:
        j=load_j();parent,coarse,review=verify_j_parent();report['j_authenticated']=parent
        h,old_windows=validate_sources(j,coarse,review);report['ancestor_sources']='SHA256_AUTHENTICATED'
        report['preflight']='PASS'
        if a.preflight:
            if a.require_parent or a.require_native or a.require_dynamic:
                raise RuntimeError('PREFLIGHT_CANNOT_SATISFY_REQUIRED_PHYSICAL_EXECUTION')
            report['scientific_decision']='PREFLIGHT_PASS__NO_NEW_IVP_NO_NEW_WIDTH_EVIDENCE'
            return 0
        L=h.loader();par=L.find_parent(a.parent_v1851s);native=a.native_npz or L.find_npz('10d',a.base)
        if par and Path(par).is_file():
            if sha(par)!=CFG['parent_v1851s_sha256']:raise RuntimeError('V1851S_PARENT_SHA_DRIFT')
            report['parent_status']='SHA256_AUTHENTICATED'
        if native and Path(native).is_file():
            if sha(native)!=CFG['native_10d_sha256']:raise RuntimeError('V1852L_NATIVE_SHA_DRIFT')
            report['native_status']='SHA256_AUTHENTICATED'
        if report['parent_status']!='SHA256_AUTHENTICATED' or report['native_status']!='SHA256_AUTHENTICATED':
            report['scientific_decision']='PARTIAL_AUTHENTIC_J_REFERENCE__NEW_PHYSICAL_IVPS_NOT_EXECUTED'
            if a.require_parent or a.require_native or a.require_dynamic:exit_code=3
            return exit_code
        temp,eng,contexts,meta=L.stage_original(par);report['engine_provenance']=meta
        results=[];refined_compare=[];ivps=[];offgrid=[];quad_checks=[];boundary=[];signs=[]
        with np.load(native,allow_pickle=False) as f:
            xi,hcount=j.xi_grid(f['common_xi'])
            if abs(xi[-1]-F['continuation_max_xi'])>1e-8:raise RuntimeError('BAD_MODEL_CONTINUATION_END')
            for seed in ('original','autonomous'):
                ctx=contexts[seed]
                reference_err,clock_err,archived=h.native_53g_reconciliation(f,xi[:hcount],seed,ctx,eng.omega)
                if clock_err>CFG['long_horizon']['material_clock_rate_rel_tolerance']:
                    raise RuntimeError('NATIVE_MATERIAL_TIME_RATE_FAIL')
                z0=np.asarray(archived[0],float)
                for kn in ('shared','original','autonomous'):
                    k=K[kn]
                    z,info=h.solve(eng,ctx,k,xi,z0);report['actual_Q550_ivps']+=1
                    if not np.all(np.isfinite(z)) or np.max(np.abs(z[:,[0,4]]))>40 or np.max(np.abs(z[:,8]))>700:
                        raise RuntimeError('BAD_NEW_INTEGRATED_Q550_STATE')
                    se=h.validated_rhs(eng,ctx,k,xi,z)
                    if se>F.get('source_rhs_rel_tolerance',CFG['three_window']['source_field_rel_tolerance']):
                        raise RuntimeError('SOURCE_RHS_DRIFT')
                    d=h.primary_metrics(xi,z,ctx,k)
                    if d['metric']['charge_range']>CFG['three_window']['charge_abs_drift_tolerance']:
                        raise RuntimeError('CHARGE_DRIFT')
                    if kn=='shared':
                        replay=float(np.max(np.abs(z[:hcount]-archived)/np.maximum(1.,np.abs(archived))))
                        if replay>CFG['three_window']['historical_replay_max_scaled_error']:
                            raise RuntimeError('HISTORICAL_SHARED_REPLAY_FAIL')
                    else:replay=None
                    acc=accumulators(xi,d['t'],d['B_U'],d['B_C'])
                    j.check_old_windows(acc,seed,kn,old_windows)
                    tau_rate=eng.omega/np.asarray(sigmoid(z[:,8]),float)
                    if tau_rate.shape!=xi.shape or not np.all(np.isfinite(tau_rate)):
                        raise RuntimeError('BAD_MATERIAL_CLOCK_RATE')
                    run=[]
                    for o in F['origins_xi']:
                        group=[]
                        for w in F['widths_xi']:
                            if o+3*w>F['continuation_max_xi']+1e-8:raise RuntimeError('OUTSIDE_PREDECLARED_DOMAIN')
                            row=make_row(seed,kn,o,w,acc,F['rank3_relative_amplitude_floor'],F['rank23_gap_floor'])
                            group.append(row)
                        transported_sign_frame(group,F['anchor_width_xi'],F['transport_min_adjacent_abs_cos'])
                        derivative_diagnostic(group,xi,d['B_U'],d['B_C'],tau_rate)
                        signs.extend(sign_candidates(group))
                        for row in group:run.append(row)
                    # Tight recomputation is performed from genuine IVP; compare several windows.
                    if kn in ('original','autonomous'):
                        z2,info2=h.solve(eng,ctx,k,xi,z0,tight=True)
                        report['actual_Q550_ivps']+=1;report['actual_refined_ivps']+=1
                        diff=float(np.max(np.abs(z-z2)/np.maximum(1.,np.abs(z2))))
                        if diff>F['integrator_refinement_max_scaled_error']:
                            raise RuntimeError('IVP_REFINEMENT_FAIL')
                        dd=h.primary_metrics(xi,z2,ctx,k);aa=accumulators(xi,dd['t'],dd['B_U'],dd['B_C'])
                        for o,w in ((12.8,12.8),(14.4,12.8),(14.4,12.775),(14.4,12.825)):
                            oldr=make_row(seed,kn,o,w,acc,F['rank3_relative_amplitude_floor'],F['rank23_gap_floor'])
                            newr=make_row(seed,kn,o,w,aa,F['rank3_relative_amplitude_floor'],F['rank23_gap_floor'])
                            metrics_diff=max(abs(oldr['m'][a]-newr['m'][a]) for a in ('transverse_norm','planar_norm','net_norm','cancellation_score'))
                            refined_compare.append(dict(seed=seed,kappa_rule=kn,origin_xi=o,width_xi=w,max_metric_abs_error=metrics_diff,
                                max_trajectory_scaled_diff=diff))
                    # Every-other sampling: independent integration grid sensitivity on parent IVP.
                    step=np.arange(0,len(xi),2)
                    if step[-1]!=len(xi)-1:step=np.r_[step,len(xi)-1]
                    aa=accumulators(xi[step],d['t'][step],d['B_U'][step],d['B_C'][step])
                    for o,w in ((12.8,12.8),(14.4,12.8),(14.4,12.775),(14.4,12.825)):
                        r0=make_row(seed,kn,o,w,acc,F['rank3_relative_amplitude_floor'],F['rank23_gap_floor'])
                        r2=make_row(seed,kn,o,w,aa,F['rank3_relative_amplitude_floor'],F['rank23_gap_floor'])
                        err=max(abs(r0['m'][v]-r2['m'][v]) for v in ('net_norm','transverse_norm','planar_norm','cancellation_score'))
                        quad_checks.append(dict(seed=seed,kappa_rule=kn,origin_xi=o,width_xi=w,max_absolute_metric_shift=err))
                    # Direct off-grid quadrature against the original underlying torque sources.
                    for o,w in ((12.8,12.8),(14.4,12.8),(14.4,12.775),(14.4,12.825)):
                        ending=o+3*w
                        U,C=interval(acc,o,ending)
                        Ui=independent_direct_window_integral(xi,d['t'],d['B_U'],o,ending)
                        Ci=independent_direct_window_integral(xi,d['t'],d['B_C'],o,ending)
                        err=float(np.linalg.norm((U+C)-(Ui+Ci))/max(1.,np.linalg.norm(U+C)))
                        if err>F['quadrature_max_scaled_error']:
                            raise RuntimeError('OFFGRID_RECONCILIATION_FAILURE')
                        offgrid.append(dict(seed=seed,kappa_rule=kn,origin_xi=o,width_xi=w,relative_offgrid_integral_error=err))
                    ivps.append(dict(seed=seed,kappa_rule=kn,kappa=k,nfev=info['nfev'],
                        historical_replay_scaled_error=replay,source_rhs_rel_error=se,
                        charge_range=d['metric']['charge_range'],material_clock_rate_check=clock_err,
                        original_53g_integral_reconciliation_error=reference_err))
                    results.extend(run)
                    print('V18_53K_DONE',seed,kn,'total_ivps',report['actual_Q550_ivps'],flush=True)
        if report['actual_Q550_ivps']!=F['expected_primary_ivps']+F['expected_refinements'] or len(results)!=2*3*len(F['origins_xi'])*len(F['widths_xi']):
            raise RuntimeError('PREDECLARED_IVP_OR_SAMPLE_COUNT_MISMATCH')
        baseline=reconcile_coarse(results,coarse)
        dump_rows('v1853j_coarse_reconciliation.csv',baseline)
        dump_rows('integrated_ivps.csv',ivps)
        dump_rows('fine_width_svd_geometry.csv',[compact_row(r) for r in results])
        dump_rows('signed_direction_change_candidates.csv',signs)
        dump_rows('integrator_refinement_reconciliation.csv',refined_compare)
        dump_rows('quadrature_decimation_sensitivity.csv',quad_checks)
        dump_rows('independent_offgrid_integral_checks.csv',offgrid)
        # Assess peak only as argmax over a finite predeclared grid, never a continuous optimizer proof.
        peaks=[];groups=defaultdict(list)
        for r in results:groups[(r['seed'],r['kn'],r['origin'])].append(r)
        for (seed,kn,o),g in groups.items():
            peak=max(g,key=lambda r:r['m']['cancellation_score']);lowest=min(g,key=lambda r:r['m']['net_norm'])
            ill=sum(not r['transport_ok'] for r in g)
            peaks.append(dict(seed=seed,kappa_rule=kn,origin_xi=o,
                 sampled_max_score_width_xi=peak['width'],sampled_max_score=peak['m']['cancellation_score'],
                 sampled_min_net_width_xi=lowest['width'],sampled_min_net_norm=lowest['m']['net_norm'],
                 undef_transported_v3_count=ill,points=len(g)))
        dump_rows('sampled_peak_and_minimum_locations.csv',peaks)
        widthshifts={str(o):[r['width'] for r in results if r['seed']=='original' and r['kn']=='shared' and r['origin']==o] for o in F['origins_xi']}
        report['summary']=dict(n_unique_physical_traces=len(ivps),n_integrations=report['actual_Q550_ivps'],
            n_fine_triplets=len(results),n_signed_coefficient_change_candidates=len(signs),
            n_transport_undefined=sum(not r['transport_ok'] for r in results),
            n_strong_near_closures=sum(r['m']['cancellation_score']>=.9 for r in results),
            peak_widths=peaks,
            max_parent_j_reconciliation_error=max(r['max_abs_difference'] for r in baseline),
            max_refinement_metric_error=max((r['max_metric_abs_error'] for r in refined_compare),default=None),
            max_offgrid_error=max(r['relative_offgrid_integral_error'] for r in offgrid),
            max_decimation_metric_difference=max(r['max_absolute_metric_shift'] for r in quad_checks),
            max_endpoint_ftc_error=max((r['boundary_derivative_error'] for r in results if r.get('boundary_derivative_error') is not None),default=None),
            finite_grid_width_scans=widthshifts)
        report['scientific_decision']=('AUTHENTIC_FINE_WIDTH_SCAN__SIGNED_TRANSVERSE_CROSSING_CANDIDATES_ONLY'
                    if signs else 'AUTHENTIC_FINE_WIDTH_SCAN__NO_RESOLVED_SIGN_CROSSING')
        report['scientific_limitations']=[
            'SVD coordinate-space v3 is not a state-space physical third axis or a conserved generator.',
            'Window boundary motion changes the integral by the fundamental theorem of calculus; FTC agreement is numerical provenance, not a dynamical mechanism.',
            'Changes in a continuously aligned signed coefficient are NOT gauge-invariant unless alignment overlaps and rank gaps pass checks.',
            'Nonzero transverse residual can arise from any near-equilateral planar triplet plus a small common out-of-plane bias; not a uniquely physical memory signature.',
            'Clock for xi>25.6 is a formal continuation of the authenticated model source, not independent historical evidence.',
            'No orbital recurrence, material temporal holonomy, zero root, or exact near-closure theorem is certified.']
        plot(results)
    except Exception as e:
        report['scientific_decision']='AUDIT_FAIL_NOT_CERTIFIED';report['failure']=str(e)
        report['traceback']=traceback.format_exc();exit_code=2
    finally:
        if temp is not None:temp.cleanup()
        dump('report.json',report)
        (OUT/'SCIENTIFIC_REPORT_SERVER.md').write_text('\n'.join([
            '# Q550 v18.53k — Fine Width Scan & Transverse SVD Direction Transport','',
            'Decision: `'+report['scientific_decision']+'`','',
            'Completed Q550 IVPs: '+str(report['actual_Q550_ivps']),'',
            'Scientific scope: Frobenius-isometric 10D bivector-coordinate SVD; the line v3 is an observer-dependent decomposition, not a new physical spatial axis.','',
            'FTC endpoint derivative: d/dw ∫_{o}^{o+3w} B(ξ)dτ = 3 B(o+3w) τ\'(o+3w). This is an exact boundary identity, not a physical conservation theorem.','',
            'Beyond ξ=25.6 only a continued external clock of the authenticated dynamical model is available.','',
            'Summary:','',json.dumps(report.get('summary',{}),indent=2),'',
            'Error (if any): '+str(report.get('failure','none')),''])+'\n')
        print('V18_53K',report['scientific_decision'],'IVPS',report['actual_Q550_ivps'],flush=True)
    if a.require_dynamic and not report['scientific_decision'].startswith('AUTHENTIC_FINE_WIDTH_SCAN'):
        return max(exit_code,3)
    if (a.require_parent and report['parent_status']!='SHA256_AUTHENTICATED') or (a.require_native and report['native_status']!='SHA256_AUTHENTICATED'):
        return max(exit_code,3)
    return exit_code

if __name__=='__main__':sys.exit(main())
