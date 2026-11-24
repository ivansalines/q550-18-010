#!/usr/bin/env python3
"""Q550 v18.52b — authenticated physical-kappa continuation.

Only the historical nonlinear model generator is used. Each extra kappa value
requires a new nonlinear trajectory and a new 10D tangent source-response IVP.
An observed sign loss is distinct from a solver failure. No finite grid proves
continuous positivity between sampled parameters or times.
"""
from __future__ import annotations
import argparse,csv,hashlib,importlib.util,io,json,os,shutil,sys,tempfile,zipfile
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'
VERSION='v18.52b'
TITLE='Signed-Margin Continuation & Physical Parameter Robustness Radius Audit'
AROOT='q550_dynamically_admissible_cross_channel_turning_robustness_v18_52a_bundle'
ANAME=AROOT+'_server_sendback.zip'
SNAME='q550_parameter_sensitivity_transport_transverse_posture_attribution_v18_51s_bundle_server_sendback.zip'
SEEDS=('original','autonomous')
ANCHORS=('autonomous','midpoint','original')
REQUIRED_NUMERIC_RADIUS=0.0


def sha_file(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(2**20),b''):h.update(b)
    return h.hexdigest()


def module_from(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    obj=importlib.util.module_from_spec(spec)
    sys.modules[name]=obj
    spec.loader.exec_module(obj)
    return obj


def find_input(filename):
    paths=[ROOT/'inputs'/filename,ROOT.parent/filename,
           Path('/root/elab1/hysteresis')/filename]
    for p in paths:
        if p.is_file():return p
    raise FileNotFoundError(f'Missing authenticated input {filename}. Place it in inputs/ or the parent workspace.')


def authenticate_parent(path):
    expected=json.loads((ROOT/'SOURCE_MANIFEST.json').read_text())['parents']['v18_52a_sha256']
    digest=sha_file(path)
    if digest!=expected:raise RuntimeError(f'Parent 52a SHA mismatch: {digest} != {expected}')
    with zipfile.ZipFile(path) as z:
        if z.testzip():raise RuntimeError('Parent 52a ZIP CRC failed')
        files=set(z.namelist())
        root=AROOT+'/'
        report=json.loads(z.read(root+'results/report.json'))
        if report.get('version')!='v18.52a' or not report.get('all_passed'):
            raise RuntimeError('52a parent results are not authenticated/passed')
        src=json.loads(z.read(root+'SOURCE_MANIFEST.json'))
        output=json.loads(z.read(root+'SHA256SUMS_OUTPUT.json'))
        for manifest,key in ((src,'files'),(output,'file_sha256')):
            for n,sha in manifest[key].items():
                if hashlib.sha256(z.read(root+n)).hexdigest()!=sha:
                    raise RuntimeError('52a internal manifest drift: '+n)
        if root+'results/physical_branch_turning_arrays.npz' not in files:
            raise RuntimeError('Missing authentic 52a turning tensors')
        arr=np.load(io.BytesIO(z.read(root+'results/physical_branch_turning_arrays.npz')))
        arrays={k:arr[k].copy() for k in arr.files};arr.close()
    return report,arrays


def anchor_rows(a):
    xi=np.asarray(a['xi'])
    midpoint=(float(a['kappa_original'])+float(a['kappa_autonomous']))/2
    half=(float(a['kappa_original'])-float(a['kappa_autonomous']))/2
    if not (half>0 and len(xi)==513 and abs(xi[-1]-.2)<1e-12):
        raise RuntimeError('Parent parameter/grid mismatch')
    rows=[];series={}
    for seed in SEEDS:
        for kind in ANCHORS:
            kval={'autonomous':midpoint-half,'midpoint':midpoint,'original':midpoint+half}[kind]
            base=f'{seed}_{kind}_'
            mixed=a[base+'rate_mixed'];self_ab=a[base+'rate_self_ab']
            other=a[base+'rate_remaining'];rate=a[base+'rate_total']
            margin=mixed+self_ab
            valid=xi>=.01
            key=f'{seed}:{kind}'
            series[key]={'xi':xi,'margin':margin,'total':rate,'mixed':mixed,'self_ab':self_ab,
                         'oriented_area':a[base+'oriented_area']}
            rows.append({'seed':seed,'label':kind,'kappa':float(kval),'normalized_radius':float(abs(kval-midpoint)/half),
                'sign_side':'center' if kind=='midpoint' else ('lower' if kind=='autonomous' else 'upper'),
                'min_margin':float(np.min(margin[valid])),'min_total':float(np.min(rate[valid])),
                'min_oriented_area':float(np.min(a[base+'oriented_area'][valid])),
                'max_self_ab':float(np.max(self_ab[valid])),
                'qualified_sections':int(valid.sum()),'status':'HISTORICAL_PASS',
                'kind':'authenticated_parent_branch','solver_nfev':None,'tangent_nfev':None,
                'max_recomposition':None,'max_budget_residual':None,
                'source_step_consistency':None,'grid_refinement_margin_max_abs':None})
    return midpoint,half,rows,series


def preflight():
    _,a=authenticate_parent(find_input(ANAME))
    midpoint,half,rows,_=anchor_rows(a)
    a1=[r for r in rows if r['seed']=='original']
    auto=[r for r in rows if r['seed']=='autonomous']
    checks={
      'three_authentic_branch_values': len(set(r['kappa'] for r in a1))==3,
      'symmetric_historical_half_width':abs((max(r['kappa'] for r in a1)-midpoint)-half)<1e-14,
      'both_seed_families_present':len(auto)==3,
      'all_historical_min_margins_positive':all(r['min_margin']>0 for r in rows),
      'all_historical_total_rates_positive':all(r['min_total']>0 for r in rows),
      'all_historical_areas_positive':all(r['min_oriented_area']>0 for r in rows),
      'parent_is_not_interpolated_extra_branch':all(r['kind']=='authenticated_parent_branch' for r in rows),
      'historical_midpoint_balanced':abs(a1[1]['kappa']-midpoint)<1e-14,
      'finite_historical_measurements':all(np.isfinite(r['min_margin']) for r in rows),
    }
    return {'version':VERSION,'mode':'lightweight_preflight','checks':checks,
            'all_passed':all(checks.values()),'n_checks':len(checks),
            'kappa_midpoint':midpoint,'half_historical_band':half,
            'min_historical_signed_margin':min(r['min_margin'] for r in rows),
            'classification':'AUTHENTIC_THREE_BRANCH_ANCHOR_ONLY__NO_EXTENDED_RADIUS_TESTED'}


def stage_parent(z,area,source_s):
    with zipfile.ZipFile(z) as arc:
        for m in arc.infolist():
            # Extraction of trusted, checksum-verified parent entries only.
            if m.is_dir():continue
            name=Path(m.filename)
            if name.is_absolute() or '..' in name.parts or name.parts[0]!=AROOT:
                raise RuntimeError('Invalid parent member path')
            dest=area/name
            dest.parent.mkdir(parents=True,exist_ok=True)
            dest.write_bytes(arc.read(m.filename))
    parent=area/AROOT
    dst=parent/'inputs'/SNAME
    dst.parent.mkdir(exist_ok=True)
    try:os.symlink(source_s.resolve(),dst)
    except OSError:shutil.copyfile(source_s,dst)
    return parent


def write_csv(path,rows):
    if not rows:raise RuntimeError('No results to write: '+path)
    with open(path,'w',encoding='utf-8',newline='') as f:
        cols=list(rows[0]);w=csv.DictWriter(f,fieldnames=cols,extrasaction='ignore')
        w.writeheader();w.writerows(rows)


def classify(summary,budget,tolerance=1e-5):
    """Classify one branch; physical sign loss is not a program failure."""
    margin=budget['rate_mixed']+budget['rate_self_ab']
    total=budget['rate_total'];xi=budget['xi']
    mask=xi>=.01
    mm=float(np.min(margin[mask]));tt=float(np.min(total[mask]));
    ar=float(np.min(budget['oriented_area'][mask]))
    neg=float(np.max(budget['rate_self_ab'][mask]))
    if not all(map(np.isfinite,[mm,tt,ar,neg])):return 'NUMERICAL_UNRESOLVED'
    if mm < -tolerance:return 'MARGIN_SIGN_LOSS'
    if tt < -tolerance:return 'TOTAL_TURNING_SIGN_LOSS'
    if ar < -1e-12:return 'ORIENTED_AREA_SIGN_LOSS'
    if neg > tolerance:return 'SELF_CHANNEL_SIGN_LOSS'
    if abs(mm)<=tolerance or abs(tt)<=tolerance or abs(ar)<=1e-12 or abs(neg)<=tolerance:
        return 'NEAR_ZERO_UNRESOLVED'
    return 'SAMPLED_PASS'


def branch_summary(seed,label,kappa,center,half,meta,budget,series):
    xi=np.asarray(budget['xi']);mask=xi>=.01
    m=budget['rate_mixed']+budget['rate_self_ab'];total=budget['rate_total']
    status=classify(meta,budget)
    row={'seed':seed,'label':label,'kappa':float(kappa),
         'normalized_radius':float(abs(kappa-center)/half),
         'sign_side':'lower' if kappa<center else 'upper',
         'min_margin':float(np.min(m[mask])),'min_total':float(np.min(total[mask])),
         'min_oriented_area':float(np.min(budget['oriented_area'][mask])),
         'max_self_ab':float(np.max(budget['rate_self_ab'][mask])),
         'qualified_sections':int(np.sum(mask)),'status':status,
         'kind':'physical_nonlinear_reintegrated',
         'solver_nfev':meta['solver_nfev'],'tangent_nfev':meta['tangent_nfev'],
         'max_recomposition':meta['max_recomposition'],
         'max_budget_residual':meta['max_budget_residual'],
         'source_step_consistency':meta['source_step_consistency'],
         'grid_refinement_margin_max_abs':meta.get('grid_refinement_margin_max_abs')}
    series[f'{seed}:{label}']={k:np.asarray(budget[k]).copy() for k in
            ('xi','rate_mixed','rate_self_ab','rate_total','oriented_area')}
    series[f'{seed}:{label}']['margin']=m.copy()
    return row


def inspect_real(seed,label,kappa,ctx,xi,a,t,q,p,v,engine,mx,center,half,series,refine=False):
    """Genuine physical IVP + recomputed model-jet source + 10D forced tangent."""
    state=q.run_one(p,v,engine,ctx,kappa,mx,xi,2e-13,2e-14,.000375)
    if not state['sol_success'] or not np.isfinite(state['Y']).all():
        raise ArithmeticError('Physical nonlinear integration failed')
    nodes=np.arange(0,len(xi),2)
    J,bs,step_check,_=a.sample_jet(xi,state['Y'],kappa,ctx,engine,nodes)
    forced,nfev=a.integrate_sources(xi,nodes,J,bs,t)
    summ,bud,plane,source_response=a.analyze_branch(seed,label,xi,state['Y'],state['Q'],(forced,nfev),bs,t)
    bud['xi']=xi
    if refine:
        coarse,_=a.integrate_sources(xi,nodes[::2],J[::2],bs[::2],t)
        # A separate solution at a coarser jet interpolation grid.
        coarse_cart=np.stack([t.j_immersion(y,q0)@resp for y,q0,resp in zip(state['Y'],state['Q'],coarse)])
        cb=a.angular_budget(xi,coarse_cart[:,2:4,:])
        cm=cb['rate_mixed']+cb['rate_self_ab'];fm=bud['rate_mixed']+bud['rate_self_ab']
        sampled=xi>=.01
        refine_err=float(np.max(np.abs(cm[sampled]-fm[sampled])))
    else:refine_err=None
    meta={'solver_nfev':int(state['nfev']),'tangent_nfev':int(nfev),
          'max_recomposition':summ['max_response_channel_additivity'],
          'max_budget_residual':summ['max_turning_budget_residual'],
          'source_step_consistency':float(step_check),
          'grid_refinement_margin_max_abs':refine_err}
    row=branch_summary(seed,label,kappa,center,half,meta,bud,series)
    if (meta['max_recomposition']>1e-8 or meta['max_budget_residual']>1e-7 or
            meta['source_step_consistency']>1e-3):
        row['status']='NUMERICAL_UNRESOLVED'
    if refine_err is not None and refine_err>max(.0001,.1*abs(row['min_margin'])):
        row['status']='NUMERICAL_UNRESOLVED'
    return row


def plot_continuation(rows,center,half):
    fig,axes=plt.subplots(2,2,figsize=(13,9),constrained_layout=True)
    styles={'original':('tab:blue','o'),'autonomous':('tab:orange','s')}
    for seed in SEEDS:
        sub=sorted([r for r in rows if r['seed']==seed and r['status'] not in ('SOLVER_UNRESOLVED','NUMERICAL_UNRESOLVED')],key=lambda r:r['kappa'])
        if not sub:continue
        x=np.array([(r['kappa']-center)/half for r in sub])
        c,marker=styles[seed]
        for ax,name in [(axes[0,0],'min_margin'),(axes[0,1],'min_total'),(axes[1,0],'min_oriented_area'),(axes[1,1],'max_self_ab')]:
            ax.plot(x,[r[name] for r in sub],label=seed,color=c,marker=marker,markersize=4)
    for ax in axes.ravel():
        ax.axvspan(-1,1,color='gray',alpha=.1,label='Historical 52a band')
        ax.axhline(0,color='black',lw=.65)
        ax.grid(alpha=.25);ax.legend(fontsize=8)
        ax.set_xlabel('(kappa0 - historical midpoint) / historical half-width')
    for ax,title in zip(axes.ravel(),['Minimum signed mixed+self margin','Minimum total turning rate',
                                    'Minimum oriented channel area','Maximum self-channel turning contribution']):
        ax.set_title(title)
    fig.suptitle('Q550 v18.52b — Dynamical continuation, sampled margins and search limits')
    fig.savefig(OUT/'physical_kappa_signed_margin_continuation.png',dpi=150)
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=TITLE)
    parser.add_argument('--preflight',action='store_true',help='Authenticate 52a anchors only; no nonlinear solver or tangent integrations')
    parser.add_argument('--levels',type=float,nargs='+',default=[2,4,8,16,32],
                        help='Progressively larger signed radii in historical half-band units; values >1 ascending')
    parser.add_argument('--bisect-steps',type=int,default=4,
                        help='New physical branch integrations to refine a real sign-loss bracket')
    parser.add_argument('--skip-refinement',action='store_true',
                        help='Skip doubled-jet-spacing refinement at outermost tested points')
    args=parser.parse_args()
    OUT.mkdir(exist_ok=True)
    if args.preflight:
        report=preflight();(OUT/'LOCAL_PREFLIGHT.json').write_text(json.dumps(report,indent=2)+'\n')
        print('Q550 v18.52b anchor preflight',sum(report['checks'].values()),'/',report['n_checks'],flush=True)
        if not report['all_passed']:raise SystemExit(2)
        return
    levels=list(args.levels)
    if len(levels)<1 or any(not np.isfinite(v) or v<=1 for v in levels) or levels!=sorted(set(levels)):
        parser.error('levels must be sorted unique finite radii >1')
    if not 0<=args.bisect_steps<=20:parser.error('bisect-steps must be between 0 and 20')
    parent_path=find_input(ANAME);source_s=find_input(SNAME)
    pre=preflight()
    if not pre['all_passed']:raise RuntimeError('Unverified inherited band')
    parent_report,anchors=authenticate_parent(parent_path)
    center,half,rows,series=anchor_rows(anchors);xi=anchors['xi']
    decisions={};failed=[];new_count=0
    with tempfile.TemporaryDirectory(prefix='q550_v1852b_continuation_') as temp:
        parent=stage_parent(parent_path,Path(temp),source_s)
        a=module_from(parent/'scripts/run_audit.py','q550_v1852b_authenticated_52a')
        t=module_from(parent/'reference/parent_v18_51t_run_audit.py','q550_v1852b_authenticated_t')
        print('AUTHENTICATED v18.52a and v18.51s parents. Staging historical 10D physical generator...',flush=True)
        _,_,rblob,srun,states,baseline=a.load_parents()
        staged=t.stage_engine(Path(temp)/'genuine_engine',rblob,srun)
        q,p,v,vr,parts,kappas,lineage,qstates=staged
        _,_,seed_api,observer,engine,unit,drpt,creport,_=parts
        if len(lineage)!=19:raise RuntimeError('Historical engine lineage count mismatch')
        if abs(kappas['midpoint']-center)>1e-12:raise RuntimeError('Kappa midpoint drift')
        mx=float(engine.params.m_x)
        for seed in SEEDS:
            ctx=seed_api.make_context(observer,seed,engine,unit,drpt,creport)
            for side in (-1,1):
                label_side='lower' if side<0 else 'upper'
                last_ok=1.
                boundary={'seed':seed,'side':label_side,'last_sampled_positive_radius':half,
                          'last_sampled_positive_kappa':center+side*half,
                          'first_detected_failure_radius':None,'first_detected_failure_kappa':None,
                          'failure_type':None,'bisection_bracket':None,'coverage':'DISCRETE_SAMPLED_ONLY'}
                previous=1.
                for level in levels:
                    kval=center+side*level*half
                    label=f'{label_side}_r{level:g}'
                    print(f'PHYSICAL BRANCH {seed} {label} kappa={kval:.12f}',flush=True)
                    if not (0<kval<1):
                        boundary['failure_type']='EXPLICIT_KAPPA_DOMAIN_GUARD';boundary['first_detected_failure_radius']=level*half
                        failed.append({'seed':seed,'side':label_side,'radius':level*half,'status':'KAPPA_DOMAIN_GUARD'});break
                    try:
                        row=inspect_real(seed,label,kval,ctx,xi,a,t,q,p,v,engine,mx,center,half,series,
                                         refine=(level==levels[-1] and not args.skip_refinement))
                    except Exception as ex:
                        boundary['failure_type']='SOLVER_UNRESOLVED'
                        boundary['first_detected_failure_radius']=level*half
                        failed.append({'seed':seed,'side':label_side,'radius':level*half,
                                       'status':'SOLVER_UNRESOLVED','error':str(ex)})
                        print('    SOLVER / MODEL UNRESOLVED:',str(ex),flush=True)
                        break
                    rows.append(row);new_count+=1
                    print(f'    signed margin minimum {row["min_margin"]:+.8f} status={row["status"]}',flush=True)
                    if row['status']=='SAMPLED_PASS':
                        last_ok=level;boundary['last_sampled_positive_radius']=level*half
                        boundary['last_sampled_positive_kappa']=kval
                        previous=level;continue
                    boundary['first_detected_failure_radius']=level*half
                    boundary['first_detected_failure_kappa']=kval
                    boundary['failure_type']=row['status']
                    if row['status'] not in ('MARGIN_SIGN_LOSS','TOTAL_TURNING_SIGN_LOSS','ORIENTED_AREA_SIGN_LOSS','SELF_CHANNEL_SIGN_LOSS'):
                        break
                    low=last_ok;high=level;best=level
                    for attempt in range(args.bisect_steps):
                        mid=(low+high)/2;k_mid=center+side*mid*half
                        probe_label=f'{label_side}_bisect{attempt+1}'
                        try:
                            probe=inspect_real(seed,probe_label,k_mid,ctx,xi,a,t,q,p,v,engine,mx,center,half,series)
                        except Exception as ex:
                            failed.append({'seed':seed,'side':label_side,'radius':mid*half,
                                           'status':'BISECTION_SOLVER_UNRESOLVED','error':str(ex)})
                            boundary['failure_type']='BISECTION_UNRESOLVED';break
                        rows.append(probe);new_count+=1
                        if probe['status']=='SAMPLED_PASS':low=mid
                        elif probe['status'] in ('MARGIN_SIGN_LOSS','TOTAL_TURNING_SIGN_LOSS','ORIENTED_AREA_SIGN_LOSS','SELF_CHANNEL_SIGN_LOSS'):
                            high=mid
                        else:
                            boundary['failure_type']='BISECTION_UNRESOLVED';break
                    boundary['bisection_bracket']=[low*half,high*half]
                    boundary['last_sampled_positive_radius']=low*half
                    boundary['last_sampled_positive_kappa']=center+side*low*half
                    break
                else:
                    boundary['failure_type']='NONE_WITHIN_TESTED_RADIUS'
                decisions[seed+'_'+label_side]=boundary
    # Numerical validity and scientific robustness are distinct. A failed physical
    # branch remains a result, while a failed solver remains unresolved.
    checks={
      'preflight_52a_verified':pre['all_passed'],
      'both_seed_families':all(any(r['seed']==seed for r in rows) for seed in SEEDS),
      'anchors_are_actual_52a_rows':sum(r['kind']=='authenticated_parent_branch' for r in rows)==6,
      'all_new_rows_finite':all(np.isfinite(r['min_margin']) for r in rows),
      'new_recompositions_small':all((r['max_recomposition'] or 0.)<1e-8 for r in rows),
      'new_turning_budgets_small':all((r['max_budget_residual'] or 0.)<1e-7 for r in rows),
      'all_new_sources_consistent':all((r['source_step_consistency'] or 0.)<1e-3 for r in rows),
      'radius_assessment_never_claims_interval_proof':True,
    }
    all_success=(all(checks.values())
        and not any(r['status']=='NUMERICAL_UNRESOLVED' for r in rows)
        and not any(e['status'] in ('SOLVER_UNRESOLVED','BISECTION_SOLVER_UNRESOLVED') for e in failed))
    if any(e['status'] in ('SOLVER_UNRESOLVED','BISECTION_SOLVER_UNRESOLVED') for e in failed):
        scientific='INCOMPLETE_DUE_TO_SOLVER_LIMIT'
    elif any(d['failure_type'] in ('MARGIN_SIGN_LOSS','TOTAL_TURNING_SIGN_LOSS','ORIENTED_AREA_SIGN_LOSS','SELF_CHANNEL_SIGN_LOSS') for d in decisions.values()):
        scientific='SAMPLED_SIGN_FAILURE_BRACKET_FOUND'
    elif any(d['failure_type'] in ('BISECTION_UNRESOLVED','NEAR_ZERO_UNRESOLVED','NUMERICAL_UNRESOLVED') for d in decisions.values()):
        scientific='SIGN_MARGIN_NEAR_THRESHOLD_UNRESOLVED'
    else:scientific='SAMPLED_STABILITY_UP_TO_TESTED_LIMIT__NO_RADIUS_BOUNDARY_FOUND'
    report={
      'version':VERSION,'title':TITLE,'all_passed':all_success,
      'numeric_checks':checks,'numerical_audit_status':'PASS' if all_success else 'PARTIAL_OR_UNRESOLVED',
      'scientific_decision':scientific,
      'classification':'PHYSICAL_PARAMETER_KAPPA_DISCRETE_CONTINUATION_NUMERICAL_DIAGNOSTIC_ONLY',
      'kappa_midpoint':center,'historical_half_width':half,'historical_interval':[center-half,center+half],
      'levels_requested':levels,'new_physical_branches_completed':new_count,
      'total_branch_records':len(rows),'boundary_search':decisions,'execution_failures':failed,
      'sources':{'v18_52a_parent_sha256':sha_file(parent_path),'v18_51s_sha256':sha_file(source_s)},
      'qualified_xi_interval':[.01,.2],
      'limitations':['A grid of sampled parameter values is not a continuous robustness certificate.',
        'No rigorous interval bounds in kappa0 or xi are computed.',
        'Solver breakdown is never classified as physical loss of dominance.',
        'Only the model kappa0 is varied; neither other physical parameters nor seed preparations are scanned.',
        'Signed channel decomposition is model-coordinate-specific.',
        'The observed radius is a tested sampling extent, not a proven maximal stability radius.']}
    write_csv(OUT/'physical_parameter_continuation_scan.csv',rows)
    if failed:write_csv(OUT/'solver_and_physical_failures.csv',failed)
    (OUT/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    np.savez_compressed(OUT/'physical_parameter_continuation_curves.npz',
        **{f'{seed}__{label}__{name}':val for pair,s in series.items()
           for seed,label in [pair.split(':',1)] for name,val in s.items()},
           kappa_midpoint=np.array(center),historical_half_width=np.array(half))
    plot_continuation(rows,center,half)
    lines=[f'# Q550 {VERSION}: {TITLE}','','## Numerical status',
        f'- {sum(checks.values())}/{len(checks)} checks passed',
        f'- Scientific decision: `{scientific}`',f'- Additional nonlinear/forced-tangent branches: {new_count}',
        '','## Boundary summary','',
        '| Seed | Direction | Largest positively sampled radius | First failure / limit | Status |',
        '|---|---|---:|---:|---|']
    for d in decisions.values():
        lines.append(f'| {d["seed"]} | {d["side"]} | {d["last_sampled_positive_radius"]:.10g} | '
          f'{d["first_detected_failure_radius"] if d["first_detected_failure_radius"] is not None else "not seen"} | {d["failure_type"]} |')
    lines+=['','## Scientific scope','',
            'The estimated extent is only a set of physically integrated parameter samples. It is not a continuous-kappa0 sign certificate.',
            'Zero sign loss means the search has not located a boundary within the tested range.',
            'Numerical solver failure is a separate unresolved endpoint, never a physical transition.']
    (OUT/'SCIENTIFIC_REPORT_SERVER.md').write_text('\n'.join(lines)+'\n')
    print('Q550 v18.52b numerical:',sum(checks.values()),'/',len(checks),flush=True)
    print('Q550 v18.52b scientific:',scientific,flush=True)
    if not all_success:print('Numerical warning; see report.json and failure ledger.',flush=True)


if __name__=='__main__':main()
