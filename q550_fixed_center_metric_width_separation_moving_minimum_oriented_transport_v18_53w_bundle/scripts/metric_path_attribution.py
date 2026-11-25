#!/usr/bin/env python3
"""Q550 v18.53w — metric, width-selection, and center-motion attribution.

No new physical IVP. Every physical trace comes from authenticated v18.53v -> u -> t.
The 45D geometric Shapley split is exact as a vector telescoping identity,
not as an additive law for angles or a causal attribution to material rotation.
All science-facing strings and comments are in English.
"""
from __future__ import annotations
import argparse,csv,hashlib,importlib.util,io,itertools,json,math,sys,zipfile
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import minimize_scalar
ROOT=Path(__file__).resolve().parents[1]
CFG=json.loads((ROOT/'inputs_config.json').read_text())
OUT=ROOT/'results'
OUT.mkdir(parents=True,exist_ok=True)
SPEC=importlib.util.spec_from_file_location('authenticated_u_metric_geometry',ROOT/'inputs/u_metric_geometry_authenticated.py')
mg=importlib.util.module_from_spec(SPEC);sys.modules[SPEC.name]=mg;SPEC.loader.exec_module(mg)
PAIRS45=mg.PAIRS45

def sha(x):return hashlib.sha256(x).hexdigest()
def verify_zip(z,prefix=None):
    error=z.testzip()
    if error is not None:raise RuntimeError('ZIP_CRC_ERROR: '+error)
    if prefix is None:
        roots={n.split('/')[0] for n in z.namelist() if '/' in n}
        if len(roots)!=1:raise RuntimeError('AMBIGUOUS_ANCESTRAL_ZIP_ROOT')
        prefix=next(iter(roots))+'/'
    counts={}
    for m in ('SOURCE_MANIFEST.json','OUTPUT_MANIFEST.json'):
        d=json.loads(z.read(prefix+m))['files']
        for name,digest in d.items():
            if sha(z.read(prefix+name))!=digest:raise RuntimeError('MANIFEST_MISMATCH '+m+' '+name)
        counts[m]=len(d)
    return counts,prefix

def authenticate(light=False):
    payload=(ROOT/'inputs'/CFG['parent_filename']).read_bytes()
    if sha(payload)!=CFG['parent_sha256']:raise RuntimeError('WRONG_V_PARENT_ZIP')
    helper=(ROOT/'inputs/u_metric_geometry_authenticated.py').read_bytes()
    if sha(helper)!=CFG['u_helper_sha256']:raise RuntimeError('WRONG_U_HELPER')
    with zipfile.ZipFile(io.BytesIO(payload)) as z:
        vcounts,p=verify_zip(z)
        vreport=json.loads(z.read(p+'results/report.json'))
        if vreport['scientific_decision']!=CFG['parent_classification'] or vreport['n_cases']!=204:
            raise RuntimeError('WRONG_V_SCIENTIFIC_CLASSIFICATION')
        rows=list(csv.DictReader(io.StringIO(z.read(p+'results/metric_lambda_continuation.csv').decode())))
        curves=list(csv.DictReader(io.StringIO(z.read(p+'results/metric_lambda_valley_curves.csv').decode())))
        oriented=np.load(io.BytesIO(z.read(p+'results/metric_lambda_oriented_45d.npz')),allow_pickle=False)
        polygons={name:np.asarray(oriented[name]).copy() for name in oriented.files};oriented.close()
        ublob=z.read(p+'inputs/v18_53u_authenticated_server_sendback.zip')
    if len(rows)!=204 or len(polygons)!=204 or len(curves)!=28764:
        raise RuntimeError('INCOMPLETE_V_PARENT_ROWS')
    archived_paths={}
    for cr in curves:
        k=(cr['seed'],cr['kappa_rule'],cr['epoch'],float(cr['lambda_value']))
        archived_paths.setdefault(k,[]).append((float(cr['origin']),float(cr['optimized_width'])))
    if any(len(v)!=141 for v in archived_paths.values()) or len(archived_paths)!=204:
        raise RuntimeError('BAD_ARCHIVED_PATH_SHAPE')
    with zipfile.ZipFile(io.BytesIO(ublob)) as z:
        ucounts,up=verify_zip(z)
        if sha(z.read(up+'scripts/metric_geometry.py'))!=CFG['u_helper_sha256']:
            raise RuntimeError('MISMATCH_WITH_AUTHENTIC_U_HELPER')
        tblob=z.read(up+'inputs/v18_53t_authenticated_server_sendback.zip')
    with zipfile.ZipFile(io.BytesIO(tblob)) as z:
        tcounts,tp=verify_zip(z)
        if light:return rows,polygons,archived_paths,{},dict(v=vcounts,u=ucounts,t=tcounts)
        with np.load(io.BytesIO(z.read(tp+'results/augmented21_primitives_and_oriented_valleys.npz')),allow_pickle=False) as f:
            traces={}
            for seed in CFG['seeds']:
                for rule in CFG['kappa_rules']:
                    tag=seed+'__'+rule
                    traces[tag]=tuple(np.asarray(f[tag+'__'+k]).copy() for k in ('xi','z21','C21','g21'))
    return rows,polygons,archived_paths,traces,dict(v=vcounts,u=ucounts,t=tcounts)

def scalar_angle(a,b):
    na=np.linalg.norm(a);nb=np.linalg.norm(b)
    if na<1e-14 or nb<1e-14:return float('nan')
    cos=float(np.clip(np.dot(a,b)/(na*nb),-1,1))
    return float(np.degrees(np.arccos(cos)))

def width_minimizer(prox,H,a,bounds,grid=21):
    lo,hi=bounds
    ws=np.linspace(lo,hi,grid)
    # Identical search logic and tolerance to the authenticated v18.53v.
    val=lambda w:float((prox(a+3*w)-prox(a))@H@(prox(a+3*w)-prox(a)))
    vals=np.array([val(w) for w in ws])
    trials=[(vals[0],ws[0]),(vals[-1],ws[-1])]
    inds=[i for i in range(1,grid-1) if vals[i]<=vals[i-1] and vals[i]<=vals[i+1]]
    if not inds:inds=[int(np.argmin(vals[1:-1]))+1]
    for i in inds:
        fit=minimize_scalar(val,bounds=(ws[i-1],ws[i+1]),method='bounded',options={'xatol':2e-11})
        trials.append((float(fit.fun),float(fit.x)))
    _,w=min(trials,key=lambda x:x[0]);return float(w)

def archived_path(prox,points,n_samples):
    # Reuse exactly the v18.53v optimized width samples on the fixed center path.
    # A 41-point smoke samples a nested subset of the 141-point archived curve.
    stride=(len(points)-1)//(n_samples-1)
    if 1+stride*(n_samples-1)!=len(points):
        # Smoke uses the closest evenly spaced indices; full uses all 141.
        ids=np.round(np.linspace(0,len(points)-1,n_samples)).astype(int)
    else:ids=np.arange(0,len(points),stride)
    a=np.array([points[i][0] for i in ids]);w=np.array([points[i][1] for i in ids])
    N=np.array([prox(x+3*y)-prox(x) for x,y in zip(a,w)])
    return a,w,N

def path_from_selection(prox,center,Hselect,w_bounds,n_samples):
    a=np.linspace(center-CFG['origin_halfspan'],center+CFG['origin_halfspan'],n_samples)
    w=np.array([width_minimizer(prox,Hselect,float(x),w_bounds,CFG['width_search_grid']) for x in a])
    N=np.array([prox(x+3*y)-prox(x) for x,y in zip(a,w)])
    if not np.isfinite(N).all():raise RuntimeError('INVALID_CURVE')
    return a,w,N

def four_interactions_shapley(vertices,factors):
    """Permutation-averaged exact telescoping in one common vector space."""
    d={f:np.zeros_like(next(iter(vertices.values()))) for f in factors}
    for order in itertools.permutations(factors):
        have=frozenset()
        for factor in order:
            after=have|{factor}
            d[factor]+=(vertices[after]-vertices[have])/math.factorial(len(factors))
            have=after
    total=vertices[frozenset(factors)]-vertices[frozenset()]
    closure=float(np.linalg.norm(total-sum(d.values())))
    return d,total,closure

def shp_metric(H):return mg.canonical_sqrt(H)

def vector_at_p(prox,p,H):
    a,w=p
    N=prox(a+3*w)-prox(a)
    return shp_metric(H)@N

def figure(rows):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for seed in CFG['seeds']:
        fig,axes=plt.subplots(2,2,figsize=(11,8),constrained_layout=True)
        for epoch in CFG['epochs']:
            tab=[r for r in rows if r['seed']==seed and r['kappa_rule']=='shared' and r['epoch']==epoch]
            lam=[r['lambda_value'] for r in tab]
            axes[0,0].plot(lam,[r['archived_fixed_center_angle_deg'] for r in tab],marker='.',label=epoch)
            axes[0,1].plot(lam,[r['new_moving_center_angle_deg'] for r in tab],marker='.',label=epoch)
            axes[1,0].plot(lam,[r['shapley_metric_fraction'] for r in tab],marker='.',label=epoch)
            axes[1,1].plot(lam,[r['shapley_center_fraction'] for r in tab],marker='.',label=epoch)
        for ax,ylabel in zip(axes.flat,(
            'Fixed-center 45D angle [deg]','Moving-center 45D angle [deg]',
            'Metric fraction of 45D vector displacement', 'Center fraction of 45D vector displacement')):
            ax.set(xlabel='Metric strength λ',ylabel=ylabel);ax.grid(alpha=.25);ax.legend(fontsize=9)
        fig.suptitle(f'Q550 v18.53w: three-factor 45D decomposition ({seed}, shared κ0)')
        fig.savefig(OUT/f'factorial_metric_path_{seed}.png',dpi=130);plt.close(fig)

def csvwrite(p,records):
    if not records:raise ValueError('EMPTY_RECORD_SET '+str(p))
    with p.open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)

def run(mode):
    parent_rows,archived,parent_paths,trace,ancestry=authenticate(light=(mode=='preflight'))
    print('V1853W_ANCESTRY_PASS',ancestry,flush=True)
    if mode=='preflight':return
    idx={(r['seed'],r['kappa_rule'],r['epoch'],float(r['lambda_value'])):r for r in parent_rows}
    chosen=([0.,.25,1.] if mode=='smoke' else CFG['lambdas'])
    seeds=(CFG['seeds'][:1] if mode=='smoke' else CFG['seeds'])
    rules=(CFG['kappa_rules'][:1] if mode=='smoke' else CFG['kappa_rules'])
    epochs=(CFG['epochs'][:1] if mode=='smoke' else CFG['epochs'])
    ns=(41 if mode=='smoke' else CFG['origin_curve_samples'])
    all_rows=[];minimal=[];store={};all_checks=[]
    for seed in seeds:
        for rule in rules:
            tag=seed+'__'+rule
            xi,z,C,g=trace[tag]
            prox=CubicHermiteSpline(xi,C,g,axis=0,extrapolate=False)
            for epoch in epochs:
                S=mg.occupation_covariance(xi,z,CFG['m_x'],'preformal' if epoch=='preformal' else 'full')
                H0=np.eye(10)
                r0=idx[(seed,rule,epoch,0.)]
                p0=np.array([float(r0['origin']),float(r0['width'])])
                w_bounds=(p0[1]-CFG['width_halfband'],p0[1]+CFG['width_halfband'])
                old_path=archived_path(prox,parent_paths[(seed,rule,epoch,0.)],ns)
                base_polygon=mg.oriented_polygon(old_path[2],H0)[1]
                for li,lam in enumerate(chosen):
                    row=idx[(seed,rule,epoch,lam)]
                    p1=np.array([float(row['origin']),float(row['width'])])
                    H1=mg.induced(mg.occupied_metric(S,lam))
                    centers={0:p0[0],1:p1[0]}
                    hs={0:H0,1:H1}
                    # 4 underlying curves, each observed through 2 metrics.
                    paths={}
                    for c in (0,1):
                        for width_rule in (0,1):
                            if c==0 and width_rule==0:paths[c,width_rule]=old_path
                            elif c==1 and abs(centers[1]-centers[0])<1e-12:paths[c,width_rule]=paths[0,width_rule]
                            elif c==0 and width_rule==1:paths[c,width_rule]=archived_path(prox,parent_paths[(seed,rule,epoch,lam)],ns)
                            else:paths[c,width_rule]=path_from_selection(prox,centers[c],hs[width_rule],w_bounds,ns)
                    vertices={};lengths={}
                    for (c,b),(_,widths,N) in paths.items():
                        for metric in (0,1):
                            arc,pol,_=mg.oriented_polygon(N,hs[metric])
                            key=frozenset(f for f,flag in [('center',c),('width',b),('metric',metric)] if flag)
                            vertices[key]=pol;lengths[key]=arc
                    factors=('center','width','metric')
                    contrib,total,closure=four_interactions_shapley(vertices,factors)
                    P000=vertices[frozenset()];P011=vertices[frozenset(('width','metric'))]
                    P111=vertices[frozenset(factors)]
                    P100=vertices[frozenset(('center',))]
                    P010=vertices[frozenset(('width',))]
                    P001=vertices[frozenset(('metric',))]
                    pname=f'{tag}__{epoch}__lam{int(round(lam/.25)):02d}__polygon45'
                    historical=archived[pname]
                    archive_err=float(np.linalg.norm(P011-historical))
                    if mode=='full' and archive_err>CFG['parent_polygon_tolerance']:
                        raise RuntimeError(f'V_PARENT_POLYGON_REPLAY_MISMATCH {archive_err} {pname}')
                    if closure>CFG['closure_tol']:raise RuntimeError('FAIL_45D_FACTORIAL_CLOSURE')
                    denom=float(np.dot(total,total))
                    fractions={k:(float(np.dot(v,total))/denom if denom>1e-20 else float('nan')) for k,v in contrib.items()}
                    # Kinetic-metric map of the 10D residual at the moving minimum.
                    q00=vector_at_p(prox,p0,H0);q10=vector_at_p(prox,p1,H0)
                    q01=vector_at_p(prox,p0,H1);q11=vector_at_p(prox,p1,H1)
                    qvertices={frozenset():q00,frozenset(('center',)):q10,
                               frozenset(('metric',)):q01,frozenset(('center','metric')):q11}
                    qcon,qtotal,qclosure=four_interactions_shapley(qvertices,('center','metric'))
                    if qclosure>CFG['closure_tol']:raise RuntimeError('FAIL_10D_FACTORIAL_CLOSURE')
                    qden=float(np.dot(qtotal,qtotal))
                    qfractions={k:(float(np.dot(v,qtotal))/qden if qden>1e-20 else float('nan')) for k,v in qcon.items()}
                    # Norm comparison includes cross evaluations, not just differing metric values.
                    nr=lambda q:float(np.linalg.norm(q))
                    n00,n10,n01,n11=map(nr,(q00,q10,q01,q11))
                    shift_effect=.5*((n10-n00)+(n11-n01))
                    metric_effect=.5*((n01-n00)+(n11-n10))
                    norm_close=abs(n11-n00-shift_effect-metric_effect)
                    if norm_close>CFG['closure_tol']:raise RuntimeError('FAIL_NORM_SHAPLEY_CLOSURE')
                    record=dict(seed=seed,kappa_rule=rule,epoch=epoch,lambda_value=float(lam),
                                baseline_origin=float(p0[0]),moving_origin=float(p1[0]),
                                baseline_width=float(p0[1]),moving_width=float(p1[1]),
                                origin_displacement=float(p1[0]-p0[0]),width_displacement=float(p1[1]-p0[1]),
                                archived_fixed_center_angle_deg=scalar_angle(P000,P011),
                                new_moving_center_angle_deg=scalar_angle(P000,P111),
                                center_only_angle_deg=scalar_angle(P000,P100),
                                width_only_angle_deg=scalar_angle(P000,P010),
                                metric_only_angle_deg=scalar_angle(P000,P001),
                                moving_vs_fixed_center_deg=scalar_angle(P011,P111),
                                fixed_center_vector_change_norm=float(np.linalg.norm(P011-P000)),
                                full_45d_vector_change_norm=float(np.linalg.norm(total)),
                                shapley_center_norm=float(np.linalg.norm(contrib['center'])),
                                shapley_width_norm=float(np.linalg.norm(contrib['width'])),
                                shapley_metric_norm=float(np.linalg.norm(contrib['metric'])),
                                shapley_center_fraction=fractions['center'],
                                shapley_width_fraction=fractions['width'],
                                shapley_metric_fraction=fractions['metric'],
                                shapley_45d_closure_error=closure,
                                archived_polygon45_replay_error=archive_err,
                                path_arc_canonical=float(lengths[frozenset()]),
                                path_arc_fixed_center_metric_and_width=float(lengths[frozenset(('width','metric'))]),
                                path_arc_moving_center=float(lengths[frozenset(factors)]),
                                q10d_full_angle_deg=scalar_angle(q00,q11),
                                q10d_moving_minimum_canonical_angle_deg=scalar_angle(q00,q10),
                                q10d_fixed_minimum_metric_angle_deg=scalar_angle(q00,q01),
                                q10d_total_change_norm=float(np.linalg.norm(qtotal)),
                                q10d_center_shapley_norm=float(np.linalg.norm(qcon['center'])),
                                q10d_metric_shapley_norm=float(np.linalg.norm(qcon['metric'])),
                                q10d_center_fraction=qfractions['center'],
                                q10d_metric_fraction=qfractions['metric'],
                                q10d_closure_error=qclosure,
                                canonical_minimum_norm_initial=n00,
                                canonical_norm_of_moved_minimum=n10,
                                weighted_norm_of_initial_minimum=n01,
                                weighted_norm_of_moved_minimum=n11,
                                norm_change_center_shapley=shift_effect,
                                norm_change_metric_shapley=metric_effect,
                                scalar_norm_closure_error=norm_close,
                                curve_center_formal=bool(p1[0]+3*p1[1]>CFG['formal_clock_boundary']))
                    all_rows.append(record)
                    label=f'{tag}__{epoch}__lam{int(round(lam/.25)):02d}'
                    if mode=='full':
                        store[label+'__P000']=P000;store[label+'__P011']=P011;store[label+'__P111']=P111
                        for k,v in contrib.items():store[label+'__shapley_'+k]=v
                        store[label+'__N_at_baseline']=prox(p0[0]+3*p0[1])-prox(p0[0])
                        store[label+'__N_at_moved_minimum']=prox(p1[0]+3*p1[1])-prox(p1[0])
                    print('METRIC_PATH_SPLIT',tag,epoch,'lambda',lam,
                          'fixed-center-angle',f'{record["archived_fixed_center_angle_deg"]:.6f}',
                          'moving-center-angle',f'{record["new_moving_center_angle_deg"]:.6f}',
                          'closure',f'{closure:.2e}',flush=True)
    if mode=='smoke':
        print('AUTHENTIC_V1853W_SMOKE_PASS',len(all_rows),'ROWS',flush=True);return
    if len(all_rows)!=204:raise RuntimeError('EXPECTED_204_NODES')
    csvwrite(OUT/'metric_path_factorial_45d.csv',all_rows)
    np.savez_compressed(OUT/'metric_path_factorial_vectors.npz',**store)
    max_align=max(r['archived_polygon45_replay_error'] for r in all_rows)
    maxerr=max(r['shapley_45d_closure_error'] for r in all_rows)
    maxq=max(r['q10d_closure_error'] for r in all_rows)
    maxh=max(r['scalar_norm_closure_error'] for r in all_rows)
    summary=[]
    for seed in CFG['seeds']:
        for rule in CFG['kappa_rules']:
            for epoch in CFG['epochs']:
                group=[r for r in all_rows if (r['seed'],r['kappa_rule'],r['epoch'])==(seed,rule,epoch)]
                subset=[r for r in group if r['lambda_value']>0]
                summary.append(dict(seed=seed,kappa_rule=rule,epoch=epoch,n_nodes=len(group),
                                    max_parent_replay_error=max(r['archived_polygon45_replay_error'] for r in group),
                                    max_fixed_center_angle_deg=max(r['archived_fixed_center_angle_deg'] for r in subset),
                                    max_moving_center_angle_deg=max(r['new_moving_center_angle_deg'] for r in subset),
                                    max_center_only_angle_deg=max(r['center_only_angle_deg'] for r in subset),
                                    max_45d_closure_error=max(r['shapley_45d_closure_error'] for r in group),
                                    max_q10d_closure_error=max(r['q10d_closure_error'] for r in group)))
    csvwrite(OUT/'metric_path_family_summary.csv',summary)
    figure(all_rows)
    decision=('AUTHENTIC_THREE_FACTOR_METRIC_PATH_DECOMPOSITION__NO_MATERIAL_HOLONOMY'
              if max_align<=CFG['parent_polygon_tolerance'] and maxerr<CFG['closure_tol'] and maxq<CFG['closure_tol']
              else 'METRIC_PATH_DECOMPOSITION_AUDIT_WARNING__NO_MATERIAL_HOLONOMY')
    report=dict(version=CFG['version'],scientific_decision=decision,
      mode='AUTHENTIC_V_U_T_POSTPROCESS_NO_NEW_PHYSICAL_IVPS',
      n_cases=len(all_rows),n_metric_families=2,n_kappa_seed_cases=6,n_lambda=17,
      parent_counts=ancestry,max_archived_fixed_center_polygon_replay_error=max_align,
      max_45d_vector_closure_error=maxerr,max_10d_vector_closure_error=maxq,
      max_scalar_norm_closure_error=maxh,
      max_fixed_center_45d_angle_deg=max(r['archived_fixed_center_angle_deg'] for r in all_rows),
      max_moving_center_45d_angle_deg=max(r['new_moving_center_angle_deg'] for r in all_rows),
      max_center_only_45d_angle_deg=max(r['center_only_angle_deg'] for r in all_rows),
      limits=['The archived v18.53v angle uses a fixed origin center; it does not include migration of the global minimum.',
              'Three-factor Shapley attribution is exact for 45D vector differences, but not additive in angles.',
              'The symmetric square-root metric identification defines the shared 45D comparison frame; results are convention-dependent.',
              'Velocity-occupation metrics are conditioned sensitivity probes, not unique canonical kinetic metrics.',
              'Current analysis is on 8961-sample 21D traces and the continuous Hermite interpolation thereof, not a rigorous continuum certificate.',
              'Window endpoints remain within formal continuation beyond xi=25.6; no material recurrence or holonomy is claimed.'])
    (OUT/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    (OUT/'SCIENTIFIC_REPORT_SERVER.md').write_text('# Q550 v18.53w — Server scientific report\n\n'+
      decision+'\n\n'+f"Rows: {len(all_rows)}. Max archived polygon replay defect: {max_align:.3e}. "+
      f"Shapley 45D closure: {maxerr:.3e}. Minimum-state 10D closure: {maxq:.3e}.\n\n"+
      'All attribution fractions are projections of vector differences onto the total difference; signed values may lie outside [0,1]. Angles are not additive.\n')
    print('V1853W',decision,'CASES',len(all_rows),'ARCHIVED_POLYGON_REPLAY_MAX',max_align,'CLOSURE_MAX',maxerr,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--preflight',action='store_true');p.add_argument('--smoke',action='store_true');p.add_argument('--full',action='store_true')
    a=p.parse_args();run('full' if a.full else 'smoke' if a.smoke else 'preflight')
