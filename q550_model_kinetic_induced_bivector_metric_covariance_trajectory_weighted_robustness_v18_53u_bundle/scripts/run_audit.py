#!/usr/bin/env python3
"""Q550 v18.53u: canonical kinetic covariance and noncanonical model-velocity metrics.

NO NEW PHYSICAL IVPs. The exact user server v18.53t checkpoint is authenticated.
No material closure or holonomy inferred. All scientific text is in English.
"""
from __future__ import annotations
import argparse, csv, hashlib, importlib.util, io, json, sys, zipfile
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import minimize, minimize_scalar
from metric_geometry import (induced,physical_velocity,kinetic_jacobian,
                             occupation_covariance,occupied_metric,metric_norm,
                             oriented_polygon,verify_pullback)
ROOT=Path(__file__).resolve().parents[1]
CFG=json.loads((ROOT/'inputs_config.json').read_text())
PARENT=ROOT/'inputs'/CFG['parent_filename']
OUT=ROOT/'results'

def sha(b):return hashlib.sha256(b).hexdigest()

def authenticate_parent():
    if not PARENT.exists() or sha(PARENT.read_bytes())!=CFG['parent_sha256']:
        raise RuntimeError('V18_53T_PARENT_SHA256_MISMATCH')
    with zipfile.ZipFile(PARENT) as z:
        if z.testzip() is not None:raise RuntimeError('V18_53T_PARENT_ZIP_CRC_FAILED')
        names=z.namelist();roots={x.split('/')[0] for x in names if '/' in x}
        if len(roots)!=1:raise RuntimeError('UNEXPECTED_PARENT_ROOT')
        prefix=next(iter(roots))+'/'
        counts={}
        for m in ('SOURCE_MANIFEST.json','OUTPUT_MANIFEST.json'):
            items=json.loads(z.read(prefix+m))['files']
            for rel,expect in items.items():
                if sha(z.read(prefix+rel))!=expect:raise RuntimeError('PARENT_MANIFEST_MISMATCH '+rel)
            counts[m]=len(items)
        report=json.loads(z.read(prefix+'results/report.json'))
        if report['scientific_decision']!=CFG['parent_decision']:
            raise RuntimeError('UNEXPECTED_PARENT_SCIENTIFIC_DECISION')
        with np.load(io.BytesIO(z.read(prefix+'results/augmented21_primitives_and_oriented_valleys.npz')),allow_pickle=False) as f:
            data={x:np.array(f[x],copy=True) for x in f.files if x.endswith(('__xi','__z21','__C21','__g21'))}
        ctext=z.read(prefix+'results/valley_geometry_21d.csv').decode()
    rows=list(csv.DictReader(io.StringIO(ctext)))
    centers={(x['seed'],x['kappa_rule']):(float(x['origin_center']),float(x['width_center']))
             for x in rows if x['variant']=='on_the_fly_coupled21'}
    if len(centers)!=6:raise RuntimeError('EXPECTED_6_AUTHENTIC_BASELINE_CASES')
    return data,centers,counts,report

def historical_embedding_validate(data):
    f=ROOT/'inputs/torque_math_v18_53k_authenticated.py'
    spec=importlib.util.spec_from_file_location('q550_v1853u_exact_historical_embed',f)
    h=importlib.util.module_from_spec(spec);sys.modules[spec.name]=h;spec.loader.exec_module(h)
    mx=float(CFG['m_x'])
    worst=0.
    for seed in CFG['seeds']:
        for rule in CFG['kappa_rules']:
            z=data[f'{seed}__{rule}__z21']
            take=np.unique(np.linspace(0,len(z)-1,17).astype(int))
            _,v,_,_=h.embed(z[take,:10],mx)
            err=np.max(np.abs(v-physical_velocity(z[take],mx)))
            worst=max(worst,float(err))
            if err>5e-12:raise RuntimeError('HISTORICAL_EMBEDDING_VELOCITY_MISMATCH')
    return worst

def residual(prox,a,w):return prox(a+3*w)-prox(a)

def best_width(prox,H,a,wlo,whi,grid=25):
    ws=np.linspace(wlo,whi,grid)
    vals=np.array([residual(prox,a,w)@H@residual(prox,a,w) for w in ws])
    candidates=[(vals[0],ws[0]),(vals[-1],ws[-1])]
    for i in range(1,len(ws)-1):
        if vals[i]<=vals[i-1] and vals[i]<=vals[i+1]:
            s=minimize_scalar(lambda w:residual(prox,a,w)@H@residual(prox,a,w),
                              bounds=(ws[i-1],ws[i+1]),method='bounded',
                              options={'xatol':1e-11,'maxiter':120})
            candidates.append((float(s.fun),float(s.x)))
    if len(candidates)==2:
        j=int(np.argmin(vals));s=minimize_scalar(lambda w:residual(prox,a,w)@H@residual(prox,a,w),
                              bounds=(ws[max(0,j-1)],ws[min(len(ws)-1,j+1)]),
                              method='bounded',options={'xatol':1e-11})
        candidates.append((float(s.fun),float(s.x)))
    val,w=min(candidates,key=lambda p:p[0])
    return w,bool(abs(w-wlo)<1e-7 or abs(w-whi)<1e-7)

def optimize_2d(prox,H,oc,wc,bounds):
    (alo,ahi),(wlo,whi)=bounds
    seeds=[(oc,wc),(alo+.2*(ahi-alo),wc),(alo+.8*(ahi-alo),wc),
           (oc,wlo+.2*(whi-wlo)),(oc,whi-.2*(whi-wlo))]
    fun=lambda p:0.5*residual(prox,p[0],p[1])@H@residual(prox,p[0],p[1])
    def jac(p):
        a,w=p;N=residual(prox,a,w);b=a+3*w
        return np.array([(prox(b,1)-prox(a,1))@H@N,3*prox(b,1)@H@N])
    candidates=[]
    for pt in seeds:
        opt=minimize(fun,pt,jac=jac,bounds=bounds,method='L-BFGS-B',
                     options={'ftol':5e-15,'gtol':1e-12,'maxiter':350,'maxls':50})
        candidates.append(opt)
    best=min(candidates,key=lambda s:s.fun)
    aa,ww=best.x
    bound=any(abs(v-lo)<1e-6 or abs(v-hi)<1e-6
              for v,(lo,hi) in zip(best.x,bounds))
    return dict(origin=float(aa),width=float(ww),min_norm=float(np.sqrt(2*best.fun)),
                gradient_norm=float(np.linalg.norm(jac(best.x))),
                bound=bool(bound),success=bool(best.success))

def normalize(x):return np.asarray(x,float)/np.linalg.norm(x)

def geometry_case(seed,rule,name,G,xi,z,prox,oc,wc,points):
    H=induced(G)
    spec=np.linalg.eigvalsh(H)
    origins=np.linspace(oc-CFG['origin_halfspan'],oc+CFG['origin_halfspan'],points)
    wlo,wup=wc-CFG['width_halfband'],wc+CFG['width_halfband']
    if origins[0]<xi[0] or origins[-1]+3*wup>xi[-1]:raise RuntimeError('OUTSIDE_TRACE')
    ww=[];bounds=[]
    for a in origins:
        w,b=best_width(prox,H,float(a),wlo,wup,grid=CFG['width_grid'])
        ww.append(w);bounds.append(b)
    ww=np.asarray(ww)
    N=np.array([residual(prox,float(a),float(w)) for a,w in zip(origins,ww)])
    norms=metric_norm(N,H)
    arclength,polygon,u=oriented_polygon(N,H)
    best=optimize_2d(prox,H,oc,wc,[(oc-CFG['origin_halfspan'],oc+CFG['origin_halfspan']),(wlo,wup)])
    i=len(origins)//2
    # At fixed origin changing the metric changes the norm; the H-induced direction
    # is represented in the principal symmetric SPD square-root frame.
    ptest=[]
    for at in np.unique(np.linspace(0,len(origins)-1,9).astype(int)):
        idx=int(np.argmin(abs(xi-origins[at])))
        ptest.append(verify_pullback(z[idx],G,N[at],CFG['m_x']))
    out=dict(seed=seed,kappa_rule=rule,metric=name,metric_5d_cond=float(np.linalg.cond(G)),
             metric_10d_cond=float(spec[-1]/spec[0]),
             metric_10d_min_eig=float(spec[0]),
             min_origin=best['origin'],min_width=best['width'],min_norm=best['min_norm'],
             min_grad=best['gradient_norm'],min_on_boundary=best['bound'],
             optimizer_success=best['success'],
             curve_center_width=float(ww[i]),curve_center_norm=float(norms[i]),
             curve_norm_min=float(norms.min()),curve_arclength=float(arclength),
             polygon45_norm=float(np.linalg.norm(polygon)),
             bound_count=int(sum(bounds)),
             pullback_max_relative_error=float(max(p['rel_error'] for p in ptest)),
             pullback_max_jac_cond=float(max(p['jac_cond'] for p in ptest)),
             formal_endpoints=bool(np.all(origins+3*ww>25.6)))
    details=[dict(seed=seed,kappa_rule=rule,metric=name,origin=float(a),width=float(w),
                  norm=float(n),x0=float(N[j,0]),x1=float(N[j,1]),
                  u0=float(u[j,0]),u1=float(u[j,1]))
             for j,(a,w,n) in enumerate(zip(origins,ww,norms))]
    return out,details,(N,ww,origins,polygon)

def csvwrite(path,records):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,list(records[0]));w.writeheader();w.writerows(records)

def run(full=False):
    data,centers,ancestry,report=authenticate_parent()
    matched=historical_embedding_validate(data)
    if not full:
        print('AUTHENTIC_V1853T_SOURCE_PASS',ancestry,'HISTORICAL_5D_EMBED_MATCH',matched)
        return
    cases=[];curves=[];store={};errors=[];reference={}
    for seed in CFG['seeds']:
        for rule in CFG['kappa_rules']:
            tag=f'{seed}__{rule}'
            xi=data[tag+'__xi'];z=data[tag+'__z21'];C=data[tag+'__C21'];g=data[tag+'__g21']
            if len(xi)!=8961 or z.shape!=(8961,11) or C.shape!=(8961,10):
                raise RuntimeError('UNEXPECTED_AUTHENTIC_4X_TRACE')
            prox=CubicHermiteSpline(xi,C,g,axis=0,extrapolate=False)
            oc,wc=centers[(seed,rule)]
            S1=occupation_covariance(xi,z,CFG['m_x'],'preformal')
            S2=occupation_covariance(xi,z,CFG['m_x'],'full')
            if abs(np.trace(S1)-1)>1e-12 or abs(np.trace(S2)-1)>1e-12:raise RuntimeError('COVARIANCE_TRACE_FAILURE')
            families=[('kinetic_exact',np.eye(5))]
            for epoch,S in [('preformal',S1),('full_formal',S2)]:
                for lam in CFG['occupation_strengths']:
                    if lam==0:continue
                    families.append((f'occupation_{epoch}_lambda{lam:g}',occupied_metric(S,lam)))
            for name,G in families:
                stat,detail,raw=geometry_case(seed,rule,name,G,xi,z,prox,oc,wc,CFG['curve_points'])
                cases.append(stat);curves.extend(detail)
                store[tag+'__'+name+'__N10']=raw[0]
                store[tag+'__'+name+'__widths']=raw[1]
                store[tag+'__'+name+'__origins']=raw[2]
                store[tag+'__'+name+'__polygon45']=raw[3]
                if name=='kinetic_exact':reference[tag]=(stat,raw)
                print('METRIC_CASE',tag,name,'N_min',f"{stat['min_norm']:.9g}",
                      'origin',f"{stat['min_origin']:.7f}",'arc',f"{stat['curve_arclength']:.8f}",flush=True)
    # Relative changes are comparisons among distinct metrics, not coordinate invariants.
    comparison=[]
    for x in cases:
        k=x['seed']+'__'+x['kappa_rule'];r=reference[k][0]
        comparison.append(dict(seed=x['seed'],kappa_rule=x['kappa_rule'],metric=x['metric'],
                               shift_origin=x['min_origin']-r['min_origin'],
                               shift_width=x['min_width']-r['min_width'],
                               min_norm_ratio=x['min_norm']/r['min_norm'],
                               arclength_delta_rad=x['curve_arclength']-r['curve_arclength'],
                               full_positive_min=x['min_norm']>CFG['zero_screen']))
    csvwrite(OUT/'metric_valley_summary.csv',cases)
    csvwrite(OUT/'metric_valley_curves.csv',curves)
    csvwrite(OUT/'metric_comparison.csv',comparison)
    np.savez_compressed(OUT/'metric_oriented_10d_data.npz',**store)
    plot_cases(cases)
    screen=all(abs(x['shift_origin'])<=0.05 and abs(x['shift_width'])<=0.01
               and abs(x['arclength_delta_rad'])<=0.10 and x['full_positive_min']
               for x in comparison) and not any(x['min_on_boundary'] for x in cases)
    cls='KINETIC_CANONICAL_COVARIANCE_AND_CONDITIONED_METRIC_SENSITIVITY__NO_MATERIAL_CLOSURE'
    result=dict(version='v18.53u',scientific_decision=cls,mode='POSTPROCESS_AUTHENTIC_21D_T_NO_NEW_IVPS',
        checked_six_physical_cases=True,metric_sensitivity_families=5,case_count=len(cases),
        descriptive_kinematic_weight_robustness_screen=bool(screen),
        historical_embedding_max_error=matched,parent_manifest_counts=ancestry,
        max_pullback_relative_error=max(x['pullback_max_relative_error'] for x in cases),
        max_minimum_origin_shift=max(abs(x['shift_origin']) for x in comparison),
        max_minimum_width_shift=max(abs(x['shift_width']) for x in comparison),
        max_arclength_change=max(abs(x['arclength_delta_rad']) for x in comparison),
        min_weighted_residual=min(x['min_norm'] for x in cases),
        all_minima_nonzero=all(x['full_positive_min'] for x in comparison),
        any_minimum_on_boundary=any(x['min_on_boundary'] for x in cases),
        limitations=['Canonical kinetic metric already Euclidean in historical mass-whitened 5D embedding.',
          'Correct log-polar pullback is an exact coordinate change and cannot independently certify robustness.',
          'Velocity occupation weighting depends on epoch and lambda, is not a uniquely imposed physical metric.',
          'The accumulated bivector is not the full physical 11D state and does not certify recurrence.',
          'All relevant far endpoints lie in formal clock continuation beyond xi=25.6.',
          'No material holonomy or exact continuum invariant is claimed.'])
    (OUT/'report.json').write_text(json.dumps(result,indent=2))
    (OUT/'SCIENTIFIC_REPORT_SERVER.md').write_text('# Q550 v18.53u server result\n\n'+cls+'\n\n'+
      f"Cases: {len(cases)}; pullback relative error max: {result['max_pullback_relative_error']:.4g}; "+
      f"minimum shift max in xi: {result['max_minimum_origin_shift']:.4g}; "+
      f"minimum residual: {result['min_weighted_residual']:.8g}.\n\n"+
      'The kinetic canonical 5D metric is Euclidean in the historical mass-weighted Cartesian embedding. '+
      'Additional frozen metrics are model-conditioned sensitivity probes only.\n')
    print('V1853U',cls,'CASES',len(cases))

def plot_cases(cases):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for seed in CFG['seeds']:
        rows=[r for r in cases if r['seed']==seed and r['kappa_rule']=='shared']
        fig,axs=plt.subplots(1,2,figsize=(11,4.2))
        for r in rows:
            axs[0].scatter(r['min_origin'],r['min_width'],label=r['metric'])
            axs[1].scatter(r['min_origin'],r['curve_arclength'],label=r['metric'])
        axs[0].set(xlabel='Origin xi',ylabel='Optimized width',title='Minimum locus (frozen metrics)')
        axs[1].set(xlabel='Origin xi',ylabel='Direction arclength [rad]',title='Weighted angular arclength')
        axs[0].legend(fontsize=7);fig.suptitle('Q550 v18.53u — '+seed+' seed; shared kappa')
        fig.tight_layout();fig.savefig(OUT/f'kinetic_metric_valley_{seed}.png',dpi=140);plt.close(fig)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--preflight',action='store_true')
    ap.add_argument('--full',action='store_true');args=ap.parse_args()
    run(full=args.full)
