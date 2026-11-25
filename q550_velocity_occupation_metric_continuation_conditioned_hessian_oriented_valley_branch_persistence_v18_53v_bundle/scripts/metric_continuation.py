#!/usr/bin/env python3
"""Q550 v18.53v: analytic metric continuation and conditioned minimum audit.

No physical IVP is run. Historical u and t ZIPs must authenticate. Mathematical
comments, labels and artifact content are intentionally written in English.
"""
from __future__ import annotations
import argparse,csv,hashlib,importlib.util,io,json,sys,zipfile
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import minimize,root as scipy_root,minimize_scalar
ROOT=Path(__file__).resolve().parents[1]
CFG=json.loads((ROOT/'inputs_config.json').read_text())
OUT=ROOT/'results'
OUT.mkdir(parents=True,exist_ok=True)
HELPER=ROOT/'inputs/u_metric_geometry_authenticated.py'
_spec=importlib.util.spec_from_file_location('u_metric_geometry_authenticated',HELPER)
mg=importlib.util.module_from_spec(_spec);sys.modules[_spec.name]=mg;_spec.loader.exec_module(mg)
PAIRS=mg.PAIRS

def sha(x:bytes)->str:return hashlib.sha256(x).hexdigest()
def verify_zip(z,prefix=None):
    if z.testzip() is not None:raise RuntimeError('BAD_ZIP_CRC')
    if prefix is None:
        roots={n.split('/')[0] for n in z.namelist() if '/' in n}
        if len(roots)!=1:raise RuntimeError('BAD_ZIP_ROOTS')
        prefix=next(iter(roots))+'/'
    counts={}
    for m in ('SOURCE_MANIFEST.json','OUTPUT_MANIFEST.json'):
        files=json.loads(z.read(prefix+m))['files']
        for n,h in files.items():
            if sha(z.read(prefix+n))!=h:raise RuntimeError('BAD_PARENT_MANIFEST '+n)
        counts[m]=len(files)
    return counts,prefix

def load_authentic_parent(light=False):
    blob=(ROOT/'inputs'/CFG['parent_filename']).read_bytes()
    if sha(blob)!=CFG['parent_sha256']:raise RuntimeError('BAD_U_PARENT_SHA')
    if sha(HELPER.read_bytes())!=CFG['u_helper_sha256']:raise RuntimeError('BAD_HISTORICAL_HELPER_SHA')
    with zipfile.ZipFile(io.BytesIO(blob)) as z:
        counts,p=verify_zip(z)
        if sha(z.read(p+'scripts/metric_geometry.py'))!=CFG['u_helper_sha256']:
            raise RuntimeError('HELPER_NOT_HISTORICAL_U')
        ureport=json.loads(z.read(p+'results/report.json'))
        if ureport['scientific_decision']!=CFG['expected_parent_classification']:
            raise RuntimeError('UNEXPECTED_U_SCIENTIFIC_RESULT')
        old=list(csv.DictReader(io.StringIO(z.read(p+'results/metric_valley_summary.csv').decode())))
        tblob=z.read(p+'inputs/v18_53t_authenticated_server_sendback.zip')
        if sha(tblob)!=CFG['expected_t_sha256']:raise RuntimeError('BAD_T_PARENT_SHA')
    with zipfile.ZipFile(io.BytesIO(tblob)) as z:
        tcounts,tp=verify_zip(z)
        treport=json.loads(z.read(tp+'results/report.json'))
        if treport['scientific_decision']!='AUTHENTIC_21D_ON_THE_FLY_PRIMITIVE__CONTROLLED_11D_STEP_EFFECT__NO_MATERIAL_HOLONOMY':
            raise RuntimeError('UNEXPECTED_T_SCIENTIFIC_RESULT')
        if not light:
            with np.load(io.BytesIO(z.read(tp+'results/augmented21_primitives_and_oriented_valleys.npz')),allow_pickle=False) as f:
                data={k:np.array(f[k],copy=True) for k in f.files if k.endswith(('__xi','__z21','__C21','__g21'))}
    source={}
    if not light:
        for seed in CFG['seeds']:
            for rule in CFG['kappa_rules']:
                tag=seed+'__'+rule
                xi=data[tag+'__xi'];z=data[tag+'__z21'];C=data[tag+'__C21'];g=data[tag+'__g21']
                if len(xi)!=8961 or xi.ndim!=1 or z.shape!=(8961,11) or C.shape!=(8961,10) or g.shape!=(8961,10):
                    raise RuntimeError('UNEXPECTED_AUTHENTIC_TRACE_SHAPE '+tag)
                if np.any(np.diff(xi)<=0) or not all(np.isfinite(x).all() for x in (z,C,g)):
                    raise RuntimeError('INVALID_AUTHENTIC_TRACE '+tag)
                source[tag]=(xi,z,C,g)
    if len(old)!=30:raise RuntimeError('UNEXPECTED_U_CASE_COUNT')
    return source,old,dict(parent_u=counts,parent_t=tcounts),ureport

def analytic_metric(S,lam):
    """G=(I+lambda*S)det(I+lambda*S)^(-1/5) and analytic derivative."""
    S=np.asarray(S,float);A=np.eye(5)+lam*S
    sign,logdet=np.linalg.slogdet(A)
    if sign<=0:raise RuntimeError('NON_SPD_5D_METRIC')
    G=A*np.exp(-logdet/5.)
    M=np.linalg.solve(A,S)
    Gp=np.exp(-logdet/5.)*(S-np.trace(M)/5.*A)
    H=mg.induced(G)
    Hp=np.array([[Gp[i,k]*G[j,l]+G[i,k]*Gp[j,l]
                    -Gp[i,l]*G[j,k]-G[i,l]*Gp[j,k] for k,l in PAIRS]
                   for i,j in PAIRS])
    Hp=(Hp+Hp.T)*.5
    if np.linalg.eigvalsh(H)[0]<=0:raise RuntimeError('NON_SPD_INDUCED_10D')
    return G,Gp,H,Hp

class ResidualObjective:
    def __init__(self,prox,H,Hp):self.prox=prox;self.H=H;self.Hp=Hp
    def terms(self,p):
        a,w=np.asarray(p,float);b=a+3*w
        C=self.prox; N=C(b)-C(a)
        ga=C(a,1);gb=C(b,1); gpa=C(a,2);gpb=C(b,2)
        J=np.stack([gb-ga,3*gb],axis=1)
        # The 2x2x10 derivative tensors reflect the dependence b=a+3*w.
        Daa=gpb-gpa; Daw=3*gpb;Dww=9*gpb
        return N,J,(Daa,Daw,Dww)
    def value_grad_hess(self,p):
        N,J,(Daa,Daw,Dww)=self.terms(p)
        f=.5*N@self.H@N
        grad=J.T@self.H@N
        Hess=J.T@self.H@J
        Hess+=np.array([[Daa@self.H@N,Daw@self.H@N],
                        [Daw@self.H@N,Dww@self.H@N]])
        return f,grad,(Hess+Hess.T)/2
    def lambda_derivatives(self,p):
        N,J,_=self.terms(p)
        f0=.5*N@self.Hp@N
        gp=J.T@self.Hp@N
        return f0,gp

def improve_root(obj,start,bounds):
    start=np.asarray(start,float)
    def inside(p):return all(lo-1e-9<=v<=hi+1e-9 for v,(lo,hi) in zip(p,bounds))
    ans=scipy_root(lambda x:obj.value_grad_hess(x)[1],start,
                   jac=lambda x:obj.value_grad_hess(x)[2],method='hybr',
                   options={'xtol':1e-11})
    if np.isfinite(ans.x).all() and inside(ans.x) and np.linalg.norm(obj.value_grad_hess(ans.x)[1])<3e-9:
        return np.asarray(ans.x,float), 'stationary_root'
    def fj(p):f,g,_=obj.value_grad_hess(p);return f,g
    fit=minimize(fj,start,jac=True,method='L-BFGS-B',bounds=bounds,
                 options={'ftol':1e-15,'gtol':2e-12,'maxiter':250,'maxls':50})
    return np.asarray(fit.x,float),'bounded_minimize'

def competitor_search(obj,start,bounds):
    (alo,ahi),(wlo,whi)=bounds
    a,w=start
    trial=[(a,w),((alo+ahi)/2,(wlo+whi)/2),
           (alo+.2*(ahi-alo),w),(alo+.8*(ahi-alo),w),
           (a,wlo+.2*(whi-wlo)),(a,whi-.2*(whi-wlo))]
    results=[]
    for x in trial:
        p,method=improve_root(obj,x,bounds)
        if not all(lo-1e-7<=v<=hi+1e-7 for v,(lo,hi) in zip(p,bounds)):continue
        f,g,h=obj.value_grad_hess(p)
        results.append((float(f),p,float(np.linalg.norm(g))))
    if not results:raise RuntimeError('NO_INDEPENDENT_COMPETITOR')
    return min(results,key=lambda x:x[0]),len(results)

def best_width(prox,H,a,bounds,grid=21):
    lo,hi=bounds
    ws=np.linspace(lo,hi,grid)
    vals=np.array([float((prox(a+3*w)-prox(a))@H@(prox(a+3*w)-prox(a))) for w in ws])
    trials=[(vals[0],ws[0]),(vals[-1],ws[-1])]
    inds=[i for i in range(1,grid-1) if vals[i]<=vals[i-1] and vals[i]<=vals[i+1]]
    if not inds:inds=[int(np.argmin(vals[1:-1]))+1]
    for i in inds:
        fit=minimize_scalar(lambda w:float((prox(a+3*w)-prox(a))@H@(prox(a+3*w)-prox(a))),
             bounds=(ws[i-1],ws[i+1]),method='bounded',options={'xatol':2e-11})
        trials.append((fit.fun,fit.x))
    _,w=min(trials,key=lambda x:x[0]);return float(w),bool(min(abs(w-lo),abs(w-hi))<1e-7)

def path_geometry(prox,H,ref,oc,wc,bounds):
    al,ah=bounds[0]; wl,wh=bounds[1]
    a=np.linspace(oc-CFG['origin_halfspan'],oc+CFG['origin_halfspan'],CFG['origin_curve_samples'])
    if a[0]<al-1e-9 or a[-1]>ah+1e-9:raise RuntimeError('GEOMETRY_OUTSIDE_ORIGIN_BOUNDS')
    w=np.array([best_width(prox,H,float(x),(wl,wh),CFG['width_search_grid'])[0] for x in a])
    N=np.array([prox(x+3*y)-prox(x) for x,y in zip(a,w)])
    arc,polygon,unit=mg.oriented_polygon(N,H)
    if ref is None:diff=float('nan');pangle=float('nan')
    else:
        # A common symmetric square-root identification is a metric-dependent convention.
        polygon0=ref['polygon'];norm0=np.linalg.norm(polygon0);norm=np.linalg.norm(polygon)
        c=np.clip(float(polygon@polygon0/max(norm*norm0,1e-300)),-1,1)
        pangle=float(np.degrees(np.arccos(c)))
    return dict(arclength=float(arc),polygon=polygon,oriented_angle_deg=pangle,
                max_width=float(w.max()),min_width=float(w.min()),
                width_total_variation=float(np.abs(np.diff(w)).sum()),
                turning_count=int(np.sum(np.diff(w)[:-1]*np.diff(w)[1:]< -1e-10)),
                endpoint_mean=float(np.mean(a+3*w)),
                boundary_hits=int(np.sum(np.minimum(np.abs(w-wl),np.abs(w-wh))<1e-7))), (a,w,N)

def csvwrite(p,rows):
    if not rows:raise RuntimeError('EMPTY_CSV '+str(p))
    with p.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,list(rows[0]));w.writeheader();w.writerows(rows)

def run(full=False,smoke=False):
    data,old,ancestry,report=load_authentic_parent(light=not(full or smoke))
    print('AUTHENTIC_U_T_PARENT_PASS',ancestry,flush=True)
    if not full and not smoke:return
    olddict={(r['seed'],r['kappa_rule'],r['metric']):r for r in old}
    lamvals=[0.,.5,1.] if smoke else CFG['lambdas']
    seeds=['original'] if smoke else CFG['seeds']
    rules=['shared'] if smoke else CFG['kappa_rules']
    epochs=['preformal'] if smoke else CFG['epochs']
    summary=[];curve_rows=[];audit_u=[];arr_store={};stats=[]
    for seed in seeds:
        for rule in rules:
            tag=seed+'__'+rule;xi,z,C,g=data[tag]
            prox=CubicHermiteSpline(xi,C,g,axis=0,extrapolate=False)
            u_base=olddict[(seed,rule,'kinetic_exact')]
            oc,wc=float(u_base['min_origin']),float(u_base['min_width'])
            bounds=[(oc-CFG['origin_halfspan'],oc+CFG['origin_halfspan']),
                    (wc-CFG['width_halfband'],wc+CFG['width_halfband'])]
            for epoch in epochs:
                S=mg.occupation_covariance(xi,z,CFG['m_x'],'preformal' if epoch=='preformal' else 'full')
                if abs(np.trace(S)-1)>1e-12:raise RuntimeError('INVALID_OCCUPATION_NORMALIZATION')
                prev_p=np.array([oc,wc]);prev_lam=0.;prev_dp=None;refpath=None;case_rows=[]
                for index,lam in enumerate(lamvals):
                    G,Gp,H,Hp=analytic_metric(S,float(lam));obj=ResidualObjective(prox,H,Hp)
                    predictor=prev_p.copy() if prev_dp is None else prev_p+(lam-prev_lam)*prev_dp
                    predictor=np.clip(predictor,[b[0]+1e-7 for b in bounds],[b[1]-1e-7 for b in bounds])
                    p,method=improve_root(obj,predictor,bounds)
                    f,grad,Hess=obj.value_grad_hess(p)
                    scaling=np.diag([CFG['origin_halfspan'],CFG['width_halfband']]);Hs=scaling@Hess@scaling
                    eig=np.linalg.eigvalsh(Hs);eval0=np.linalg.eigvalsh(Hess)
                    pos=bool(eig[0]>CFG['min_hessian_scaled_eigen_screen'])
                    if pos:
                        _,gl=obj.lambda_derivatives(p)
                        dp=np.linalg.solve(Hess,-gl)
                    else:dp=np.array([float('nan')]*2)
                    dfn,_=obj.lambda_derivatives(p)
                    norm=np.sqrt(max(2*f,0))
                    normlam=dfn/max(norm,1e-30)
                    # Every other node, plus lambda = 1 and 3, gets independent multistart.
                    do_multi=bool(index%2==0 or abs(lam-1)<1e-9 or abs(lam-3)<1e-9 or index==len(lamvals)-1)
                    if do_multi:
                        (fbest,pbest,gbest),num=competitor_search(obj,p,bounds)
                        branch_dist=float(np.linalg.norm(p-pbest));cost_gap=float(f-fbest)
                    else:
                        branch_dist=float('nan');cost_gap=float('nan');num=0
                    boundary=any(min(abs(v-lo),abs(v-hi))<1e-6 for v,(lo,hi) in zip(p,bounds))
                    ref=(refpath if refpath is not None else None)
                    geom,raw=path_geometry(prox,H,ref,oc,wc,bounds)
                    if refpath is None:refpath=geom
                    origin_shift=float(p[0]-oc);width_shift=float(p[1]-wc)
                    Udiag=None
                    if lam in (0.,1.,3.):
                        metric='kinetic_exact' if lam==0 else f'occupation_{epoch}_lambda{lam:g}'
                        uref=olddict.get((seed,rule,metric))
                        if uref is None:raise RuntimeError('MISSING_U_BASELINE '+metric)
                        Udiag=dict(seed=seed,kappa_rule=rule,epoch=epoch,lam=lam,
                             origin_delta=p[0]-float(uref['min_origin']),
                             width_delta=p[1]-float(uref['min_width']),
                             norm_delta=norm-float(uref['min_norm']),
                             arclength_delta=geom['arclength']-float(uref['curve_arclength']))
                        audit_u.append(Udiag)
                    row=dict(seed=seed,kappa_rule=rule,epoch=epoch,lambda_value=float(lam),
                      origin=float(p[0]),width=float(p[1]),minimum_norm=float(norm),
                      minimum_energy=float(f),gradient_norm=float(np.linalg.norm(grad)),
                      hessian_min_eig=float(eval0[0]),hessian_max_eig=float(eval0[-1]),
                      hessian_scaled_min_eig=float(eig[0]),
                      hessian_scaled_cond=float(eig[-1]/eig[0]) if eig[0]>0 else float('inf'),
                      metric_5d_cond=float(np.linalg.cond(G)),metric_10d_cond=float(np.linalg.cond(H)),
                      d_origin_d_lambda=float(dp[0]),d_width_d_lambda=float(dp[1]),
                      d_minimum_norm_d_lambda=float(normlam),
                      origin_shift_from_canonical=origin_shift,width_shift_from_canonical=width_shift,
                      is_internal=not boundary,positive_hessian=pos,stationarity_ok=bool(np.linalg.norm(grad)<CFG['minimum_root_tol']),
                      method=method,independent_search=do_multi,independent_candidates=num,
                      competing_distance=branch_dist,competing_cost_gap=cost_gap,
                      arclength=geom['arclength'],oriented_polygon_angle_from_zero_deg=geom['oriented_angle_deg'],
                      curve_width_total_variation=geom['width_total_variation'],
                      curve_width_turning_count=geom['turning_count'],curve_boundary_hits=geom['boundary_hits'],
                      min_residual_nonzero=bool(norm>CFG['zero_screen']),
                      far_endpoint_formal=bool(p[0]+3*p[1]>CFG['formal_clock_boundary']))
                    summary.append(row);case_rows.append(row)
                    for a,w,N in zip(*raw):
                        curve_rows.append(dict(seed=seed,kappa_rule=rule,epoch=epoch,lambda_value=float(lam),
                                origin=float(a),optimized_width=float(w),
                                norm10=float(np.sqrt(N@H@N))))
                    arr_store[f'{tag}__{epoch}__lam{index:02d}__polygon45']=geom['polygon']
                    prev_p=np.asarray(p);prev_lam=float(lam);prev_dp=dp if pos else None
                    if not smoke:
                        print('METRIC_CONTINUE',tag,epoch,lam,'origin',f'{p[0]:.8f}',
                              'norm',f'{norm:.9f}','hessian_scaled_min',f'{eig[0]:.4g}',flush=True)
                for ia in range(1,len(case_rows)-1):
                    l=case_rows[ia];before=case_rows[ia-1];after=case_rows[ia+1]
                    dl=after['lambda_value']-before['lambda_value']
                    for key,dkey in (('origin','d_origin_d_lambda'),('width','d_width_d_lambda'),('minimum_norm','d_minimum_norm_d_lambda')):
                        l['fd_'+dkey]=(after[key]-before[key])/dl
                        l['fd_error_'+dkey]=l['fd_'+dkey]-l[dkey]
                for ix in (0,len(case_rows)-1):
                    for k in ('d_origin_d_lambda','d_width_d_lambda','d_minimum_norm_d_lambda'):
                        case_rows[ix]['fd_'+k]=float('nan');case_rows[ix]['fd_error_'+k]=float('nan')
                # row insertion already points at the same mutable dicts
                stats.append(dict(seed=seed,kappa_rule=rule,epoch=epoch,
                    origin_range=max(x['origin'] for x in case_rows)-min(x['origin'] for x in case_rows),
                    max_hessian_cond=max(x['hessian_scaled_cond'] for x in case_rows),
                    min_hessian_scaled_eig=min(x['hessian_scaled_min_eig'] for x in case_rows),
                    max_multistart_branch_distance=max((x['competing_distance'] for x in case_rows if x['independent_search']),default=0),
                    min_norm=min(x['minimum_norm'] for x in case_rows),
                    max_orientation_angle_deg=max(x['oriented_polygon_angle_from_zero_deg'] for x in case_rows)))
    if smoke:
        print('AUTHENTIC_LIGHT_METRIC_CONTINUATION_PASS',len(summary),'nodes');return
    csvwrite(OUT/'metric_lambda_continuation.csv',summary)
    csvwrite(OUT/'metric_lambda_valley_curves.csv',curve_rows)
    csvwrite(OUT/'metric_lambda_parent_u_replay.csv',audit_u)
    csvwrite(OUT/'metric_lambda_family_summary.csv',stats)
    np.savez_compressed(OUT/'metric_lambda_oriented_45d.npz',**arr_store)
    make_figures(summary)
    gap=max(abs(x[k]) for x in audit_u for k in ('origin_delta','width_delta','norm_delta','arclength_delta'))
    if gap>CFG['relative_u_audit_tolerance']:
        raise RuntimeError(f'PARENT_U_REGRESSION_BREACH {gap:.6g}')
    all_positive=all(x['positive_hessian'] for x in summary)
    all_internal=all(x['is_internal'] for x in summary)
    stat=all(x['stationarity_ok'] for x in summary)
    multi=max((x['competing_distance'] for x in summary if x['independent_search']),default=0.)
    branch=multi<=CFG['max_branch_vs_global_param_diff']
    category=('REGULAR_CONDITIONED_NUMERIC_METRIC_CONTINUATION__NO_MATERIAL_HOLONOMY'
              if all_positive and all_internal and stat and branch else
              'METRIC_CONTINUATION_NUMERIC_BRANCH_AMBIGUITY_OR_CONDITIONING_WARNING__NO_MATERIAL_HOLONOMY')
    report=dict(version='v18.53v',scientific_decision=category,
        mode='AUTHENTIC_U_T_POSTPROCESS_NO_NEW_PHYSICAL_IVPS',
        n_cases=len(summary),n_metric_families=2,n_physical_cases=6,
        n_lambda_values=len(CFG['lambdas']),n_geometry_rows=len(curve_rows),
        n_ancestral_u_replays=len(audit_u),parent_manifest_counts=ancestry,
        max_u_replay_absolute_difference=gap,all_interior_minima=all_internal,
        all_positive_scaled_hessians=all_positive,all_stationary_minima=stat,
        independent_branch_screen_pass=branch,
        max_multistart_param_displacement=multi,
        min_scaled_hessian_eigenvalue=min(x['hessian_scaled_min_eig'] for x in summary),
        max_scaled_hessian_condition=max(x['hessian_scaled_cond'] for x in summary),
        max_origin_shift=max(abs(x['origin_shift_from_canonical']) for x in summary),
        max_oriented_polygon_angle_deg=max(x['oriented_polygon_angle_from_zero_deg'] for x in summary),
        limitations=['The canonical kinetic metric is Euclidean in the historical mass-whitened Cartesian embedding.',
        'Velocity-occupation-weighted metrics form model-conditioned sensitivity probes, not new canonical physical metrics.',
        'Implicit derivatives and Hessians are exact for the archived CubicHermiteSpline interpolant, not certified for the continuous flow.',
        'The 45D oriented polygon uses the principal symmetric square-root metric identification, which is convention-dependent.',
        'Far window endpoints remain in formal clock continuation beyond xi=25.6.',
        'Neither a nonzero accumulated 10D bivector minimum nor this continuation proves 11D recurrence, material memory or material holonomy.'])
    (OUT/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    (OUT/'SCIENTIFIC_REPORT_SERVER.md').write_text(
        '# Q550 v18.53v — Server scientific result\n\n'+category+'\n\n'+
        f"Nodes: {len(summary)}; metric families: 2; physical cases: 6; parent u replays: {len(audit_u)}.\n\n"+
        f"Minimum scaled Hessian eigenvalue: {report['min_scaled_hessian_eigenvalue']:.8g}; "+
        f"maximum normalized Hessian condition: {report['max_scaled_hessian_condition']:.8g}.\n\n"+
        f"Maximum origin shift: {report['max_origin_shift']:.8g}; maximum oriented polygon angle: "+
        f"{report['max_oriented_polygon_angle_deg']:.8g} deg.\n\n"+
        'All stated invariances and thresholds apply to the interpolated bivector and selected metric family, not to material closure.\n')
    print('V1853V',category,'CASES',len(summary),'U_REPLAY_MAX',gap)

def make_figures(summary):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for seed in CFG['seeds']:
        fig,axs=plt.subplots(2,2,figsize=(11,8),constrained_layout=True)
        for epoch in CFG['epochs']:
            s=[x for x in summary if x['seed']==seed and x['kappa_rule']=='shared' and x['epoch']==epoch]
            L=[x['lambda_value'] for x in s]
            axs[0,0].plot(L,[x['origin'] for x in s],marker='.',label=epoch)
            axs[0,1].plot(L,[x['minimum_norm'] for x in s],marker='.',label=epoch)
            axs[1,0].plot(L,[x['hessian_scaled_min_eig'] for x in s],marker='.',label=epoch)
            axs[1,1].plot(L,[x['oriented_polygon_angle_from_zero_deg'] for x in s],marker='.',label=epoch)
        for ax,yl in zip(axs.flat,['Minimum origin ξ','10D weighted minimum norm','Min scaled Hessian eigenvalue','Oriented 45D angle from λ=0 [deg]']):
            ax.set(xlabel='λ',ylabel=yl);ax.grid(alpha=.25);ax.legend(fontsize=8)
        fig.suptitle(f'Q550 v18.53v: metric continuation ({seed}, shared κ₀)')
        fig.savefig(OUT/f'metric_continuation_{seed}.png',dpi=150);plt.close(fig)

if __name__=='__main__':
    a=argparse.ArgumentParser()
    a.add_argument('--preflight',action='store_true')
    a.add_argument('--smoke',action='store_true')
    a.add_argument('--full',action='store_true')
    args=a.parse_args()
    run(full=args.full,smoke=args.smoke)
