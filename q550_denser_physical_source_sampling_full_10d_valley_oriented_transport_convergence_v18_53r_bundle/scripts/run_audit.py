#!/usr/bin/env python3
"""Q550 v18.53r: denser *physical* 11D source sampling and oriented 10D valley convergence.

All physical model calls are imported from the authenticated ancestor v18.53o.
Do not run full integrations on a laptop or replace the Q550 equations with synthetic data.
"""
from __future__ import annotations
import argparse, csv, hashlib, importlib.util, io, json, sys, tempfile, traceback, zipfile
from pathlib import Path
from contextlib import ExitStack
import numpy as np
from scipy.integrate import solve_ivp
from source_valley_math import Primitive, optimized_width, curve

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'
CFG=json.loads((ROOT/'inputs_config.json').read_text())
PARENT=ROOT/'inputs'/CFG['parent_filename']
S=lambda b:hashlib.sha256(b).hexdigest()

def save_csv(name, rows):
    if not rows:raise ValueError('REFUSE_EMPTY_'+name)
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with (OUT/name).open('w',newline='') as fd:
        w=csv.DictWriter(fd,fieldnames=fields);w.writeheader();w.writerows(rows)

def save_json(name, obj):
    (OUT/name).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')

def check_zip(z,required_decision,label):
    if z.testzip() is not None:raise ValueError(label+'_ZIP_CRC_FAIL')
    names=z.namelist()
    if any(Path(x).is_absolute() or '..' in Path(x).parts for x in names):raise ValueError(label+'_ZIP_UNSAFE_MEMBER')
    roots={x.split('/')[0] for x in names if '/' in x}
    if len(roots)!=1:raise ValueError(label+'_INVALID_ROOTS')
    prefix=next(iter(roots))+'/'
    counts={}
    for mname in ('SOURCE_MANIFEST.json','OUTPUT_MANIFEST.json'):
        entries=json.loads(z.read(prefix+mname))['files']
        for rel,digest in entries.items():
            if S(z.read(prefix+rel))!=digest:raise ValueError(label+'_MANIFEST_SHA_FAIL '+rel)
        counts[mname]=len(entries)
    report=json.loads(z.read(prefix+'results/report.json'))
    if report['scientific_decision']!=required_decision:
        raise ValueError(label+'_UNEXPECTED_PARENT_CLASSIFICATION')
    return prefix,counts,report

def authenticated_parents():
    if not PARENT.is_file() or S(PARENT.read_bytes())!=CFG['parent_sha256']:
        raise RuntimeError('Q_PARENT_SHA256_OR_FILE_MISMATCH')
    with zipfile.ZipFile(PARENT) as q:
        qp,qcounts,qr=check_zip(q,CFG['parent_q_decision'],'Q')
        if qr.get('new_physical_ivps')!=0:raise RuntimeError('Q_MUST_BE_POSTPROCESS_ONLY')
        qmath=q.read(qp+'scripts/source_valley_math.py')
        if S(qmath)!=CFG['parent_source_math_sha256']:raise RuntimeError('Q_GEOMETRY_MATH_MISMATCH')
        qconfig=json.loads(q.read(qp+'inputs_config.json'))
        pbytes=q.read(qp+'inputs/'+qconfig['parent_filename'])
        if S(pbytes)!=qconfig['parent_sha256']:raise RuntimeError('Q_P_ANCESTOR_SHA_MISMATCH')
        qsummary=list(csv.DictReader(io.StringIO(q.read(qp+'results/convergence_summary.csv').decode())))
    with zipfile.ZipFile(io.BytesIO(pbytes)) as p:
        pp,pcounts,pr=check_zip(p,CFG['parent_p_decision'],'P')
        pconfig=json.loads(p.read(pp+'inputs_config.json'))
        oname='inputs/v18_53o_authenticated_server_sendback.zip'
        obytes=p.read(pp+oname)
        if S(obytes)!=pconfig['parent_v18_53o_sha256']:raise RuntimeError('P_O_ANCESTOR_SHA_MISMATCH')
    with zipfile.ZipFile(io.BytesIO(obytes)) as o:
        op,ocounts,orr=check_zip(o,CFG['parent_o_decision'],'O')
        if orr.get('actual_Q550_ivps')!=10 or orr.get('refined_Q550_ivps')!=4:
            raise RuntimeError('PARENT_O_NOT_GENUINELY_EXECUTED')
        with np.load(io.BytesIO(o.read(op+'results/derived_primitives.npz')),allow_pickle=False) as data:
            archived={k:np.array(data[k],copy=True) for k in data.files}
        oldfixed=list(csv.DictReader(io.StringIO(o.read(op+'results/fixed_width_oriented_transport.csv').decode())))
    return obytes,archived,oldfixed,{'q':qcounts,'p':pcounts,'o':ocounts,'sha_q':CFG['parent_sha256'],'sha_o':S(obytes)},qsummary

def import_file(name,path):
    path=Path(path)
    if str(path.parent) not in sys.path:sys.path.insert(0,str(path.parent))
    spec=importlib.util.spec_from_file_location(name,path)
    mod=importlib.util.module_from_spec(spec);sys.modules[name]=mod;spec.loader.exec_module(mod)
    return mod

def require_ancestry(oarchive,stack):
    d=Path(stack.enter_context(tempfile.TemporaryDirectory(prefix='q550_v1853r_o_')))
    with zipfile.ZipFile(io.BytesIO(oarchive)) as z:
        if any(Path(x).is_absolute() or '..' in Path(x).parts for x in z.namelist()):raise RuntimeError('UNSAFE_O_ARCHIVE')
        z.extractall(d)
        r={x.split('/')[0] for x in z.namelist() if '/' in x};oroot=d/next(iter(r))
    o=import_file('q550_v1853r_authenticated_o',oroot/'scripts/run_audit.py')
    nroot,best,ncounts,nreport=o.authenticate_parent_n(stack.enter_context(tempfile.TemporaryDirectory(prefix='q550_v1853r_n_')))
    n=o.import_n(nroot)
    mroot,prior,mcounts,mreport=n.authenticate_parent_m(stack.enter_context(tempfile.TemporaryDirectory(prefix='q550_v1853r_m_')))
    m=n.import_m(mroot)
    lroot,oldl,lc=m.authenticate_l(stack.enter_context(tempfile.TemporaryDirectory(prefix='q550_v1853r_l_')))
    l=m.import_l(lroot)
    kroot,krows,kreport=l.verified_parent(stack.enter_context(tempfile.TemporaryDirectory(prefix='q550_v1853r_k_')))
    k=l.import_k(kroot)
    j=k.load_j();k.verify_j_parent()
    h=j.hist_module();h.provenance()
    for fun in ('loader','solve','rhs','validated_rhs','primary_metrics','native_53g_reconciliation'):
        if not callable(getattr(h,fun,None)):raise RuntimeError('MISSING_AUTHENTIC_PHYSICAL_FUNCTION_'+fun)
    return o,n,m,l,k,j,h,ncounts,mcounts

def refine(x,factor):
    x=np.asarray(x,float)
    if factor==1:return x.copy()
    if factor<1 or int(factor)!=factor or np.any(np.diff(x)<=0):raise ValueError('BAD_REFINEMENT')
    y=(x[:-1,None]+np.diff(x)[:,None]*np.arange(factor)[None,:]/factor).ravel()
    y=np.r_[y,x[-1]]
    if not np.allclose(y[::factor],x,rtol=0,atol=3e-13):raise ValueError('GRID_NOT_NESTED')
    return y

def archived_case(archived,seed,rule):
    tag=seed+'__'+rule
    return dict(xi=archived['common_xi'],C=archived[tag+'__primitive'],
        g=archived[tag+'__sampled_source_10d']*archived[tag+'__clock_rate'][:,None])

def find_centers(oldfixed):
    centers={}
    for seed in ('original','autonomous'):
        for rule in ('shared','original','autonomous'):
            rr=sorted((v for v in oldfixed if v['seed']==seed and v['kappa_rule']==rule),key=lambda r:int(r['index']))
            if len(rr)!=141:raise RuntimeError('INVALID_O_BASELINE_ORIGIN_COUNT '+seed+' '+rule)
            centers[seed,rule]=(float(rr[70]['origin_xi']),float(rr[70]['width_xi']))
    return centers

def curve_case(seed,rule,level,history,center,allrows,curve_rows,checkpoints):
    xi=history['xi'];C=history['C'];g=history['g']
    oc,wc=center; ox=np.linspace(oc-CFG['origin_halfspan'],oc+CFG['origin_halfspan'],CFG['curve_points'])
    a,b=wc-CFG['width_halfband'],wc+CFG['width_halfband']
    if ox[0]<xi[0] or ox[-1]+3*b>xi[-1]:raise RuntimeError('CURVE_OUT_OF_TRACE')
    for method in ('hermite','pchip'):
        p=Primitive(xi,C,g,method)
        mins=[optimized_width(p,float(o),a,b,nsearch=CFG['search_brackets']) for o in ox]
        widths=np.array([v['w'] for v in mins]);der=np.array([v['implicit_wprime'] for v in mins])
        cr=curve(p,ox,widths,der)
        center_index=(len(ox)-1)//2
        ends=ox+3*widths
        N=p.eval(ends)-p.eval(ox)
        rnorm=np.linalg.norm(N,axis=1)
        if np.min(rnorm)<=1e-12:raise RuntimeError('UNDEFINED_DIRECTION_AT_NEAR_ZERO')
        vecs=cr['u'];center_u=vecs[center_index]
        baseline=dict(seed=seed,kappa_rule=rule,density_factor=level,method=method,source_knots=len(xi),
             width_center=float(widths[center_index]),residual_norm_center=float(rnorm[center_index]),
             direction_speed_center=float(cr['directional_speed'][center_index]),
             min_residual_norm=float(np.min(rnorm)),max_direction_speed=float(cr['max_speed']),
             arclength_direct=float(cr['direct']),arclength_tangent=float(cr['chain']),
             arclength_difference=float(cr['arc_diff']),
             width_total_variation=float(cr['width_total_variation']),
             width_sign_changes=int(cr['width_sign_changes']),
             width_candidate_switches=int(sum(x['label']!=y['label'] for x,y in zip(mins[1:],mins[:-1]))),
             num_width_bounds=sum(v['bound'] for v in mins),
             tangent_width_derivative_masks=sum(not np.isfinite(v['implicit_wprime']) for v in mins),
             oriented_polygon_norm=float(np.linalg.norm(cr['oriented_polygon'])),
             oriented_tangent_norm=float(np.linalg.norm(cr['oriented_tangent'])),
             origin_center=oc,width_reference=wc,
             start_xi=float(ox[0]),end_xi=float(ends[-1]),
             endpoints_are_formal_clock_continuation=bool(np.all(ends>CFG['formal_clock_boundary'])))
        allrows.append(baseline)
        prefix=f'{seed}__{rule}__x{level}__{method}'
        checkpoints[prefix+'__origins']=ox
        checkpoints[prefix+'__widths']=widths
        checkpoints[prefix+'__u10d']=vecs
        checkpoints[prefix+'__oriented_polygon45d']=cr['oriented_polygon']
        checkpoints[prefix+'__oriented_tangent45d']=cr['oriented_tangent']
        for j,m in enumerate(mins):
            curve_rows.append(dict(seed=seed,kappa_rule=rule,density_factor=level,method=method,
                sample=j,origin_xi=float(ox[j]),width_xi=float(widths[j]),
                end_xi=float(ends[j]),net_norm=float(rnorm[j]),
                direction_speed=float(cr['directional_speed'][j]),
                width_implicit_derivative=float(m['implicit_wprime']),stationarity=float(m['stationarity']),
                stationary_curvature=float(m['curvature']),width_candidate_label=m['label'],
                width_candidate_count=m['num_candidates'],
                width_bounds=m['bound']))
    return

def write_plots(rows,curve_rows):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        for seed in ('original','autonomous'):
            fig,ax=plt.subplots(2,1,figsize=(10,7),sharex=True)
            for lev in (1,2,4):
                sub=[x for x in curve_rows if x['seed']==seed and x['kappa_rule']=='original' and x['density_factor']==lev and x['method']=='hermite']
                if not sub:continue
                x=np.array([r['origin_xi'] for r in sub])
                ax[0].plot(x,[r['width_xi'] for r in sub],label=f'{lev}x physical source')
                ax[1].plot(x,[r['direction_speed'] for r in sub],label=f'{lev}x physical source')
            ax[0].set_ylabel('Optimized window width')
            ax[1].set_ylabel('Unit 10D directional speed / xi')
            ax[1].set_xlabel('Window origin xi')
            ax[0].set_title('v18.53r — genuine 11D source-density convergence: '+seed)
            for a in ax:a.legend();a.grid(alpha=.25)
            fig.tight_layout();fig.savefig(OUT/f'dense_physical_source_valley_{seed}.png',dpi=140);plt.close(fig)
    except Exception as e:(OUT/'PLOT_WARNING.txt').write_text(str(e)+'\n')

def run_full(args,archive,old,centers,pins,stack,report):
    o,n,m,l,k,j,h,ncount,mcount=require_ancestry(archive,stack)
    loader=h.loader()
    physical_parent=loader.find_parent(args.parent_v1851s)
    native=args.native_npz or loader.find_npz('10d',args.base)
    sha_file=lambda p:S(Path(p).read_bytes())
    if not physical_parent or sha_file(physical_parent)!=o.CFG['physical_parent_source_sha256']:
        raise RuntimeError('HISTORICAL_V1851S_PARENT_MISSING_OR_SHA_DRIFT')
    if not native or sha_file(native)!=o.CFG['physical_native_source_sha256']:
        raise RuntimeError('HISTORICAL_10D_NATIVE_SOURCE_MISSING_OR_SHA_DRIFT')
    engine_stage,engine,contexts,meta=loader.stage_original(physical_parent)
    try:
        with np.load(native,allow_pickle=False) as f:
            xi,nh=j.xi_grid(f['common_xi'])
            if not np.array_equal(xi,old['common_xi']):
                if np.max(np.abs(xi-old['common_xi']))>1e-12:raise RuntimeError('BASELINE_GRID_MISMATCH')
            checkpoints={}
            geomrows=[];samples=[];phys=[];matches=[]
            baseline_by_case={}
            states={}
            for seed in ('original','autonomous'):
                ctx=contexts[seed]
                rerr,cerr,hist=h.native_53g_reconciliation(f,xi[:nh],seed,ctx,engine.omega)
                initial=np.asarray(hist[0],float)
                for rule in ('shared','original','autonomous'):
                    tag=seed+'__'+rule
                    baseline=archived_case(old,seed,rule)
                    baseline_by_case[tag]=baseline
                    if baseline['xi'].shape!=xi.shape:raise RuntimeError('PARENT_O_GRID_SHAPE_MISMATCH')
                    curve_case(seed,rule,1,baseline,centers[(seed,rule)],geomrows,samples,checkpoints)
                    for factor in CFG['refinement_factors']:
                        x=refine(xi,factor)
                        rt=min(float(h.C['refined_rtol']),float(h.C['integrator_rtol']))
                        at=min(float(h.C['refined_atol']),float(h.C['integrator_atol']))
                        maxstep=float(h.C['integrator_max_step_xi'])/factor
                        kp=float(k.K[rule])
                        sol=solve_ivp(lambda t,z:h.rhs(engine,ctx,kp,t,z),(float(x[0]),float(x[-1])),initial,
                               t_eval=x,method='DOP853',rtol=rt,atol=at,max_step=maxstep)
                        if not sol.success or not np.all(np.isfinite(sol.y)) or sol.y.shape!=(11,len(x)):
                            raise RuntimeError('DENSE_PHYSICAL_IVP_FAILED '+tag+' x'+str(factor)+' '+str(sol.message))
                        z=sol.y.T
                        maxrhs=float(h.validated_rhs(engine,ctx,kp,x,z))
                        d=h.primary_metrics(x,z,ctx,kp)
                        charge=float(d['metric']['charge_range'])
                        limit=float(k.CFG['three_window']['charge_abs_drift_tolerance'])
                        if charge>limit:raise RuntimeError('DENSE_PHYSICAL_CHARGE_DRIFT '+tag+' '+str(charge))
                        from three_window_geometry import accumulators
                        from torque_math import sigmoid
                        a=accumulators(x,d['t'],d['B_U'],d['B_C'])
                        _,C=m.primitive_from_acc(a)
                        B=np.stack([m.to10(v) for v in d['B_U']+d['B_C']])
                        rate=engine.omega/np.asarray(sigmoid(z[:,8]),float)
                        g=B*rate[:,None]
                        if C.shape!=(len(x),10) or not np.all(np.isfinite(g)) or not np.all(np.isfinite(C)):
                            raise RuntimeError('BAD_DENSE_10D_SOURCE')
                        data=dict(xi=x,C=C,g=g)
                        curve_case(seed,rule,factor,data,centers[seed,rule],geomrows,samples,checkpoints)
                        base=baseline_by_case[tag]
                        matchedC=C[::factor];matchedg=g[::factor]
                        if matchedC.shape!=base['C'].shape:raise RuntimeError('BAD_NESTED_C_MATCH')
                        deltaC=np.linalg.norm(matchedC-base['C'],axis=1)
                        deltag=np.linalg.norm(matchedg-base['g'],axis=1)
                        source_primitive=np.vstack([C[0],C[0]+np.cumsum((g[:-1]+g[1:])*.5*np.diff(x)[:,None],axis=0)])
                        drift=float(np.max(np.linalg.norm(source_primitive-C,axis=1)))
                        matches.append(dict(seed=seed,kappa_rule=rule,density_factor=factor,
                             max_primitive_difference_from_archived_1x=float(np.max(deltaC)),
                             rms_primitive_difference_from_archived_1x=float(np.sqrt(np.mean(deltaC**2))),
                             max_physical_source_difference_from_archived_1x=float(np.max(deltag)),
                             rms_physical_source_difference_from_archived_1x=float(np.sqrt(np.mean(deltag**2))),
                             max_sampled_source_integral_nodal_drift=drift,
                             source_max_norm=float(np.max(np.linalg.norm(g,axis=1)))))
                        phys.append(dict(seed=seed,kappa_rule=rule,density_factor=factor,
                             actual_ivp=True,n_state_samples=len(x),physical_model='AUTHENTIC_11D_DOP853',
                             rtol=rt,atol=at,max_step_xi=maxstep,nfev=int(sol.nfev),
                             max_rhs_error=maxrhs,charge_drift=charge,
                             baseline_referr=rerr,baseline_clockerr=cerr))
                        report['new_physical_ivps']+=1
                        if factor==2:states[tag]=z.copy()
                        else:
                            if tag not in states:raise RuntimeError('RUN_2X_BEFORE_4X')
                            z2=states.pop(tag)
                            if z[::2].shape!=z2.shape:raise RuntimeError('BAD_DENSE_STATE_INDEX')
                            phys[-1]['max_scaled_11d_difference_from_2x']=float(np.max(abs(z[::2]-z2)/np.maximum(1,abs(z[::2]))))
                        checkpoints[tag+'__x'+str(factor)+'__xi']=x
                        checkpoints[tag+'__x'+str(factor)+'__C10d']=C
                        checkpoints[tag+'__x'+str(factor)+'__g10d']=g
                        checkpoints[tag+'__x'+str(factor)+'__z11d']=z
                        print('V1853R_GENUINE_PHYSICAL_CASE',seed,rule,'factor',factor,'IVPS',report['new_physical_ivps'], 'nfev',sol.nfev,flush=True)
            if len(phys)!=CFG['require_converged_cases']*len(CFG['refinement_factors']):
                raise RuntimeError('INCOMPLETE_DENSE_PHYSICAL_INTEGRATIONS')
            # Cross-level metrics preserve case, observer convention, and parameter domain.
            cross=[]
            for seed in ('original','autonomous'):
                for rule in ('shared','original','autonomous'):
                    for method in ('hermite','pchip'):
                        subset={(r['density_factor']):r for r in geomrows if r['seed']==seed and r['kappa_rule']==rule and r['method']==method}
                        for a,b in ((1,2),(2,4),(1,4)):
                            aa=subset[a];bb=subset[b]
                            ta=f'{seed}__{rule}__x{a}__{method}'
                            tb=f'{seed}__{rule}__x{b}__{method}'
                            u=checkpoints[ta+'__u10d'];v=checkpoints[tb+'__u10d']
                            dot=np.clip(np.sum(u*v,axis=1),-1,1)
                            xarea=checkpoints[ta+'__oriented_polygon45d'];yarea=checkpoints[tb+'__oriented_polygon45d']
                            cos=float(np.clip(xarea@yarea/max(np.linalg.norm(xarea)*np.linalg.norm(yarea),1e-100),-1,1))
                            cross.append(dict(seed=seed,kappa_rule=rule,method=method,level_a=a,level_b=b,
                                max_width_absolute_difference=float(np.max(abs(checkpoints[ta+'__widths']-checkpoints[tb+'__widths']))),
                                max_unit_direction_angle_deg=float(np.rad2deg(np.max(np.arccos(dot)))),
                                oriented_area_angle_deg=float(np.rad2deg(np.arccos(cos))),
                                direct_arclength_absolute_difference=abs(aa['arclength_direct']-bb['arclength_direct']),
                                center_direction_speed_absolute_difference=abs(aa['direction_speed_center']-bb['direction_speed_center']),
                                min_residual_norm_absolute_difference=abs(aa['min_residual_norm']-bb['min_residual_norm'])))
            save_csv('physical_11d_integrations.csv',phys)
            save_csv('source_knot_convergence.csv',matches)
            save_csv('full10d_valley_convergence.csv',geomrows)
            save_csv('full10d_valley_samples.csv',samples)
            save_csv('cross_density_comparison.csv',cross)
            np.savez_compressed(OUT/'dense_11d_10d_primitives.npz',**checkpoints)
            write_plots(geomrows,samples)
            report['source_reintegration_distinguished_from_postprocessing']=True
            report['comparisons']=len(cross)
            report['densification_factors']=CFG['refinement_factors']
            report['max_source_difference_2x4x']=max(r.get('max_scaled_11d_difference_from_2x',0) for r in phys)
            report['max_source_integral_nodal_drift']=max(r['max_sampled_source_integral_nodal_drift'] for r in matches)
            report['max_2x_to_4x_width_difference']=max(r['max_width_absolute_difference'] for r in cross if r['level_a']==2 and r['level_b']==4)
            report['max_2x_to_4x_oriented_area_angle_deg']=max(r['oriented_area_angle_deg'] for r in cross if r['level_a']==2 and r['level_b']==4)
            report['max_2x_to_4x_direct_arclength_difference']=max(r['direct_arclength_absolute_difference'] for r in cross if r['level_a']==2 and r['level_b']==4)
            report['minimum_10d_norm_all_runs']=min(r['min_residual_norm'] for r in geomrows)
            report['formal_clock_extension_beyond_25p6']=True
            screen=CFG['near_zero_screen']
            if report['minimum_10d_norm_all_runs']<screen:
                report['scientific_decision']='GENUINE_DENSER_11D_SOURCE__NEAR_ZERO_SCREEN_ONLY__NO_MATERIAL_CLOSURE'
            elif (report['max_2x_to_4x_width_difference']<1e-3 and
                  report['max_2x_to_4x_oriented_area_angle_deg']<.2 and
                  report['max_2x_to_4x_direct_arclength_difference']<.02):
                report['scientific_decision']='GENUINE_DENSER_11D_SOURCE__NUMERICALLY_ROBUST_ORIENTED_VALLEY__NO_MATERIAL_HOLONOMY'
            else:
                report['scientific_decision']='GENUINE_DENSER_11D_SOURCE__DENSITY_SENSITIVITY_REMAINS__NO_MATERIAL_HOLONOMY'
    finally:
        engine_stage.cleanup()

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--preflight',action='store_true')
    p.add_argument('--smoke',action='store_true')
    p.add_argument('--base',type=Path,default=ROOT.parent)
    p.add_argument('--parent-v1851s',type=Path)
    p.add_argument('--native-npz',type=Path)
    p.add_argument('--require-parent',action='store_true')
    p.add_argument('--require-native',action='store_true')
    p.add_argument('--require-dynamic',action='store_true')
    args=p.parse_args()
    OUT.mkdir(exist_ok=True)
    report=dict(version='v18.53r',scientific_decision='NOT_EVALUATED',new_physical_ivps=0,inherited_physical_ivps=10)
    try:
        if args.preflight and (args.require_parent or args.require_native or args.require_dynamic):
            raise RuntimeError('PREFLIGHT_WITH_PHYSICAL_FLAGS_IS_FORBIDDEN')
        obytes,old,oldrows,pins,qsummary=authenticated_parents()
        report['authenticated_parent_manifest_counts']=pins
        centers=find_centers(oldrows)
        if args.smoke:
            h=archived_case(old,'original','shared')
            p=Primitive(h['xi'],h['C'],h['g'],'hermite')
            x,w=centers['original','shared']
            a=optimized_width(p,x,w-.05,w+.05)
            if not np.isfinite(a['norm']):raise RuntimeError('SMOKE_FAILED')
            print('V1853R_HISTORICAL_DATA_SMOKE_PASS',a['norm'],flush=True)
            return 0
        with ExitStack() as stack:
            if args.preflight:
                _,_,_,_,_,_,_,nc,mc=require_ancestry(obytes,stack)
                print('V1853R_PREFLIGHT_PASS_Q_P_O_MODEL_ANCESTRY',pins,nc,mc,flush=True)
                return 0
            run_full(args,obytes,old,centers,pins,stack,report)
        report['limitations']=[
          'Denser same-model 11D numerical integrations are not independent observed material trajectories.',
          'The clock/force source beyond xi=25.6 uses formal model continuation.',
          'The objective is the Euclidean norm of 10 bivector components, not a full-state metric.',
          'Finite 1x-2x-4x differences are numerical evidence, not certified continuous-time error bounds.',
          'No material holonomy, physical spin or true return is certified.' ]
        save_json('report.json',report)
        (OUT/'SCIENTIFIC_REPORT_SERVER.md').write_text('# Q550 — v18.53r physical source-density audit\n\n'
            +'**Scientific decision:** `'+report['scientific_decision']+'`\n\n'
            +'Actual new Q550 11D IVPs: '+str(report['new_physical_ivps'])+'.\n\n'
            +'```json\n'+json.dumps({k:v for k,v in report.items() if k.startswith(('max_','min_'))},indent=2)+'\n```\n\n'
            +'The tail is a formal clock continuation, with no physical closure or material holonomy certificate.\n')
        print('V1853R',report['scientific_decision'],'IVPS',report['new_physical_ivps'],flush=True)
        return 0
    except Exception as exc:
        report.update(scientific_decision='FAILED_OR_INCOMPLETE__NO_SCIENTIFIC_CERTIFICATE',error=str(exc),traceback=traceback.format_exc())
        save_json('report.json',report)
        print('V1853R_FAILED',str(exc),flush=True)
        return 2

if __name__=='__main__':sys.exit(main())
