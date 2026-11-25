#!/usr/bin/env python3
"""Q550 v18.53s - authentic nested quadrature/physical-step separation.

Default execution on user server performs six genuinely tighter 11D IVPs at the
same 4x evaluation nodes. Preflight and quadrature-only modes never run IVPs.
All scientific conclusions are restricted to the formal Q550 clock model.
"""
from __future__ import annotations
import argparse,csv,hashlib,importlib.util,io,json,sys,zipfile,tempfile,traceback
from contextlib import ExitStack
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from quadrature_math import (cumulative_trapezoid,baseline_quadrature,source_pchip_integral,
                             decompose_nested_2x4x,max_norm,norm_rows)

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'
CFG=json.loads((ROOT/'inputs_config.json').read_text())
PARENT=ROOT/'inputs'/CFG['parent_filename']
SHA=lambda b:hashlib.sha256(b).hexdigest()


def csv_write(name,rows):
    if not rows:raise RuntimeError('EMPTY_OUTPUT_'+name)
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with (OUT/name).open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=keys)
        writer.writeheader();writer.writerows(rows)


def module_from(path,name):
    path=Path(path)
    if str(path.parent) not in sys.path:sys.path.insert(0,str(path.parent))
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec)
    sys.modules[name]=m
    spec.loader.exec_module(m)
    return m


def secure_zip(z,label):
    if z.testzip() is not None:raise RuntimeError(label+'_CRC_FAILURE')
    names=z.namelist()
    if any(Path(p).is_absolute() or '..' in Path(p).parts for p in names):raise RuntimeError(label+'_UNSAFE_PATH')
    roots={p.split('/')[0] for p in names if '/' in p}
    if len(roots)!=1:raise RuntimeError(label+'_BAD_ROOT')
    prefix=next(iter(roots))+'/'
    counts={}
    for manifest in ('SOURCE_MANIFEST.json','OUTPUT_MANIFEST.json'):
        d=json.loads(z.read(prefix+manifest))['files']
        for rel,s in d.items():
            if SHA(z.read(prefix+rel))!=s:raise RuntimeError(label+'_HASH_FAIL_'+rel)
        counts[manifest]=len(d)
    return prefix,counts


def authenticated_parent(stack):
    if not PARENT.is_file() or SHA(PARENT.read_bytes())!=CFG['parent_sha256']:
        raise RuntimeError('V1853R_PARENT_SHA256_MISMATCH')
    with zipfile.ZipFile(PARENT) as z:
        prefix,counts=secure_zip(z,'R')
        report=json.loads(z.read(prefix+'results/report.json'))
        if report['scientific_decision']!=CFG['parent_r_decision'] or report.get('new_physical_ivps')!=12:
            raise RuntimeError('WRONG_PARENT_R_PHYSICAL_PROVENANCE')
        work=Path(stack.enter_context(tempfile.TemporaryDirectory(prefix='q550_v1853s_r_')))
        z.extractall(work)
    rr=work/prefix.rstrip('/')
    r=module_from(rr/'scripts/run_audit.py','q550_v1853s_authentic_r')
    archive,old,old_fixed,pins,_=r.authenticated_parents()
    centers=r.find_centers(old_fixed)
    return r,rr,old,centers,pins,counts,report,archive


def read_r_npz(rr):
    path=rr/'results/dense_11d_10d_primitives.npz'
    return np.load(path,allow_pickle=False)


def histories_for(data,old,r,seed,rule):
    tag=seed+'__'+rule
    base=r.archived_case(old,seed,rule)
    histories={1:(base['xi'],base['C'],base['g'])}
    for k in (2,4):
        histories[k]=tuple(data[tag+f'__x{k}__'+ending] for ending in ('xi','C10d','g10d'))
    x1=histories[1][0]
    if len(histories[2][0])!=2*(len(x1)-1)+1 or len(histories[4][0])!=4*(len(x1)-1)+1:
        raise RuntimeError('INVALID_NESTED_GRID_LENGTH')
    for k in (1,2):
        if not np.allclose(histories[4][0][::4//k],histories[k][0],atol=1e-12,rtol=0):
            raise RuntimeError('INVALID_NESTED_GRID_COORDINATES')
    return histories


def quadrature_case(hs,seed,rule,records,candidates):
    x1,C1,g1=hs[1];x2,C2,g2=hs[2];x4,C4,g4=hs[4]
    c0=np.asarray(C4[0],float)
    # Quadrature-only hierarchy: identical g4 at all coarser knots.
    q1=baseline_quadrature(x4[::4],g4[::4],1,'trapezoid',c0)
    q2=baseline_quadrature(x4[::2],g4[::2],2,'trapezoid',c0)
    q4=baseline_quadrature(x4,g4,4,'trapezoid',c0)
    s2=baseline_quadrature(x4[::2],g4[::2],2,'simpson',c0)
    s4=baseline_quadrature(x4,g4,4,'simpson',c0)
    rich_t=(4*q4-q2)/3  # heuristic truncation cancellation for smooth sources
    rich_s=(16*s4-s2)/15
    for a,b,left,right in ((1,2,q1,q2),(2,4,q2,q4)):
        d=left-right
        val,rms=norm_rows(d)
        records.append(dict(seed=seed,kappa_rule=rule,contrast=f'T{a}_vs_T{b}_same_g4',
                            max_cumulative_10d_difference=val,rms_cumulative_10d_difference=rms,
                            change_in_source_samples=False,method='trapezoid'))
    for label,mat in [('S2_vs_S4_same_g4',s2-s4),('S4_vs_Richardson_S',s4-rich_s),
                      ('T4_vs_S4_same_g4',q4-s4),('T4_vs_Richardson_T',q4-rich_t),
                      ('S4_vs_Richardson_T',s4-rich_t),('archived_C4_vs_T4',C4[::4]-q4),
                      ('archived_C2_vs_T2',C2[::2]-baseline_quadrature(x2,g2,2,'trapezoid',c0)),
                      ('archived_C1_vs_T1',C1-baseline_quadrature(x1,g1,1,'trapezoid',c0))]:
        value,rms=norm_rows(mat)
        records.append(dict(seed=seed,kappa_rule=rule,contrast=label,
                            max_cumulative_10d_difference=value,rms_cumulative_10d_difference=rms,
                            change_in_source_samples=False,method='mixed'))
    q1own=baseline_quadrature(x1,g1,1,'trapezoid',c0)
    q2own=baseline_quadrature(x2,g2,2,'trapezoid',c0)
    q4own=q4
    for label,mat in [('T1_own_vs_T2_own',q1own-q2own),('T2_own_vs_T4_own',q2own-q4own)]:
        value,rms=norm_rows(mat)
        records.append(dict(seed=seed,kappa_rule=rule,contrast=label,
                            max_cumulative_10d_difference=value,rms_cumulative_10d_difference=rms,
                            change_in_source_samples=True,method='trapezoid'))
    split=decompose_nested_2x4x(x2,g2,x4,g4,c0)
    if split['identity_max_error']>5e-13:raise RuntimeError('TELESCOPING_IDENTITY_FAILURE')
    for key in ('total','quadrature_on_common_g4','physical_sample_change'):
        val,rms=norm_rows(split[key][::2])
        records.append(dict(seed=seed,kappa_rule=rule,contrast='2x_4x_decomposition_'+key,
                            max_cumulative_10d_difference=val,rms_cumulative_10d_difference=rms,
                            change_in_source_samples=key=='physical_sample_change',method='exact_telescope'))
    sourceP=source_pchip_integral(x4,g4,c0)
    val,rms=norm_rows(sourceP[::4]-s4)
    records.append(dict(seed=seed,kappa_rule=rule,contrast='source_PCHIP_integral_vs_S4',
                        max_cumulative_10d_difference=val,rms_cumulative_10d_difference=rms,
                        change_in_source_samples=False,method='source_interpolation'))
    candidates.update({'archived':(x4,C4,g4),'trapezoid4x':(x4,cumulative_trapezoid(x4,g4,c0),g4),
                       'source_pchip_integral':(x4,sourceP,g4)})
    return split['identity_max_error']


def geometry_case(r,case,center,label,history,records,store):
    x,C,g=history
    p=r.Primitive(x,C,g,'hermite')
    oc,wc=center
    origins=np.linspace(oc-CFG['origin_halfspan'],oc+CFG['origin_halfspan'],CFG['curve_points'])
    wlo,wup=wc-CFG['width_halfband'],wc+CFG['width_halfband']
    if origins[0]<x[0] or origins[-1]+3*wup>x[-1]:raise RuntimeError('GEOMETRY_EXCEEDS_AUTHENTIC_GRID')
    minima=[r.optimized_width(p,float(o),wlo,wup,nsearch=CFG['search_brackets']) for o in origins]
    w=np.array([item['w'] for item in minima])
    implicit=np.array([item['implicit_wprime'] for item in minima])
    cr=r.curve(p,origins,w,implicit)
    mid=(len(origins)-1)//2
    record=dict(seed=case[0],kappa_rule=case[1],variant=label,
                origin_center=float(oc),width_center=float(w[mid]),
                min_norm_10d=float(np.min(cr['norm'])),norm_center=float(cr['norm'][mid]),
                direction_speed_center=float(cr['directional_speed'][mid]),
                arclength_direct=float(cr['direct']),arclength_tangent=float(cr['chain']),
                arclength_direct_minus_tangent=float(cr['arc_diff']),
                width_total_variation=float(cr['width_total_variation']),
                width_sign_changes=int(cr['width_sign_changes']),
                num_boundary_minima=int(sum(m['bound'] for m in minima)),
                oriented_polygon_norm=float(np.linalg.norm(cr['oriented_polygon'])),
                max_angular_step=float(cr['max_angular_step']),
                all_formal_clock_endpoints=bool(np.all(origins+3*w>CFG['formal_clock_boundary'])))
    records.append(record)
    key=case[0]+'__'+case[1]+'__'+label
    store[key+'__origins']=origins
    store[key+'__widths']=w
    store[key+'__u10d']=cr['u']
    store[key+'__oriented45d']=cr['oriented_polygon']
    return record


def generate_plots(rows):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        for seed in CFG['seeds']:
            rr=[r for r in rows if r['seed']==seed and r['kappa_rule']=='original']
            if not rr:continue
            fig,ax=plt.subplots(1,2,figsize=(12,4.4))
            ax[0].bar([r['variant'] for r in rr],[r['min_norm_10d'] for r in rr])
            ax[0].set_ylabel('Minimum Euclidean 10D residual norm')
            ax[1].bar([r['variant'] for r in rr],[r['arclength_direct'] for r in rr])
            ax[1].set_ylabel('Direct angular arclength (radians)')
            for a in ax:
                a.tick_params(axis='x',rotation=60)
                a.grid(axis='y',alpha=.25)
            fig.suptitle('v18.53s: quadrature/step sensitivity — '+seed)
            fig.tight_layout()
            fig.savefig(OUT/f'quadrature_geometry_separation_{seed}.png',dpi=140)
            plt.close(fig)
    except Exception as exc:
        (OUT/'PLOT_WARNING.txt').write_text(str(exc)+'\n')


def full_physical_reference(r,archive,stack,refs,args):
    # Only this function executes new 11D physics.
    o,n,m,l,k,j,h,nc,mc=r.require_ancestry(archive,stack)
    loader=h.loader()
    physical_parent=loader.find_parent(args.parent_v1851s)
    native=args.native_npz or loader.find_npz('10d',args.base)
    if not physical_parent or SHA(Path(physical_parent).read_bytes())!=o.CFG['physical_parent_source_sha256']:
        raise RuntimeError('HISTORICAL_PHYSICAL_PARENT_MISSING_OR_CHANGED')
    if not native or SHA(Path(native).read_bytes())!=o.CFG['physical_native_source_sha256']:
        raise RuntimeError('HISTORICAL_10D_NATIVE_MISSING_OR_CHANGED')
    stage,engine,contexts,meta=loader.stage_original(physical_parent)
    try:
        from torque_math import sigmoid
        with np.load(native,allow_pickle=False) as f:
            records=[];runs={}; nfevs=[]
            for seed in CFG['seeds']:
                ctx=contexts[seed]
                for rule in CFG['kappa_rules']:
                    tag=seed+'__'+rule
                    xi=refs[tag]['xi4']
                    xx,nh=j.xi_grid(f['common_xi'])
                    if not np.allclose(xi[::4],xx,rtol=0,atol=1e-12):raise RuntimeError('CLOCK_BASE_GRID_MISMATCH')
                    rerr,cerr,hist=h.native_53g_reconciliation(f,xx[:nh],seed,ctx,engine.omega)
                    initial=np.asarray(hist[0],float)
                    kp=float(k.K[rule])
                    base_rt=min(float(h.C['refined_rtol']),float(h.C['integrator_rtol']))
                    base_at=min(float(h.C['refined_atol']),float(h.C['integrator_atol']))
                    rt=base_rt*CFG['physical_reference_rtol_ratio']
                    at=base_at*CFG['physical_reference_atol_ratio']
                    hmax=float(h.C['integrator_max_step_xi'])/4*CFG['physical_reference_max_step_ratio']
                    sol=solve_ivp(lambda t,z:h.rhs(engine,ctx,kp,t,z),
                                  (float(xi[0]),float(xi[-1])),initial,
                                  t_eval=xi,method='DOP853',rtol=rt,atol=at,max_step=hmax)
                    if not sol.success or sol.y.shape!=(11,len(xi)) or not np.all(np.isfinite(sol.y)):
                        raise RuntimeError('TIGHT_11D_IVP_FAILURE_'+tag+' '+str(sol.message))
                    z=sol.y.T
                    rhs_err=float(h.validated_rhs(engine,ctx,kp,xi,z))
                    data=h.primary_metrics(xi,z,ctx,kp)
                    charge=float(data['metric']['charge_range'])
                    if charge>float(k.CFG['three_window']['charge_abs_drift_tolerance']):
                        raise RuntimeError('TIGHT_11D_CHARGE_DRIFT_'+tag)
                    g=np.stack([m.to10(v) for v in data['B_U']+data['B_C']])*(engine.omega/np.asarray(sigmoid(z[:,8]),float))[:,None]
                    zold=refs[tag]['z4'];gold=refs[tag]['g4']
                    if zold.shape!=z.shape or gold.shape!=g.shape:raise RuntimeError('TIGHT_SAME_GRID_MISMATCH')
                    cdiff=max_norm(cumulative_trapezoid(xi,g)-cumulative_trapezoid(xi,gold))
                    source_diff=max_norm(g-gold)
                    state_diff=max_norm(z-zold)
                    records.append(dict(seed=seed,kappa_rule=rule,nfev=int(sol.nfev),
                                        same_output_knots=len(xi),rtol=rt,atol=at,max_step_xi=hmax,
                                        max_state_11d_difference_from_r4x=state_diff,
                                        max_source_10d_difference_from_r4x=source_diff,
                                        max_trapezoid_primitive_difference_from_r4x=cdiff,
                                        rhs_residual=rhs_err,charge_range=charge,
                                        model='AUTHENTIC_11D_DOP853_SAME_EVALUATION_GRID'))
                    runs[tag]=(xi,g,z)
                    nfevs.append(sol.nfev)
                    print('V1853S_AUTHENTIC_TIGHT_PHYSICAL_CASE',seed,rule,'IVPS',len(records),'NFEV',sol.nfev,flush=True)
            return runs,records
    finally:
        stage.cleanup()


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--preflight',action='store_true')
    parser.add_argument('--quadrature-only',action='store_true',help='No new physical IVPs; diagnostic only, no full physical certificate')
    parser.add_argument('--smoke',action='store_true')
    parser.add_argument('--base',type=Path,default=ROOT.parent)
    parser.add_argument('--parent-v1851s',type=Path)
    parser.add_argument('--native-npz',type=Path)
    parser.add_argument('--require-parent',action='store_true')
    parser.add_argument('--require-native',action='store_true')
    parser.add_argument('--require-dynamic',action='store_true')
    args=parser.parse_args()
    OUT.mkdir(exist_ok=True)
    report=dict(version='v18.53s',scientific_decision='NOT_EVALUATED',
                new_physical_ivps=0,inherited_r_physical_ivps=12)
    try:
        if args.preflight and (args.require_parent or args.require_native or args.require_dynamic):
            raise RuntimeError('PREFLIGHT_PHYSICAL_FLAGS_NOT_ALLOWED')
        with ExitStack() as stack:
            r,rr,old,centers,pins,counts,parent_report,archive=authenticated_parent(stack)
            report['r_parent_manifests']=counts
            report['r_parent_sha256']=CFG['parent_sha256']
            if args.preflight:
                # Verify complete nested model ancestry, without physical IVP calls.
                _,_,_,_,_,_,_,nc,mc=r.require_ancestry(archive,stack)
                print('V1853S_PREFLIGHT_PASS_R_Q_P_O_N_M_L_K_J_H',counts,pins,nc,mc,flush=True)
                return 0
            with read_r_npz(rr) as data:
                quad=[];geo=[];splits=[];diag=[];store={};refs={}
                for seed in CFG['seeds']:
                    for rule in CFG['kappa_rules']:
                        tag=seed+'__'+rule
                        hs=histories_for(data,old,r,seed,rule)
                        history={}
                        epsilon=quadrature_case(hs,seed,rule,quad,history)
                        splits.append(epsilon)
                        refs[tag]=dict(xi4=hs[4][0],g4=hs[4][2],z4=np.array(data[tag+'__x4__z11d']))
                        if args.smoke:
                            p=r.Primitive(*history['source_pchip_integral'],'hermite')
                            val=r.optimized_width(p,centers[seed,rule][0],centers[seed,rule][1]-.05,centers[seed,rule][1]+.05)
                            print('V1853S_SMOKE_AUTHENTIC_QUADRATURE_PASS',tag,round(val['norm'],9),flush=True)
                            continue
                        for label,hist in history.items():
                            geometry_case(r,(seed,rule),centers[seed,rule],label,hist,geo,store)
                if args.smoke:return 0
                if not args.quadrature_only:
                    runs,physical=full_physical_reference(r,archive,stack,refs,args)
                    for seed in CFG['seeds']:
                        for rule in CFG['kappa_rules']:
                            tag=seed+'__'+rule
                            xi,g,z=runs[tag]
                            C=source_pchip_integral(xi,g)
                            geometry_case(r,(seed,rule),centers[seed,rule],
                                          'tight_11d_source_pchip_integral',(xi,C,g),geo,store)
                            store[tag+'__tight_xi']=xi
                            store[tag+'__tight_g10d']=g
                            store[tag+'__tight_z11d']=z
                            del refs[tag]['z4']
                            # 4x step effect on a *fixed* xi grid and identical quadrature.
                            srcdiff=max_norm(g-refs[tag]['g4'])
                            tadiff=max_norm(cumulative_trapezoid(xi,g)-cumulative_trapezoid(xi,refs[tag]['g4']))
                            diag.append(dict(seed=seed,kappa_rule=rule,contrast='tight_vs_r4x_same_4x_grid',
                                max_source_10d_difference=srcdiff,max_cumulative_trapezoid_difference=tadiff))
                    report['new_physical_ivps']=len(physical)
                    csv_write('physical_step_reference.csv',physical)
                    csv_write('same_grid_physical_source_comparison.csv',diag)
                else:
                    report['new_physical_ivps']=0
                csv_write('quadrature_order_bias.csv',quad)
                csv_write('valley_geometry_quadrature_bias.csv',geo)
                np.savez_compressed(OUT/'oriented_valley_quadrature_comparison.npz',**store)
                generate_plots(geo)
                report['max_telescoping_identity_error']=max(splits)
                report['number_quadrature_comparisons']=len(quad)
                report['number_valley_geometries']=len(geo)
                report['minimum_full_10d_residual']=float(min(v['min_norm_10d'] for v in geo))
                report['max_width_sign_changes']=int(max(v['width_sign_changes'] for v in geo))
                report['max_physically_constrained_quadrature_difference']=float(max(v['max_cumulative_10d_difference'] for v in quad if v['contrast']=='T4_vs_S4_same_g4'))
                report['formal_clock_continuation_beyond_25p6']=True
                report['scientific_decision']=(
                    'AUTHENTIC_QUADRATURE_ONLY__NO_NEW_PHYSICAL_STEP_CERTIFICATE'
                    if args.quadrature_only else
                    'AUTHENTIC_FIXED_GRID_QUADRATURE_VS_TIGHT_11D_STEP_SEPARATION__NO_MATERIAL_HOLONOMY')
                report['limitations']=[
                    'Nested trapezoid/Simpson differences are empirical error indicators, not rigorous upper bounds.',
                    'Source PCHIP integral is a reconstructed source, not the exact continuum model integral.',
                    'The same-grid tighter DOP853 replay isolates a numerical step/tolerance sensitivity, not a unique exact ODE error.',
                    'The Euclidean 10D bivector metric is not a full physical-state distance.',
                    'All relevant late endpoints belong to formal clock continuation beyond xi=25.6.',
                    'No material closure, true state recurrence, or material holonomy is certified.' ]
                (OUT/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
                (OUT/'SCIENTIFIC_REPORT_SERVER.md').write_text(
                    '# Q550 v18.53s — Numerical quadrature/ODE separation\n\n'
                    +'Decision: `'+report['scientific_decision']+'`\n\n'
                    +'New 11D integrations: '+str(report['new_physical_ivps'])+'\n\n'
                    +'See `quadrature_order_bias.csv`, `valley_geometry_quadrature_bias.csv`, '
                    +'and, for full runs, `physical_step_reference.csv`.\n\n'
                    +'No continuous-time or material holonomy certificate is implied.\n')
                print('V1853S',report['scientific_decision'],'NEW_IVPS',report['new_physical_ivps'],flush=True)
                return 0
    except Exception as exc:
        report.update(scientific_decision='FAILED_OR_INCOMPLETE__NO_SCIENTIFIC_CERTIFICATE',
                      error=str(exc),traceback=traceback.format_exc())
        (OUT/'report.json').write_text(json.dumps(report,indent=2)+'\n')
        print('V1853S_FAILED',exc,flush=True)
        return 2

if __name__=='__main__':sys.exit(main())
