#!/usr/bin/env python3
"""v18.53o: 10D oriented near-return valley and two-ended source audit.

All source and report text is English. Never substitute synthetic dynamics.
Heavy Q550 IVP execution is reserved for the user's server.
"""
from __future__ import annotations
import argparse,csv,hashlib,importlib.util,io,json,sys,tempfile,traceback,zipfile
from pathlib import Path
import numpy as np
from oriented_valley_geometry import derive_curves,net_geometry,curve_stats,fd_derivative_check,angle_deg
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results';OUT.mkdir(exist_ok=True)
CFG=json.loads((ROOT/'inputs_config.json').read_text())
PARENT=ROOT/'inputs/v18_53n_fix1_authenticated_server_sendback.zip'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
shab=lambda b:hashlib.sha256(b).hexdigest()

def output_csv(name,rows):
    if not rows:return
    fields=list(dict.fromkeys(f for r in rows for f in r))
    with (OUT/name).open('w',newline='') as fd:
        w=csv.DictWriter(fd,fieldnames=fields);w.writeheader();w.writerows(rows)

def output_json(name,rows):
    (OUT/name).write_text(json.dumps(rows,indent=2,allow_nan=False)+'\n')

def authenticate_parent_n(td):
    if sha(PARENT)!=CFG['parent_v1853n_fix1_sha256']:
        raise RuntimeError('PARENT_N_FIX1_SHA256_MISMATCH')
    with zipfile.ZipFile(PARENT) as z:
        if z.testzip() is not None:raise RuntimeError('PARENT_N_CRC_FAIL')
        members=z.namelist()
        roots={p.split('/')[0] for p in members if '/' in p}
        if len(roots)!=1:raise RuntimeError('INVALID_PARENT_ROOT_COUNT')
        root=next(iter(roots));prefix=root+'/'
        if any(Path(s).is_absolute() or '..' in Path(s).parts for s in members):
            raise RuntimeError('UNSAFE_PARENT_MEMBER')
        counted={}
        for key in ('SOURCE_MANIFEST.json','OUTPUT_MANIFEST.json'):
            man=json.loads(z.read(prefix+key))
            for rel,digest in man['files'].items():
                if shab(z.read(prefix+rel))!=digest:
                    raise RuntimeError('PARENT_N_MANIFEST_HASH_FAIL '+rel)
            counted[key]=len(man['files'])
        report=json.loads(z.read(prefix+'results/report.json'))
        if (report.get('scientific_decision')!=CFG['parent_required_decision'] or
            report.get('actual_Q550_ivps')!=10 or report.get('refined_Q550_ivps')!=4):
            raise RuntimeError('PARENT_N_NOT_AUTHENTIC_PHYSICAL_RUN')
        rows=list(csv.DictReader(io.StringIO(z.read(prefix+'results/best_full_10d_minimum_per_case.csv').decode())))
        if len(rows)!=CFG['parent_minima_count']:
            raise RuntimeError('PARENT_N_MINIMA_COUNT_MISMATCH')
        z.extractall(td)
    return Path(td)/root,rows,counted,report

def import_n(root):
    sys.path.insert(0,str(root/'scripts'))
    s=importlib.util.spec_from_file_location('q550_parent_v1853nfix1',root/'scripts/run_audit.py')
    mod=importlib.util.module_from_spec(s);sys.modules[s.name]=mod;s.loader.exec_module(mod)
    return mod

def scalar_source_at(xi,b10,rate):
    xi=np.asarray(xi,float);v=np.asarray(b10,float)*np.asarray(rate,float)[:,None]
    def interp(t):
        if t<xi[0]-1e-10 or t>xi[-1]+1e-10:raise ValueError('SOURCE_ENDPOINT_OUTSIDE_TRACE')
        return np.array([np.interp(t,xi,v[:,j]) for j in range(10)])
    return interp

def case_derived(seed,rule,xi,primitive,b10,rate,parent_best,n,checkpoint,results):
    from full_net_geometry import make_eval
    at,slope=make_eval(xi,primitive)
    po=float(parent_best['origin_xi']);pw=float(parent_best['width_xi'])
    r0=net_geometry(po,pw,at,slope)
    err=abs(r0['r']-float(parent_best['net_norm']))
    if err>CFG['minimum_replay_norm_abs_tol']:
        raise RuntimeError(f'PARENT_N_MINIMUM_REPLAY_MISMATCH: {seed} {rule} {err}')
    checker=scalar_source_at(xi,b10,rate)
    curves=derive_curves(xi,at,slope,po,pw,CFG['valley_halfwidth_origin_xi'],
        CFG['n_curve_points'],n.CFG['width_range_xi'],checker)
    for row in curves['fixed_rows']+curves['valley_rows']:
        row.update(seed=seed,kappa_rule=rule)
    results['fixed'].extend(curves['fixed_rows'])
    results['valley'].extend(curves['valley_rows'])
    def short_stats(kind):
        return dict(seed=seed,kappa_rule=rule,**curve_stats(curves,kind))
    results['summary'].extend([short_stats('fixed_width'),short_stats('valley')])
    testpoints=np.linspace(po-.33,po+.33,17)
    chk=[fd_derivative_check(p,pw,at,slope,CFG['derivative_fd_step_xi'],xi) for p in testpoints]
    tests=[t for t in chk if t['evaluated']]
    maxdef=max((t['relative_defect'] for t in tests),default=0.)
    if maxdef>CFG['derivative_check_tol_relative']:
        raise RuntimeError(f'ENDPOINT_SOURCE_DIFFERENTIAL_MISMATCH {seed} {rule}: {maxdef}')
    project_grad=np.array([np.dot(r0['net'],r0['jo']),np.dot(r0['net'],r0['jw'])])
    src0phys=checker(po);src1phys=checker(po+3*pw)
    row=dict(seed=seed,kappa_rule=rule,origin_xi=po,width_xi=pw,end_xi=r0['e'],
        parent_minimum_norm=float(parent_best['net_norm']),replayed_minimum_norm=r0['r'],
        parent_replay_abs_error=err,
        projected_gradient_origin=float(project_grad[0]),projected_gradient_width=float(project_grad[1]),
        gradient_norm=float(np.linalg.norm(project_grad)),
        residue_dot_start_source=float(r0['net']@r0['g0']),
        residue_dot_end_source=float(r0['net']@r0['g1']),
        residue_cos_start_source=r0['g0_cos'],residue_cos_end_source=r0['g1_cos'],
        residual_direction_turning_speed=r0['ang_speed'],
        residual_scalar_slope_fixed_width=r0['radial'],
        signed_turning_form_norm=float(np.linalg.norm(r0['form'])),
        physical_start_source_relative_defect=float(np.linalg.norm(r0['g0']-src0phys)/max(np.linalg.norm(src0phys),1e-30)),
        physical_end_source_relative_defect=float(np.linalg.norm(r0['g1']-src1phys)/max(np.linalg.norm(src1phys),1e-30)),
        numerical_fd_samples_evaluated=len(tests),numerical_fd_samples_skipped=len(chk)-len(tests),
        max_relative_fd_tangent_defect=maxdef,
        is_formal_continuation_past_25p6=bool(r0['e']>CFG['formal_continuation_boundary_xi']))
    results['orthogonality'].append(row)
    # Repeat evaluation using the exact same primitive under two interpolation conventions.
    p_at,p_slope=make_eval(xi,primitive,'pchip')
    p_geom=net_geometry(po,pw,p_at,p_slope)
    ss=dict(seed=seed,kappa_rule=rule,test='pchip_same_physical_primitive',
        norm_at_parent_min=p_geom['r'],absolute_norm_change=abs(p_geom['r']-r0['r']),
        full10d_unit_angle_deg=angle_deg(r0['u'],p_geom['u']),
        direction_speed_change=abs(p_geom['ang_speed']-r0['ang_speed']))
    results['sensitivity'].append(ss)
    # Every-other-sample source integration, not a new physical model.
    from three_window_geometry import accumulators
    select=np.arange(0,len(xi),2)
    if select[-1]!=len(xi)-1:select=np.r_[select,len(xi)-1]
    # The parent accumulators path is saved outside this function; caller provides coarse primitive.
    prefix=seed+'__'+rule
    checkpoint[prefix+'__primitive']=np.asarray(primitive,float)
    checkpoint[prefix+'__sampled_source_10d']=np.asarray(b10,float)
    checkpoint[prefix+'__clock_rate']=np.asarray(rate,float)
    return r0,curves

def plot_curves(results):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        for seed in ('original','autonomous'):
            fig,ax=plt.subplots(2,1,figsize=(8,7),sharex=True)
            for rule in ('original','autonomous','shared'):
                rows=[r for r in results['fixed'] if r['seed']==seed and r['kappa_rule']==rule]
                xx=np.array([r['origin_xi'] for r in rows])
                ax[0].plot(xx,[r['norm_10d'] for r in rows],label=rule)
                ax[1].plot(xx,[r['directional_angular_speed_per_xi'] for r in rows],label=rule)
            ax[0].set_ylabel('Full 10D residual norm')
            ax[1].set_ylabel('10D unit-vector angular speed / xi')
            ax[1].set_xlabel('Window origin xi')
            ax[0].set_title('v18.53o — unprojected residual through the near-return valley: '+seed)
            for a in ax:a.grid(alpha=.25);a.legend()
            fig.tight_layout();fig.savefig(OUT/f'full10d_residual_valley_oriented_turning_{seed}.png',dpi=145)
            plt.close(fig)
    except Exception as e:
        (OUT/'PLOT_WARNING.txt').write_text(str(e)+'\n')

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--preflight',action='store_true')
    p.add_argument('--base',type=Path,default=ROOT.parent)
    p.add_argument('--parent-v1851s',type=Path)
    p.add_argument('--native-npz',type=Path)
    p.add_argument('--require-parent',action='store_true')
    p.add_argument('--require-native',action='store_true')
    p.add_argument('--require-dynamic',action='store_true')
    args=p.parse_args()
    report={'version':'v18.53o','actual_Q550_ivps':0,'refined_Q550_ivps':0,
            'scientific_decision':'NOT_EVALUATED',
            'formal_clock_extension_after_xi_25p6':True,
            'objective':'UNPROJECTED_FULL_10D_RESIDUAL_ORIENTED_TRANSPORT'}
    failure=0
    try:
        with tempfile.TemporaryDirectory(prefix='q550_v1853o_parentn_') as nd:
            nroot,best,parent_pins,parent_report=authenticate_parent_n(nd)
            report['parent_n_manifests']=parent_pins
            best_by_case={(b['seed'],b['kappa_rule']):b for b in best}
            if len(best_by_case)!=6:raise RuntimeError('PARENT_BEST_CASES_NOT_UNIQUE')
            n=import_n(nroot)
            with tempfile.TemporaryDirectory(prefix='q550_v1853o_parentm_') as md:
                mroot,prior,pins,pr=n.authenticate_parent_m(md)
                m=n.import_m(mroot)
                with tempfile.TemporaryDirectory(prefix='q550_v1853o_parentl_') as ld:
                    lroot,old_l,lc=m.authenticate_l(ld)
                    l=m.import_l(lroot)
                    with tempfile.TemporaryDirectory(prefix='q550_v1853o_parentk_') as kd:
                        kroot,krows,krep=l.verified_parent(kd)
                        k=l.import_k(kroot)
                        j=k.load_j();k.verify_j_parent()
                        h=j.hist_module()
                        if any(not callable(getattr(h,fn,None)) for fn in
                            ('loader','solve','validated_rhs','primary_metrics','native_53g_reconciliation')):
                            raise RuntimeError('PHYSICAL_ANCESTRY_INCOMPLETE')
                        if args.preflight:
                            if args.require_parent or args.require_native or args.require_dynamic:
                                raise RuntimeError('PHYSICAL_FLAGS_FORBIDDEN_IN_PREFLIGHT')
                            print('V1853O_GENEALOGY_PREFLIGHT_PASS',parent_pins,flush=True)
                            return 0
                        loader=h.loader()
                        parent=loader.find_parent(args.parent_v1851s)
                        native=args.native_npz or loader.find_npz('10d',args.base)
                        if not parent or not Path(parent).is_file() or sha(parent)!=CFG['physical_parent_source_sha256']:
                            raise RuntimeError('ORIGINAL_PHYSICAL_PARENT_MISSING_OR_HASH_MISMATCH')
                        if not native or not Path(native).is_file() or sha(native)!=CFG['physical_native_source_sha256']:
                            raise RuntimeError('NATIVE_Q550_SOURCE_MISSING_OR_HASH_MISMATCH')
                        from three_window_geometry import accumulators
                        from torque_math import sigmoid
                        from full_net_geometry import make_eval
                        staged,engine,contexts,metadata=loader.stage_original(parent)
                        try:
                            results=dict(fixed=[],valley=[],summary=[],orthogonality=[],sensitivity=[],ivps=[])
                            checkpoints={}
                            with np.load(native,allow_pickle=False) as f:
                                xi,nh=j.xi_grid(f['common_xi'])
                                checkpoints['common_xi']=xi
                                n.check_box(xi,n.CFG['origin_range_xi'],n.CFG['width_range_xi'])
                                for seed in ('original','autonomous'):
                                    ctx=contexts[seed]
                                    referr,clockerr,history=h.native_53g_reconciliation(f,xi[:nh],seed,ctx,engine.omega)
                                    z0=np.asarray(history[0],float)
                                    for rule in ('shared','original','autonomous'):
                                        kp=float(k.K[rule]);z,info=h.solve(engine,ctx,kp,xi,z0)
                                        report['actual_Q550_ivps']+=1
                                        rhs=h.validated_rhs(engine,ctx,kp,xi,z)
                                        d=h.primary_metrics(xi,z,ctx,kp)
                                        if d['metric']['charge_range']>k.CFG['three_window']['charge_abs_drift_tolerance']:
                                            raise RuntimeError('PHYSICAL_CHARGE_DRIFT')
                                        acc=accumulators(xi,d['t'],d['B_U'],d['B_C'])
                                        _,c=m.primitive_from_acc(acc)
                                        b10=np.stack([m.to10(v) for v in d['B_U']+d['B_C']])
                                        rate=engine.omega/np.asarray(sigmoid(z[:,8]),float)
                                        target=best_by_case[(seed,rule)]
                                        geom,curve=case_derived(seed,rule,xi,c,b10,rate,target,n,checkpoints,results)
                                        # Coarser quadrature recomputed from actual 11D solver output.
                                        sel=np.arange(0,len(xi),2)
                                        if sel[-1]!=len(xi)-1:sel=np.r_[sel,len(xi)-1]
                                        ac=accumulators(xi[sel],d['t'][sel],d['B_U'][sel],d['B_C'][sel])
                                        _,c2=m.primitive_from_acc(ac)
                                        c_at,c_slope=make_eval(xi[sel],c2)
                                        cc=net_geometry(geom['o'],geom['w'],c_at,c_slope)
                                        results['sensitivity'].append(dict(seed=seed,kappa_rule=rule,
                                            test='every_other_sample_source_quadrature',
                                            norm_at_parent_min=cc['r'],absolute_norm_change=abs(cc['r']-geom['r']),
                                            full10d_unit_angle_deg=angle_deg(geom['u'],cc['u']),
                                            direction_speed_change=abs(cc['ang_speed']-geom['ang_speed'])))
                                        if rule in ('original','autonomous'):
                                            ztight,_=h.solve(engine,ctx,kp,xi,z0,tight=True)
                                            report['actual_Q550_ivps']+=1
                                            report['refined_Q550_ivps']+=1
                                            dev=float(np.max(np.abs(z-ztight)/np.maximum(1.,np.abs(ztight))))
                                            if dev>CFG['max_refined_scaled_trajectory_deviation']:
                                                raise RuntimeError('TIGHT_Q550_TRAJECTORY_DEVIATION')
                                            dr=h.primary_metrics(xi,ztight,ctx,kp)
                                            ar=accumulators(xi,dr['t'],dr['B_U'],dr['B_C'])
                                            _,cr=m.primitive_from_acc(ar)
                                            tt,dd=make_eval(xi,cr)
                                            rr=net_geometry(geom['o'],geom['w'],tt,dd)
                                            results['sensitivity'].append(dict(seed=seed,kappa_rule=rule,
                                                test='tight_11d_IVP_same_parent_minimum',
                                                norm_at_parent_min=rr['r'],absolute_norm_change=abs(rr['r']-geom['r']),
                                                full10d_unit_angle_deg=angle_deg(geom['u'],rr['u']),
                                                direction_speed_change=abs(rr['ang_speed']-geom['ang_speed']),
                                                max_scaled_11d_trajectory_difference=dev))
                                        results['ivps'].append(dict(seed=seed,kappa_rule=rule,kappa0=kp,
                                            nfev=info['nfev'],max_rhs_error=rhs,
                                            charge_drift=d['metric']['charge_range'],
                                            native_reference_error=referr,native_clock_error=clockerr))
                                        print('V1853O_PHYSICAL_CASE',seed,rule,'IVPS',report['actual_Q550_ivps'],flush=True)
                            if (report['actual_Q550_ivps']!=CFG['expected_actual_ivps'] or
                                report['refined_Q550_ivps']!=CFG['expected_refined_ivps'] or
                                len(results['orthogonality'])!=6):
                                raise RuntimeError('INCOMPLETE_Q550_PHYSICAL_EXPERIMENT')
                            for name,key in [('fixed_width_oriented_transport.csv','fixed'),
                                ('local_valley_oriented_transport.csv','valley'),
                                ('minimum_two_endpoint_orthogonality.csv','orthogonality'),
                                ('valley_curve_summary.csv','summary'),
                                ('quadrature_pchip_refinement_sensitivity.csv','sensitivity'),
                                ('genuine_physical_ivps.csv','ivps')]:output_csv(name,results[key])
                            np.savez_compressed(OUT/'derived_primitives.npz',**checkpoints)
                            plot_curves(results)
                            rows=results['orthogonality']
                            report['summary']={
                                'n_fixed_width_curve_rows':len(results['fixed']),
                                'n_valley_follow_curve_rows':len(results['valley']),
                                'max_parent_minimum_replay_abs_error':max(x['parent_replay_abs_error'] for x in rows),
                                'max_fd_tangent_relative_error':max(x['max_relative_fd_tangent_defect'] for x in rows),
                                'maximum_minimum_full10d_norm':max(x['replayed_minimum_norm'] for x in rows),
                                'minimum_minimum_full10d_norm':min(x['replayed_minimum_norm'] for x in rows),
                                'minimum_turning_speed_at_parent_min':min(x['residual_direction_turning_speed'] for x in rows),
                                'maximum_turning_speed_at_parent_min':max(x['residual_direction_turning_speed'] for x in rows),
                                'max_end_orthogonality_cos_abs':max(abs(x['residue_cos_end_source']) for x in rows),
                                'max_start_orthogonality_cos_abs':max(abs(x['residue_cos_start_source']) for x in rows),
                                'max_interpolant_vs_sampled_physical_source_relative_defect':max(max(x['physical_start_source_relative_defect'],x['physical_end_source_relative_defect']) for x in rows),
                                'any_width_optimizer_active_bound':any(x['num_width_bounds']>0 for x in results['summary']),
                                'all_minimum_endpoints_in_formal_continuation':all(x['is_formal_continuation_past_25p6'] for x in rows),
                            }
                            any_zero=any(x['replayed_minimum_norm']<CFG['numerical_closure_screen'] for x in rows)
                            report['scientific_decision']=('NUMERICAL_CLOSURE_SCREEN__PHYSICAL_RETURN_NOT_CERTIFIED'
                                if any_zero else 'AUTHENTIC_ORIENTED_10D_VALLEY_AND_BOUNDARY_GEOMETRY__NO_MATERIAL_HOLONOMY')
                            report['limits']=[
                                'The 10D bivector norm is not a complete 11D state-distance diagnostic.',
                                'A turning orientation in the 10D residual coordinate space is not physical angular spin or material holonomy.',
                                'The full 10D data after xi=25.6 uses a formal model clock continuation.',
                                'Near-valley boundary source differences are exact for the chosen cumulative interpolant, not a theorem of continuous-source derivatives.',
                                'The local width optimizer may be constrained at domain edges; report all such cases.',
                            ]
                        finally:
                            staged.cleanup()
    except Exception as exc:
        failure=2;report['scientific_decision']='FAILED_OR_INCOMPLETE__NO_SCIENTIFIC_CERTIFICATE'
        report['error']=str(exc);report['traceback']=traceback.format_exc()
    finally:
        output_json('report.json',report)
        (OUT/'SCIENTIFIC_REPORT_SERVER.md').write_text('# Q550 v18.53o — Full-10D Oriented Valley Audit\n\n'
            '**Scientific decision:** `'+report['scientific_decision']+'`\n\n'
            'Actual physical Q550 11D integrations: '+str(report['actual_Q550_ivps'])+'.\n\n'
            'This study does not certify full material closure, physical spin, or holonomy.\n\n'
            '```json\n'+json.dumps(report.get('summary',{}),indent=2)+'\n```\n\n'
            'Error: '+str(report.get('error','none'))+'\n')
        print('V1853O',report['scientific_decision'],'IVPS',report['actual_Q550_ivps'],flush=True)
    if args.require_dynamic and not report['scientific_decision'].startswith(('AUTHENTIC_ORIENTED_10D','NUMERICAL_CLOSURE_SCREEN')):
        return max(3,failure)
    return failure

if __name__=='__main__':sys.exit(main())
