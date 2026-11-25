#!/usr/bin/env python3
"""Q550 v18.53l genuine Q550 IVP replay + projector-geometric discrimination.

Authenticated v18.53k is an immutable input; the original 10D Q550 engine
is replayed on the research server. No synthetic or cached CSV is called a
new physical integration. All CSV projections are in the declared kinetic
Frobenius embedding; transported comparisons use a specified radial lift.
"""
from __future__ import annotations
import argparse,csv,hashlib,importlib.util,io,json,sys,tempfile,traceback,zipfile
from pathlib import Path
from collections import defaultdict
import numpy as np
from frame_projector_geometry import (to10,svd_rank2,cosine,projector_angle,
    transport_curve,dynamic_projector,secant_product_decomposition,
    derivative_sample)
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results';OUT.mkdir(exist_ok=True)
CFG=json.loads((ROOT/'inputs_config.json').read_text())
PARENT=ROOT/'inputs/v18_53k_authenticated_server_sendback.zip'

def sha_bytes(data):return hashlib.sha256(data).hexdigest()
def sha(p):return sha_bytes(Path(p).read_bytes())
def dump(name,obj):(OUT/name).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')
def write_csv(name,rows):
    if not rows:return
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with (OUT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)

def verified_parent(temp_dir):
    if sha(PARENT)!=CFG['parent_k_sha256']:raise RuntimeError('PARENT_K_PIN_FAIL')
    with zipfile.ZipFile(PARENT) as z:
        if z.testzip() is not None:raise RuntimeError('PARENT_K_CRC_FAIL')
        names=z.namelist();roots={n.split('/')[0] for n in names if '/' in n}
        if len(roots)!=1:raise RuntimeError('PARENT_K_MANY_ROOTS')
        prefix=next(iter(roots))+'/'
        for lname in ('SOURCE_MANIFEST.json','OUTPUT_MANIFEST.json'):
            d=json.loads(z.read(prefix+lname));m=d['files']
            for rel,digest in m.items():
                if sha_bytes(z.read(prefix+rel))!=digest:raise RuntimeError('PARENT_K_'+lname+'_FAIL_'+rel)
        report=json.loads(z.read(prefix+'results/report.json'))
        if report['scientific_decision']!=CFG['parent_k_report_decision'] or report['actual_Q550_ivps']!=10:
            raise RuntimeError('PARENT_K_PHYSICAL_AUDIT_MISSING')
        rows=list(csv.DictReader(io.StringIO(z.read(prefix+'results/fine_width_svd_geometry.csv').decode())))
        if len(rows)!=CFG['parent_k_expected_rows']:raise RuntimeError('PARENT_K_GRID_SIZE_WRONG')
        for name in names:
            path=Path(name)
            if path.is_absolute() or '..' in path.parts:raise RuntimeError('PARENT_K_ZIP_PATH_UNSAFE')
        z.extractall(temp_dir)
        root=Path(temp_dir)/next(iter(roots))
    return root,rows,report

def import_k(kroot):
    sys.path.insert(0,str(kroot/'scripts'))
    path=kroot/'scripts/run_audit.py'
    spec=importlib.util.spec_from_file_location('q550_v1853l_parent_k',path)
    mod=importlib.util.module_from_spec(spec);sys.modules[spec.name]=mod;spec.loader.exec_module(mod)
    return mod

def interp_source(xi,values,x):
    xi=np.asarray(xi);values=np.asarray(values)
    j=int(np.searchsorted(xi,x,side='right'))-1
    j=max(0,min(len(xi)-2,j))
    a=float(np.clip((x-xi[j])/(xi[j+1]-xi[j]),0.,1.))
    return (1-a)*values[j]+a*values[j+1]

def actual_endpoint_derivative(xi,B,tau_rate,origin,width):
    x=origin+3.*width
    return 3.*to10(interp_source(xi,B,x))*float(interp_source(xi,tau_rate,x))

def rows_for_group(seed,kn,origin,widths,acc,xi,B,tau_rate,C,parent_index):
    from three_window_geometry import interval
    anchor=float(CFG['anchor_width_xi']);wa=float(origin+3*anchor)
    M_anchor=[]
    for aa,bb in zip([origin,origin+anchor,origin+2*anchor],[origin+anchor,origin+2*anchor,origin+3*anchor]):
        U,V=interval(acc,aa,bb);M_anchor.append(to10(U+V))
    P0,Sa=svd_rank2(np.array(M_anchor))
    results=[];frames=[];max_ref=0.
    for w in widths:
        bounds=[origin,origin+w,origin+2*w,origin+3*w]
        parts=[]
        for aa,bb in zip(bounds[:-1],bounds[1:]):
            U,V=interval(acc,aa,bb);parts.append(to10(U+V))
        M=np.stack(parts);N=M.sum(axis=0)
        Pm,S=svd_rank2(M)
        Pd,L,T,maperr=dynamic_projector(P0,C,wa,float(bounds[-1]))
        # A second genuinely fixed-frame comparison, not a fixed re-expression
        # of the same per-width SVD plane.
        vf=P0@N;vm=Pm@N;vd=Pd@N;vpull=L.T@vd
        old=parent_index[(seed,kn,round(origin,6),round(w,6))]
        match=max(abs(float(old[f'net_{i:02d}'])-N[i]) for i in range(10))
        match=max(match,abs(float(old['rank2_net_norm'])-np.linalg.norm(vm)))
        if match>CFG['parent_net_absolute_tolerance']:raise RuntimeError('V1853K_REPLAY_DRIFT '+str((seed,kn,origin,w,match)))
        gap=float((S[1]-S[2])/max(S[1],1.e-14))
        idem=max(np.linalg.norm(Pm@Pm-Pm),np.linalg.norm(Pd@Pd-Pd),np.linalg.norm(P0@P0-P0))
        if idem>CFG['max_projector_idempotence_error']:raise RuntimeError('PROJECTOR_NOT_IDEMPOTENT')
        ortho=np.linalg.norm(T.T@T-np.eye(5))
        if ortho>CFG['max_model_rotation_orthogonality_error']:
            raise RuntimeError('MODEL_CONNECTION_NOT_ORTHOGONAL')
        if maperr>CFG['max_kinematic_mapping_error']:raise RuntimeError('MODEL_CONNECTION_U_MISMATCH')
        angle_fixed,_=projector_angle(P0,Pm)
        angle_model,_=projector_angle(Pd,Pm)
        r=dict(seed=seed,kappa_rule=kn,origin_xi=origin,width_xi=w,
            xi_end=bounds[-1],cancellation_score=old['cancellation_score'],
            rank2_gap_relative=gap,rank2_is_conditioned=bool(gap>=CFG['projector_gap_floor']),
            plane_angle_moving_to_fixed_deg=angle_fixed,
            plane_angle_moving_to_model_transport_deg=angle_model,
            fixed_norm=float(np.linalg.norm(vf)),moving_norm=float(np.linalg.norm(vm)),
            model_transport_norm=float(np.linalg.norm(vd)),
            fixed_complement_norm=float(np.linalg.norm(N-vf)),
            moving_complement_norm=float(np.linalg.norm(N-vm)),
            model_complement_norm=float(np.linalg.norm(N-vd)),
            parent_k_replay_abs_error=float(match),
            projector_idempotence_error=float(idem),model_rotation_orthogonality_error=float(ortho),
            model_unit_radial_mapping_error=float(maperr),
            source_ftc_norm=float(np.linalg.norm(actual_endpoint_derivative(xi,B,tau_rate,origin,w))))
        for prefix,value in (('net',N),('fixed_planar',vf),('moving_planar',vm),('model_planar',vd),('model_pullback',vpull)):
            for j,a in enumerate(value):r[f'{prefix}_{j:02d}']=float(a)
        results.append(r)
        frames.append(dict(N=N,Pfixed=P0,Pmoving=Pm,Pmodel=Pd,L=L,
                           vf=vf,vm=vm,vd=vd,vpull=vpull))
        max_ref=max(max_ref,match)
    diag=[];secants=[]
    arr=np.stack([x['N'] for x in frames]);xx=np.asarray(widths,float)
    for i,(r,f) in enumerate(zip(results,frames)):
        actual=actual_endpoint_derivative(xi,B,tau_rate,origin,xx[i])
        for count in (3,5):
            if i in (0,len(xx)-1) and count==3:continue
            fd=derivative_sample(xx,arr,i,count)
            err=float(np.linalg.norm(fd-actual)/max(CFG['boundary_error_scales_floor'],np.linalg.norm(fd),np.linalg.norm(actual)))
            d=dict(seed=seed,kappa_rule=kn,origin_xi=origin,width_xi=xx[i],
               stencil_nodes=count,nonuniform_corrected=True,
               scaled_derivative_error=err,derivative_norm=float(np.linalg.norm(fd)),
               ftc_endpoint_derivative_norm=float(np.linalg.norm(actual)),
               is_fine_band=bool(CFG['boundary_fine_band_xi'][0]<=xx[i]<=CFG['boundary_fine_band_xi'][1]))
            diag.append(d)
        if i==0:continue
        fa=frames[i-1];fb=f
        left,right=float(xx[i-1]),float(xx[i])
        exact_net=to10(sum(interval(acc,origin+3*left,origin+3*right)))
        integral_check=np.linalg.norm(fb['N']-fa['N']-exact_net)
        if integral_check>CFG['max_exact_secant_reconstruction_error']:
            raise RuntimeError('EXACT_BOUNDARY_INTEGRAL_SECANT_DRIFT')
        for typ, key in (('fixed','Pfixed'),('svd_moving','Pmoving'),('model_transported','Pmodel')):
            total,source,drift,ident=secant_product_decomposition(fa[key],fb[key],fa['N'],fb['N'])
            if ident>1e-10:raise RuntimeError('SECANT_PRODUCT_IDENTITY_FAILURE')
            secants.append(dict(seed=seed,kappa_rule=kn,origin_xi=origin,
                width_left_xi=left,width_right_xi=right,plane_choice=typ,
                delta_projected_norm=float(np.linalg.norm(total)),
                true_boundary_source_projected_norm=float(np.linalg.norm(source)),
                geometric_projector_drift_norm=float(np.linalg.norm(drift)),
                boundary_integral_reconciliation_error=float(integral_check),
                exact_product_rule_error=ident,
                drift_to_source_ratio=float(np.linalg.norm(drift)/max(np.linalg.norm(source),1e-15))))
    pivot={w:r for w,r in zip(widths,results)}
    a,b=[float(w) for w in (12.775,12.825)]
    fa,fb=frames[list(widths).index(a)],frames[list(widths).index(b)]
    orient=dict(seed=seed,kappa_rule=kn,origin_xi=origin,
       fixed_planar_cosine=cosine(fa['vf'],fb['vf']),
       moving_planar_cosine=cosine(fa['vm'],fb['vm']),
       model_planar_ambient_cosine=cosine(fa['vd'],fb['vd']),
       model_planar_pullback_cosine=cosine(fa['vpull'],fb['vpull']),
       fixed_min_planar_norm=float(min(np.linalg.norm(fa['vf']),np.linalg.norm(fb['vf']))),
       model_min_planar_norm=float(min(np.linalg.norm(fa['vd']),np.linalg.norm(fb['vd']))),
       moving_min_planar_norm=float(min(np.linalg.norm(fa['vm']),np.linalg.norm(fb['vm']))),
       pair_rank2_gap_min_relative=float(min(next(r['rank2_gap_relative'] for r in results if r['width_xi']==a),next(r['rank2_gap_relative'] for r in results if r['width_xi']==b))),
       pair_rank2_conditioned=bool(min(next(r['rank2_gap_relative'] for r in results if r['width_xi']==a),next(r['rank2_gap_relative'] for r in results if r['width_xi']==b))>=CFG['projector_gap_floor']),
       model_rotation_relative_angle_to_anchor_left_deg=next(r['plane_angle_moving_to_model_transport_deg'] for r in results if r['width_xi']==a),
       model_rotation_relative_angle_to_anchor_right_deg=next(r['plane_angle_moving_to_model_transport_deg'] for r in results if r['width_xi']==b))
    return results,diag,secants,orient,frames

def plots(rows):
    try:
        import matplotlib;matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:return
    for seed in ('original','autonomous'):
        group=sorted([r for r in rows if r['seed']==seed and r['kappa_rule']=='shared' and r['origin_xi']==14.4],key=lambda r:r['width_xi'])
        w=[r['width_xi'] for r in group]
        fig,ax=plt.subplots(figsize=(9,5.2),layout='constrained')
        for y,lab in [('fixed_norm','fixed 2-plane'),('moving_norm','per-width SVD'),('model_transport_norm','model-induced radial lift')]:
            ax.plot(w,[r[y] for r in group],marker='.',label=lab)
        ax.axvline(12.8,ls='--',color='gray')
        ax.grid(alpha=.25);ax.legend();ax.set_xlabel('Equal three-window width (xi)');ax.set_ylabel('Projected net norm (Frobenius)')
        ax.set_title('Q550 v18.53l | Fixed vs moving vs model-transported plane ('+seed+')')
        fig.savefig(OUT/f'plane_comparison_{seed}.png',dpi=140);plt.close(fig)
        fig,ax=plt.subplots(figsize=(9,5.2),layout='constrained')
        ax.plot(w,[r['plane_angle_moving_to_fixed_deg'] for r in group],marker='.',label='moving vs fixed')
        ax.plot(w,[r['plane_angle_moving_to_model_transport_deg'] for r in group],marker='.',label='moving vs model radial transport')
        ax.grid(alpha=.25);ax.legend();ax.set_xlabel('Equal three-window width (xi)');ax.set_ylabel('Maximum principal angle (deg)')
        fig.savefig(OUT/f'projector_drift_{seed}.png',dpi=140);plt.close(fig)

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--preflight',action='store_true')
    p.add_argument('--base',type=Path,default=ROOT.parent)
    p.add_argument('--parent-v1851s',type=Path)
    p.add_argument('--native-npz',type=Path)
    p.add_argument('--require-parent',action='store_true')
    p.add_argument('--require-native',action='store_true')
    p.add_argument('--require-dynamic',action='store_true')
    a=p.parse_args()
    report={'version':'v18.53l','scientific_decision':'NOT_EVALUATED',
       'new_actual_Q550_ivps':0,'new_refined_Q550_ivps':0,
       'physical_claim_of_unique_connection':False,
       'model_connection':'minimal radial SO(5) kinetic-frame lift induced on exterior-square (10D) bivectors',
       'continued_clock_beyond_xi25_6':True}
    exitcode=0
    try:
        with tempfile.TemporaryDirectory(prefix='q550_v1853l_stage_') as td:
            parentroot,parentrows,oldreport=verified_parent(td)
            report['parent_k_authentic']=True
            report['parent_k_pin_sha256']=CFG['parent_k_sha256']
            report['parent_k_recorded_ivps']=oldreport['actual_Q550_ivps']
            k=import_k(parentroot)
            j=k.load_j();k.verify_j_parent()
            # k's verifier covers pinned j and the independent j SVD review.
            hprobe=j.hist_module()
            if not all(callable(getattr(hprobe,name,None)) for name in ('loader','solve','validated_rhs','primary_metrics','native_53g_reconciliation')):
                raise RuntimeError('PARENT_J_H_DYNAMIC_FUNCTIONS_UNAVAILABLE')
            if a.preflight:
                if a.require_dynamic or a.require_parent or a.require_native:
                    raise RuntimeError('PREFLIGHT_CANNOT_CLAIM_PHYSICAL_EXECUTION')
                report['scientific_decision']='PREFLIGHT_PASS__NO_NEW_Q550_IVPS'
                return 0
            # Re-use the authenticated original source-staging methods.
            h=hprobe
            L=h.loader()
            par=L.find_parent(a.parent_v1851s)
            native=a.native_npz or L.find_npz('10d',a.base)
            if par and Path(par).is_file() and sha(par)==CFG['original_source_parent_sha256']:
                report['original_source_parent']='SHA256_AUTHENTICATED'
            else:report['original_source_parent']='MISSING_OR_MISMATCHED'
            if native and Path(native).is_file() and sha(native)==CFG['native_10d_sha256']:
                report['native_10d']='SHA256_AUTHENTICATED'
            else:report['native_10d']='MISSING_OR_MISMATCHED'
            if report['original_source_parent']!='SHA256_AUTHENTICATED' or report['native_10d']!='SHA256_AUTHENTICATED':
                report['scientific_decision']='AUTHENTIC_PARENT_K_ONLY__NEW_IVPS_NOT_EXECUTED'
                exitcode=3
                return exitcode
            parent_index={(r['seed'],r['kappa_rule'],round(float(r['origin_xi']),6),round(float(r['width_xi']),6)):r for r in parentrows}
            scan=k.F
            widths=[float(x) for x in scan['widths_xi']]
            origins=[float(x) for x in scan['origins_xi']]
            from three_window_geometry import accumulators
            from torque_math import embed,sigmoid
            contexttemp,engine,contexts,eng_meta=L.stage_original(par)
            try:
                allrows=[];boundary=[];secants=[];pairs=[];ivps=[];refs=[]
                maxref=0.
                with np.load(native,allow_pickle=False) as f:
                    xi,hcount=j.xi_grid(f['common_xi'])
                    for seed in ('original','autonomous'):
                        ctx=contexts[seed]
                        referr,clockerr,archived=h.native_53g_reconciliation(f,xi[:hcount],seed,ctx,engine.omega)
                        z0=np.asarray(archived[0],float)
                        for kn in ('shared','original','autonomous'):
                            kp=float(k.K[kn])
                            z,info=h.solve(engine,ctx,kp,xi,z0)
                            report['new_actual_Q550_ivps']+=1
                            source_err=h.validated_rhs(engine,ctx,kp,xi,z)
                            dd=h.primary_metrics(xi,z,ctx,kp)
                            if dd['metric']['charge_range']>k.CFG['three_window']['charge_abs_drift_tolerance']:
                                raise RuntimeError('CHARGE_DRIFT')
                            acc=accumulators(xi,dd['t'],dd['B_U'],dd['B_C'])
                            B=dd['B_U']+dd['B_C']
                            tau_rate=engine.omega/np.asarray(sigmoid(z[:,8]),float)
                            Q,V,_,_=embed(z,k.CFG['params']['m_x'])
                            C=transport_curve(xi,Q)
                            # Derived numeric radial connection is genuinely tied to Q(xi), not to SVD.
                            for origin in origins:
                                row,diag,sec,orient,_=rows_for_group(seed,kn,origin,widths,acc,xi,B,tau_rate,C,parent_index)
                                allrows.extend(row);boundary.extend(diag);secants.extend(sec);pairs.append(orient)
                                maxref=max(maxref,max(r['parent_k_replay_abs_error'] for r in row))
                            if kn in ('original','autonomous'):
                                tight_z,_=h.solve(engine,ctx,kp,xi,z0,tight=True)
                                report['new_actual_Q550_ivps']+=1;report['new_refined_Q550_ivps']+=1
                                scaled=float(np.max(np.abs(z-tight_z)/np.maximum(1.,np.abs(tight_z))))
                                refs.append({'seed':seed,'kappa_rule':kn,'trajectory_scaled_difference':scaled})
                                if scaled>CFG['refinement_scaled_tolerance']:
                                    raise RuntimeError('NEW_Q550_SOLVER_REFINEMENT_FAIL')
                            ivps.append({'seed':seed,'kappa_rule':kn,'kappa':kp,'nfev':info['nfev'],
                                'rhs_relative_error':source_err,'charge_drift':dd['metric']['charge_range'],
                                'parent53g_reconciliation_error':referr,'material_clock_error':clockerr,
                                'kinematic_radial_transport_endpoint_norm':float(np.linalg.norm(C['U'][-1]))})
                            print('V1853L_PHYSICAL_DONE',seed,kn,'ivps',report['new_actual_Q550_ivps'],flush=True)
                if report['new_actual_Q550_ivps']!=10 or len(allrows)!=756 or len(pairs)!=36:
                    raise RuntimeError('PHYSICAL_SAMPLE_COUNT_MISMATCH')
                if maxref>CFG['parent_net_absolute_tolerance']:raise RuntimeError('PARENT_REPLAY_MAX_DRIFT')
                write_csv('three_plane_comparison.csv',allrows)
                write_csv('nonuniform_endpoint_derivative.csv',boundary)
                write_csv('exact_secant_boundary_vs_projector_drift.csv',secants)
                write_csv('fixed_vs_moving_vs_model_opposition.csv',pairs)
                write_csv('genuine_integrated_ivps.csv',ivps)
                write_csv('refined_integration.csv',refs)
                fopp=[r for r in pairs if r['fixed_planar_cosine'] is not None and
                  r['fixed_planar_cosine']<CFG['fixed_opposition_cosine_threshold'] and
                  r['pair_rank2_conditioned'] and r['fixed_min_planar_norm']>=CFG['minimum_projection_norm']]
                mopp=[r for r in pairs if r['moving_planar_cosine'] is not None and
                  r['moving_planar_cosine']<CFG['fixed_opposition_cosine_threshold'] and r['pair_rank2_conditioned']]
                topp=[r for r in pairs if r['model_planar_pullback_cosine'] is not None and
                  r['model_planar_pullback_cosine']<CFG['fixed_opposition_cosine_threshold'] and
                  r['pair_rank2_conditioned'] and r['model_min_planar_norm']>=CFG['minimum_projection_norm']]
                deriv_summary={}
                for n in (3,5):
                    ss=[r for r in boundary if r['stencil_nodes']==n]
                    deriv_summary[str(n)]={'global_max_scaled_error':max(x['scaled_derivative_error'] for x in ss),
                       'fine_band_max_scaled_error':max(x['scaled_derivative_error'] for x in ss if x['is_fine_band']),
                       'global_median_scaled_error':float(np.median([x['scaled_derivative_error'] for x in ss]))}
                report['summary']={'n_rows':len(allrows),'n_groups':len(pairs),
                   'n_secants':len(secants),'n_boundary_stencils':len(boundary),
                   'fixed_opposition_groups':len(fopp),'moving_opposition_groups':len(mopp),
                   'model_pullback_opposition_groups':len(topp),
                   'max_k_csv_reconciliation_abs_error':maxref,
                   'max_exact_secant_identity_error':max(r['exact_product_rule_error'] for r in secants),
                   'max_ftc_integral_secant_error':max(r['boundary_integral_reconciliation_error'] for r in secants),
                   'max_model_u_mapping_error':max(r['model_unit_radial_mapping_error'] for r in allrows),
                   'boundary_derivative_errors':deriv_summary,
                   'max_projector_drift_over_source_ratio':max(r['drift_to_source_ratio'] for r in secants),
                   'maximum_moving_fixed_plane_angle_deg':max(r['plane_angle_moving_to_fixed_deg'] for r in allrows),
                   'maximum_moving_model_plane_angle_deg':max(r['plane_angle_moving_to_model_transport_deg'] for r in allrows)}
                report['scientific_decision']='AUTHENTIC_FIXED_MOVING_MODEL_PROJECTOR_SEPARATION__NO_UNIQUE_CONNECTION_CLAIM'
                report['scientific_limits']=[
                   'The original physical Q550 model is replayed; changes in widths are POSTPROCESSING of the same six traces.',
                   'The model-induced plane uses a DECLARED minimal radial kinetic-frame SO(5) transport; it is not a unique dynamical material connection or a demonstrated physical holonomy.',
                   'Fixed-vs-moving planar opposition in bivector-coordinate space does not establish a separate hidden degree of freedom.',
                   'Nonuniform polynomial stencil mismatch reflects truncation and/or interpolation; exact secant FTC attribution is more fundamental.',
                   'The clock at xi>25.6 is a continued historical model clock, not independently measured material recurrence.',
                   'The two-parameter SVD rank-2 plane is conditioned using sigma2-sigma3 gap; near degeneracy is reported, not smoothed into a spurious axis.']
                plots(allrows)
            finally:contexttemp.cleanup()
    except Exception as exc:
        report['scientific_decision']='AUDIT_FAILED_OR_INCOMPLETE__NO_CERTIFICATE'
        report['failure']=str(exc);report['traceback']=traceback.format_exc();exitcode=2
    finally:
        dump('report.json',report)
        (OUT/'SCIENTIFIC_REPORT_SERVER.md').write_text('# Q550 v18.53l — Fixed-Frame Opposition and Model-Transport Audit\n\n'
            'Decision: `'+report['scientific_decision']+'`\n\n'
            'Actual new Q550 integrations: '+str(report['new_actual_Q550_ivps'])+'\n\n'
            'Conclusions cannot be read as a proof of a unique material-frame connection, a physical flip, or temporal holonomy.\n\n'
            'Detailed summary:\n\n```json\n'+json.dumps(report.get('summary',{}),indent=2)+'\n```\n\n'
            'Failure, if present: '+str(report.get('failure','none'))+'\n')
        print('V1853L',report['scientific_decision'],'ACTUAL_IVPS',report['new_actual_Q550_ivps'],flush=True)
    if a.require_dynamic and report['scientific_decision'].startswith('AUTHENTIC_') is False:
        return max(3,exitcode)
    return exitcode

if __name__=='__main__':sys.exit(main())
