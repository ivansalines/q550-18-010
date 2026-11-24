#!/usr/bin/env python3
"""Q550 v18.52e: signed relative acceleration work and physical turning.

Preflight reuses authenticated prior kick norms; server reproduces historical
10D nonlinear trajectories and measures signed vectors at fixed material time.
All expensive authentic integrations are reserved for the user's server.
"""
from __future__ import annotations
import argparse,csv,hashlib,importlib.util,io,json,os,platform,sys,tempfile,zipfile
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicSpline
import scipy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'
DR='q550_perturbation_to_lateral_acceleration_transfer_recovery_v18_52d_bundle'
DZIP=ROOT/'inputs'/(DR+'_server_sendback.zip')
SNAME='q550_parameter_sensitivity_transport_transverse_posture_attribution_v18_51s_bundle_server_sendback.zip'
SR='q550_parameter_sensitivity_transport_transverse_posture_attribution_v18_51s_bundle'
RR='q550_long_horizon_paired_kappa0_transverse_transport_return_section_search_v18_51r_bundle'
SEEDS=('original','autonomous')
UVEL=np.array([1,3,5,7,9],int)
QPOS=np.array([0,2,4,6,8],int)

def sha(blob):return hashlib.sha256(blob).hexdigest()
def file_hash(path):
    with open(path,'rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def mag(x):return np.linalg.norm(x,axis=-1)
def imp(path,alias):
    spec=importlib.util.spec_from_file_location(alias,path)
    m=importlib.util.module_from_spec(spec);sys.modules[alias]=m;spec.loader.exec_module(m)
    return m

def verify_sources():
    mf=json.loads((ROOT/'SOURCE_MANIFEST.json').read_text())
    for p,h in mf['file_sha256'].items():
        loc=ROOT/p
        if not loc.is_file() or file_hash(loc)!=h:raise RuntimeError('Source digest mismatch: '+p)
    with zipfile.ZipFile(DZIP) as z:
        if z.testzip() is not None:raise RuntimeError('Parent d ZIP CRC failure')
        manifest=json.loads(z.read(DR+'/SHA256SUMS_OUTPUT.json'))
        for f,h in manifest['file_sha256'].items():
            if sha(z.read(DR+'/'+f))!=h:raise RuntimeError('Parent d output hash failure: '+f)
        report=json.loads(z.read(DR+'/results/report.json'))
        if report['version']!='v18.52d' or not report['all_numerical_checks_passed'] or report['transient_kicks_completed']!=8:
            raise RuntimeError('Authentic d results incomplete')
        for orig,ours in [
            ('scripts/run_audit.py','reference/parent_v18_52d_run_audit.py'),
            ('reference/parent_v18_51t_run_audit.py','reference/parent_v18_51t_run_audit.py'),
            ('reference/kinetic_metric_v18_49h.py','reference/kinetic_metric_v18_49h.py')]:
            if z.read(DR+'/'+orig)!=(ROOT/ours).read_bytes():raise RuntimeError('Ancestor source not byte-identical: '+ours)
        state=z.read(DR+'/inputs/historical_v18_51r_physical_states.npz')
        origkick=z.read(DR+'/results/physical_kick_trajectories_and_diagnostics.npz')
    return report,state,origkick,mf

def prior_preflight(origkick):
    traces=[];plot=[]
    with np.load(io.BytesIO(origkick),allow_pickle=False) as a:
        labels=sorted(k.removesuffix('__tau') for k in a.files if k.endswith('__tau'))
        if len(labels)!=8:raise RuntimeError('Unexpected number of historical physical kicks')
        for lab in labels:
            tau=np.asarray(a[lab+'__tau']);dv=np.asarray(a[lab+'__dv']);dq=np.asarray(a[lab+'__dq'])
            kinetic=.5*dv**2
            p=CubicSpline(tau,kinetic).derivative()(tau)
            g=CubicSpline(tau,.5*dq**2).derivative()(tau)
            qualified=slice(2,-2)
            r=dict(label=lab,velocity_ratio_terminal_initial=float(dv[-1]/dv[0]),
                position_distance_terminal=float(dq[-1]),
                negative_relative_power_fraction=float(np.mean(p[qualified]<0)),
                positive_position_separation_rate_fraction=float(np.mean(g[qualified]>0)),
                preflight_only=True)
            traces.append(r);plot.append((lab,tau,dv/dv[0],dq,p))
    return traces,plot

def get_s_archive(mf):
    candidates=[ROOT/'inputs'/SNAME,ROOT.parent/SNAME,Path('/root/elab1/hysteresis')/SNAME]
    found=next((p for p in candidates if p.is_file()),None)
    if found is None:raise FileNotFoundError('Place exact v18.51s sendback at workspace root or inputs/')
    if file_hash(found)!=mf['parent_archives']['v18_51s_sha256']:
        raise RuntimeError('v18.51s sendback SHA-256 mismatch')
    with zipfile.ZipFile(found) as z:
        if z.testzip() is not None:raise RuntimeError('v18.51s CRC mismatch')
        nested=z.read(SR+'/inputs/'+RR+'_server_sendback.zip')
        sscript=z.read(SR+'/scripts/run_audit.py')
    return nested,sscript

def base_history(states,seed,engine,kin,context,mx,kappa):
    """Read archival genuine trajectory and independently compute acceleration."""
    pre=seed+'_midpoint';xi=states['xi']
    Y=states[pre+'_Y'];Q=states[pre+'_Q'];V=states[pre+'_V']
    clock=states[pre+'_clock'];tau=states[pre+'_tau']; aa=[];clockdef=[];embdef=[]
    for i,x in enumerate(xi):
        yy=Y[i]; t=float(context['event']+x)
        z=np.r_[yy,float(context['theta'](t))]
        f=np.asarray(engine.jet.structural_jet(z,kappa,engine.params,engine.slope)[0][:10],float)
        h,j=kin.embedding_and_jac(yy,mx,1.)
        c=float((yy[UVEL]@f[QPOS])/(yy[UVEL]@yy[UVEL]))
        aa.append(j[5:]@f/c)
        clockdef.append(abs(c-clock[i]));embdef.append(float(np.max(abs(h-np.c_[Q,V][i]))))
    return dict(Y=Y,Q=Q,V=V,acc=np.asarray(aa),tau=tau,clock=clock,xi=xi,
        max_clock_defect=float(max(clockdef)),max_embedding_defect=float(max(embdef)))

def full_geometry(base,pert,label,eps,parent_norms):
    end=min(float(base['tau'][-1]),float(pert['tau'][-1]))
    tau=np.linspace(0,end,513)
    u={};w={}
    for values,source in ((u,base),(w,pert)):
        for key in ('Q','V','acc'):
            values[key]=CubicSpline(source['tau'],source[key],axis=0)(tau)
        values['xi']=CubicSpline(source['tau'],source['xi'])(tau)
    dq=w['Q']-u['Q'];dv=w['V']-u['V'];da=w['acc']-u['acc']
    qn=mag(dq);vn=mag(dv);an=mag(da)
    e=u['V']/mag(u['V'])[:,None]
    par=np.sum(da*e,axis=1)
    da_lat=da-par[:,None]*e
    p=np.einsum('ij,ij->i',dv,da)
    g=np.einsum('ij,ij->i',dq,dv)
    inward=np.einsum('ij,ij->i',dq,da)
    second=vn**2+inward
    E=.5*vn**2;D=.5*qn**2
    # Two independent numerical derivative routes: derivative of a cubic
    # spline of the norm, and derivative of the full vector components.
    p_spline=CubicSpline(tau,E).derivative()(tau)
    g_spline=CubicSpline(tau,D).derivative()(tau)
    dv_spline=CubicSpline(tau,dv,axis=0).derivative()(tau)
    dq_spline=CubicSpline(tau,dq,axis=0).derivative()(tau)
    p_from_v_spline=np.einsum('ij,ij->i',dv,dv_spline)
    g_from_q_spline=np.einsum('ij,ij->i',dq,dq_spline)
    cp=[];speeds=[];turns=[]
    for branch in (u,w):
        vel=branch['V'];acc=branch['acc'];speed=mag(vel)
        tangent=vel/speed[:,None]
        aperp=acc-np.einsum('ij,ij->i',acc,tangent)[:,None]*tangent
        cp.append(mag(aperp)/speed**2)
        turns.append(mag(aperp)/speed)
        speeds.append(speed)
    curvature_delta=cp[1]-cp[0]
    # Compare with norm-only parent archive; this is a sanity check, not an
    # independent dynamics audit (the new IVPs are independently rerun).
    nerr={}
    for name,value in [('tau',tau),('dq',qn),('dv',vn),('da',an),('da_lat',mag(da_lat)),('da_par',par)]:
        key=label+'__'+name
        if key not in parent_norms:raise RuntimeError('Parent missing trace '+key)
        expected=np.asarray(parent_norms[key]);nerr[name]=float(np.max(abs(value-expected)))
    qmask=qn>max(1e-10,1e-6*float(qn[-1]));
    valid=np.arange(len(tau));mid=(valid>=2)&(valid<len(tau)-2)
    c0=float(vn[-1]/vn[0]); pneg=float(np.mean(p[mid]<0));gpos=float(np.mean(g[mid]>0))
    inward_fraction=float(np.mean(inward[qmask]<0)) if qmask.any() else None
    slope=np.diff(qn);negative_windows=int(np.sum(slope<0))
    status=('RELATIVE_VELOCITY_CONTRACTION_WITH_POSITION_DIVERGENCE'
       if c0<1 and pneg>=.99 and gpos>=.99 else 'MIXED_OR_NO_VELOCITY_CONTRACTION')
    frac=mag(da_lat)/np.maximum(an,1e-25)
    summary=dict(label=label,seed=label.split('_')[0],relative_kick=float(eps),
       elapsed_material_tau=float(tau[-1]),
       velocity_difference_terminal_over_initial=c0,
       relative_speed_squared_terminal_over_initial=float(E[-1]/E[0]),
       position_difference_final=float(qn[-1]),
       position_distance_increasing_sample_fraction=float(np.mean(np.diff(qn)>0)),
       negative_signed_relative_power_fraction=pneg,
       positive_signed_position_separation_fraction=gpos,
       negative_position_acceleration_dot_fraction_qualified=inward_fraction,
       position_restoration_sample_count=int(np.sum(g[mid]<0)),
       qualified_inward_acceleration_sample_count=int(qmask.sum()),
       position_distance_decrease_steps=negative_windows,
       integrated_signed_relative_power=float(np.trapezoid(p,tau)),
       exact_terminal_relative_kinetic_change=float(E[-1]-E[0]),
       power_integral_closure_abs=float(abs(np.trapezoid(p,tau)-(E[-1]-E[0]))),
       max_power_vs_energy_spline_abs=float(np.max(abs(p[2:-2]-p_spline[2:-2]))),
       max_power_vs_vector_spline_abs=float(np.max(abs(p[2:-2]-p_from_v_spline[2:-2]))),
       max_sep_rate_vs_distance_spline_abs=float(np.max(abs(g[2:-2]-g_spline[2:-2]))),
       max_sep_rate_vs_vector_spline_abs=float(np.max(abs(g[2:-2]-g_from_q_spline[2:-2]))),
       terminal_curvature_reference=float(cp[0][-1]),
       terminal_curvature_kick=float(cp[1][-1]),
       terminal_curvature_change=float(curvature_delta[-1]),
       curvature_sign_change_count=int(np.sum(np.sign(curvature_delta[2:-1])*np.sign(curvature_delta[1:-2])<0)),
       terminal_delta_acceleration_lateral_fraction=float(frac[-1]),
       max_abs_external_phase_difference_at_fixed_tau=float(np.max(abs(w['xi']-u['xi']))),
       max_abs_norm_trace_replay_difference=max(nerr.values()),
       classification=status,physical_position_return_certified=False,
       physical_asymptotic_stability_certified=False)
    vectors=dict(tau=tau,dQ=dq,dV=dv,da=da,
       distance_Q=qn,distance_V=vn,relative_power=p,
       relative_position_rate=g,position_acceleration_dot=inward,
       second_position_separation_derivative=second,
       curvature_reference=cp[0],curvature_kick=cp[1],
       curvature_delta=curvature_delta,lateral_acceleration_fraction=frac,
       delta_external_phase_xi=w['xi']-u['xi'],
       relative_energy=E,relative_distance_squared=D,
       power_spline=p_spline,sep_rate_spline=g_spline)
    return summary,vectors

def save_rows(path,rows):
    with path.open('w',newline='',encoding='utf8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def plot_preflight(traces):
    fig,axs=plt.subplots(1,2,figsize=(11,4.1),constrained_layout=True)
    for name,tau,relative,dq,p in traces:
        axs[0].plot(tau,relative,label=name,lw=1.3)
        axs[1].plot(tau,p/abs(p[2]) if abs(p[2])>1e-20 else p,lw=1.3)
    axs[0].axhline(1,color='k',lw=.7,ls=':');axs[0].set_title('Relative velocity magnitude / initial')
    axs[1].axhline(0,color='k',lw=.7,ls=':');axs[1].set_title('Signed relative power: archive-norm derivative')
    for a in axs:a.grid(alpha=.2);a.set_xlabel('Material tau')
    axs[0].set_ylabel('Ratio');axs[1].set_ylabel('Scaled derivative (diagnostic only)')
    fig.savefig(OUT/'PRELIMINARY_parent_52d_velocity_contracts.png',dpi=140)
    plt.close(fig)

def plot_server(results):
    fig,axs=plt.subplots(2,2,figsize=(12,8),constrained_layout=True)
    for label,v in results.items():
        t=v['tau'];s=1e4 if '1.0e-04' in label else 1e4/3
        axs[0,0].plot(t,v['relative_power']*s*s,label=label,lw=1.15)
        axs[0,1].plot(t,v['relative_position_rate']*s*s,lw=1.15)
        axs[1,0].plot(t,v['distance_V']/v['distance_V'][0],lw=1.15)
        axs[1,1].plot(t,v['distance_Q']*s,lw=1.15)
    axs[0,0].axhline(0,color='k',ls=':',lw=.7);axs[0,0].set_title('Signed relative acceleration work  ΔV · Δa')
    axs[0,1].axhline(0,color='k',ls=':',lw=.7);axs[0,1].set_title('Position separation rate  ΔQ · ΔV')
    axs[1,0].axhline(1,color='k',ls=':',lw=.7);axs[1,0].set_title('Relative speed contraction')
    axs[1,1].set_title('Position separation accumulation')
    axs[0,0].set_ylabel('Scaled power');axs[0,1].set_ylabel('Scaled separation rate')
    axs[1,0].set_ylabel('|ΔV| / |ΔV(0)|');axs[1,1].set_ylabel('Scaled |ΔQ|')
    for ax in axs.ravel():ax.set_xlabel('Common material tau');ax.grid(alpha=.2)
    axs[0,0].legend(fontsize=6,ncol=2,loc='best')
    fig.savefig(OUT/'relative_work_vs_position_separation.png',dpi=140);plt.close(fig)
    fig,axs=plt.subplots(1,3,figsize=(14,4.4),constrained_layout=True)
    for label,v in results.items():
        t=v['tau'];axs[0].plot(t,v['curvature_delta'],lw=1.2,label=label)
        axs[1].plot(t,v['position_acceleration_dot'],lw=1.2)
        axs[2].plot(t,v['lateral_acceleration_fraction'],lw=1.2)
    axs[0].axhline(0,color='k',ls=':',lw=.7);axs[1].axhline(0,color='k',ls=':',lw=.7)
    axs[0].set_title('Kick minus baseline curvature, 5D Cartesian')
    axs[1].set_title('Position-oriented acceleration  ΔQ · Δa')
    axs[2].set_title('Lateral fraction of Δa')
    for ax in axs:ax.grid(alpha=.2);ax.set_xlabel('Common material tau')
    axs[0].legend(fontsize=6,ncol=1,loc='best')
    fig.savefig(OUT/'lateral_curvature_and_restoring_acceleration.png',dpi=140);plt.close(fig)

def main():
    arg=argparse.ArgumentParser(description=__doc__)
    arg.add_argument('--mode',choices=('preflight','server'),default='server')
    args=arg.parse_args();OUT.mkdir(exist_ok=True,parents=True)
    dreport,state_blob,kickblob,mf=verify_sources()
    pre,plot=prior_preflight(kickblob)
    (OUT/'LOCAL_PREFLIGHT.json').write_text(json.dumps({
        'version':'v18.52e','stage':'ARCHIVED_NORM_DERIVATIVE_ONLY',
        'source':'authenticated v18.52d genuine transient kick traces',
        'results':pre,'direct_vector_signed_work_verified':False,
        'position_return_proven':False},indent=2)+'\n')
    plot_preflight(plot)
    print('Q550 v18.52e authenticated preflight: 8 kick norm traces verified',flush=True)
    if args.mode=='preflight':
        print('Signed vector work, curvature and restoring acceleration: SERVER PENDING',flush=True)
        return
    nested,srun=get_s_archive(mf)
    parent=imp(ROOT/'reference/parent_v18_52d_run_audit.py','v1852e_parent_d')
    parentt=imp(ROOT/'reference/parent_v18_51t_run_audit.py','v1852e_parent_t')
    kin=imp(ROOT/'reference/kinetic_metric_v18_49h.py','v1852e_kin')
    with np.load(io.BytesIO(state_blob),allow_pickle=False) as tmp:
        states={key:tmp[key] for key in tmp.files}
    xi=states['xi']
    with np.load(io.BytesIO(kickblob),allow_pickle=False) as tmp:
        historical={key:tmp[key] for key in tmp.files}
    with tempfile.TemporaryDirectory(prefix='q550_v1852e_historical_engine_') as tmp:
        print('Recover authenticated original nonlinear model from s/r lineage',flush=True)
        parts=parentt.stage_engine(Path(tmp),nested,srun)[4]
        _,_,seed_api,observer,engine,unit,drpt,creport,_=parts
        mx=float(engine.params.m_x)
        km=.5*(float(states['kappa_original'])+float(states['kappa_autonomous']))
        allsumm=[];allvec={};baselines=[]
        for seed in SEEDS:
            ctx=seed_api.make_context(observer,seed,engine,unit,drpt,creport)
            base=base_history(states,seed,engine,kin,ctx,mx,km)
            baselines.append({'seed':seed,'max_clock_defect':base['max_clock_defect'],
               'max_embedding_defect':base['max_embedding_defect']})
            print('Authentic midpoint baseline acceleration:',seed,flush=True)
            y0=states[seed+'_midpoint_Y'][0]
            for eps in (1e-4,3e-4):
                for sign in (-1,+1):
                    lab=f'{seed}_{"positive" if sign>0 else "negative"}_{eps:.1e}'
                    changed,dv0=parent.kick_vector(y0,kin,mx,sign*eps)
                    kick=parent.integrate_kick(ctx,engine,km,changed,xi,kin,mx,
                                              rtol=2e-10,atol=2e-11,max_step=.00125)
                    kick['xi']=xi
                    summ,vec=full_geometry(base,kick,lab,sign*eps,historical)
                    summ['nfev']=int(kick['nfev']);summ['initial_kick_lateral_residual']=float(abs(np.dot(dv0,base['V'][0])))
                    allsumm.append(summ);allvec[lab]=vec
                    print('Kick',lab,'P negative fraction',f'{summ["negative_signed_relative_power_fraction"]:.6f}',
                        'position rate positive fraction',f'{summ["positive_signed_position_separation_fraction"]:.6f}',flush=True)
        checks={
          'parent_v18_52d_crc_source_and_output_hashes':True,
          'authentic_v18_51s_sha256':True,
          'all_eight_exact_nonlinear_kick_integrations':len(allsumm)==8,
          'all_physical_clock_and_embedding_residuals_small':all(
            x['max_clock_defect']<2e-8 and x['max_embedding_defect']<2e-8 for x in baselines),
          'all_signed_kick_vectors_finite':all(all(np.isfinite(v).all() for v in group.values()) for group in allvec.values()),
          'initial_kick_lateral':all(x['initial_kick_lateral_residual']<1e-10 for x in allsumm),
          'full_state_kick_traces_reproduce_parent_norms':all(x['max_abs_norm_trace_replay_difference']<2e-7 for x in allsumm),
          'relative_work_power_vs_energy_cubic_spline':all(x['max_power_vs_energy_spline_abs']<1e-9 for x in allsumm),
          'relative_work_power_vs_vector_spline':all(x['max_power_vs_vector_spline_abs']<1e-9 for x in allsumm),
          'position_separation_rate_vs_distance_spline':all(x['max_sep_rate_vs_distance_spline_abs']<1e-9 for x in allsumm),
          'position_separation_rate_vs_vector_spline':all(x['max_sep_rate_vs_vector_spline_abs']<1e-9 for x in allsumm),
          'relative_power_integral_closes':all(x['power_integral_closure_abs']<1e-9 for x in allsumm),
          'positive_finite_speed_curvature':all(all(np.isfinite(group[k]).all() for k in ('curvature_reference','curvature_kick')) for group in allvec.values()),
        }
        report={
          'version':'v18.52e','title':'Lateral Curvature, Relative-Acceleration Work & Restoring-Tendency Audit',
          'stage':'SERVER_AUTHENTIC_NONLINEAR_SIGNED_VECTOR_DIAGNOSTIC',
          'all_numerical_checks_passed':all(checks.values()),'checks':checks,
          'physical_kicks':allsumm,'baseline_reconstruction':baselines,
          'physical_recovery_certified':False,
          'asymptotic_stability_or_dissipation_certified':False,
          'qualified_interpretation':'Negative dV·da is contraction of squared relative velocity, NOT physical energy dissipation. Positive dQ·dV means position distance still increasing. Signed deltaQ·deltaa describes inward acceleration projection only.',
          'limitations':['finite xi <= 0.2 and fixed-kappa physical kicks',
              'common tau can differ in external drive phase',
              'Cartesian kinetic coordinate chart and its Euclidean norm are specified',
              'impulses are velocity-state perturbations; global charge invariance not proven',
              'nonlinear fixed-time sample diagnostics, not interval bounds or asymptotic stability'],
          'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,
        }
        (OUT/'report.json').write_text(json.dumps(report,indent=2)+'\n')
        save_rows(OUT/'kick_classifications.csv',allsumm)
        records=[];savearrays={}
        for lab,v in allvec.items():
            for k,arr in v.items():savearrays[lab+'__'+k]=arr
            for i in range(len(v['tau'])):
                records.append(dict(label=lab,index=i,tau=float(v['tau'][i]),
                    relative_work_power=float(v['relative_power'][i]),
                    relative_position_rate=float(v['relative_position_rate'][i]),
                    position_acceleration_dot=float(v['position_acceleration_dot'][i]),
                    delta_Q_norm=float(v['distance_Q'][i]),delta_V_norm=float(v['distance_V'][i]),
                    relative_curvature_difference=float(v['curvature_delta'][i]),
                    acceleration_lateral_fraction=float(v['lateral_acceleration_fraction'][i]),
                    delta_external_phase_xi=float(v['delta_external_phase_xi'][i])))
        save_rows(OUT/'relative_work_curvature_ledger.csv',records)
        np.savez_compressed(OUT/'physical_relative_vector_diagnostics.npz',**savearrays)
        plot_server(allvec)
        (OUT/'SCIENTIFIC_REPORT_SERVER.md').write_text('\n'.join([
            '# Q550 v18.52e — Server numerical results','',
            'Numerical checks passed: **'+str(sum(checks.values()))+'/'+str(len(checks))+'**.','',
            'Eight genuine 10D nonlinear kicks were integrated at fixed kappa and independently compared with the archived v18.52d norm traces.','',
            '| Seed and kick | Terminal relative speed ratio | P_rel < 0 fraction | d(|delta Q|²)/dtau > 0 fraction | Qualified inward acceleration fraction |',
            '|---|---:|---:|---:|---:|',
            *[f'| {r["label"]} | {r["velocity_difference_terminal_over_initial"]:.9f} | {r["negative_signed_relative_power_fraction"]:.6f} | {r["positive_signed_position_separation_fraction"]:.6f} | {r["negative_position_acceleration_dot_fraction_qualified"]:.6f} |' for r in allsumm],
            '', '**Do not infer recovery solely from negative relative work.**',
            'The power is a Euclidean relative-kinematics proxy, not physical dissipated energy; the full position separation and clock must also be examined.',
            '', 'Scientific classification and numerical tolerances are in `report.json`.'
        ])+'\n')
        print('Q550 v18.52e checks',sum(checks.values()),'/',len(checks),flush=True)
        if not all(checks.values()):raise SystemExit('NUMERICAL VALIDATION FAILURE — inspect results/report.json')

if __name__=='__main__':main()
