#!/usr/bin/env python3
"""Q550 v18.52d: perturbation-to-lateral-acceleration and finite-kick recovery.

All numerical work is derived from authenticated historical 10D nonlinear states.
The local preflight differentiates independently archived physical velocities in
material time tau. Server mode additionally evaluates the historical nonlinear
field (no surrogate dynamics), then evolves finite transverse velocity kicks
at fixed kappa0 and fixed drive phase.

A permanent change of kappa0 is NOT a transient kick; its trajectories cannot
by themselves certify relaxation after removal of a perturbation.
"""
from __future__ import annotations
import argparse,csv,hashlib,importlib.util,io,json,os,platform,sys,tempfile,zipfile
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
import scipy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'
VERSION='v18.52d'
SEEDS=('original','autonomous')
BRANCHES=('autonomous','midpoint','original')
QPOS=np.array([0,2,4,6,8],dtype=int)
UVEL=np.array([1,3,5,7,9],dtype=int)
SNAME='q550_parameter_sensitivity_transport_transverse_posture_attribution_v18_51s_bundle_server_sendback.zip'
SR='q550_parameter_sensitivity_transport_transverse_posture_attribution_v18_51s_bundle'
RR='q550_long_horizon_paired_kappa0_transverse_transport_return_section_search_v18_51r_bundle'
CR='q550_signed_margin_parameter_derivative_continuous_robustness_bound_v18_52c_bundle'


def digest(data):return hashlib.sha256(data).hexdigest()
def sha_file(path):
    with open(path,'rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def norm(v):return np.linalg.norm(v,axis=-1)

def import_file(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec)
    sys.modules[name]=m
    spec.loader.exec_module(m)
    return m


def source_manifest():return json.loads((ROOT/'SOURCE_MANIFEST.json').read_text())

def verify_sources(server=False):
    manifest=source_manifest()
    for name,h in manifest['file_sha256'].items():
        p=ROOT/name
        if not p.is_file() or sha_file(p)!=h:
            raise RuntimeError('Source file SHA256 mismatch: '+name)
    with zipfile.ZipFile(ROOT/'inputs'/('q550_signed_margin_parameter_derivative_continuous_robustness_bound_v18_52c_bundle_server_sendback.zip')) as z:
        if z.testzip() is not None:raise RuntimeError('v18.52c input CRC failure')
        cr=json.loads(z.read(CR+'/results/report.json'))
        if cr['version']!='v18.52c' or not cr['all_numerical_checks_passed'] or cr['new_physical_midpoints_completed']!=24:
            raise RuntimeError('v18.52c parent did not complete 24 real midpoint probes')
    ap=None
    if server:
        candidates=(ROOT/'inputs'/SNAME,ROOT.parent/SNAME,Path('/root/elab1/hysteresis')/SNAME)
        ap=next((x for x in candidates if x.is_file()),None)
        if ap is None:raise FileNotFoundError('Place original v18.51s server sendback in inputs/ or parent workspace')
        if sha_file(ap)!=manifest['parent_archives']['v18_51s_sha256']:
            raise RuntimeError('Authentic v18.51s SHA256 mismatch')
        with zipfile.ZipFile(ap) as z:
            if z.testzip() is not None:raise RuntimeError('v18.51s ZIP CRC failure')
            rr=z.read(SR+'/inputs/'+RR+'_server_sendback.zip')
            sr=z.read(SR+'/scripts/run_audit.py')
        if digest(rr)!=manifest['parent_archives']['v18_51r_sha256']:
            raise RuntimeError('Nested v18.51r SHA256 mismatch')
        with zipfile.ZipFile(io.BytesIO(rr)) as z:
            if z.testzip() is not None:raise RuntimeError('Nested v18.51r CRC failure')
            raw=z.read(RR+'/results/extended_paired_kappa_nonlinear_states.npz')
        if digest(raw)!=manifest['file_sha256']['inputs/historical_v18_51r_physical_states.npz']:
            raise RuntimeError('The cached r states differ from the exact nested v18.51s provenance')
        return ap,rr,sr
    return None,None,None


def load_states():
    p=ROOT/'inputs/historical_v18_51r_physical_states.npz'
    with np.load(p,allow_pickle=False) as arr:
        a={key:arr[key].copy() for key in arr.files}
    xi=a['xi']
    if len(xi)!=513 or abs(xi[-1]-.2)>1e-14:raise RuntimeError('Unexpected r physical span')
    for seed in SEEDS:
        for branch in BRANCHES:
            tag=seed+'_'+branch
            for k,shape in (('Y',(513,10)),('Q',(513,5)),('V',(513,5)),('tau',(513,)),('clock',(513,))):
                if a[tag+'_'+k].shape!=shape:raise RuntimeError('Unexpected genuine state '+tag+'_'+k)
            if not np.all(np.diff(a[tag+'_tau'])>0):raise RuntimeError('Non-monotone physical clock')
        assert np.max(abs(a[seed+'_original_Y'][0]-a[seed+'_autonomous_Y'][0]))<1e-14
    return a


def sample_at_common_tau(a,seed,kind0,kind1,source='kinematic',model_a=None,n_samples=513):
    """A same-SEED fixed-material-time paired comparison, NOT fixed-drive time."""
    key0=seed+'_'+kind0;key1=seed+'_'+kind1
    max_tau=min(float(a[key0+'_tau'][-1]),float(a[key1+'_tau'][-1]))
    tau=np.linspace(0,max_tau,n_samples)
    interp={}
    for idx,key in enumerate((key0,key1)):
        t=a[key+'_tau']
        for quantity in ('Q','V'):
            interp[(idx,quantity)]=CubicSpline(t,a[key+'_'+quantity],axis=0)(tau)
        if model_a is None:
            interp[(idx,'acc')]=CubicSpline(t,a[key+'_V'],axis=0).derivative()(tau)
        else:
            interp[(idx,'acc')]=CubicSpline(t,model_a[key],axis=0)(tau)
        interp[(idx,'xi')]=CubicSpline(t,a['xi'])(tau)
    vel0=interp[(0,'V')];vv=norm(vel0)
    if np.min(vv)<1e-12:raise RuntimeError('Cannot define material tangent at near-zero speed')
    e=vel0/vv[:,None]
    da=interp[(1,'acc')]-interp[(0,'acc')]
    dQ=interp[(1,'Q')]-interp[(0,'Q')]
    dV=interp[(1,'V')]-interp[(0,'V')]
    da_par=np.sum(da*e,axis=1)
    da_lat=da-da_par[:,None]*e
    # Do NOT equate delta(P_perp a) with P_perp delta(a).
    a0=interp[(0,'acc')];a1=interp[(1,'acc')]
    vel1=interp[(1,'V')]; e1=vel1/norm(vel1)[:,None]
    native0=a0-np.sum(a0*e,axis=1)[:,None]*e
    native1=a1-np.sum(a1*e1,axis=1)[:,None]*e1
    native_delta=native1-native0
    return {
      'tau':tau,'delta_a':da,'delta_a_norm':norm(da),'delta_a_parallel':da_par,
      'delta_a_lateral_norm':norm(da_lat),'delta_a_lateral_vector':da_lat,'baseline_velocity_unit':e,
      'delta_a_lateral_fraction':norm(da_lat)/np.maximum(norm(da),1e-25),
      'delta_Q_norm':norm(dQ),'delta_V_norm':norm(dV),
      'delta_speed':norm(vel1)-norm(vel0),
      'native_lateral_acceleration_difference_norm':norm(native_delta),
      'frame_rotation_correction_norm':norm(native_delta-da_lat),
      'external_drive_phase_difference_xi':interp[(1,'xi')]-interp[(0,'xi')],
    }


def spline_decimation_difference(a,seed):
    s=sample_at_common_tau(a,seed,'original','autonomous',n_samples=257)
    dt=s['tau']
    # Compare a derivative produced from half of the raw input knots, not merely a resampled curve.
    vv=[]
    for tag in (seed+'_original',seed+'_autonomous'):
        vv.append(CubicSpline(a[tag+'_tau'][::2],a[tag+'_V'][::2],axis=0).derivative()(dt))
    return float(np.max(norm((vv[1]-vv[0])-s['delta_a'])))


def summarize_pair(c,seed,mode):
    n=len(c['tau']);start=max(2,int(.05/.2*(n-1)))
    frac=c['delta_a_lateral_fraction'][start:]
    q=c['delta_Q_norm'];v=c['delta_V_norm']
    return dict(seed=seed,kind='fixed_kappa_difference',acceleration_source=mode,
      terminal_tau=float(c['tau'][-1]),
      max_abs_drive_clock_delta_at_matched_tau=float(np.max(abs(c['external_drive_phase_difference_xi']))),
      terminal_delta_acceleration=float(c['delta_a_norm'][-1]),
      terminal_lateral_acceleration=float(c['delta_a_lateral_norm'][-1]),
      terminal_parallel_acceleration=float(c['delta_a_parallel'][-1]),
      terminal_lateral_fraction=float(c['delta_a_lateral_fraction'][-1]),
      median_lateral_fraction_after_xi_approx_005=float(np.median(frac)),
      terminal_delta_position=float(q[-1]),terminal_delta_velocity=float(v[-1]),
      max_delta_position=float(q.max()),max_delta_velocity=float(v.max()),
      last_quarter_deltaQ_increasing_fraction=float(np.mean(np.diff(q[int(.75*n):])>0)),
      last_quarter_deltaV_increasing_fraction=float(np.mean(np.diff(v[int(.75*n):])>0)),
      recovery_status='NOT_APPLICABLE_TO_PERMANENT_KAPPA_SHIFT')


def model_acceleration_field(a,parts,kin):
    """Exact original 10D RHS -> physical acceleration dV/dtau.

    No numerical differentiation of a fitted trajectory is used for this branch.
    """
    _,_,seed_api,observer,engine,unit,drpt,creport,_=parts
    mx=float(engine.params.m_x)
    kvals={
      'original':float(a['kappa_original']),
      'autonomous':float(a['kappa_autonomous']),
    }
    kvals['midpoint']=.5*(kvals['original']+kvals['autonomous'])
    out={};model_residuals=[]
    for seed in SEEDS:
        ctx=seed_api.make_context(observer,seed,engine,unit,drpt,creport)
        for br in BRANCHES:
            tag=seed+'_'+br;acc=[];clockerr=[];embederr=[]
            for i,x in enumerate(a['xi']):
                Y=a[tag+'_Y'][i]
                tim=float(ctx['event']+x)
                z=np.r_[Y,float(ctx['theta'](tim))]
                f=np.asarray(engine.jet.structural_jet(z,kvals[br],engine.params,engine.slope)[0][:10],float)
                h,j=kin.embedding_and_jac(Y,mx,1.)
                c=float((Y[UVEL]@f[QPOS])/(Y[UVEL]@Y[UVEL]))
                if c<=0:raise RuntimeError('Nonpositive true material clock')
                acc.append((j[5:]@f)/c)
                clockerr.append(abs(c-a[tag+'_clock'][i]))
                embederr.append(max(np.max(abs(h[:5]-a[tag+'_Q'][i])),np.max(abs(h[5:]-a[tag+'_V'][i]))))
            out[tag]=np.asarray(acc)
            model_residuals.append(dict(seed=seed,branch=br,
                max_archived_clock_difference=float(max(clockerr)),
                max_archived_cartesian_embedding_difference=float(max(embederr))))
            print('Exact historical acceleration:',tag,flush=True)
    return out,model_residuals,(seed_api,observer,engine,unit,drpt,creport,mx,kvals)


def kick_vector(y,kin,mx,relative_amplitude,axis=3):
    h,j=kin.embedding_and_jac(y,mx,1.)
    v=h[5:];speed=float(np.linalg.norm(v))
    e=v/speed
    target=np.eye(5)[axis]
    transverse=target-(target@e)*e
    transverse/=np.linalg.norm(transverse)
    dv=float(relative_amplitude)*speed*transverse
    du=np.linalg.solve(j[5:,UVEL],dv)
    yy=np.array(y,copy=True);yy[UVEL]+=du
    actual,jj=kin.embedding_and_jac(yy,mx,1.)
    if np.linalg.norm(actual[:5]-h[:5])>1e-12:
        raise RuntimeError('Position drift during same-position velocity kick')
    if np.linalg.norm(actual[5:]-h[5:]-dv)>1e-11:
        raise RuntimeError('Velocity kick did not reproduce requested Cartesian delta')
    return yy,dv


def integrate_kick(ctx,engine,kappa,initial,xi,kin,mx,rtol,atol,max_step):
    from scipy.interpolate import CubicSpline
    start=float(ctx['event']);end=start+float(xi[-1])
    def rhs(t,y):
        z=np.r_[np.asarray(y,float),float(ctx['theta'](t))]
        return np.asarray(engine.jet.structural_jet(z,kappa,engine.params,engine.slope)[0][:10],float)
    sol=solve_ivp(rhs,(start,end),np.asarray(initial,float),method='DOP853',t_eval=start+xi,
      rtol=rtol,atol=atol,max_step=max_step)
    if not sol.success:raise RuntimeError('Historical nonlinear kick IVP failed: '+sol.message)
    ys=sol.y.T; qs=[];vs=[];cs=[]; aa=[]
    for x,y in zip(xi,ys):
        t=start+float(x);f=rhs(t,y)
        h,j=kin.embedding_and_jac(y,mx,1.)
        c=float((y[UVEL]@f[QPOS])/(y[UVEL]@y[UVEL]))
        if c<=0:raise RuntimeError('Nonpositive material clock in kick branch')
        qs.append(h[:5]);vs.append(h[5:]);cs.append(c)
        aa.append(j[5:]@f/c)
    cs=np.asarray(cs)
    return dict(Y=ys,Q=np.asarray(qs),V=np.asarray(vs),acc=np.asarray(aa),
                tau=np.asarray(CubicSpline(xi,cs).antiderivative()(xi)),clock=cs,
                nfev=int(sol.nfev))


def compare_kick_at_fixed_tau(base,pert,relative_amplitude):
    max_tau=min(base['tau'][-1],pert['tau'][-1]);tau=np.linspace(0,max_tau,513)
    vals=[]
    for r in (base,pert):
        vals.append({k:CubicSpline(r['tau'],r[k],axis=0)(tau) for k in ('Q','V','acc')})
    u,v=vals
    dq=v['Q']-u['Q'];dv=v['V']-u['V'];da=v['acc']-u['acc']
    tangent=u['V']/norm(u['V'])[:,None]
    along=np.sum(da*tangent,axis=1)
    lat=da-along[:,None]*tangent
    nv0=float(np.linalg.norm(u['V'][0]));length=nv0*max_tau
    phase=np.sqrt((norm(dq)/max(length,1e-12))**2+(norm(dv)/max(nv0,1e-12))**2)
    peak=float(max(phase));initial=float(phase[0]);terminal=float(phase[-1])
    # Diagnostic only. There is no theorem guaranteeing asymptotic recovery from this window.
    has_phase_return=bool(terminal<.5*initial and terminal<.5*peak)
    return dict(tau=tau,dq=norm(dq),dv=norm(dv),da=norm(da),da_lat=norm(lat),da_par=along,
         lateral_fraction=norm(lat)/np.maximum(norm(da),1e-25),phase=phase,
         relative_amplitude=relative_amplitude,
         phase_initial=initial,phase_peak=peak,phase_final=terminal,
         ratio_final_to_initial=terminal/max(initial,1e-30),
         ratio_final_to_peak=terminal/max(peak,1e-30),
         finite_window_phase_return_candidate=has_phase_return,
         rigorous_asymptotic_recovery_certified=False)


def write_csv(path,rows):
    if not rows:return
    with path.open('w',newline='',encoding='utf8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def plot_pairs(pairs,kicks):
    fig,axs=plt.subplots(2,2,figsize=(12,8),constrained_layout=True)
    for seed,color in (('original','tab:blue'),('autonomous','tab:orange')):
        p=pairs[seed]
        t=p['tau'];a=p['delta_a_norm'];ll=p['delta_a_lateral_norm']
        axs[0,0].plot(t,ll,label=seed+' lateral',color=color)
        axs[0,0].plot(t,abs(p['delta_a_parallel']),linestyle=':',color=color,label=seed+' |longitudinal|')
        axs[0,1].plot(t,p['delta_a_lateral_fraction'],label=seed,color=color)
        axs[1,0].plot(t,p['delta_Q_norm'],color=color,label=seed+' position')
        axs[1,0].plot(t,p['delta_V_norm'],color=color,linestyle='--',label=seed+' velocity')
    axs[0,0].set_title('Fixed physical time: acceleration change decomposition')
    axs[0,1].set_title('Lateral share of |delta acceleration|')
    axs[1,0].set_title('Permanent kappa change: separation is NOT recovery')
    axs[1,1].set_title('Finite-duration response of fixed-parameter velocity kicks')
    if kicks:
        for seed,color in (('original','tab:blue'),('autonomous','tab:orange')):
            for key,v in kicks.items():
                if not key.startswith(seed+'_'):continue
                ls='-' if '_positive' in key else ':'
                axs[1,1].plot(v['tau'],v['phase']/max(v['phase_initial'],1e-30),
                    color=color,ls=ls,alpha=.75,label=key)
    else:
        axs[1,1].text(.5,.5,'Server kick integrations pending',ha='center',va='center',transform=axs[1,1].transAxes)
    for ax in axs.flat:
        ax.set_xlabel('Elapsed material tau');ax.grid(alpha=.25)
        ax.legend(fontsize=7,loc='best') if ax.has_data() else None
    axs[0,0].set_ylabel('Acceleration units / tau');axs[0,1].set_ylabel('fraction')
    axs[1,0].set_ylabel('Cartesian norm');axs[1,1].set_ylabel('kick deviation / initial')
    fig.suptitle('Q550 v18.52d — Lateral acceleration vs genuine recovery',fontsize=13)
    fig.savefig(OUT/'lateral_acceleration_vs_recovery.png',dpi=145)
    plt.close(fig)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--mode',choices=('preflight','server'),default='server')
    ap.add_argument('--kick-relative-amplitudes',nargs='*',type=float,default=[1e-4,3e-4])
    ap.add_argument('--skip-kicks',action='store_true',help='Run only exact nonlinear RHS validation; no new trajectories')
    args=ap.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    parent,rr,srun=verify_sources(args.mode=='server')
    states=load_states();xi=states['xi']
    pre={}
    for seed in SEEDS:
        pre[seed]=sample_at_common_tau(states,seed,'original','autonomous')
    pre_summ={seed:summarize_pair(pre[seed],seed,'archived_velocity_spline_d_dtau') for seed in SEEDS}
    for seed in SEEDS:
        pre_summ[seed]['stride_two_vs_one_derivative_discrepancy_max']=spline_decimation_difference(states,seed)
    (OUT/'LOCAL_PREFLIGHT.json').write_text(json.dumps({'version':VERSION,'source':'v18.51r physical states',
        'prior_52c_midpoints_used_as_robustness_context_only':True,'summary':pre_summ,
        'recovery_claim_for_permanent_kappa_change':False},indent=2)+'\n')
    if args.mode=='preflight':
        plot_pairs(pre,{})
        print('Q550 v18.52d PRELIGHT: archived tau-kinematics only. Exact generator and finite kick pending.',flush=True)
        for seed,p in pre_summ.items():
            print(seed,p['terminal_lateral_fraction'],p['terminal_delta_acceleration'],flush=True)
        return
    kin=import_file(ROOT/'reference/kinetic_metric_v18_49h.py','q550_v1852d_kin')
    t=import_file(ROOT/'reference/parent_v18_51t_run_audit.py','q550_v1852d_t')
    with tempfile.TemporaryDirectory(prefix='q550_v1852d_auth_model_') as tmp:
        print('Recovering authenticated physical generator from 19-ancestor lineage',flush=True)
        parts=t.stage_engine(Path(tmp),rr,srun)[4]
        model,model_checks,ctxpack=model_acceleration_field(states,parts,kin)
        seed_api,observer,engine,unit,drpt,creport,mx,kappas=ctxpack
        physical={seed:sample_at_common_tau(states,seed,'original','autonomous',model_a=model) for seed in SEEDS}
        physical_summ={seed:summarize_pair(physical[seed],seed,'authentic_10D_rhs_to_dV_dtau') for seed in SEEDS}
        discrepancies={seed:float(np.max(norm(physical[seed]['delta_a']-pre[seed]['delta_a']))) for seed in SEEDS}
        relative_discrepancies={seed:float(np.max(
            norm(physical[seed]['delta_a'][26:]-pre[seed]['delta_a'][26:]) /
            np.maximum(1e-7,physical[seed]['delta_a_norm'][26:]))) for seed in SEEDS}
        fixed_xi={};
        for seed in SEEDS:
            v=states[seed+'_original_V'];e=v/norm(v)[:,None]
            da=model[seed+'_autonomous']-model[seed+'_original']
            pa=np.sum(da*e,axis=1);ll=norm(da-pa[:,None]*e);ss=norm(da)
            fixed_xi[seed]={'terminal_delta_acceleration':float(ss[-1]),
                'terminal_lateral_fraction':float(ll[-1]/max(ss[-1],1e-25))}
        # Only model-based acceleration is eligible for source attribution.
        kicks={};kick_rows=[];failures=[]; baseline_errors=[]
        if not args.skip_kicks:
            for seed in SEEDS:
                ctx=seed_api.make_context(observer,seed,engine,unit,drpt,creport)
                refkey=seed+'_midpoint'
                y0=states[refkey+'_Y'][0]
                reference=integrate_kick(ctx,engine,kappas['midpoint'],y0,xi,kin,mx,2e-10,2e-11,.00125)
                baseline_errors.append(dict(seed=seed,
                    max_archived_position_difference=float(np.max(norm(reference['Q']-states[refkey+'_Q']))),
                    max_archived_velocity_difference=float(np.max(norm(reference['V']-states[refkey+'_V']))),
                    nfev=int(reference['nfev'])))
                print('Historical fixed-kappa baseline checked:',seed,flush=True)
                for eps in args.kick_relative_amplitudes:
                    for sgn in (-1,+1):
                        lab=seed+'_'+('positive' if sgn>0 else 'negative')+'_'+format(eps,'.1e')
                        try:
                            changed,dv=kick_vector(y0,kin,mx,sgn*eps)
                            pert=integrate_kick(ctx,engine,kappas['midpoint'],changed,xi,kin,mx,2e-10,2e-11,.00125)
                            vv=compare_kick_at_fixed_tau(reference,pert,sgn*eps)
                            kicks[lab]=vv
                            rec=dict(seed=seed,label=lab,relative_kick=sgn*eps,initial_cartesian_velocity_kick_norm=float(norm(dv)),
                               initial_lateral_kick_parallel_component=float(np.dot(dv,reference['V'][0])/norm(reference['V'][0])),
                               **{k:vv[k] for k in ('phase_initial','phase_peak','phase_final',
                                  'ratio_final_to_initial','ratio_final_to_peak',
                                  'finite_window_phase_return_candidate','rigorous_asymptotic_recovery_certified')},
                               terminal_lateral_fraction=float(vv['lateral_fraction'][-1]),
                               nfev=int(pert['nfev']))
                            kick_rows.append(rec)
                            print('Physical kick',lab,'final/initial',rec['ratio_final_to_initial'],flush=True)
                        except Exception as exc:
                            failures.append({'seed':seed,'label':lab,'reason':str(exc)})
        checks={
         'parent_source_hash_verified':True,
         'all_model_acceleration_arrays_finite':all(np.isfinite(v).all() for v in model.values()),
         'cartesian_embedding_reproduces_archived_states':max(x['max_archived_cartesian_embedding_difference'] for x in model_checks)<2e-8,
         'physical_clock_reproduces_archive':max(x['max_archived_clock_difference'] for x in model_checks)<2e-8,
         'independent_velocity_derivative_approximately_agrees':max(discrepancies.values())<1e-6 and max(relative_discrepancies.values())<.03,
         'time_aligned_lateral_projection_orthogonality':all(
             np.max(abs(np.einsum('ij,ij->i',physical[s]['delta_a_lateral_vector'],
                        physical[s]['baseline_velocity_unit'])))<1e-11 for s in SEEDS),
         'kick_baselines_match_authentic_history':all(x['max_archived_position_difference']<1e-6 and x['max_archived_velocity_difference']<1e-6 for x in baseline_errors),
         'all_requested_kick_ivps_complete':args.skip_kicks or (not failures and len(kick_rows)==4*len(args.kick_relative_amplitudes)),
         'kick_initial_velocity_lateral':all(abs(x['initial_lateral_kick_parallel_component'])<1e-10 for x in kick_rows),
        }
        # Check exact pointwise projection, with the matched-time reference vector.
        for s in SEEDS:
            p=physical[s]
            checks[s+'_pointwise_parallel_plus_perpendicular_closure']=bool(
                np.max(abs(norm(p['delta_a'])**2-(p['delta_a_parallel']**2+norm(p['delta_a_lateral_vector'])**2)))<1e-12)
        classifications={s:('SAMPLED_MODEL_LATERAL_DOMINANT' if physical_summ[s]['median_lateral_fraction_after_xi_approx_005']>.5
                    else 'SAMPLED_MODEL_NOT_LATERAL_DOMINANT') for s in SEEDS}
        positive_kick_return=sum(x['finite_window_phase_return_candidate'] for x in kick_rows)
        report={
          'version':VERSION,'title':'Perturbation-to-Lateral-Acceleration Transfer & Recovery Audit',
          'mode':'authentic_10D_nonlinear_field_plus_transient_velocity_kicks',
          'all_numerical_checks_passed':bool(all(checks.values())),
          'checks':checks,'scientific_classification':classifications,
          'parent_v18_52c_midpoints_completed':24,
          'historical_fixed_kappa_pairing':physical_summ,'archived_kinematic_pairing':pre_summ,
          'model_acceleration_vs_velocity_spline_absolute_discrepancy':discrepancies,
          'model_acceleration_vs_velocity_spline_relative_discrepancy_qualified':relative_discrepancies,
          'same_external_phase_xi_comparison':fixed_xi,
          'model_rhs_immersion_residuals':model_checks,'kick_baseline_reproduction':baseline_errors,
          'transient_kicks_completed':len(kick_rows),'transient_kick_failures':failures,
          'finite_window_phase_return_candidates':int(positive_kick_return),
          'physical_recovery_certified':False,
          'permanent_parameter_shift_is_recovery_test':False,
          'limits':['Fixed tau is not fixed external drive phase; both clocks are recorded.',
              'Lateral component depends on the kinetic Cartesian chart and reference tangent.',
              'The kick changes one admissible state coordinate, without proving global charge-constraint preservation.',
              'Finite-window closeness or contraction is not asymptotic stability or dissipation.',
              'This is a finite sampled numerical experiment, not interval arithmetic.'],
          'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,
        }
        (OUT/'report.json').write_text(json.dumps(report,indent=2)+'\n')
        write_csv(OUT/'fixed_tau_lateral_acceleration_ledger.csv',[
            dict(seed=seed,index=i,tau=float(p['tau'][i]),
                 delta_acceleration=float(p['delta_a_norm'][i]),
                 lateral_acceleration=float(p['delta_a_lateral_norm'][i]),
                 signed_parallel_acceleration=float(p['delta_a_parallel'][i]),
                 lateral_share=float(p['delta_a_lateral_fraction'][i]),
                 delta_position=float(p['delta_Q_norm'][i]),
                 delta_velocity=float(p['delta_V_norm'][i]),
                 drive_phase_difference_xi=float(p['external_drive_phase_difference_xi'][i]),
                 native_lateral_difference=float(p['native_lateral_acceleration_difference_norm'][i]),
                 moving_tangent_correction=float(p['frame_rotation_correction_norm'][i]))
            for seed,p in physical.items() for i in range(len(p['tau']))])
        write_csv(OUT/'transverse_kick_recovery_ledger.csv',kick_rows)
        write_csv(OUT/'transient_kick_failures.csv',failures)
        if kicks:
            np.savez_compressed(OUT/'physical_kick_trajectories_and_diagnostics.npz',**{
              key+'__'+name:v for key,d in kicks.items() for name,v in d.items() if isinstance(v,np.ndarray)})
        plot_pairs(physical,kicks)
        (OUT/'SCIENTIFIC_REPORT_SERVER.md').write_text('\n'.join([
          '# Q550 v18.52d — Server scientific audit','',
          f'Numerical checks passed: **{all(checks.values())}**.',
          f'Physical kicks completed: **{len(kick_rows)}**; unresolved: **{len(failures)}**.',
          '', 'The fixed-parameter pair is **not** a transient recovery experiment.',
          'Finite kicks are performed at fixed physical model parameters and common initial drive phase.',
          'All recovered trajectories are compared again at the same elapsed material time.',
          'No finite-horizon test alone certifies asymptotic material recovery.',
          '', '## Pairwise acceleration diagnostics','',
          *[f'- {seed}: final lateral fraction {physical_summ[seed]["terminal_lateral_fraction"]:.9f}, '
             f'final |delta a| {physical_summ[seed]["terminal_delta_acceleration"]:.9g}, '
             f'generator-vs-spline maximum {discrepancies[seed]:.4g}.' for seed in SEEDS],
          '', '## Conclusions','',
          'See `report.json` and both CSV ledgers. Physical interpretation is limited by the stated metrics and horizons.',
          '']))
        print('Q550 v18.52d numerical checks',sum(checks.values()),'/',len(checks),flush=True)
        print('Q550 v18.52d scientific classification',classifications,flush=True)
        if not all(checks.values()):raise SystemExit('NUMERICAL CHECK FAILURE — inspect report.json')


if __name__=='__main__':main()
