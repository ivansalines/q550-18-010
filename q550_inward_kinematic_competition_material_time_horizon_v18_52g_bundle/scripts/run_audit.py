#!/usr/bin/env python3
"""Q550 v18.52g: authentic nonlinear horizon extension and event competition.

No surrogate force, extrapolated spline, or altered kappa. All event roots are
located on matched-material-time dense DOP853 trajectories, subject to finite
precision and the event-search grid. Historical lineage is SHA-verified.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import importlib.util
import io
import json
import platform
import sys
import tempfile
import zipfile
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import scipy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results'
FROOT = 'q550_relative_separation_convexity_inward_acceleration_sufficiency_v18_52f_bundle'
EROOT = 'q550_lateral_curvature_relative_acceleration_work_restoring_tendency_v18_52e_bundle'
DROOT = 'q550_perturbation_to_lateral_acceleration_transfer_recovery_v18_52d_bundle'
SROOT = 'q550_parameter_sensitivity_transport_transverse_posture_attribution_v18_51s_bundle'
RROOT = 'q550_long_horizon_paired_kappa0_transverse_transport_return_section_search_v18_51r_bundle'
FNAME = FROOT + '_server_sendback.zip'
SNAME = SROOT + '_server_sendback.zip'
S_SHA256 = '553775d6f46c613b6be384fb39abef81a5ab429a3188f691c125200a185d2c8b'
SEEDS = ('original', 'autonomous')
QPOS = np.array([0, 2, 4, 6, 8])
UVEL = np.array([1, 3, 5, 7, 9])


def sha(data):
    return hashlib.sha256(data).hexdigest()


def file_sha(path):
    with open(path, 'rb') as file:
        return hashlib.file_digest(file, 'sha256').hexdigest()


def module_from(path, tag):
    spec = importlib.util.spec_from_file_location(tag, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[tag] = mod
    spec.loader.exec_module(mod)
    return mod


def checked_zip(blob):
    z = zipfile.ZipFile(io.BytesIO(blob))
    if z.testzip() is not None:
        raise RuntimeError('Ancestor archive CRC failure')
    return z


def checked_entries(z, folder, manifest_name, file_key):
    manifest = json.loads(z.read(folder + '/' + manifest_name))
    for name, digest in manifest[file_key].items():
        if sha(z.read(folder + '/' + name)) != digest:
            raise RuntimeError(f'Ancestor digest mismatch: {folder}/{name}')
    return manifest


def preflight_parent():
    """Validate local bundle and available inherited f/e/d source+results."""
    manifest = json.loads((ROOT/'SOURCE_MANIFEST.json').read_text())
    for name, expected in manifest['sha256'].items():
        if file_sha(ROOT/name) != expected:
            raise RuntimeError('v18.52g source SHA mismatch: ' + name)
    fblob = (ROOT/'inputs'/FNAME).read_bytes()
    with checked_zip(fblob) as f:
        fm = checked_entries(f, FROOT, 'SOURCE_MANIFEST.json', 'files')
        checked_entries(f, FROOT, 'OUTPUT_MANIFEST.json', 'files')
        freport = json.loads(f.read(FROOT + '/results/report.json'))
        assert freport['version']=='v18.52f' and freport['all_numerical_checks_passed']
        eblob = f.read(FROOT + '/inputs/' + EROOT + '_server_sendback.zip')
        if sha(eblob) != fm['parent_archive_sha256']:
            raise RuntimeError('f->e archive digest mismatch')
    with checked_zip(eblob) as e:
        em = checked_entries(e, EROOT, 'SOURCE_MANIFEST.json', 'file_sha256')
        checked_entries(e, EROOT, 'OUTPUT_MANIFEST.json', 'file_sha256')
        er = json.loads(e.read(EROOT+'/results/report.json'))
        assert er['version']=='v18.52e' and er['all_numerical_checks_passed']
        dblob = e.read(EROOT + '/inputs/' + DROOT + '_server_sendback.zip')
        if sha(dblob) != em['file_sha256']['inputs/' + DROOT + '_server_sendback.zip']:
            raise RuntimeError('e->d digest mismatch')
        for rel, reference in [('scripts/run_audit.py','parent_v18_52e_run_audit.py'),
                              ('reference/parent_v18_51t_run_audit.py','parent_v18_51t_run_audit.py'),
                              ('reference/parent_v18_52d_run_audit.py','parent_v18_52d_run_audit.py'),
                              ('reference/kinetic_metric_v18_49h.py','kinetic_metric_v18_49h.py')]:
            if e.read(EROOT+'/'+rel) != (ROOT/'reference'/reference).read_bytes():
                raise RuntimeError('Genealogy mismatch: ' + reference)
        old_vectors_blob = e.read(EROOT+'/results/physical_relative_vector_diagnostics.npz')
    with checked_zip(dblob) as d:
        dm = checked_entries(d, DROOT, 'SOURCE_MANIFEST.json', 'file_sha256')
        checked_entries(d, DROOT, 'SHA256SUMS_OUTPUT.json', 'file_sha256')
        dr = json.loads(d.read(DROOT+'/results/report.json'))
        assert dr['version']=='v18.52d' and dr['all_numerical_checks_passed']
        states_blob = d.read(DROOT+'/inputs/historical_v18_51r_physical_states.npz')
        assert sha(states_blob)==dm['file_sha256']['inputs/historical_v18_51r_physical_states.npz']
    with np.load(io.BytesIO(states_blob),allow_pickle=False) as a:
        states = {name: a[name].copy() for name in a.files}
    with np.load(io.BytesIO(old_vectors_blob),allow_pickle=False) as a:
        old_vectors = {name: a[name].copy() for name in a.files}
    assert states['xi'].shape == (513,) and abs(states['xi'][-1] - .2) < 1e-13
    for seed in SEEDS:
        assert states[seed+'_midpoint_Y'].shape == (513,10)
        assert states[seed+'_midpoint_tau'].shape == (513,)
    return states, old_vectors, {'f_sha256':sha(fblob), 'e_sha256':sha(eblob), 'd_sha256':sha(dblob),
                                  'f_checks':len(freport['checks']), 'historical_xi_end':.2,
                                  'historical_state_sha256':sha(states_blob),
                                  'parent_r_archive_sha256':dm['parent_archives']['v18_51r_sha256']}


def get_model_source(s_override=None):
    candidates = ([Path(s_override)] if s_override else []) + [ROOT/'inputs'/SNAME, ROOT.parent/SNAME,
                 Path('/root/elab1/hysteresis')/SNAME]
    source = next((x for x in candidates if x.is_file()), None)
    if source is None:
        raise FileNotFoundError('The exact v18.51s server sendback is required. Set --source-s-path or place '
                                +SNAME+' next to the extracted bundle in /root/elab1/hysteresis.')
    if file_sha(source) != S_SHA256:
        raise RuntimeError('v18.51s source SHA-256 mismatch. Refusing to change the underlying model.')
    with zipfile.ZipFile(source) as s:
        if s.testzip() is not None: raise RuntimeError('v18.51s CRC mismatch')
        nested = s.read(SROOT+'/inputs/'+RROOT+'_server_sendback.zip')
        sscript = s.read(SROOT+'/scripts/run_audit.py')
    return nested, sscript, {'path':str(source), 'sha256':S_SHA256}


def model_rhs(ctx,engine,km,kin,mx):
    start = float(ctx['event'])
    def f(x, state):
        y = state[:10]
        t = start + float(x)
        z = np.r_[y, float(ctx['theta'](t))]
        field = np.asarray(engine.jet.structural_jet(z,km,engine.params,engine.slope)[0][:10],float)
        if not np.isfinite(field).all():raise RuntimeError('Nonfinite original 10D nonlinear RHS')
        vv = y[UVEL]
        den = float(vv@vv)
        if den <= 1e-20:raise RuntimeError('Physical material clock is degenerate')
        clock = float(vv@field[QPOS]/den)
        if not np.isfinite(clock) or clock <= 0:
            raise RuntimeError('Material clock became nonpositive: cannot compare matched tau')
        return np.r_[field,clock]
    return f


def integrate(ctx,engine,km,kin,mx,initial,xi_end,rtol,atol,max_step):
    f = model_rhs(ctx,engine,km,kin,mx)
    z0=np.r_[np.asarray(initial,float),0.]
    sol=solve_ivp(f,(0.,float(xi_end)),z0,method='DOP853',rtol=rtol,atol=atol,
                  max_step=max_step,dense_output=True)
    if not sol.success or not np.isfinite(sol.y).all():
        raise RuntimeError('Genuine nonlinear IVP failed: '+sol.message)
    if min(np.diff(sol.y[10])) <= 0:raise RuntimeError('Integrated material clock not monotonic')
    return sol


def inverse_tau(sol, tau):
    """Recover xi using the actual augmented dense IVP, not an extrapolation."""
    xs=sol.t; ts=sol.y[10]
    if tau <= 0: return 0.
    if tau>=float(ts[-1]):return float(xs[-1])
    ix=int(np.searchsorted(ts,tau,side='right'))
    lo=float(xs[ix-1]); hi=float(xs[ix])
    return float(brentq(lambda x:float(sol.sol(x)[10]-tau),lo,hi,
                        xtol=5e-14,rtol=9e-15))


def physical_at_tau(sol, tau,ctx,engine,km,kin,mx):
    xi=inverse_tau(sol,tau)
    Y=np.asarray(sol.sol(xi)[:10],float)
    field=model_rhs(ctx,engine,km,kin,mx)(xi,np.r_[Y,tau])
    H,J=kin.embedding_and_jac(Y,mx,1.)
    acceleration=J[5:]@field[:10]/field[10]
    kinematic_residual=float(np.linalg.norm((J[:5]@field[:10])/field[10]-H[5:]))
    return {'xi':xi,'Y':Y,'Q':H[:5], 'V':H[5:], 'a':acceleration,'clock':field[10],
            'Qdot_residual':kinematic_residual}


def compare_at_tau(tau,base,kick,ctx,engine,km,kin,mx):
    a=physical_at_tau(base,tau,ctx,engine,km,kin,mx)
    b=physical_at_tau(kick,tau,ctx,engine,km,kin,mx)
    q=b['Q']-a['Q'];v=b['V']-a['V'];acc=b['a']-a['a']
    P=float(v@v);I=float(-q@acc);G=float(q@v)
    D=float(.5*(q@q));C=P-I
    return dict(tau=float(tau),D=D,G=G,C=C,P=P,I=I,R=I/P if P>1e-28 else None,
                distance=float(np.linalg.norm(q)),distance_v=float(np.linalg.norm(v)),
                xi_base=a['xi'],xi_kick=b['xi'],delta_xi=b['xi']-a['xi'],
                dQ=q,dV=v,da=acc,
                material_identity_position_error=max(a['Qdot_residual'],b['Qdot_residual']))


# The current evaluator is bound per trace so event roots evaluate both actual
# nonlinear dense solutions at the same material time.
def find_events(times, vals, evaluator):
    def roots(key):
        arr=np.array([r[key] for r in vals],float)
        out=[]
        for i in range(len(times)-1):
            left,right=float(times[i]),float(times[i+1])
            u,v=float(arr[i]),float(arr[i+1])
            if u*v<0:
                t=float(brentq(lambda t:evaluator(t)[key],left,right,xtol=1e-12))
            elif v==0 and u!=0:
                t=right
            else:
                continue
            if not t>0 or any(abs(t-r['tau'])<1e-10 for r in out): continue
            width=max(1e-8, (right-left)*1e-3)
            va=evaluator(max(0.,t-width))[key];vb=evaluator(min(times[-1],t+width))[key]
            direction=('positive_to_negative' if va>0 and vb<0 else
                       'negative_to_positive' if va<0 and vb>0 else 'tangent_or_unresolved')
            p=evaluator(t)
            out.append(dict(tau=t,crossing=direction,bracket=[left,right],
                            xi_reference=p['xi_base'],xi_kick=p['xi_kick'],
                            drive_phase_delta_xi=p['delta_xi'],D=p['D'],G=p['G'],C=p['C'],
                            R=p['R']))
        return out
    c=roots('C'); g=roots('G')
    # A genuine maximum of D requires a positive-to-negative G crossing.
    dmax=[v for v in g if v['crossing']=='positive_to_negative']
    return dict(convexity_zero=c,rate_zero=g,position_separation_maxima=dmax,
                n_grid_brackets=len(times)-1,continuous_sign_certified=False,
                root_method='Brent on two true DOP853 dense flows inverted to common material tau')


def historical_replay(seed, histories, old_vectors, model, ctx, engine,km,kin,mx):
    """Compare recovered nonlinear physical flow with exact archived 0-.2 trace."""
    xi_history=histories['xi'];Y_history=histories[seed+'_midpoint_Y']
    qerrs=[];verrs=[];tauerrs=[];yerrs=[]
    sol=model
    for i in (0,16,64,128,256,384,512):
        xi=float(xi_history[i]);z=sol.sol(xi)
        h,_=kin.embedding_and_jac(z[:10],mx,1.)
        yerrs.append(float(np.max(np.abs(z[:10]-Y_history[i]))))
        qerrs.append(float(np.max(np.abs(h[:5]-histories[seed+'_midpoint_Q'][i]))))
        verrs.append(float(np.max(np.abs(h[5:]-histories[seed+'_midpoint_V'][i]))))
        tauerrs.append(float(abs(z[10]-histories[seed+'_midpoint_tau'][i])))
    return dict(max_Y_abs=max(yerrs),max_Q_abs=max(qerrs),max_V_abs=max(verrs),max_tau_abs=max(tauerrs))


def historical_kick_replay(label, old_vectors, evaluator):
    tt=old_vectors[label+'__tau'];indices=(0,32,128,256,384,512)
    discrepancies={key:[] for key in ('D','G','C','dQ','dV','da')}
    for i in indices:
        t=float(tt[i]);p=evaluator(t)
        expected={'D':old_vectors[label+'__relative_distance_squared'][i],
                  'G':old_vectors[label+'__relative_position_rate'][i],
                  'C':old_vectors[label+'__second_position_separation_derivative'][i],
                  'dQ':old_vectors[label+'__dQ'][i],
                  'dV':old_vectors[label+'__dV'][i],
                  'da':old_vectors[label+'__da'][i]}
        for key in discrepancies:
            discrepancies[key].append(float(np.linalg.norm(np.atleast_1d(p[key])-np.atleast_1d(expected[key]))))
    return {k+'_max_abs':max(v) for k,v in discrepancies.items()}


def check_convergence(ref_f,ref_c,kick_f,kick_c,tt,ctx,engine,km,kin,mx):
    diffs=[]
    for tau in np.linspace(0,float(tt),17):
        f=compare_at_tau(tau,ref_f,kick_f,ctx,engine,km,kin,mx)
        c=compare_at_tau(tau,ref_c,kick_c,ctx,engine,km,kin,mx)
        diffs.append(dict(tau=tau,deltaQ_norm=float(np.linalg.norm(f['dQ']-c['dQ'])),
                          deltaV_norm=float(np.linalg.norm(f['dV']-c['dV'])),
                          deltaa_norm=float(np.linalg.norm(f['da']-c['da'])),
                          D_abs=abs(f['D']-c['D']),G_abs=abs(f['G']-c['G']),
                          C_abs=abs(f['C']-c['C'])))
    return {k+'_max':float(max(v[k] for v in diffs)) for k in diffs[0] if k!='tau'}


def save_csv(path, rows):
    if not rows:
        # An empty event ledger is an explicit valid outcome, not a missing file.
        if path.name=='material_horizon_events.csv':
            path.write_text('label,kind,tau,crossing,bracket,xi_reference,xi_kick,drive_phase_delta_xi,D,G,C,R\n')
        return
    with open(path,'w',newline='',encoding='utf8') as f:
        wr=csv.DictWriter(f,fieldnames=list(rows[0].keys()))
        wr.writeheader();wr.writerows(rows)


def make_plot(traces, summaries, events, historical_tau):
    fig,axs=plt.subplots(2,2,figsize=(14,9),constrained_layout=True)
    for label, data in traces.items():
        t=np.asarray(data['tau']);code='-' if label.startswith('original') else '--'
        axs[0,0].plot(t,data['P'],code,lw=1.1,label=label+' P')
        axs[0,0].plot(t,data['I'],code,lw=.7,alpha=.5,label=label+' I')
        axs[0,1].plot(t,np.asarray(data['C'])/np.maximum(np.asarray(data['P']),1e-30),code,lw=1.3,label=label)
        axs[1,0].plot(t,data['G'],code,lw=1.3,label=label)
        axs[1,1].plot(t,data['D'],code,lw=1.3,label=label)
        for typ,ax in [('convexity_zero',axs[0,1]),('rate_zero',axs[1,0]),('position_separation_maxima',axs[1,1])]:
            for evt in events[label][typ]:ax.axvline(evt['tau'],color='gray',lw=.4,alpha=.3)
    axs[0,0].set(title='P=||delta V||² and I=-delta Q·delta a',ylabel='Physical Cartesian units')
    axs[0,1].set(title='Relative convexity C/P = 1-I/P',ylabel='C/P')
    axs[1,0].set(title='Separation rate G=D\u2032',ylabel='G')
    axs[1,1].set(title='Squared-position separation D=||delta Q||²/2',ylabel='D')
    for ax in axs.ravel():
        ax.axvline(historical_tau,color='black',ls=':',alpha=.7,lw=1)
        ax.axhline(0,color='gray',lw=.75)
        ax.grid(alpha=.2);ax.set_xlabel('Matched material time tau')
    axs[0,1].legend(fontsize=6,ncol=2)
    fig.suptitle('Q550 v18.52g — authentic horizon extension: not an extrapolation')
    fig.savefig(OUT/'inward_kinematic_material_horizon.png',dpi=155)
    plt.close(fig)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--mode',choices=('preflight','server'),default='server')
    ap.add_argument('--source-s-path',default=None)
    ap.add_argument('--xi-end',type=float,default=.8,help='Extended event-clock horizon; historical = 0.2')
    ap.add_argument('--samples',type=int,default=2049,help='Uniform material-time event-bracketing grid')
    ap.add_argument('--max-step',type=float,default=.00125)
    ap.add_argument('--rtol',type=float,default=2e-10)
    ap.add_argument('--atol',type=float,default=2e-11)
    ap.add_argument('--refine',action=argparse.BooleanOptionalAction,default=True,help='Independent full fine step/tolerance IVPs (default: on)')
    args=ap.parse_args()
    if args.xi_end<=.2 or args.xi_end>3.: ap.error('--xi-end must be in (0.2,3.0]')
    if args.samples<129 or args.samples>20001:ap.error('samples out of supported range 129..20001')
    OUT.mkdir(exist_ok=True,parents=True)
    states,old_vectors,source=preflight_parent()
    print('v18.52g f/e/d lineage source and result SHA/CRC authenticated.',flush=True)
    if args.mode=='preflight':
        (OUT/'PREFLIGHT.json').write_text(json.dumps(dict(stage='PREFLIGHT_AUTHENTICATED_ONLY',source=source,
                              server_nonlinear_integration_completed=False),indent=2)+'\n')
        print('PREFLIGHT PASS — no nonlinear integration, no horizon result claimed',flush=True)
        return
    nested,srun,ex_source=get_model_source(args.source_s_path)
    if sha(nested)!=source['parent_r_archive_sha256']:
        raise RuntimeError('v18.51s nested v18.51r archive hash differs from parent d lineage')
    with checked_zip(nested) as r_ancestor:
        expected=r_ancestor.read(RROOT+'/results/extended_paired_kappa_nonlinear_states.npz')
        if sha(expected)!=source['historical_state_sha256']:
            raise RuntimeError('v18.51s v18.51r historical states differ from parent d saved seed history')
    source['v18_51s']=ex_source
    parent=module_from(ROOT/'reference/parent_v18_52d_run_audit.py','g52_parent_d')
    parentt=module_from(ROOT/'reference/parent_v18_51t_run_audit.py','g52_parent_t')
    kin=module_from(ROOT/'reference/kinetic_metric_v18_49h.py','g52_kin')
    results={}; events={}; traces={}; eventrows=[];ledger=[];checks={}; raw_states={}
    with tempfile.TemporaryDirectory(prefix='q550_v1852g_engine_') as temp:
        print('Staging the authenticated v18.51s->r->q->19 ancestors field...',flush=True)
        parts=parentt.stage_engine(Path(temp),nested,srun)[4]
        _,_,seed_api,observer,engine,unit,drpt,creport,_=parts
        mx=float(engine.params.m_x)
        km=.5*(float(states['kappa_original'])+float(states['kappa_autonomous']))
        xi_save=np.linspace(0,args.xi_end,args.samples)
        raw_states['common_xi']=xi_save
        for seed in SEEDS:
            ctx=seed_api.make_context(observer,seed,engine,unit,drpt,creport)
            y0=states[seed+'_midpoint_Y'][0]
            print('Integrating genuine base nonlinear 10D, seed',seed,flush=True)
            base=integrate(ctx,engine,km,kin,mx,y0,args.xi_end,args.rtol,args.atol,args.max_step)
            raw_states[seed+'__reference__Y_tau']=base.sol(xi_save).T
            old_check=historical_replay(seed,states,old_vectors,base,ctx,engine,km,kin,mx)
            checks[seed+'.historical_reference_replay']=bool(max(old_check.values()) < 3e-5)
            print(' Historical 0.2 replay:',old_check,flush=True)
            if args.refine:
                base2=integrate(ctx,engine,km,kin,mx,y0,args.xi_end,
                                args.rtol/4,args.atol/4,args.max_step/2)
                raw_states[seed+'__reference_refined__Y_tau']=base2.sol(xi_save).T
            for eps in (1e-4,3e-4):
                for signum in (-1,+1):
                    label=f'{seed}_{"positive" if signum>0 else "negative"}_{eps:.1e}'
                    ykick,dv0=parent.kick_vector(y0,kin,mx,signum*eps)
                    print(' Integrating 10D kick',label,flush=True)
                    kick=integrate(ctx,engine,km,kin,mx,ykick,args.xi_end,
                                   args.rtol,args.atol,args.max_step)
                    raw_states[label+'__kick__Y_tau']=kick.sol(xi_save).T
                    end=min(float(base.y[10,-1]),float(kick.y[10,-1]))
                    time=np.linspace(0.,end,args.samples)
                    cache={}
                    def evaluate(t):
                        key=float(t)
                        if key not in cache:
                            cache[key]=compare_at_tau(key,base,kick,ctx,engine,km,kin,mx)
                        return cache[key]
                    data=[evaluate(t) for t in time]
                    ev=find_events(time,data,evaluate)
                    hist=historical_kick_replay(label,old_vectors,evaluate)
                    checks[label+'.historical_full_kick_replay']=bool(
                        hist['dQ_max_abs']<5e-5 and hist['dV_max_abs']<5e-5 and hist['da_max_abs']<1e-3)
                    checks[label+'.identities_and_embedding']=bool(
                        max(z['material_identity_position_error'] for z in data)<1e-8 and
                        max(abs(z['C']-(z['P']-z['I'])) for z in data)<1e-15 and
                        min(z['P'] for z in data)>=0)
                    refinement=None
                    if args.refine:
                        k2=integrate(ctx,engine,km,kin,mx,ykick,args.xi_end,
                                     args.rtol/4,args.atol/4,args.max_step/2)
                        raw_states[label+'__kick_refined__Y_tau']=k2.sol(xi_save).T
                        common=min(end,float(base2.y[10,-1]),float(k2.y[10,-1]))
                        refinement=check_convergence(base,base2,kick,k2,common,
                                                     ctx,engine,km,kin,mx)
                        checks[label+'.refinement_QV']=bool(refinement['deltaQ_norm_max']<1e-5
                                                           and refinement['deltaV_norm_max']<1e-5)
                        # A candidate crossing must survive the independently refined ODE.
                        # Rebracket on the fine flow; this is numerical convergence,
                        # not an interval-arithmetic existence certificate.
                        def fine_evaluate(t):
                            return compare_at_tau(t,base2,k2,ctx,engine,km,kin,mx)
                        refined_events=[]
                        for category,field in [('convexity_zero','C'),('rate_zero','G')]:
                            for evt in ev[category]:
                                lo,hi=evt['bracket']
                                hi=min(hi,common)
                                if hi<=lo:
                                    refined_events.append(dict(category=category,coarse_tau=evt['tau'],fine_tau=None,
                                                               passed=False,reason='outside refined common horizon'))
                                    continue
                                fa,fb=fine_evaluate(lo)[field],fine_evaluate(hi)[field]
                                if fa*fb<0:
                                    te=float(brentq(lambda x:fine_evaluate(x)[field],lo,hi,xtol=1e-12))
                                    passed=abs(te-evt['tau'])<max(1e-5,2*(hi-lo))
                                    refined_events.append(dict(category=category,coarse_tau=evt['tau'],
                                                               fine_tau=te,passed=passed,abs_difference=abs(te-evt['tau'])))
                                else:
                                    refined_events.append(dict(category=category,coarse_tau=evt['tau'],fine_tau=None,
                                                               passed=False,reason='no fine-flow opposite signs'))
                        refinement['events']=refined_events
                        checks[label+'.event_roots_survive_refinement']=all(x['passed'] for x in refined_events)
                    R=np.array([z['R'] if z['R'] is not None else np.nan for z in data])
                    C=np.array([z['C'] for z in data]);G=np.array([z['G'] for z in data]);D=np.array([z['D'] for z in data])
                    if any(not np.isfinite(z) for z in C) or any(not np.isfinite(z) for z in G):
                        raise RuntimeError('Nonfinite physical material diagnostic')
                    has_c=bool(ev['convexity_zero'])
                    has_d=bool(ev['position_separation_maxima'])
                    has_rate=bool(ev['rate_zero'])
                    label_summary=dict(label=label,seed=seed,relative_kick=signum*eps,
                        xi_end_requested=float(args.xi_end),tau_end=end,n_samples=args.samples,
                        nfev_reference=int(base.nfev),nfev_kick=int(kick.nfev),
                        terminal_R=(float(R[-1]) if np.isfinite(R[-1]) else None),
                        maximum_R_sampled=(float(np.nanmax(R)) if np.isfinite(R).any() else None),
                        minimum_C_over_P_sampled=float(np.min(C/np.maximum(np.array([z['P'] for z in data]),1e-28))),
                        minimum_C_sampled=float(np.min(C)),maximum_D_sampled=float(np.max(D)),
                        terminal_G=float(G[-1]),terminal_D=float(D[-1]),
                        terminal_distance=float(data[-1]['distance']),
                        terminal_drive_phase_mismatch=float(data[-1]['delta_xi']),
                        convexity_zero_detected=has_c,rate_zero_detected=has_rate,
                        genuine_position_separation_maximum_detected=has_d,
                        D_sample_steps_decreasing=int(np.sum(np.diff(D)<0)),
                        historical_comparison=hist,refinement=refinement,
                        event_search_certified_exhaustive=False,
                        classification=('POSITION_SEPARATION_MAXIMUM_DETECTED' if has_d else
                                        'CONVEXITY_REVERSAL_DETECTED_WITHOUT_RETURN' if has_c else
                                        'RATE_ZERO_WITHOUT_CONFIRMED_MAXIMUM' if has_rate else
                                        'NO_DETECTED_CONVEXITY_OR_RATE_REVERSAL_IN_EXTENDED_HORIZON'))
                    results[label]=label_summary;events[label]=ev
                    traces[label]={k:np.array([z[k] for z in data],float).tolist()
                                   for k in ('tau','D','G','C','P','I','distance','distance_v','delta_xi')}
                    for typ,vv in ev.items():
                        if not isinstance(vv,list):continue
                        for e in vv:eventrows.append({'label':label,'kind':typ,**e})
                    for idx,z in enumerate(data):
                        ledger.append(dict(label=label,index=idx,tau=z['tau'],
                           xi_reference=z['xi_base'],xi_kick=z['xi_kick'],
                           delta_drive_phase_xi=z['delta_xi'],D=z['D'],G=z['G'],P=z['P'],
                           I=z['I'],C=z['C'],R=z['R'],distance=z['distance'],
                           relative_velocity_norm=z['distance_v']))
                    print('  Event counts C/G/Dmax:',len(ev['convexity_zero']),len(ev['rate_zero']),
                           len(ev['position_separation_maxima']),'R final',f'{R[-1]:.7f}',flush=True)
    source['model_kappa_midpoint']=km
    good=bool(all(checks.values()))
    report=dict(version='v18.52g',title='Inward–Kinematic Competition & Material-Time Horizon Audit',
       stage='SERVER_AUTHENTIC_NONLINEAR_EXTENSION',settings=vars(args),
       source_lineage=source,results=results,events=events,
       all_numerical_checks_passed=good,refinement_performed=bool(args.refine),numerical_checks=checks,
       n_nontrivial_separation_maxima=sum(len(v['position_separation_maxima']) for v in events.values()),
       n_convexity_zeros=sum(len(v['convexity_zero']) for v in events.values()),
       n_rate_zeros=sum(len(v['rate_zero']) for v in events.values()),
       mathematical_interpretation='D=0.5|dq|^2; G=dD/dtau=dq.dv; C=dG/dtau=|dv|^2+dq.da=P-I. '
        'C=0 is not return; true maximum in position separation requires G positive-to-negative. '
        'I>0 alone is inward acceleration, not sufficiency for convergence.',
       limitations=['Numerical dense ODE root finding, not certified continuous-time enclosure',
         'Sign-sampling may miss even-multiplicity or rapidly oscillating roots between search sections',
         'Model parameter kappa fixed at historical midpoint; no change in force',
         'Physical Euclidean five-dimensional Cartesian kinetic embedding specified, not gauge-independent',
         'Material-clock coincidence does not force common external drive phase',
         'No conclusion about long-time asymptotic recovery or physical dissipation',
         'Extension beyond previously authenticated xi<=0.2 is newly simulated, not historically validated'],
       python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__)
    (OUT/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    save_csv(OUT/'material_horizon_ledger.csv',ledger)
    save_csv(OUT/'material_horizon_events.csv',eventrows)
    save_csv(OUT/'kick_horizon_summary.csv',[{k:v for k,v in x.items() if isinstance(v,(str,bool,int,float))} for x in results.values()])
    # Retain physical scalar diagnostics, and the exact source/solver recipe for replay.
    data_arrays={f'{label}__{k}':np.asarray(v) for label,series in traces.items() for k,v in series.items()}
    np.savez_compressed(OUT/'material_horizon_traces.npz',**data_arrays)
    np.savez_compressed(OUT/'extended_10d_nonlinear_states.npz',**raw_states)
    make_plot(traces,results,events,min(float(old_vectors[f'{seed}_negative_1.0e-04__tau'][-1]) for seed in SEEDS))
    lines=['# Q550 v18.52g — Server scientific report','',
      f'Authenticated lineage: f/e/d and historical v18.51s. Numerical checks: {sum(checks.values())}/{len(checks)}.',
      f'New 10D nonlinear drive-phase horizon: xi = [0, {args.xi_end}].',
      'All events are numerically bracketed on dense original-model DOP853 solutions at common material time.',
      'This is not a continuous rigorous enclosure or a global proof of return.','',
      '| Label | Final tau | Final I/P | C=0 roots | G=0 roots | D maxima | Classification |',
      '|---|---:|---:|---:|---:|---:|---|']
    for label,r in results.items():
        v=events[label]
        rr='undefined' if r['terminal_R'] is None else f'{r["terminal_R"]:.7g}'
        lines.append(f'| {label} | {r["tau_end"]:.7g} | {rr} | '
            f'{len(v["convexity_zero"])} | {len(v["rate_zero"])} | {len(v["position_separation_maxima"])} | {r["classification"]} |')
    lines.extend(['','An inward acceleration projection I>0 is not sufficient for C<=0; C<=0 is not sufficient for G<0.',
       'If numerical checks fail, the scientific summary is preliminary and must not be promoted.',
       'Inspect event brackets and optional --refine reruns before interpreting any new crossing.'])
    (OUT/'SCIENTIFIC_REPORT_SERVER.md').write_text('\n'.join(lines)+'\n')
    print('NUMERICAL VALIDATION',sum(checks.values()),'/',len(checks), 'PASS' if good else 'FAIL',flush=True)
    if not good:raise SystemExit('Numerical validation failed — retain diagnostics, do not interpret scientifically')


if __name__=='__main__':main()
