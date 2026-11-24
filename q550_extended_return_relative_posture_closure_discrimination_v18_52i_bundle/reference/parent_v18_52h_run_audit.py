#!/usr/bin/env python3
"""Q550 v18.52h: genuine post-crest continuation from SHA-verified v18.52g 10D endpoints.

Eight paired trajectories are continued under the unchanged nonlinear engine,
fixed midpoint kappa and independently integrated physical material clocks.
No extrapolated force, interpolated physical tail or promised return is allowed.
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
import scipy
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results'
PARENT_ROOT = 'q550_inward_kinematic_competition_material_time_horizon_v18_52g_bundle'
PARENT_NAME = PARENT_ROOT + '_server_sendback.zip'
PARENT = ROOT / 'inputs' / PARENT_NAME
XI_SEAM = .8
SOURCE_S_NAME = 'q550_parameter_sensitivity_transport_transverse_posture_attribution_v18_51s_bundle_server_sendback.zip'
SEEDS = ('original', 'autonomous')


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def file_sha(path: Path) -> str:
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def import_script(path: Path, tag: str):
    spec = importlib.util.spec_from_file_location(tag, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[tag] = mod
    spec.loader.exec_module(mod)
    return mod


def validate_parent():
    """No writes, no integration. Verify g archive, results, and all source entries."""
    with zipfile.ZipFile(PARENT) as z:
        if z.testzip() is not None:
            raise RuntimeError('Parent g ZIP CRC failure')
        pref = PARENT_ROOT + '/'
        source = json.loads(z.read(pref+'SOURCE_MANIFEST.json'))
        output = json.loads(z.read(pref+'OUTPUT_MANIFEST.json'))
        if source.get('version') != 'v18.52g' or output.get('version') != 'v18.52g':
            raise RuntimeError('Parent manifests are not v18.52g')
        for rel, digest in source['sha256'].items():
            if sha(z.read(pref+rel)) != digest:
                raise RuntimeError('Parent g source SHA mismatch: '+rel)
        for rel, digest in output['file_sha256'].items():
            if sha(z.read(pref+rel)) != digest:
                raise RuntimeError('Parent g output SHA mismatch: '+rel)
        if z.read(pref+'scripts/run_audit.py') != (ROOT/'reference/parent_v18_52g_run_audit.py').read_bytes():
            raise RuntimeError('Local g parent generator is not the archived generator')
        report = json.loads(z.read(pref+'results/report.json'))
        if (report.get('version') != 'v18.52g' or
                report.get('stage') != 'SERVER_AUTHENTIC_NONLINEAR_EXTENSION' or
                not report.get('all_numerical_checks_passed') or
                not report.get('refinement_performed')):
            raise RuntimeError('Parent g is not an authenticated successful refined nonlinear run')
        if (abs(report['settings']['xi_end'] - XI_SEAM)>1e-13 or
                len(report['results']) != 8 or
                report['n_convexity_zeros'] != 8):
            raise RuntimeError('Unexpected parent g horizon or event topology')
        if any(item['rate_zero_detected'] or item['genuine_position_separation_maximum_detected']
               for item in report['results'].values()):
            raise RuntimeError('Parent g has a return; h post-crest assumption must be revised')
        with np.load(io.BytesIO(z.read(pref+'results/extended_10d_nonlinear_states.npz')),
                     allow_pickle=False) as data:
            state = {key:data[key].copy() for key in data.files}
        with np.load(io.BytesIO(z.read(pref+'results/material_horizon_traces.npz')),
                     allow_pickle=False) as data:
            traces = {key:data[key].copy() for key in data.files}
        if abs(float(state['common_xi'][-1])-XI_SEAM)>1e-13:
            raise RuntimeError('Parent 10D state grid has wrong endpoint')
        for seed in SEEDS:
            for suffix in ('__reference__Y_tau', '__reference_refined__Y_tau'):
                assert state[seed+suffix].shape[1]==11
                assert np.isfinite(state[seed+suffix]).all()
            for label, res in report['results'].items():
                if res['seed'] != seed:
                    continue
                for suffix in ('__kick__Y_tau','__kick_refined__Y_tau'):
                    assert state[label+suffix].shape==state[seed+'__reference__Y_tau'].shape
                    assert np.isfinite(state[label+suffix]).all()
                gtau = min(state[seed+'__reference__Y_tau'][-1,10],
                           state[label+'__kick__Y_tau'][-1,10])
                if abs(float(gtau)-res['tau_end'])>5e-10:
                    raise RuntimeError('Parent G terminal common-material time mismatch')
                for scalar in ('tau','D','G','C','P','I'):
                    if label+'__'+scalar not in traces:
                        raise RuntimeError('Missing parent trace '+label+'__'+scalar)
        return {'report':report,'states':state,'traces':traces,
                'parent_zip_sha256':file_sha(PARENT),
                'parent_source_entries':len(source['sha256']),
                'parent_output_entries':len(output['file_sha256']),
                'parent_v18_51s_sha256':source['required_external_v18_51s_archive_sha256']}


def stage_parent_module(temp: Path, g):
    """Stage only g resources used by the exact parent preflight/10D engine."""
    staged = temp / PARENT_ROOT
    with zipfile.ZipFile(PARENT) as z:
        source=json.loads(z.read(PARENT_ROOT+'/SOURCE_MANIFEST.json'))
        needed=['SOURCE_MANIFEST.json',*source['sha256'].keys()]
        for rel in needed:
            p=staged/rel
            p.parent.mkdir(parents=True,exist_ok=True)
            p.write_bytes(z.read(PARENT_ROOT+'/'+rel))
    g.ROOT=staged
    g.OUT=staged/'results'
    return staged


def resume_ivp(g, ctx, engine, km, kin, mx, zseam, xi_end, rtol, atol, max_step):
    zseam=np.asarray(zseam,float)
    if zseam.shape != (11,) or not np.isfinite(zseam).all():
        raise RuntimeError('Invalid archived 10D+clock endpoint')
    if zseam[10] <= 0:
        raise RuntimeError('Nonpositive archived material-time endpoint')
    # Crucial: xi is GLOBAL (0.8 to xi_end). The exact parent model evaluates
    # theta(ctx.event + xi), thus the drive phase is never reset at restart.
    sol=solve_ivp(g.model_rhs(ctx,engine,km,kin,mx),(XI_SEAM,float(xi_end)),zseam,
                  method='DOP853',rtol=rtol,atol=atol,max_step=max_step,
                  dense_output=True)
    if not sol.success or not np.isfinite(sol.y).all():
        raise RuntimeError('Authentic 10D continuation failed: '+sol.message)
    if np.min(np.diff(sol.y[10]))<=0:
        raise RuntimeError('Material clock became non-monotone on continuation')
    if abs(float(sol.sol(XI_SEAM)[10])-float(zseam[10]))>2e-12:
        raise RuntimeError('Lost clock continuity at g->h seam')
    return sol


def scalar_trace(g,sol_b,sol_k,lo,hi,count,ctx,engine,km,kin,mx):
    ts=np.linspace(lo,hi,count)
    cache={}
    def eval_at(t):
        t=float(t)
        if t < lo - 1e-12 or t > hi+1e-12:
            raise RuntimeError('Attempt to evaluate outside shared material-time horizon')
        if t not in cache:
            cache[t]=g.compare_at_tau(t,sol_b,sol_k,ctx,engine,km,kin,mx)
        return cache[t]
    values=[eval_at(t) for t in ts]
    return ts,values,eval_at


def crossing_events(ts,values,evaluator,key,tag):
    """Detect sign-changing roots; tangencies outside brackets are not certified."""
    arr=np.array([v[key] for v in values],float)
    results=[]
    for j in range(len(ts)-1):
        a,b=float(ts[j]),float(ts[j+1])
        fa,fb=float(arr[j]),float(arr[j+1])
        if fa*fb<0:
            root=float(brentq(lambda t:evaluator(t)[key],a,b,
                             xtol=1e-12,rtol=9e-15))
        elif fb==0 and fa!=0:
            root=b
        else:
            continue
        if any(abs(root-e['tau'])<1e-10 for e in results):
            continue
        pad=min(0.02*(b-a),root-ts[0],ts[-1]-root)
        if pad<=0:
            direction='boundary_or_unresolved'
        else:
            x=evaluator(root-pad)[key];y=evaluator(root+pad)[key]
            direction=('positive_to_negative' if x>0 and y<0 else
                       'negative_to_positive' if x<0 and y>0 else 'tangent_or_unresolved')
        p=evaluator(root)
        results.append(dict(kind=tag,tau=root,crossing=direction,bracket=[a,b],
                            xi_reference=p['xi_base'],xi_kick=p['xi_kick'],
                            phase_mismatch_xi=p['delta_xi'],
                            D=p['D'],G=p['G'],C=p['C'],P=p['P'],I=p['I'],R=p['R']))
    return results


def ratio_crest_events(ts,values,evaluator,ratio_floor):
    """Locate *sample-bracketed* local peaks of I/P, refine on original IVP dense output."""
    R=np.array([float(v['I']/v['P']) if v['P']>ratio_floor else np.nan for v in values])
    out=[]
    for i in range(1,len(ts)-1):
        if not (np.isfinite(R[i-1:i+2]).all() and R[i-1]<R[i] and R[i]>R[i+1]):
            continue
        a,b=float(ts[i-1]),float(ts[i+1])
        obj=lambda tau: -(lambda p:p['I']/p['P'] if p['P']>ratio_floor else -np.inf)(evaluator(tau))
        r=minimize_scalar(obj,bounds=(a,b),method='bounded',options={'xatol':1e-11})
        if not (r.success and np.isfinite(r.fun)):
            continue
        t=float(r.x);v=evaluator(t)
        if v['P']<=ratio_floor or t<=a or t>=b:
            continue
        out.append(dict(kind='ratio_local_crest',tau=t,crossing='sample_bracketed_local_maximum',
                        bracket=[a,b],R=v['I']/v['P'],G=v['G'],D=v['D'],C=v['C'],
                        xi_reference=v['xi_base'],xi_kick=v['xi_kick'],
                        phase_mismatch_xi=v['delta_xi']))
    return out


def parent_crest(parent,label):
    v=parent['traces']
    I=np.asarray(v[label+'__I']);P=np.asarray(v[label+'__P']);T=np.asarray(v[label+'__tau'])
    good=P>max(1e-28,float(P.max())*1e-12)
    R=np.divide(I,P,out=np.full_like(P,np.nan),where=good)
    ix=int(np.nanargmax(R))
    return dict(tau=float(T[ix]),R=float(R[ix]),G=float(v[label+'__G'][ix]),
                xi_max_parent=XI_SEAM,method='sampled_from_SHA_verified_g_trace')


def event_convergence(events, fine_times, fine_evaluator, floor):
    """Independently re-bracket every coarse sign-changing root on refined 10D flow."""
    results=[]
    for ev in events:
        if ev['kind'] not in ('convexity_zero','rate_zero'):
            continue
        lo,hi=ev['bracket'];field='C' if ev['kind']=='convexity_zero' else 'G'
        width=hi-lo
        # Allow small movement of bracket under re-integration.
        a=max(float(fine_times[0]),lo-2*width)
        b=min(float(fine_times[1]),hi+2*width)
        if b<=a:
            results.append(dict(kind=ev['kind'],coarse_tau=ev['tau'],fine_tau=None,passed=False,
                                reason='event beyond refined material horizon'))
            continue
        q=np.linspace(a,b,17)
        vv=[fine_evaluator(t)[field] for t in q]
        bs=[(q[i],q[i+1]) for i in range(16) if vv[i]*vv[i+1]<0]
        if not bs:
            results.append(dict(kind=ev['kind'],coarse_tau=ev['tau'],fine_tau=None,passed=False,
                                reason='no fine-flow sign-changing bracket'))
            continue
        root=min((float(brentq(lambda x:fine_evaluator(x)[field],float(l),float(r),xtol=1e-12))
                  for l,r in bs),key=lambda t:abs(t-ev['tau']))
        difference=abs(root-ev['tau'])
        results.append(dict(kind=ev['kind'],coarse_tau=ev['tau'],fine_tau=root,
                            difference=difference,passed=difference<max(2e-5,3*width)))
    return results


def write_csv(path,rows,cols):
    with path.open('w',newline='',encoding='utf-8') as stream:
        writer=csv.DictWriter(stream,fieldnames=cols,extrasaction='ignore')
        writer.writeheader();writer.writerows(rows)


def render(traces,reports,events,parent):
    fig,axs=plt.subplots(2,2,figsize=(15,9),constrained_layout=True)
    seam_taus=[]
    for label,dat in traces.items():
        ts=np.asarray(dat['tau']);p=np.asarray(dat['P']);inn=np.asarray(dat['I'])
        r=np.divide(inn,p,out=np.full_like(p,np.nan),where=p>max(1e-28,float(p.max())*1e-12))
        n='--' if label.startswith('autonomous') else '-'
        lines=axs[0,0].plot(ts,r,n,lw=1.0,label=label)
        color=lines[0].get_color()
        axs[0,1].plot(ts,1-r,n,lw=1.0,color=color)
        crest=reports[label]['g_parent_ratio_crest']
        seam_taus.append(reports[label]['parent_tau_end'])
        gt=parent['traces']; old_t=np.asarray(gt[label+'__tau']); old_P=np.asarray(gt[label+'__P']);old_I=np.asarray(gt[label+'__I'])
        old_R=np.divide(old_I,old_P,out=np.full_like(old_P,np.nan),where=old_P>max(1e-28,float(old_P.max())*1e-12))
        old_mask=old_t>=max(0.0,crest['tau']-0.035)
        axs[0,0].plot(old_t[old_mask],old_R[old_mask],n,lw=.65,color=color,alpha=.4)
        axs[0,1].plot(old_t[old_mask],1-old_R[old_mask],n,lw=.65,color=color,alpha=.4)
        axs[0,0].scatter([crest['tau']],[crest['R']],color=color,s=10)
        d=np.asarray(dat['D']);gg=np.asarray(dat['G'])
        norm=max(reports[label]['parent_terminal_G'],1e-30)
        axs[1,0].plot(ts,gg/norm,n,color=color,lw=1)
        axs[1,1].plot(ts,d/max(reports[label]['parent_terminal_D'],1e-30),n,color=color,lw=1)
        for e in events[label]['rate_zero']:
            axs[1,0].scatter([e['tau']],[0],color='red',marker='x',s=15)
        for e in events[label]['convexity_zero']:
            axs[0,1].scatter([e['tau']],[0],color='black',s=9)
    axs[0,0].set(title='Post-crest competition I/P (g crest marked)',ylabel='I/P')
    axs[0,1].set(title='Relative curvature C/P = 1 - I/P',ylabel='C/P')
    axs[1,0].set(title='Separation rate G relative to g-end G',ylabel='G/G(g endpoint)')
    axs[1,1].set(title='Squared separation D relative to g-end D',ylabel='D/D(g endpoint)')
    seam_tau=float(np.median(seam_taus))
    for ax in axs.flat:
        ax.axvline(seam_tau,color='gray',lw=.7,ls=':',alpha=.65)
        ax.grid(alpha=.22);ax.set_xlabel('Common physical material time tau')
    axs[0,0].axhline(1,color='black',ls='--',lw=.7)
    axs[0,1].axhline(0,color='black',ls='--',lw=.7)
    axs[1,0].axhline(0,color='black',ls='--',lw=.7)
    axs[0,0].legend(fontsize=6,ncol=2)
    fig.suptitle('Q550 v18.52h — authentic nonlinear post-crest continuation (unchanged model)')
    fig.savefig(OUT/'post_crest_competition_rate_arrest.png',dpi=150)
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=['preflight','server'],default='server')
    parser.add_argument('--source-s-path',default=None)
    parser.add_argument('--xi-end',type=float,default=1.6,help='Global event xi, starts at g endpoint 0.8')
    parser.add_argument('--samples',type=int,default=4097,help='Common material-time root bracketing samples')
    parser.add_argument('--max-step',type=float,default=.00125)
    parser.add_argument('--rtol',type=float,default=2e-10)
    parser.add_argument('--atol',type=float,default=2e-11)
    parser.add_argument('--refine',action=argparse.BooleanOptionalAction,default=True)
    args=parser.parse_args()
    if not (XI_SEAM<args.xi_end<=3.):parser.error('--xi-end must be >0.8 and <=3')
    if not 129<=args.samples<=30001:parser.error('--samples must be 129..30001')
    if not (args.max_step>0 and args.rtol>0 and args.atol>0):parser.error('positive tolerances/step required')
    OUT.mkdir(parents=True,exist_ok=True)
    parent=validate_parent()
    print('Authenticated v18.52g source+results',parent['parent_source_entries'],
          parent['parent_output_entries'],'files',flush=True)
    if args.mode=='preflight':
        (OUT/'PREFLIGHT.json').write_text(json.dumps(dict(version='v18.52h',
            stage='AUTHENTICATED_G_PARENT_ONLY',nonlinear_extension_executed=False,
            parent_zip_sha256=parent['parent_zip_sha256'],parent_g_passed=True),indent=2)+'\n')
        print('PREFLIGHT PASS; no physical horizon extension calculated',flush=True)
        return
    reports={};events={};traces={};ledger=[];event_rows=[];checks={};saved={}
    with tempfile.TemporaryDirectory(prefix='q550_v1852h_') as td:
        tmp=Path(td)
        g=import_script(ROOT/'reference/parent_v18_52g_run_audit.py','h52_ancestor_g')
        stage_parent_module(tmp,g)
        states,old,deep=g.preflight_parent()
        if deep['f_sha256']!=parent['report']['source_lineage']['f_sha256']:
            raise RuntimeError('Parent g->f deep audit divergence')
        nested,srun,source_s=g.get_model_source(args.source_s_path)
        if sha(nested)!=deep['parent_r_archive_sha256']:
            raise RuntimeError('External s->r archive differs from authenticated parent lineage')
        with g.checked_zip(nested) as parent_r:
            archived_history=parent_r.read(g.RROOT+'/results/extended_paired_kappa_nonlinear_states.npz')
            if sha(archived_history)!=deep['historical_state_sha256']:
                raise RuntimeError('External model ancestry and historical initial physical states diverged')
        parts=g.module_from(g.ROOT/'reference/parent_v18_51t_run_audit.py','h52_ancestor_t')
        parent_d=g.module_from(g.ROOT/'reference/parent_v18_52d_run_audit.py','h52_ancestor_d')
        kin=g.module_from(g.ROOT/'reference/kinetic_metric_v18_49h.py','h52_kinetic')
        engineroot=tmp/'engine'
        engineroot.mkdir(parents=True,exist_ok=True)
        print('Restaging authentic v18.51s->r->q->19 nonlinear engine...',flush=True)
        _,_,seed_api,observer,engine,unit,drpt,creport,_=parts.stage_engine(engineroot,nested,srun)[4]
        km=.5*(float(states['kappa_original'])+float(states['kappa_autonomous']))
        mx=float(engine.params.m_x)
        expected_km=parent['report']['source_lineage']['model_kappa_midpoint']
        if abs(km-expected_km)>1e-15:
            raise RuntimeError('Midpoint kappa inconsistent with parent g')
        xs=np.linspace(XI_SEAM,args.xi_end,args.samples)
        saved['common_xi']=xs
        for seed in SEEDS:
            ctx=seed_api.make_context(observer,seed,engine,unit,drpt,creport)
            base0=parent['states'][seed+'__reference__Y_tau'][-1]
            print('Continuing reference',seed,'from xi=0.8 with original clock',flush=True)
            base=resume_ivp(g,ctx,engine,km,kin,mx,base0,args.xi_end,
                            args.rtol,args.atol,args.max_step)
            saved[seed+'__reference__Y_tau']=base.sol(xs).T
            if args.refine:
                base20=parent['states'][seed+'__reference_refined__Y_tau'][-1]
                base2=resume_ivp(g,ctx,engine,km,kin,mx,base20,args.xi_end,
                                 args.rtol/4,args.atol/4,args.max_step/2)
                saved[seed+'__reference_refined__Y_tau']=base2.sol(xs).T
            for label, pr in parent['report']['results'].items():
                if pr['seed']!=seed:continue
                kick0=parent['states'][label+'__kick__Y_tau'][-1]
                print(' Continuing kick',label,flush=True)
                kick=resume_ivp(g,ctx,engine,km,kin,mx,kick0,args.xi_end,
                                args.rtol,args.atol,args.max_step)
                saved[label+'__kick__Y_tau']=kick.sol(xs).T
                lo=max(float(base.y[10,0]),float(kick.y[10,0]))
                hi=min(float(base.y[10,-1]),float(kick.y[10,-1]))
                if hi<=lo:raise RuntimeError('No shared physical material horizon after g seam')
                ts,data,eval_at=scalar_trace(g,base,kick,lo,hi,args.samples,
                                             ctx,engine,km,kin,mx)
                coarse_c=crossing_events(ts,data,eval_at,'C','convexity_zero')
                coarse_g=crossing_events(ts,data,eval_at,'G','rate_zero')
                pvals=np.array([x['P'] for x in data]);floor=max(1e-28,float(pvals.max())*1e-12)
                peaks=ratio_crest_events(ts,data,eval_at,floor)
                dmax=[x for x in coarse_g if x['crossing']=='positive_to_negative']
                dmin=[x for x in coarse_g if x['crossing']=='negative_to_positive']
                ev=dict(convexity_zero=coarse_c,rate_zero=coarse_g,
                        genuine_separation_maximum=dmax,local_separation_minimum=dmin,
                        ratio_local_crest=peaks,
                        ratio_return_to_below_one=[x for x in coarse_c if x['crossing']=='negative_to_positive'])
                events[label]=ev
                gparent=parent_crest(parent,label)
                gg=np.array([x['G'] for x in data]);dd=np.array([x['D'] for x in data]);
                c=np.array([x['C'] for x in data]);ii=np.array([x['I'] for x in data]);
                ratio=np.divide(ii,pvals,out=np.full_like(pvals,np.nan),where=pvals>floor)
                pg=parent['traces'];g_last_G=float(pg[label+'__G'][-1]);g_last_D=float(pg[label+'__D'][-1])
                # g last trace is at min(parent base clock, parent kick clock),
                # while the h common start is at max(...). The offset is ~1--3 microseconds.
                seam_dt=lo-float(pg[label+'__tau'][-1])
                seam_G_rel=abs(float(gg[0])-g_last_G)/max(abs(g_last_G),1e-30)
                seam_D_rel=abs(float(dd[0])-g_last_D)/max(abs(g_last_D),1e-30)
                checks[label+'.parent_seam_continuity']=bool(0<=seam_dt<5e-5 and
                                    seam_G_rel<.01 and seam_D_rel<.01)
                checks[label+'.physical_identity']=bool(np.isfinite(c).all() and np.isfinite(gg).all()
                     and np.max(abs(c-(pvals-ii)))<1e-14 and
                     max(x['material_identity_position_error'] for x in data)<1e-8)
                # Two independent numerical derivative identities checked via small
                # symmetric differences of actual dense nonlinear trajectories.
                derivative_err=[]
                half=min(1e-4,.15*(hi-lo)/(args.samples-1))
                for ix in (max(2,args.samples//6),args.samples//2,5*args.samples//6):
                    t=float(ts[ix]);aa=eval_at(t-half);bb=eval_at(t+half)
                    center=eval_at(t)
                    dD=(bb['D']-aa['D'])/(2*half)
                    dG=(bb['G']-aa['G'])/(2*half)
                    derivative_err.append((abs(dD-center['G'])/max(abs(center['G']),1e-16),
                                           abs(dG-center['C'])/max(abs(center['C']),1e-16)))
                deriv_max=max(max(a,b) for a,b in derivative_err)
                checks[label+'.independent_D_G_C_derivatives']=bool(deriv_max<.025)
                refined=None
                if args.refine:
                    k20=parent['states'][label+'__kick_refined__Y_tau'][-1]
                    k2=resume_ivp(g,ctx,engine,km,kin,mx,k20,args.xi_end,
                                  args.rtol/4,args.atol/4,args.max_step/2)
                    saved[label+'__kick_refined__Y_tau']=k2.sol(xs).T
                    lo2=max(float(base2.y[10,0]),float(k2.y[10,0]))
                    hi2=min(float(base2.y[10,-1]),float(k2.y[10,-1]))
                    common_a=max(lo,lo2);common_b=min(hi,hi2)
                    if common_b<=common_a:raise RuntimeError('No refined shared material-time horizon')
                    def fine_eval(t):
                        return g.compare_at_tau(t,base2,k2,ctx,engine,km,kin,mx)
                    times=np.linspace(common_a,common_b,19)
                    errs=[]
                    for t in times:
                        a=eval_at(float(t));b=fine_eval(float(t))
                        errs.append((np.linalg.norm(a['dQ']-b['dQ']),
                                     np.linalg.norm(a['dV']-b['dV']),
                                     np.linalg.norm(a['da']-b['da'])))
                    max_err=np.max(errs,axis=0)
                    checks[label+'.refined_physical_vectors']=bool(max_err[0]<2e-5 and
                                                                    max_err[1]<2e-5 and max_err[2]<2e-3)
                    fine_events=event_convergence(coarse_c+coarse_g,(common_a,common_b),fine_eval,floor)
                    checks[label+'.crossings_survive_refinement']=bool(all(e['passed'] for e in fine_events))
                    refined=dict(dQ_dV_da_max_abs=max_err.tolist(),events=fine_events,
                                 shared_tau_interval=[common_a,common_b])
                finite_ratio=ratio[np.isfinite(ratio)]
                min_g_idx=int(np.argmin(gg))
                initial_normalizer=max(g_last_G,1e-30)
                status=('CONFIRMED_NUMERICAL_D_MAXIMUM' if dmax else
                        'RATE_ZERO_NO_POSITIVE_TO_NEGATIVE_CROSSING' if coarse_g else
                        'SECOND_CURVATURE_REVERSAL_WITHOUT_DETECTED_ARREST' if ev['ratio_return_to_below_one'] else
                        'POST_CREST_CONTINUATION_NO_RATE_ZERO_DETECTED')
                out=dict(label=label,seed=seed,kick=float(pr['relative_kick']),
                    classification=status,xi_start=XI_SEAM,xi_end=args.xi_end,
                    tau_start=lo,tau_end=hi,parent_tau_end=float(pg[label+'__tau'][-1]),
                    parent_terminal_G=g_last_G,parent_terminal_D=g_last_D,
                    seam_time_gap=seam_dt,seam_G_relative_change=seam_G_rel,
                    seam_D_relative_change=seam_D_rel,
                    g_parent_ratio_crest=gparent,
                    terminal_G=float(gg[-1]),terminal_D=float(dd[-1]),
                    terminal_C=float(c[-1]),terminal_R=(float(ratio[-1]) if np.isfinite(ratio[-1]) else None),
                    maximum_R_sampled=(float(np.max(finite_ratio)) if finite_ratio.size else None),
                    minimum_R_sampled=(float(np.min(finite_ratio)) if finite_ratio.size else None),
                    minimum_G_sampled=float(np.min(gg)),
                    minimum_G_over_parent_terminal_G=float(gg[min_g_idx]/initial_normalizer),
                    G_minimum_sample_time=float(ts[min_g_idx]),
                    n_C_zeros=len(coarse_c),n_G_zeros=len(coarse_g),
                    n_D_maxima=len(dmax),n_R_secondary_crests=len(peaks),
                    n_C_recoveries=len(ev['ratio_return_to_below_one']),
                    D_decrease_samples=int(np.sum(np.diff(dd)<0)),
                    ratio_floor_P=floor,ratio_undefined_samples=int(np.sum(~np.isfinite(ratio))),
                    material_phase_mismatch_end=float(data[-1]['delta_xi']),
                    derivative_identity_max_relative_defect=deriv_max,
                    refinement=refined,bracket_search_not_exhaustive=True)
                reports[label]=out
                traces[label]={k:np.asarray([x[k] for x in data],float).tolist()
                               for k in ('tau','D','G','C','P','I','delta_xi','distance','distance_v')}
                for event_class,ee in ev.items():
                    for x in ee:
                        event_rows.append({'label':label,'event_class':event_class,**x})
                for ix,x in enumerate(data):
                    ledger.append(dict(label=label,index=ix,tau=x['tau'],xi_reference=x['xi_base'],
                         xi_kick=x['xi_kick'],delta_xi=x['delta_xi'],D=x['D'],G=x['G'],
                         C=x['C'],P=x['P'],I=x['I'],R=(x['I']/x['P'] if x['P']>floor else ''),
                         separation=x['distance'],relative_speed=x['distance_v']))
                print(' ',label,'C/G/Dmax/R-crest',len(coarse_c),len(coarse_g),
                      len(dmax),len(peaks),'terminal G',f'{gg[-1]:.8e}',flush=True)
    all_checks=bool(checks) and all(checks.values())
    source=dict(parent_g_sendback_sha256=parent['parent_zip_sha256'],
                ancestral_f_e_d_sha256=deep,
                original_v18_51s=source_s,
                fixed_model_kappa=km,
                physical_metric='Original 5D Cartesian kinetic embedding v18.49h',
                global_event_clock_seam=XI_SEAM)
    report=dict(version='v18.52h',title='Post-Crest Competition & Separation-Rate Arrest Audit',
          stage='SERVER_AUTHENTIC_POST_CREST_10D_CONTINUATION',settings=vars(args),source_lineage=source,
          results=reports,events=events,numerical_checks=checks,
          all_numerical_checks_passed=all_checks,refinement_performed=bool(args.refine),
          n_rate_zeros=sum(len(v['rate_zero']) for v in events.values()),
          n_genuine_position_maxima=sum(len(v['genuine_separation_maximum']) for v in events.values()),
          n_C_zeros=sum(len(v['convexity_zero']) for v in events.values()),
          n_return_to_R_below_one=sum(len(v['ratio_return_to_below_one']) for v in events.values()),
          interpretation='C=D\'\'=P-I, G=D\', D=0.5||delta Q||^2. '
             'A new C zero does not imply return; only G>0 crossing to G<0 establishes a local D maximum.',
          limitations=['Authentic 10D continuation, but not an interval-certified continuous-time proof',
                       'Crossings bracketed by a finite material-time grid; tangencies may be missed',
                       'Parent g 0.8 endpoint and independent material clocks must remain valid',
                       'Ratio I/P undefined when kinetic P is too small',
                       'Fixed historical midpoint kappa and original 5D Cartesian metric',
                       'Coincident material time is not equal drive phase; monitor delta xi',
                       'No claim of asymptotic relaxation or material-memory certification'],
          python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__)
    (OUT/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    write_csv(OUT/'post_crest_ledger.csv',ledger,
        ['label','index','tau','xi_reference','xi_kick','delta_xi','D','G','C','P','I','R','separation','relative_speed'])
    write_csv(OUT/'post_crest_events.csv',event_rows,
        ['label','event_class','kind','tau','crossing','bracket','xi_reference','xi_kick','phase_mismatch_xi','D','G','C','P','I','R'])
    write_csv(OUT/'post_crest_summary.csv',reports.values(),
        ['label','seed','kick','classification','xi_start','xi_end','tau_start','tau_end',
         'parent_terminal_G','terminal_G','terminal_C','terminal_R','terminal_D','minimum_G_sampled',
         'minimum_G_over_parent_terminal_G','n_C_zeros','n_G_zeros','n_D_maxima','n_R_secondary_crests',
         'n_C_recoveries','D_decrease_samples','seam_time_gap','seam_G_relative_change',
         'seam_D_relative_change','ratio_undefined_samples'])
    np.savez_compressed(OUT/'post_crest_traces.npz',
       **{label+'__'+key:np.asarray(vals) for label,item in traces.items() for key,vals in item.items()})
    np.savez_compressed(OUT/'post_crest_continued_10d_states.npz',**saved)
    render(traces,reports,events,parent)
    lines=['# Q550 v18.52h — Scientific server sendback','',
       f'Exact parent g archive SHA-256: `{parent["parent_zip_sha256"]}`.',
       f'Original nonlinear model resumed at xi={XI_SEAM} to xi={args.xi_end}.',
       f'Numerical checks: {sum(checks.values())}/{len(checks)}; all passed: **{all_checks}**.','',
       '| Kick | New C zeros | G zeros | D maxima | Return R<1 | Terminal G | Classification |',
       '|---|---:|---:|---:|---:|---:|---|']
    for lbl,r in reports.items():
        lines.append(f'| {lbl} | {r["n_C_zeros"]} | {r["n_G_zeros"]} | '
                     f'{r["n_D_maxima"]} | {r["n_C_recoveries"]} | '
                     f'{r["terminal_G"]:.9g} | {r["classification"]} |')
    lines.extend(['','A local maximum of squared separation D requires a '+
         'positive-to-negative zero of G. A new C zero, R crest or reduced positive G alone does not establish return.',
         'The run is a numerical finite-horizon audit with independent refinements, not a rigorous interval enclosure.',
         'If any numerical check fails, do not promote the physical interpretation.',''])
    (OUT/'SCIENTIFIC_REPORT_SERVER.md').write_text('\n'.join(lines),encoding='utf-8')
    print('Q550 v18.52h completed; checks',sum(checks.values()),'/',len(checks),
          'event counts C/G/Dmax',report['n_C_zeros'],report['n_rate_zeros'],
          report['n_genuine_position_maxima'],flush=True)
    if not all_checks:
        raise RuntimeError('Numerical checks failed; results saved for diagnosis, not scientific promotion')


if __name__=='__main__':
    main()
