#!/usr/bin/env python3
"""v18.52k: authentic nonlinear continuation from v18.52j endpoint with radial drift and 5D posture transport.

References, comments and report output are in English. Historical 10D engine, kappa,
material clock and five-dimensional physical embedding are never modified.
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
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'
SEAM=6.4
PARENT_ROOT='q550_long_horizon_relative_breathing_posture_recurrence_v18_52j_bundle'
PARENT_NAME=PARENT_ROOT+'_server_sendback.zip'
PARENT=ROOT/'inputs'/PARENT_NAME
H_ROOT='q550_post_crest_competition_separation_rate_arrest_v18_52h_bundle'
G_ROOT='q550_inward_kinematic_competition_material_time_horizon_v18_52g_bundle'
SEEDS=('original','autonomous')


def sha_bytes(x):
    return hashlib.sha256(x).hexdigest()


def sha_file(path):
    with open(path,'rb') as f: return hashlib.file_digest(f,'sha256').hexdigest()


def load_script(path,tag):
    spec=importlib.util.spec_from_file_location(tag,path)
    mod=importlib.util.module_from_spec(spec)
    sys.modules[tag]=mod
    spec.loader.exec_module(mod)
    return mod


def validate_j():
    """Authenticate complete j source/results and all physical restart arrays."""
    if not PARENT.exists():
        raise FileNotFoundError('Missing embedded parent v18.52j sendback '+str(PARENT))
    with zipfile.ZipFile(PARENT) as z:
        if z.testzip() is not None:raise RuntimeError('Parent j CRC failure')
        pre=PARENT_ROOT+'/'
        sm=json.loads(z.read(pre+'SOURCE_MANIFEST.json'))
        om=json.loads(z.read(pre+'OUTPUT_MANIFEST.json'))
        if sm.get('version')!='v18.52j' or om.get('version')!='v18.52j':
            raise RuntimeError('Unexpected parent j manifest version')
        for rel,digest in sm['sha256'].items():
            if sha_bytes(z.read(pre+rel))!=digest:raise RuntimeError('J source digest mismatch '+rel)
        for rel,digest in om['file_sha256'].items():
            if sha_bytes(z.read(pre+rel))!=digest:raise RuntimeError('J output digest mismatch '+rel)
        if z.read(pre+'scripts/run_audit.py')!=(ROOT/'reference/parent_v18_52j_run_audit.py').read_bytes():
            raise RuntimeError('Pinned j generator mismatch')
        report=json.loads(z.read(pre+'results/report.json'))
        if (report.get('version')!='v18.52j' or
            report.get('stage')!='SERVER_AUTHENTIC_POST_I_10D_CONTINUATION' or
            not report.get('all_numerical_checks_passed') or
            not report.get('refinement_performed') or
            abs(report['settings']['xi_end']-SEAM)>1e-12 or
            len(report['results'])!=8 or report['n_local_D_minima']!=16 or
            report['n_local_D_maxima']!=16):
            raise RuntimeError('Parent j failed authentic refined 8-kick validation')
        i_rel='inputs/q550_extended_return_relative_posture_closure_discrimination_v18_52i_bundle_server_sendback.zip'
        i_bytes=z.read(pre+i_rel)
        if sha_bytes(i_bytes)!=sm['parent_i_archive_sha256'] or sha_bytes(i_bytes)!=report['source_lineage']['parent_i_sha256']:
            raise RuntimeError('Embedded i genealogy mismatch')
        with np.load(io.BytesIO(z.read(pre+'results/long_horizon_10d_states.npz')),allow_pickle=False) as f:
            state={k:np.array(f[k],copy=True) for k in f.files}
        with np.load(io.BytesIO(z.read(pre+'results/long_horizon_5d_vectors_and_traces.npz')),allow_pickle=False) as f:
            traces={k:np.array(f[k],copy=True) for k in f.files}
        if abs(float(state['common_xi'][-1])-SEAM)>1e-13 or not np.all(np.diff(state['common_xi'])>0):
            raise RuntimeError('Parent j grid endpoint or monotonicity mismatch')
        for seed in SEEDS:
            for suffix in ('__reference__Y_tau','__reference_refined__Y_tau'):
                key=seed+suffix
                if state[key].shape!=(len(state['common_xi']),11) or not np.isfinite(state[key]).all():
                    raise RuntimeError('Invalid j reference '+key)
            for lbl,item in report['results'].items():
                if item['seed']!=seed:continue
                for suffix in ('__kick__Y_tau','__kick_refined__Y_tau'):
                    key=lbl+suffix
                    if state[key].shape!=(len(state['common_xi']),11) or not np.isfinite(state[key]).all():
                        raise RuntimeError('Invalid j kick '+key)
                if item['n_D_minima']!=2 or item['n_D_maxima']!=2:
                    raise RuntimeError('Parent j extrema topology mismatch')
                for key in ('tau','D','G','C','P','I','dQ','dV','da','angle_from_h_seam_deg'):
                    if lbl+'__'+key not in traces:raise RuntimeError('Missing j trace '+lbl+'__'+key)
                tau_end=min(state[seed+'__reference__Y_tau'][-1,10],state[lbl+'__kick__Y_tau'][-1,10])
                if abs(float(tau_end)-item['tau_end'])>5e-9:
                    raise RuntimeError('Parent j final material clock mismatch')
        return dict(report=report,states=state,traces=traces,
                    parent_sha256=sha_file(PARENT),parent_i_sha256=sha_bytes(i_bytes),
                    source_entries=len(sm['sha256']),output_entries=len(om['file_sha256']))


def stage_j(temp:Path):
    """Stage byte-authenticated j source with its complete i/h/g ancestry."""
    folder=temp/PARENT_ROOT
    with zipfile.ZipFile(PARENT) as z:
        pre=PARENT_ROOT+'/'
        m=json.loads(z.read(pre+'SOURCE_MANIFEST.json'))
        for rel in ['SOURCE_MANIFEST.json',*m['sha256'].keys()]:
            p=folder/rel
            p.parent.mkdir(parents=True,exist_ok=True)
            p.write_bytes(z.read(pre+rel))
    return folder


def resume(g,ctx,engine,km,kin,mx,z0,xi_end,rtol,atol,max_step):
    z0=np.asarray(z0,dtype=float)
    if z0.shape!=(11,) or not np.isfinite(z0).all():raise RuntimeError('Invalid endpoint state')
    f=g.model_rhs(ctx,engine,km,kin,mx)
    # xi remains global: the drive is theta(ctx.event + xi), NOT theta(xi-SEAM).
    sol=solve_ivp(f,(SEAM,float(xi_end)),z0,method='DOP853',rtol=rtol,
                  atol=atol,max_step=max_step,dense_output=True)
    if not sol.success or not np.isfinite(sol.y).all():
        raise RuntimeError('10D physical extension failed: '+sol.message)
    if not np.all(np.diff(sol.y[10])>0):raise RuntimeError('Clock not monotone')
    if np.max(np.abs(sol.sol(SEAM)-z0))>3e-12:raise RuntimeError('Loss of 11D seam state')
    return sol


def vector_geometry(p,anchor,small):
    """Physical 5D vector geometry, independent of an arbitrary 2D projection."""
    q=np.asarray(p['dQ'],float);v=np.asarray(p['dV'],float)
    r=float(np.linalg.norm(q));vnorm=float(np.linalg.norm(v))
    if r<=small:
        return dict(radius=r,radial_rate=None,tangential_speed=None,
                    angular_speed=None,angle_from_h_seam_deg=None,orientation_defined=False)
    unit=q/r
    vr=float(unit@v)
    vt=float(np.linalg.norm(v-vr*unit))
    anchor_norm=float(np.linalg.norm(anchor))
    angle=float(np.degrees(np.arccos(np.clip(float(unit@(anchor/anchor_norm)),-1.,1.)))) if anchor_norm>small else None
    return dict(radius=r,radial_rate=vr,tangential_speed=vt,angular_speed=vt/r,
                angle_from_h_seam_deg=angle,orientation_defined=True)


# The external ledger column keeps the explicit h-seam name, while the in-memory
# trace retains the historical plotting name. This is a naming alias only.
def ledger_value(data, key, idx):
    return data[key][idx]


def compare_posture_vectors(prior_vec, current_vec):
    """Direction, norm and residual of two physical 5D relative positions; no closure inference."""
    p=np.asarray(prior_vec,float);q=np.asarray(current_vec,float)
    if p.shape!=(5,) or q.shape!=(5,) or not np.isfinite(p).all() or not np.isfinite(q).all():
        raise ValueError('Non-finite or non-5D physical vector')
    pr=float(np.linalg.norm(p));cr=float(np.linalg.norm(q))
    angle=float(np.degrees(np.arccos(np.clip(np.dot(p,q)/(pr*cr),-1.,1.)))) if min(pr,cr)>1e-14 else None
    return dict(angle_between_vectors_deg=angle,
                relative_norm_ratio=cr/pr if pr>1e-14 else None,
                distance_in_5d=float(np.linalg.norm(q-p)))


def write_csv(path,rows,keys):
    with open(path,'w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=keys,extrasaction='ignore')
        w.writeheader();w.writerows(rows)


def event_geometry(evaluator,events,anchor,small):
    out=[]
    for category,seq in events.items():
        for item in seq:
            p=evaluator(item['tau'])
            d=vector_geometry(p,anchor,small)
            out.append(dict(event_class=category,**item,**d))
    return out


def event_discrimination(events):
    maxes=[e for e in events['rate_zero'] if e['crossing']=='positive_to_negative']
    mins=[e for e in events['rate_zero'] if e['crossing']=='negative_to_positive']
    if mins:return 'NUMERICAL_LOCAL_D_MINIMUM_DETECTED_NOT_CONFIGURATION_CLOSURE'
    if maxes:return 'SECOND_D_MAXIMUM_WITHOUT_DETECTED_D_MINIMUM'
    if events['convexity_zero']:return 'CURVATURE_TURN_WITHOUT_DETECTED_SEPARATION_RATE_ARREST'
    return 'MONOTONE_SEGMENT_OR_UNRESOLVED_EVENTS_NO_CLOSURE_CLAIM'


def plot_all(traces,results,events,parent):
    fig,axes=plt.subplots(2,3,figsize=(17,9),constrained_layout=True)
    keys=['D_normalized','G_normalized','C_normalized','R','angle_from_h_seam_deg','tangent_over_speed']
    names=['Separation D / D(h start)','Separation rate G / D(h start)',
           'Convexity C / D(h start)','Inward-kinematic ratio I/P',
           '5D separation direction vs fixed h seam (degrees)',
           'Tangential fraction of relative speed']
    for ax,key,title in zip(axes.flat,keys,names):
        for label,d in traces.items():
            style='-' if label.startswith('original') else '--'
            ax.plot(d['tau'],d[key],style,lw=.9,label=label)
            if key=='G_normalized':
                for e in events[label]['rate_zero']:
                    ax.plot(e['tau'],0,marker='x',markersize=4,color='black')
        if key in ('G_normalized','C_normalized'):ax.axhline(0,lw=.8,c='gray')
        if key=='R':ax.axhline(1,lw=.8,c='gray',ls='--')
        ax.set(title=title,xlabel='Matched material time tau');ax.grid(alpha=.22)
    axes[0,0].legend(fontsize=6,ncol=2)
    fig.suptitle('Q550 v18.52k: radial breathing drift and 5D posture transport')
    fig.savefig(OUT/'extended_horizon_relative_breathing.png',dpi=150)
    plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,6),constrained_layout=True)
    for label,d in traces.items():
        ax.plot(d['dQ'][:,0],d['dQ'][:,1],ls='-' if label.startswith('original') else '--',lw=1,label=label)
        ax.scatter(d['dQ'][0,0],d['dQ'][0,1],s=12)
    ax.set(xlabel='Relative physical position component Q1',ylabel='Relative physical position component Q2',
           title='Projection of 5D relative motion — NOT a closure certificate')
    ax.grid(alpha=.2);ax.legend(fontsize=6,ncol=2)
    fig.savefig(OUT/'extended_horizon_projected_5d.png',dpi=140)
    plt.close(fig)


def plot_minimum_drift(results,events,parent):
    """Report dimensionless successive D-minimum drift and full 5D angle changes."""
    fig,axes=plt.subplots(2,2,figsize=(13,9),constrained_layout=True)
    for label,r in results.items():
        prior=parent['report']['events'][label]['local_D_minimum'][-1]
        seq=[prior]+events[label]['local_D_minimum']
        if len(seq)<2:continue
        ts=np.array([x['tau'] for x in seq],float)
        ds=np.array([x['D'] for x in seq],float)
        style='-' if label.startswith('original') else '--'
        baseline=max(float(prior['D']),1e-30)
        axes[0,0].plot(ts,ds/baseline,style,marker='.',ms=3,lw=.8,label=label)
        drift=r['consecutive_minimum_drift']
        axes[0,1].plot([x['to_tau'] for x in drift],
                       [x['delta_tau'] for x in drift],style,marker='.',ms=3,lw=.8)
        axes[1,0].plot([x['to_tau'] for x in drift],
                       [x['angle_between_vectors_deg'] for x in drift],style,marker='.',ms=3,lw=.8)
        axes[1,1].plot([x['to_tau'] for x in drift],
                       [x['D_ratio'] for x in drift],style,marker='.',ms=3,lw=.8)
    for ax,title in zip(axes.flat,[
      'Successive D minima normalized to last j minimum',
      'Material-time interval between successive minima',
      'Angle between successive physical 5D q vectors',
      'Consecutive D-minimum ratio']):
        ax.set(title=title,xlabel='Matched material time tau')
        ax.grid(alpha=.2)
    axes[0,0].legend(ncol=2,fontsize=6)
    axes[1,1].axhline(1,color='gray',lw=.9,ls='--')
    fig.suptitle('Q550 v18.52k: breathing-envelope drift and 5D posture comparison')
    fig.savefig(OUT/'minimum_drift_and_posture.png',dpi=145)
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=('preflight','deep-preflight','server'),default='server')
    parser.add_argument('--source-s-path',default=None)
    parser.add_argument('--xi-end',type=float,default=12.8,help='Global phase coordinate, continuing after j at xi=6.4')
    parser.add_argument('--samples',type=int,default=32769,help='Matched-material-time evaluation and event bracketing grid')
    parser.add_argument('--max-step',type=float,default=.00125)
    parser.add_argument('--rtol',type=float,default=2e-10)
    parser.add_argument('--atol',type=float,default=2e-11)
    parser.add_argument('--refine',action=argparse.BooleanOptionalAction,default=True)
    args=parser.parse_args()
    if not SEAM<args.xi_end<=19.2:parser.error('--xi-end must be >6.4 and <=19.2')
    if not 257<=args.samples<=131073:parser.error('--samples must be 257..131073')
    if not min(args.max_step,args.rtol,args.atol)>0:parser.error('Positive step and tolerances required')
    OUT.mkdir(exist_ok=True,parents=True)
    parent=validate_j()
    print('Parent j authenticated',parent['source_entries'],'source and',parent['output_entries'],
          'output entries, SHA',parent['parent_sha256'],flush=True)
    if args.mode=='preflight':
        (OUT/'PREFLIGHT.json').write_text(json.dumps(dict(version='v18.52k',
            stage='AUTHENTICATED_J_PARENT_ONLY',physical_extension_executed=False,
            parent_j_sha256=parent['parent_sha256'],parent_i_sha256=parent['parent_i_sha256'],
            parent_j_checks_passed=True),indent=2)+'\n')
        print('PREFLIGHT PASS: genuine physical extension not executed',flush=True)
        return
    rows=[];all_vectors={};traces={};events={};results={};checks={};saved={}
    with tempfile.TemporaryDirectory(prefix='q550_v1852k_') as td:
        troot=Path(td)
        staged_j=stage_j(troot)
        j=load_script(staged_j/'scripts/run_audit.py','k52_ancestor_j')
        authenticated_i=j.validate_i()
        if authenticated_i['parent_sha256']!=parent['parent_i_sha256']:
            raise RuntimeError('Deep j->i parent lineage mismatch')
        staged_i=j.stage_i(troot)
        i=load_script(staged_i/'scripts/run_audit.py','k52_ancestor_i')
        authenticated_h=i.validate_h()
        if authenticated_h['parent_sha256']!=parent['report']['source_lineage']['parent_h_sha256']:
            raise RuntimeError('Deep i->h parent lineage mismatch')
        h=load_script(staged_i/'reference/parent_v18_52h_run_audit.py','k52_ancestor_h')
        staged_h=i.stage_h(troot,h)
        authenticated_g=h.validate_parent()
        if authenticated_g['parent_zip_sha256']!=authenticated_i['report']['source_lineage']['embedded_parent_g_sha256']:
            raise RuntimeError('Deep h->g parent lineage mismatch')
        g=load_script(staged_h/'reference/parent_v18_52g_run_audit.py','k52_ancestor_g')
        h.stage_parent_module(troot,g)
        states,old,deep=g.preflight_parent()
        if deep['f_sha256']!=parent['report']['source_lineage']['parent_f_e_d']['f_sha256']:
            raise RuntimeError('Ancestral f source mismatch')
        if args.mode=='deep-preflight':
            (OUT/'PREFLIGHT.json').write_text(json.dumps(dict(
              version='v18.52k',stage='AUTHENTICATED_J_TO_F_ANCESTRY',
              physical_extension_executed=False,
              parent_j_sha256=parent['parent_sha256'],
              parent_i_sha256=parent['parent_i_sha256'],
              parent_h_sha256=authenticated_h['parent_sha256'],
              parent_g_sha256=authenticated_g['parent_zip_sha256'],
              historical_f_sha256=deep['f_sha256']),indent=2)+'\n')
            print('DEEP PREFLIGHT PASS: j->i->h->g->f ancestry; no nonlinear extension',flush=True)
            return
        nested,srun,source_s=g.get_model_source(args.source_s_path)
        if sha_bytes(nested)!=deep['parent_r_archive_sha256']:
            raise RuntimeError('External s->r model source mismatch')
        with g.checked_zip(nested) as zr:
            raw=zr.read(g.RROOT+'/results/extended_paired_kappa_nonlinear_states.npz')
            if sha_bytes(raw)!=deep['historical_state_sha256']:
                raise RuntimeError('Historical nonlinear source mismatch')
        parts=g.module_from(g.ROOT/'reference/parent_v18_51t_run_audit.py','k52_ancestor_t')
        kin=g.module_from(g.ROOT/'reference/kinetic_metric_v18_49h.py','k52_kinetic')
        engine_dir=troot/'engine';engine_dir.mkdir()
        print('Restaging archived nonlinear 10D model and kinetic 5D embedding...',flush=True)
        _,_,seed_api,observer,engine,unit,drpt,creport,_=parts.stage_engine(engine_dir,nested,srun)[4]
        km=.5*(float(states['kappa_original'])+float(states['kappa_autonomous']))
        if abs(km-parent['report']['source_lineage']['fixed_model_kappa'])>1e-15:
            raise RuntimeError('Historical midpoint kappa changed')
        mx=float(engine.params.m_x)
        xs=np.linspace(SEAM,args.xi_end,args.samples)
        saved['common_xi']=xs
        for seed in SEEDS:
            ctx=seed_api.make_context(observer,seed,engine,unit,drpt,creport)
            print('Continuing physical reference',seed,'global xi',SEAM,'to',args.xi_end,flush=True)
            a0=parent['states'][seed+'__reference__Y_tau'][-1]
            ref=resume(g,ctx,engine,km,kin,mx,a0,args.xi_end,args.rtol,args.atol,args.max_step)
            saved[seed+'__reference__Y_tau']=ref.sol(xs).T
            ref2=None
            if args.refine:
                r20=parent['states'][seed+'__reference_refined__Y_tau'][-1]
                ref2=resume(g,ctx,engine,km,kin,mx,r20,args.xi_end,args.rtol/4,args.atol/4,args.max_step/2)
                saved[seed+'__reference_refined__Y_tau']=ref2.sol(xs).T
            for label,original in parent['report']['results'].items():
                if original['seed']!=seed: continue
                print(' Continuing kick',label,flush=True)
                b0=parent['states'][label+'__kick__Y_tau'][-1]
                kick=resume(g,ctx,engine,km,kin,mx,b0,args.xi_end,args.rtol,args.atol,args.max_step)
                saved[label+'__kick__Y_tau']=kick.sol(xs).T
                start=max(float(ref.y[10,0]),float(kick.y[10,0]))
                end=min(float(ref.y[10,-1]),float(kick.y[10,-1]))
                if end<=start:raise RuntimeError('Nonoverlapping material-time extension')
                tt,vals,evaluate=h.scalar_trace(g,ref,kick,start,end,args.samples,ctx,engine,km,kin,mx)
                first=vals[0]
                q0=np.array(first['dQ'],copy=True)
                h_anchor=np.array(authenticated_i['traces'][label+'__dQ'][0],copy=True)
                radii=np.array([float(np.linalg.norm(p['dQ'])) for p in vals])
                small=max(1e-14,1e-8*float(np.max(radii)))
                geo=[vector_geometry(p,h_anchor,small) for p in vals]
                geo_k=[vector_geometry(p,q0,small) for p in vals]
                orientation_valid=all(v['orientation_defined'] for v in geo+geo_k)
                pvals=np.array([v['P'] for v in vals]);ivals=np.array([v['I'] for v in vals])
                floor=max(1e-28,float(pvals.max())*1e-12)
                ratio=np.divide(ivals,pvals,out=np.full_like(pvals,np.nan),where=pvals>floor)
                c_events=h.crossing_events(tt,vals,evaluate,'C','convexity_zero')
                g_events=h.crossing_events(tt,vals,evaluate,'G','rate_zero')
                ev={'convexity_zero':c_events,'rate_zero':g_events,
                    'local_D_minimum':[e for e in g_events if e['crossing']=='negative_to_positive'],
                    'local_D_maximum':[e for e in g_events if e['crossing']=='positive_to_negative']}
                events[label]=ev
                d=np.array([p['D'] for p in vals]);rate=np.array([p['G'] for p in vals]);curv=np.array([p['C'] for p in vals]);
                speed=np.array([p['distance_v'] for p in vals]);
                ang=np.array([np.nan if q['angular_speed'] is None else q['angular_speed'] for q in geo])
                tangent=np.array([np.nan if q['tangential_speed'] is None else q['tangential_speed'] for q in geo])
                theta=np.array([np.nan if q['angle_from_h_seam_deg'] is None else q['angle_from_h_seam_deg'] for q in geo])
                tangfrac=np.divide(tangent,speed,out=np.full_like(speed,np.nan),where=speed>1e-14)
                spin_path=float(np.trapezoid(ang,tt)) if np.isfinite(ang).all() else None
                normalized_d=d/max(d[0],1e-30)
                data=dict(tau=tt,D=d,G=rate,C=curv,P=pvals,I=ivals,R=ratio,
                     D_normalized=normalized_d,G_normalized=rate/max(d[0],1e-30),
                     C_normalized=curv/max(d[0],1e-30),
                     angle_from_seam_deg=np.array([np.nan if x['angle_from_h_seam_deg'] is None else x['angle_from_h_seam_deg'] for x in geo_k]),
                     angle_from_h_seam_deg=theta,tangent_over_speed=tangfrac,
                     radius=radii,radial_rate=np.array([np.nan if x['radial_rate'] is None else x['radial_rate'] for x in geo]),
                     tangential_speed=tangent,angular_speed=ang,
                     xi_reference=np.array([p['xi_base'] for p in vals]),
                     xi_kick=np.array([p['xi_kick'] for p in vals]),
                     delta_xi=np.array([p['delta_xi'] for p in vals]),
                     dQ=np.stack([p['dQ'] for p in vals]),
                     dV=np.stack([p['dV'] for p in vals]),
                     da=np.stack([p['da'] for p in vals]))
                traces[label]=data
                h_t=parent['traces'][label+'__tau'];h_d=parent['traces'][label+'__D'];h_g=parent['traces'][label+'__G'];
                seam_dt=float(tt[0]-h_t[-1])
                seam_d_error=float(abs(d[0]-h_d[-1])/max(abs(h_d[-1]),1e-30))
                seam_g_error=float(abs(rate[0]-h_g[-1])/max(abs(h_g[-1]),1e-30))
                checks[label+'.parent_seam_continuity']=bool(0<=seam_dt<8e-5 and seam_d_error<.002 and seam_g_error<.002)
                checks[label+'.vector_identities']=bool(
                   np.max(np.abs(curv-(pvals-ivals)))<1e-14 and
                   np.max(np.abs(d-.5*radii*radii))<1e-14 and
                   max(x['material_identity_position_error'] for x in vals)<1e-8)
                checks[label+'.physical_geometry']=bool(
                   np.isfinite(data['dQ']).all() and np.isfinite(data['dV']).all() and
                   orientation_valid and np.all(tangent>=-1e-12) and
                   np.max(np.abs((data['radial_rate']*radii)-rate))<1e-12 and
                   np.max(np.abs((data['radial_rate']**2+tangent**2)-speed**2))<1e-12)
                derivative_defects=[]
                hstep=min(1e-4,.15*(end-start)/(args.samples-1))
                for ix in (max(2,args.samples//6),args.samples//2,5*args.samples//6):
                    t=float(tt[ix]);left=evaluate(t-hstep);right=evaluate(t+hstep);mid=evaluate(t)
                    errD=abs((right['D']-left['D'])/(2*hstep)-mid['G'])/max(abs(mid['G']),1e-16)
                    errG=abs((right['G']-left['G'])/(2*hstep)-mid['C'])/max(abs(mid['C']),1e-16)
                    derivative_defects.extend([errD,errG])
                checks[label+'.material_derivative_consistency']=bool(max(derivative_defects)<.04)
                fine_details=None
                if args.refine:
                    b20=parent['states'][label+'__kick_refined__Y_tau'][-1]
                    fine=resume(g,ctx,engine,km,kin,mx,b20,args.xi_end,args.rtol/4,args.atol/4,args.max_step/2)
                    saved[label+'__kick_refined__Y_tau']=fine.sol(xs).T
                    ta=max(start,float(ref2.y[10,0]),float(fine.y[10,0]))
                    tb=min(end,float(ref2.y[10,-1]),float(fine.y[10,-1]))
                    if tb<=ta:raise RuntimeError('No refined common material-time overlap')
                    def e2(t):return g.compare_at_tau(t,ref2,fine,ctx,engine,km,kin,mx)
                    dx=[];dv=[];da=[]
                    for t in np.linspace(ta,tb,21):
                        u=evaluate(float(t));v=e2(float(t))
                        dx.append(np.linalg.norm(u['dQ']-v['dQ']))
                        dv.append(np.linalg.norm(u['dV']-v['dV']))
                        da.append(np.linalg.norm(u['da']-v['da']))
                    err=[float(np.max(dx)),float(np.max(dv)),float(np.max(da))]
                    checks[label+'.refined_physical_vectors']=bool(err[0]<3e-5 and err[1]<3e-5 and err[2]<3e-3)
                    root_check=h.event_convergence(c_events+g_events,(ta,tb),e2,floor)
                    checks[label+'.event_refinement']=bool(all(x['passed'] for x in root_check))
                    fine_details=dict(max_absolute_dQ_dV_da=err,event_convergence=root_check,
                                      common_tau=[ta,tb])
                min_index=int(np.argmin(d))
                min_event=ev['local_D_minimum'][0] if ev['local_D_minimum'] else None
                qat_min=np.array(evaluate(min_event['tau'])['dQ']) if min_event else None
                separation_reduction=float((d[0]-min(e['D'] for e in ev['local_D_minimum']))/d[0]) if ev['local_D_minimum'] else None
                label_status=event_discrimination(ev)
                prior_min=[e for e in parent['report']['events'][label]['local_D_minimum']]
                if len(prior_min)!=2:raise RuntimeError('Parent j expected two D minima')
                prior_min_tau=float(prior_min[-1]['tau'])
                prior_tt=parent['traces'][label+'__tau']
                prior_q=parent['traces'][label+'__dQ']
                prior_vec=np.array([np.interp(prior_min_tau,prior_tt,prior_q[:,k]) for k in range(5)])
                recurrence=[]
                for e in ev['local_D_minimum']:
                    current_q=np.asarray(evaluate(e['tau'])['dQ'],float)
                    recurrence.append({'tau':float(e['tau']), 'xi_reference':float(e['xi_reference']),
                                       **compare_posture_vectors(prior_vec,current_q)})
                minima_all=[dict(source='j_last',tau=float(prior_min_tau),D=float(prior_min[-1]['D']),
                                  q=np.array(prior_vec,copy=True))]
                minima_all += [dict(source='k',tau=float(e['tau']),D=float(e['D']),
                                     q=np.array(evaluate(e['tau'])['dQ'],copy=True))
                               for e in ev['local_D_minimum']]
                min_drift=[]
                for prev,curr in zip(minima_all[:-1],minima_all[1:]):
                    geo_change=compare_posture_vectors(prev['q'],curr['q'])
                    min_drift.append(dict(from_tau=prev['tau'],to_tau=curr['tau'],
                         delta_tau=curr['tau']-prev['tau'],D_from=prev['D'],D_to=curr['D'],
                         D_ratio=curr['D']/prev['D'] if prev['D']>0 else None,
                         D_difference=curr['D']-prev['D'],**geo_change))
                summary=dict(label=label,seed=seed,kick=original['kick'],classification=label_status,
                    xi_start=SEAM,xi_end=args.xi_end,tau_start=start,tau_end=end,
                    j_parent_end_tau=float(h_t[-1]),seam_material_time_gap=seam_dt,
                    seam_D_relative_error=seam_d_error,seam_G_relative_error=seam_g_error,
                    j_parent_terminal_D=float(h_d[-1]),j_parent_terminal_G=float(h_g[-1]),
                    initial_D=float(d[0]),initial_G=float(rate[0]),terminal_D=float(d[-1]),
                    terminal_G=float(rate[-1]),terminal_C=float(curv[-1]),
                    terminal_R=float(ratio[-1]) if np.isfinite(ratio[-1]) else None,
                    n_convexity_zero=len(c_events),n_rate_zero=len(g_events),
                    n_D_minima=len(ev['local_D_minimum']),n_D_maxima=len(ev['local_D_maximum']),
                    min_sample_D=float(d[min_index]),min_sample_tau=float(tt[min_index]),
                    first_D_minimum_tau=(min_event['tau'] if min_event else None),
                    first_D_minimum=(min_event['D'] if min_event else None),
                    D_reduction_after_h_seam=separation_reduction,
                    min_q_distance=float(np.min(radii)),max_angle_from_h_seam_deg=(float(np.nanmax(theta)) if orientation_valid else None),
                    terminal_angle_from_k_seam_deg=(float(data['angle_from_seam_deg'][-1]) if orientation_valid else None),
                    terminal_angle_from_h_seam_deg=(float(theta[-1]) if orientation_valid else None),
                    total_5d_directional_path_rad=spin_path,
                    max_tangent_fraction=(float(np.nanmax(tangfrac)) if np.isfinite(tangfrac).any() else None),
                    terminal_delta_xi=float(data['delta_xi'][-1]),
                    n_ratio_undefined=int(np.sum(~np.isfinite(ratio))),
                    parent_j_prior_min_tau=prior_min_tau,min_event_posture_comparisons=recurrence,
                    consecutive_minimum_drift=min_drift,
                    orientation_defined_at_all_samples=orientation_valid,
                    derivative_identity_max_relative_defect=float(max(derivative_defects)),
                    refinement=fine_details,finite_grid_event_search=True,
                    no_global_configuration_closure_certificate=True)
                results[label]=summary
                for item in event_geometry(evaluate,ev,h_anchor,small):
                    rows.append(dict(label=label,**item))
                print('  ',label,'C zeros',len(c_events),'G zeros',len(g_events),
                      'D minima',len(ev['local_D_minimum']),
                      'terminal D',f'{d[-1]:.6e}', 'angle',f'{theta[-1]:.2f}deg',flush=True)
    all_pass=bool(checks) and all(checks.values())
    np.savez_compressed(OUT/'extended_horizon_10d_states.npz',**saved)
    np.savez_compressed(OUT/'extended_horizon_5d_vectors_and_traces.npz',**{
        label+'__'+k:v for label,t in traces.items() for k,v in t.items()})
    # CSV readback for quick server-side inspection, without losing full 5D vectors in NPZ.
    ledger_keys=['label','sample','tau','xi_reference','xi_kick','delta_xi',
                 'D','G','C','P','I','R','radius','radial_rate','tangential_speed',
                 'angular_speed','angle_from_h_seam_deg','angle_from_seam_deg','tangent_over_speed']
    with (OUT/'extended_horizon_ledger.csv').open('w',newline='',encoding='utf-8') as stream:
        writer=csv.DictWriter(stream,fieldnames=ledger_keys)
        writer.writeheader()
        for label,data in traces.items():
            for idx in range(len(data['tau'])):
                line={'label':label,'sample':idx}
                for k in ledger_keys[2:]:
                    value=ledger_value(data,k,idx)
                    line[k]=float(value) if np.isfinite(value) else ''
                writer.writerow(line)
    event_cols=['label','event_class','kind','tau','crossing','bracket','xi_reference','xi_kick',
       'phase_mismatch_xi','D','G','C','P','I','R','radius','radial_rate','tangential_speed',
       'angular_speed','angle_from_h_seam_deg','orientation_defined']
    write_csv(OUT/'extended_horizon_events.csv',rows,event_cols)
    cols=['label','seed','kick','classification','xi_start','xi_end','tau_start','tau_end',
          'n_convexity_zero','n_rate_zero','n_D_minima','n_D_maxima','first_D_minimum_tau',
          'first_D_minimum','initial_D','terminal_D','terminal_G','terminal_C','terminal_R',
          'min_q_distance','max_angle_from_h_seam_deg','terminal_angle_from_h_seam_deg','terminal_angle_from_k_seam_deg',
          'total_5d_directional_path_rad','max_tangent_fraction','terminal_delta_xi']
    write_csv(OUT/'extended_horizon_summary.csv',results.values(),cols)
    plot_all(traces,results,events,parent)
    plot_minimum_drift(results,events,parent)
    source=dict(parent_j_sha256=parent['parent_sha256'],parent_i_sha256=parent['parent_i_sha256'],
      parent_h_sha256=parent['report']['source_lineage']['parent_h_sha256'],
      parent_f_e_d=parent['report']['source_lineage']['parent_f_e_d'],
      external_v18_51s=source_s,fixed_model_kappa=km,physical_metric='5D Cartesian v18.49h',
      global_event_coordinate_seam=SEAM)
    report=dict(version='v18.52k',title='Radial Breathing Drift & Five-Dimensional Posture Transport Audit',
       stage='SERVER_AUTHENTIC_POST_J_10D_CONTINUATION',settings=vars(args),
       source_lineage=source,results=results,events=events,checks=checks,
       all_numerical_checks_passed=all_pass,refinement_performed=bool(args.refine),
       n_local_D_minima=sum(x['n_D_minima'] for x in results.values()),
       n_local_D_maxima=sum(x['n_D_maxima'] for x in results.values()),
       n_C_zeros=sum(x['n_convexity_zero'] for x in results.values()),
       n_G_zeros=sum(x['n_rate_zero'] for x in results.values()),
       interpretation='D=0.5||dQ||^2; G=dQ.dV; C=||dV||^2-dQ.da. A local D minimum needs G crossing negative-to-positive. '
          'Angle compares actual 5D dQ directions; equal radius is not equal configuration.',
       limitations=['Finite-grid bracketing, no interval-certified exhaustive event search',
                    'D minimum alone does not prove full 10D return, orbit closure, memory, or periodicity',
                    'Matched material time is not identical drive xi: monitor delta_xi',
                    '5D orientation undefined at exact dQ=0; threshold diagnostics are reported',
                    'Physical 5D Euclidean embedding inherited from v18.49h',
                    'External v18.51s ancestry must match recorded SHA',
                    'No parameter, forcing or kappa modifications',
                    'Scalar D drift is not sufficient to label physical stability or instability',
                    'Consecutive minima compare full 5D vectors, not only scalar radii'],
       runtime={'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__})
    (OUT/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    lines=['# Q550 v18.52k — Scientific server report','',
      f'Parent j SHA-256: `{parent["parent_sha256"]}`',
      f'Authentic continuation from global xi={SEAM} to xi={args.xi_end}; matched material time.',
      f'Numerical checks: {sum(checks.values())}/{len(checks)}; all passed: **{all_pass}**.','',
      '| Kick | C zeros | G zeros | New D minima | Terminal D | Direction angle at end | Classification |',
      '|---|---:|---:|---:|---:|---:|---|']
    for label,r in results.items():
        lines.append(f'| {label} | {r["n_convexity_zero"]} | {r["n_rate_zero"]} | {r["n_D_minima"]} | '
                     f'{r["terminal_D"]:.8g} | {r["terminal_angle_from_h_seam_deg"]:.5g}° | {r["classification"]} |')
    lines.extend(['','Recurrent minima of D are **not** automatically recurrence of physical 5D posture, 10D state or phase.',
      'Numerical zeros were bracketed on a finite grid and checked with a refined nonlinear DOP853 flow.',
      'If any check fails, these are diagnostics only; do not promote physical claims.',''])
    (OUT/'SCIENTIFIC_REPORT_SERVER.md').write_text('\n'.join(lines),encoding='utf-8')
    print('Q550 v18.52k completed:',sum(checks.values()),'/',len(checks),
          'checks; D minima',report['n_local_D_minima'],flush=True)
    if not all_pass: raise RuntimeError('Numerical checks failed; preserve files for diagnosis only')

if __name__=='__main__':main()
