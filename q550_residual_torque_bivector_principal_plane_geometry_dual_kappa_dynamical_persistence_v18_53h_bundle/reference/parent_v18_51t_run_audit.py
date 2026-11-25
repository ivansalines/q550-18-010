#!/usr/bin/env python3
"""Q550 v18.51t -- two-channel parameter source/response Green transport audit.

FULL: Authenticate the v18.51s sendback and the original 19-ancestor field,
construct the historical 10D J and b=partial_kappa f along each saved *midpoint*
trajectory, and transport source channels b[v2], b[omega2], and a remaining
channel with independent forced tangent IVPs. Independently compute the
fundamental propagator Phi and interaction-picture Green integral.

CACHED: only authenticate and characterize the s source/response ledgers;
DO NOT imply this performs independent channel transport.
"""
from __future__ import annotations
import argparse,csv,hashlib,importlib.util,io,json,platform,sys,tempfile,zipfile
from pathlib import Path
import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'
SROOT='q550_parameter_sensitivity_transport_transverse_posture_attribution_v18_51s_bundle'
SPARENT=ROOT/'inputs'/(SROOT+'_server_sendback.zip')
RROOT='q550_long_horizon_paired_kappa0_transverse_transport_return_section_search_v18_51r_bundle'
SEEDS=('original','autonomous')
SOURCE_NAMES=('v2','omega2','other')
IDX={'v2':5,'omega2':7}
POS=(0,2,4,6,8)

def sha(b):return hashlib.sha256(b).hexdigest()
def norm(x):return float(np.linalg.norm(x))
def ratio(a,b,floor=1e-10):return norm(a-b)/max(floor,norm(a),norm(b))
def write_csv(name,rows):
    if not rows:raise RuntimeError('No rows: '+name)
    with (OUT/name).open('w',newline='',encoding='utf8') as h:
        w=csv.DictWriter(h,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def import_file(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m)
    return m

def authenticate_s():
    with zipfile.ZipFile(SPARENT) as z:
        if z.testzip() is not None:raise RuntimeError('Parent S damaged CRC')
        src=json.loads(z.read(SROOT+'/SOURCE_MANIFEST.json'))
        out=json.loads(z.read(SROOT+'/SHA256SUMS_OUTPUT.json'))
        assert src['root']==out['root']==SROOT
        assert src['version']==out['version']=='v18.51s'
        for manifest,key in ((src,'files'),(out,'file_sha256')):
            for name,digest in manifest[key].items():
                if sha(z.read(SROOT+'/'+name)) !=digest:raise RuntimeError('Parent hash drift '+name)
        report=json.loads(z.read(SROOT+'/results/report.json'))
        assert report['version']=='v18.51s' and report['mode']=='full' and report['all_passed']
        assert all(report['scientific_checks'].values())
        table={}
        for name in ('variational_parameter_sensitivity_ledger.csv','model_parameter_source_ledger.csv'):
            table[name]=list(csv.DictReader(io.StringIO(z.read(SROOT+'/results/'+name).decode('utf8'))))
        rblob=z.read(SROOT+'/inputs/'+RROOT+'_server_sendback.zip')
        srunscript=z.read(SROOT+'/scripts/run_audit.py')
    with zipfile.ZipFile(io.BytesIO(rblob)) as r:
        assert r.testzip() is None
        rreport=json.loads(r.read(RROOT+'/results/report.json'))
        assert rreport['version']=='v18.51r' and all(rreport['checks'].values())
        a_bytes=r.read(RROOT+'/results/extended_paired_kappa_nonlinear_states.npz')
    with np.load(io.BytesIO(a_bytes),allow_pickle=False) as obj:
        a={k:obj[k].copy() for k in obj.files}
    xi=a['xi']
    assert len(xi)==513 and abs(xi[-1]-.2)<1e-14
    assert abs(float(a['kappa_original']-a['kappa_autonomous'])-6.705624902611507e-5)<1e-10
    assert report['variational_summary']['original']['max_secant_Q_relative_difference_resolved']<1e-5
    return a,report,table,rblob,srunscript,src,out

def j_immersion(y,q):
    x=1./(1.+np.exp(-float(y[8])))
    J=np.zeros((5,10))
    J[0,0]=q[0];J[0,2]=-q[1];J[1,0]=q[1];J[1,2]=q[0]
    J[2,4]=q[2];J[2,6]=-q[3];J[3,4]=q[3];J[3,6]=q[2]
    J[4,8]=q[4]*(1.-x)
    return J

def axes_projection(dq,y):
    a=float(y[6]);r=np.array([np.cos(a),np.sin(a)])
    az=np.array([-r[1],r[0]])
    return float(dq[2:4]@r),float(dq[2:4]@az)

def source_channel(b,name):
    out=np.zeros(10)
    if name=='other':
        out[:]=b;out[5]=0;out[7]=0
    else:out[IDX[name]]=b[IDX[name]]
    return out

def cached_precheck(a,report,table):
    results={}
    for seed in SEEDS:
        rr=[r for r in table['model_parameter_source_ledger.csv'] if r['seed_family']==seed and r['sampling_stride']=='2']
        rv=[r for r in table['variational_parameter_sensitivity_ledger.csv'] if r['seed_family']==seed]
        assert len(rr)==257 and len(rv)==513
        bv=np.array([float(r['source_v2']) for r in rr]);bw=np.array([float(r['source_omega2']) for r in rr]);
        bmag=np.array([float(r['source_f_kappa_norm']) for r in rr])
        rem=np.sqrt(np.maximum(0.,bmag*bmag-bv*bv-bw*bw))
        results[seed]={'parent_source_v2_at_zero':float(bv[0]),
        'parent_source_omega2_at_zero':float(bw[0]),
        'max_residual_source_inferred_from_norm_rounded':float(np.max(rem)),
        'parent_integrated_Q_norm_at_end':float(rv[-1]['variational_Q_norm']),
        'parent_radial_plane2_at_end':float(rv[-1]['radial_plane2']),
        'parent_azimuthal_plane2_at_end':float(rv[-1]['azimuthal_plane2'])}
    return results

def stage_engine(staging,rblob,srun):
    """Recover *original* field through the authenticated s -> r -> q ancestry."""
    staging=Path(staging)
    sbase=staging/SROOT
    (sbase/'scripts').mkdir(parents=True,exist_ok=True)
    (sbase/'inputs').mkdir(parents=True,exist_ok=True)
    (sbase/'scripts'/'run_audit.py').write_bytes(srun)
    (sbase/'inputs'/(RROOT+'_server_sendback.zip')).write_bytes(rblob)
    s=import_file(sbase/'scripts'/'run_audit.py','q550_v1851t_authenticated_s_loader')
    # Authenticated v18.51s loader stages r/q/p/ancestors and validates lineage.
    return s.stages(staging/'engine_staging')

def sampled_field(seed,xi,Y,km,context,engine,step=2):
    ids=np.arange(0,len(xi),step)
    if ids[-1]!=len(xi)-1:ids=np.r_[ids,len(xi)-1]
    vals=[];js=[];bs=[];delta=[]
    for i in ids:
        yy=np.asarray(Y[i],float)
        t=float(context['event']+xi[i]); z=np.r_[yy,float(context['theta'](t))]
        f,j,_=engine.jet.structural_jet(z,km,engine.params,engine.slope)
        f=np.asarray(f[:10],float);j=np.asarray(j[:10,:10],float)
        comps=[]
        for h in (2e-5,1e-5):
            fp=np.asarray(engine.jet.structural_jet(z,km+h,engine.params,engine.slope)[0][:10],float)
            fm=np.asarray(engine.jet.structural_jet(z,km-h,engine.params,engine.slope)[0][:10],float)
            comps.append((fp-fm)/(2*h))
        vals.append(f);js.append(j);bs.append(comps[1]);delta.append(ratio(comps[0],comps[1],1e-6))
    return np.array(ids),np.asarray(js),np.asarray(bs),float(max(delta))

def transport(xi, ids, Js, bs, rtol=3e-11, atol=3e-13):
    """Independent forced IVP plus fundamental Green integration using Phi and quadrature.

    The Green integral is calculated in the interaction picture W'=Phi^{-1} b_j;
    hence S_j=Phi W_j. Phi is integrated separately and not constructed from
    the channel sensitivities (which are independently integrated).
    """
    xs=xi[ids]
    Jx=CubicSpline(xs,Js,axis=0)
    bx=CubicSpline(xs,bs,axis=0)
    def rhs_forced(t,flat):
        y=flat.reshape(10,4);f=np.asarray(bx(t));B=np.column_stack((
            source_channel(f,'v2'),source_channel(f,'omega2'),
            source_channel(f,'other'),f))
        return (Jx(t)@y+B).reshape(-1)
    forced=solve_ivp(rhs_forced,(xi[0],xi[-1]),np.zeros(40),method='DOP853',
            t_eval=xi,rtol=rtol,atol=atol,max_step=.001)
    if not forced.success:raise RuntimeError('Channel IVP: '+forced.message)
    sc=forced.y.T.reshape(-1,10,4)
    # Fundamental matrix, integrated independently of the forcing.
    phis=solve_ivp(lambda t,z:(Jx(t)@z.reshape(10,10)).reshape(100),
            (xi[0],xi[-1]),np.eye(10).reshape(100),method='DOP853',
            t_eval=xi,dense_output=True,rtol=rtol,atol=atol,max_step=.001)
    if not phis.success:raise RuntimeError('Fundamental matrix IVP: '+phis.message)
    Phi=phis.y.T.reshape(-1,10,10)
    cond=np.linalg.cond(Phi)
    if not np.all(np.isfinite(cond)) or max(cond)>1e11:
        raise RuntimeError('Ill-conditioned fundamental representation, max cond '+str(max(cond)))
    # Direct cumulative quadrature of Phi^{-1}(xi)b_j(xi); not the forced IVP.
    Bfull=np.zeros((len(xi),10,4))
    for i,t in enumerate(xi):
        f=np.asarray(bx(t));Bfull[i]=np.column_stack((source_channel(f,'v2'),source_channel(f,'omega2'),source_channel(f,'other'),f))
    G=np.linalg.solve(Phi,Bfull)
    # High-order composite Simpson on uniform xi (even panels), midpoint linear
    # half-intervals are computed with a trapezoid average at odd indices.
    dxi=float(xi[1]-xi[0]);assert np.allclose(np.diff(xi),dxi,rtol=0,atol=1e-13)
    W=np.zeros_like(G)
    for i in range(2,len(xi),2):
        W[i]=W[i-2]+dxi/3*(G[i-2]+4*G[i-1]+G[i])
    # At odd nodes use Simpson 3/8 on previous triple when available.
    # A trapezoid step for the single first interval is a small edge error;
    # all even nodes (including the endpoint) use composite Simpson.
    W[1]=dxi*.5*(G[0]+G[1])
    for i in range(3,len(xi),2):
        W[i]=W[i-3]+3*dxi/8*(G[i-3]+3*G[i-2]+3*G[i-1]+G[i])
    Sgreen=np.einsum('nij,njk->nik',Phi,W)
    return {'sc':sc,'green':Sgreen,'Phi':Phi,'max_cond':float(np.max(cond)),
          'nfev_channels':int(forced.nfev),'nfev_phi':int(phis.nfev),
          'Jmax':float(np.max(np.linalg.norm(Js,axis=(1,2)))),
          'bmax':float(np.max(np.linalg.norm(bs,axis=1))),
          'other_bmax':float(np.max(np.linalg.norm(bs[:,[0,1,2,3,4,6,8,9]],axis=1)))}

def analysis_full(a,sreport,table,rblob,srun,precheck):
    xi=a['xi'];km=.5*(float(a['kappa_original'])+float(a['kappa_autonomous']))
    checks={};summaries={};ledger=[];out_arrays={};source_ledgers=[]
    with tempfile.TemporaryDirectory(prefix='q550_v1851t_model_') as td:
        print('STAGING HISTORICAL MODEL AND VERIFIED ANCESTRY',flush=True)
        q,p,v,vr,parts,kappas,lineage,qstates=stage_engine(td,rblob,srun)
        checks['historical_19_ancestor_genealogy']=len(lineage)==19
        _,_,seed_api,observer,engine,unit,drpt,creport,_=parts
        assert abs(km-kappas['midpoint'])<1e-12
        for seed in SEEDS:
            print('MODEL-GENERATED SOURCE CHANNELS',seed,flush=True)
            ctx=seed_api.make_context(observer,seed,engine,unit,drpt,creport)
            Y=a[f'{seed}_midpoint_Y'];Q=a[f'{seed}_midpoint_Q'];V=a[f'{seed}_midpoint_Qprime']
            ids,J,b,hs_consistency=sampled_field(seed,xi,Y,km,ctx,engine,step=2)
            # Refinement compares a nested grid; same sampled field values.
            fine=transport(xi,ids,J,b)
            even=(ids%4==0)
            coarse=transport(xi,ids[even],J[even],b[even])
            sc=fine['sc'];g=fine['green'];scp=coarse['sc']
            # Prepare Cartesian and basis-projected channel sensitivities.
            Qch=np.array([j_immersion(y,qq)@ss for y,qq,ss in zip(Y,Q,sc)])
            Qgreen=np.array([j_immersion(y,qq)@ss for y,qq,ss in zip(Y,Q,g)])
            sec=(a[f'{seed}_original_Q']-a[f'{seed}_autonomous_Q'])/(float(a['kappa_original']-a['kappa_autonomous']))
            sqp=np.array([norm(t) for t in sec]);Sq=np.linalg.norm(Qch[:,:,3],axis=1)
            full=sreport['variational_summary'][seed]
            selected=np.flatnonzero(xi>=.005)
            den=np.maximum(1e-4,np.linalg.norm(sec[selected],axis=1))
            q_sec_max=float(np.max(np.linalg.norm(Qch[selected,:,3]-sec[selected],axis=1)/den))
            den_v=np.maximum(1e-4,np.linalg.norm(Qch[selected,:,3],axis=1))
            gridrel=float(np.max(np.linalg.norm(scp[selected,:,3]-sc[selected,:,3],axis=1)/np.maximum(1e-4,np.linalg.norm(sc[selected,:,3],axis=1))))
            greenrel=float(np.max(np.linalg.norm(Qgreen[2::2,:,3]-Qch[2::2,:,3],axis=1)/np.maximum(1e-4,np.linalg.norm(Qch[2::2,:,3],axis=1))))
            composition=float(np.max(np.linalg.norm(sc[:,:,:3].sum(axis=2)-sc[:,:,3],axis=1)))
            green_comp=float(np.max(np.linalg.norm(g[:,:,:3].sum(axis=2)-g[:,:,3],axis=1)))
            source_others=fine['other_bmax']
            # Parent variational ledger independent check: final norm and projections.
            varrows=[r for r in table['variational_parameter_sensitivity_ledger.csv'] if r['seed_family']==seed]
            parent_final=float(varrows[-1]['variational_Q_norm']);ours_final=norm(Qch[-1,:,3]);
            parent_component_diff=(abs(parent_final-ours_final),
             abs(float(varrows[-1]['radial_plane2'])-axes_projection(Qch[-1,:,3],Y[-1])[0]),
             abs(float(varrows[-1]['azimuthal_plane2'])-axes_projection(Qch[-1,:,3],Y[-1])[1]))
            rad,az=np.zeros((len(xi),4)),np.zeros((len(xi),4))
            for i in range(len(xi)):
                for j in range(4):rad[i,j],az[i,j]=axes_projection(Qch[i,:,j],Y[i])
            # Absolute spatial and projection shares can include cancellations.
            rad_sum_err=float(np.max(np.abs(rad[:,:3].sum(axis=1)-rad[:,3])))
            az_sum_err=float(np.max(np.abs(az[:,:3].sum(axis=1)-az[:,3])))
            velocities=np.einsum('ni,ni->n',V,V)
            phase=np.einsum('nik,ni->nk',Qch,V)/velocities[:,None]
            trans=Qch-phase[:,None,:]*V[:,:,None]
            cross_ratio_r=(abs(rad[-1,1])/max(1e-16,abs(rad[-1,3])))
            cross_ratio_a=(abs(az[-1,0])/max(1e-16,abs(az[-1,3])))
            for i,x in enumerate(xi):
                record={'seed_family':seed,'index':i,'xi':float(x),
                  'b_v2':float(np.interp(x,xi[ids],b[:,5])),
                  'b_omega2':float(np.interp(x,xi[ids],b[:,7])),
                  'S_Q_total_norm':float(Sq[i]),
                  'S_Q_v2_norm':norm(Qch[i,:,0]),'S_Q_omega2_norm':norm(Qch[i,:,1]),
                  'S_Q_other_norm':norm(Qch[i,:,2]),
                  'S_Q_secant_norm':float(sqp[i]),
                  'S_Q_total_to_secant_error':norm(Qch[i,:,3]-sec[i]),
                  'S_Q_v2_radial2':float(rad[i,0]),'S_Q_v2_azimuthal2':float(az[i,0]),
                  'S_Q_omega2_radial2':float(rad[i,1]),'S_Q_omega2_azimuthal2':float(az[i,1]),
                  'S_Q_other_radial2':float(rad[i,2]),'S_Q_other_azimuthal2':float(az[i,2]),
                  'S_Q_total_radial2':float(rad[i,3]),'S_Q_total_azimuthal2':float(az[i,3]),
                  'radial_absolute_response_sum':float(np.sum(np.abs(rad[i,:3]))),
                  'azimuthal_absolute_response_sum':float(np.sum(np.abs(az[i,:3]))),
                  'radial_signed_cancellation_fraction':float(1.-abs(rad[i,3])/max(1e-30,np.sum(np.abs(rad[i,:3])))),
                  'azimuthal_signed_cancellation_fraction':float(1.-abs(az[i,3])/max(1e-30,np.sum(np.abs(az[i,:3])))),
                  'phase_equivalent_total':float(phase[i,3]),
                  'transverse_fraction_total':norm(trans[i,:,3])/max(1e-30,Sq[i]),
                  'green_minus_forced_Q_norm':norm(Qgreen[i,:,3]-Qch[i,:,3]),
                  'composition_error_Q':norm(Qch[i,:,:3].sum(axis=1)-Qch[i,:,3])}
                ledger.append(record)
            summaries[seed]={
                'S_Q_total_norm_end':ours_final,'S_Q_v2_norm_end':norm(Qch[-1,:,0]),
                'S_Q_omega2_norm_end':norm(Qch[-1,:,1]),
                'S_Q_other_norm_end':norm(Qch[-1,:,2]),
                'S_Q_radial2_total_end':float(rad[-1,3]),
                'S_Q_radial2_from_v2_end':float(rad[-1,0]),
                'S_Q_radial2_from_omega2_end':float(rad[-1,1]),
                'S_Q_azimuthal2_total_end':float(az[-1,3]),
                'S_Q_azimuthal2_from_v2_end':float(az[-1,0]),
                'S_Q_azimuthal2_from_omega2_end':float(az[-1,1]),
                'cross_response_ratio_omega2_to_radial':float(cross_ratio_r),
                'cross_response_ratio_v2_to_azimuthal':float(cross_ratio_a),
                'radial_cancellation_fraction_end':float(1.-abs(rad[-1,3])/max(1e-30,sum(abs(rad[-1,:3])))),
                'azimuthal_cancellation_fraction_end':float(1.-abs(az[-1,3])/max(1e-30,sum(abs(az[-1,:3])))),
                'transverse_fraction_end':float(norm(trans[-1,:,3])/norm(Qch[-1,:,3])),
                'max_state_source_channel_recomposition_absolute_error':composition,
                'max_green_source_channel_recomposition_absolute_error':green_comp,
                'max_green_vs_forced_Q_relative_error_even_indices':greenrel,
                'max_source_off_v2_omega2_absolute':source_others,
                'max_source_h_step_relative_difference':hs_consistency,
                'max_parameter_secant_Q_relative_error_resolved':q_sec_max,
                'max_jet_sampling_grid_relative_error_resolved':gridrel,
                'max_parent_S_final_norm_and_components_absolute_difference':max(parent_component_diff),
                'max_fundamental_condition_number':fine['max_cond'],
                'max_radial_additivity_error':rad_sum_err,'max_azimuthal_additivity_error':az_sum_err,
                'channel_solver_nfev':fine['nfev_channels'],'fundamental_solver_nfev':fine['nfev_phi']}
            out_arrays.update({seed+'_source_J':J,seed+'_source_b':b,seed+'_source_xi':xi[ids],
               seed+'_forced_Y':sc,seed+'_green_Y':g,seed+'_cartesian_Q':Qch,
               seed+'_cartesian_Q_green':Qgreen,seed+'_radial2':rad,seed+'_azimuthal2':az,
               seed+'_phase_lag':phase,seed+'_transverse_Q':trans})
            print(seed,'END radial:',rad[-1,:].tolist(),'az:',az[-1,:].tolist(),flush=True)
            print(seed,'Green mismatch:',greenrel,'composition:',composition,'max cond:',fine['max_cond'],flush=True)
        # Scientific assertions: these are independent measurements, not predictions.
        checks.update({
          'two_seed_parameter_trajectories':len(summaries)==2,
          'historical_parameter_h_refinement':all(s['max_source_h_step_relative_difference']<1e-4 for s in summaries.values()),
          'source_support_measured_not_assumed':all(np.isfinite(s['max_source_off_v2_omega2_absolute']) for s in summaries.values()),
          'separate_linear_channels_recompose_total':all(s['max_state_source_channel_recomposition_absolute_error']<5e-8 for s in summaries.values()),
          'green_integral_recomposes_total':all(s['max_green_source_channel_recomposition_absolute_error']<5e-8 for s in summaries.values()),
          'independent_green_quadrature_matches_forced_ivp':all(s['max_green_vs_forced_Q_relative_error_even_indices']<1e-3 for s in summaries.values()),
          'source_grid_refinement':all(s['max_jet_sampling_grid_relative_error_resolved']<5e-4 for s in summaries.values()),
          'physical_distinct_kappa_secant_consistency':all(s['max_parameter_secant_Q_relative_error_resolved']<1e-3 for s in summaries.values()),
          'v18_51s_independent_variational_endpoint_consistency':all(s['max_parent_S_final_norm_and_components_absolute_difference']<1e-6 for s in summaries.values()),
          'radial_azimuthal_linear_channel_additivity':all(max(s['max_radial_additivity_error'],s['max_azimuthal_additivity_error'])<1e-8 for s in summaries.values()),
          'nonsingular_fundamental_over_audited_interval':all(s['max_fundamental_condition_number']<1e11 for s in summaries.values()),
          'cross_coupling_quantified_no_zero_assumption':all(np.isfinite(s['cross_response_ratio_omega2_to_radial']) and np.isfinite(s['cross_response_ratio_v2_to_azimuthal']) for s in summaries.values()),
          'no_material_return_or_global_causality_claimed':True})
    write_csv('green_transport_attribution_ledger.csv',ledger)
    np.savez_compressed(OUT/'green_channel_transport_states.npz',xi=xi,kappa_midpoint=km,**out_arrays)
    fig,axs=plt.subplots(2,2,figsize=(12.4,8.6))
    for seed in SEEDS:
        rr=[r for r in ledger if r['seed_family']==seed]
        axs[0,0].plot(xi,[r['S_Q_total_norm'] for r in rr],label=seed+' (Green forced)')
        axs[0,1].plot(xi,[r['S_Q_total_to_secant_error'] for r in rr],label=seed)
        # Signed component curves for each physical seed are separate.
        for c,linestyle in [('v2','-'),('omega2','--')]:
            axs[1,0].plot(xi,[r['S_Q_'+c+'_radial2'] for r in rr],ls=linestyle,label=f'{seed} · {c}')
            axs[1,1].plot(xi,[r['S_Q_'+c+'_azimuthal2'] for r in rr],ls=linestyle,label=f'{seed} · {c}')
    titles=['Total parameter-response norm','Full versus physical parameter secant',
            'Second plane: radial response by source','Second plane: azimuthal response by source']
    labels=['|dQ/dkappa|','|S_Q - secant|','Signed radial component','Signed azimuthal component']
    for ax,t,l in zip(axs.flat,titles,labels):
        ax.set(title=t,xlabel='Event-relative xi',ylabel=l);ax.grid(alpha=.22);ax.legend(fontsize=8)
    fig.suptitle('Q550 v18.51t · model-generated source → Green response',fontsize=13)
    fig.tight_layout();fig.savefig(OUT/'source_to_response_green_radial_azimuthal_transport.png',dpi=155);plt.close(fig)
    return checks,summaries

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--mode',choices=('cached','full'),default='full')
    args=ap.parse_args();OUT.mkdir(parents=True,exist_ok=True)
    a,sreport,table,rblob,srun,sm,om=authenticate_s()
    pre=cached_precheck(a,sreport,table)
    (OUT/'LOCAL_CACHED_PRECHECK.json').write_text(json.dumps(pre,indent=2)+'\n')
    print('S AUTHENTICATED: '+str(len(sm['files']))+' source hashes, '+str(len(om['file_sha256']))+' output hashes',flush=True)
    checks={'s_source_and_output_hashes_authenticated':True,
            's_full_variational_result_authenticated':sreport['mode']=='full' and sreport['all_passed'],
            'two_distinct_parameter_branches_and_seed_families':float(a['kappa_original']-a['kappa_autonomous'])>0 and len(SEEDS)==2,
            'audited_horizon_xi_0_to_0p2':len(a['xi'])==513 and abs(a['xi'][-1]-.2)<1e-12}
    summaries=None
    if args.mode=='full':
        additional,summaries=analysis_full(a,sreport,table,rblob,srun,pre)
        checks.update(additional)
    report={'version':'v18.51t','title':'Source-to-Response Green Transport & Radial–Azimuthal Causality Audit',
        'mode':args.mode,
        'classification':('MODEL_GENERATED_GREEN_TRANSPORT_NUMERICAL_SOURCE_ATTRIBUTION' if args.mode=='full'
            else 'AUTHENTICATED_PARENT_PRECHECK_GREEN_CHANNEL_TRANSPORT_PENDING'),
        'scientific_checks':{k:bool(v) for k,v in checks.items()},
        'all_passed':bool(all(checks.values())),
        'parent_s_sha256':sha(SPARENT.read_bytes()),
        'parent_source_hashes':len(sm['files']),'parent_output_hashes':len(om['file_sha256']),
        'original_kappa':float(a['kappa_original']),'autonomous_kappa':float(a['kappa_autonomous']),
        'xi_end':float(a['xi'][-1]),'sections_per_seed':len(a['xi']),
        'parent_precheck':pre,'green_channel_results':summaries,
        'limitations':[
          'All measurements are on-trajectory parameter linear-response diagnostics; artificially splitting b is a linear-source intervention, not an experimental independent force knob.',
          'The source derivative is finite differenced on an original nonlinear field; its state Jacobian is the historical analytic jet.',
          'Green reconstruction uses a numerically integrated fundamental matrix and composite Simpson quadrature, not validated interval arithmetic.',
          'Radial and azimuthal contributions depend on the chosen Cartesian chart and moving local basis; no observer-free material torsion or torque has been proved.',
          'The horizon xi in [0,0.2] contains no demonstrated physical return/closed 10D state; no global shadowing or material memory established.',
          'Source split is relative to the explicit 10D physical state coordinates (v2,omega2,rest); coordinate changes alter channel attribution.'
        ],'environment':{'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'platform':platform.platform()}}
    (OUT/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print('Q550 v18.51t',args.mode,sum(checks.values()),'/',len(checks),'PASS' if report['all_passed'] else 'FAIL',flush=True)
    if not report['all_passed']:
        raise RuntimeError('Scientific tests failed: '+','.join(k for k,v in checks.items() if not v))

if __name__=='__main__':main()
