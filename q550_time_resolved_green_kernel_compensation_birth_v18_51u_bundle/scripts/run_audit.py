#!/usr/bin/env python3
"""Q550 v18.51u: time-resolved Green kernel and compensation-birth audit.

Uses the complete authenticated v18.51t server sendback without reintegrating
its expensive original 10D nonlinear model. Reintegrates only its *stored*
model-generated 10x10 tangent generator J along two midpoint trajectories.

The triangular K(observation xi, source u) is conditional on the original
explicit 10D coordinates, the v2/omega2/rest source split, and the moving
second-plane radial/azimuthal projection at observation time xi.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,platform,sys,zipfile
from pathlib import Path
import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'
TROOT='q550_source_to_response_green_transport_radial_azimuthal_causality_v18_51t_bundle'
SROOT='q550_parameter_sensitivity_transport_transverse_posture_attribution_v18_51s_bundle'
RROOT='q550_long_horizon_paired_kappa0_transverse_transport_return_section_search_v18_51r_bundle'
TPARENT=ROOT/'inputs'/(TROOT+'_server_sendback.zip')
SEEDS=('original','autonomous')
CHANNELS=('v2','omega2','other','total')
PROJECTIONS=('radial','azimuthal')


def digest_bytes(data):return hashlib.sha256(data).hexdigest()
def safe_ratio(a,b,floor=1e-10):
    return float(np.linalg.norm(a-b)/max(floor,np.linalg.norm(a),np.linalg.norm(b)))
def write_csv(path,rows):
    if not rows:raise ValueError('No records: '+str(path))
    with path.open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)

def authenticate_parent():
    """Authenticate both T file manifests and nested S/R provenance; no solvers."""
    with zipfile.ZipFile(TPARENT) as tz:
        bad=tz.testzip()
        if bad:raise RuntimeError('Damaged T archive at '+bad)
        src=json.loads(tz.read(TROOT+'/SOURCE_MANIFEST.json'))
        out=json.loads(tz.read(TROOT+'/SHA256SUMS_OUTPUT.json'))
        assert src['root']==out['root']==TROOT
        assert src['version']==out['version']=='v18.51t'
        for manifest,key in ((src,'files'),(out,'file_sha256')):
            for relative,expected in manifest[key].items():
                actual=digest_bytes(tz.read(TROOT+'/'+relative))
                if actual!=expected:raise RuntimeError('T hash mismatch: '+relative)
        report=json.loads(tz.read(TROOT+'/results/report.json'))
        assert report['version']=='v18.51t' and report['mode']=='full'
        assert report['all_passed'] and all(report['scientific_checks'].values())
        raw_data=tz.read(TROOT+'/results/green_channel_transport_states.npz')
        nested_s=tz.read(TROOT+'/inputs/'+SROOT+'_server_sendback.zip')
    with np.load(io.BytesIO(raw_data),allow_pickle=False) as obj:
        a={k:obj[k].copy() for k in obj.files}
    with zipfile.ZipFile(io.BytesIO(nested_s)) as sz:
        assert sz.testzip() is None,'Invalid nested S archive'
        sr=json.loads(sz.read(SROOT+'/results/report.json'))
        assert sr['version']=='v18.51s' and sr['mode']=='full' and sr['all_passed']
        nested_r=sz.read(SROOT+'/inputs/'+RROOT+'_server_sendback.zip')
    with zipfile.ZipFile(io.BytesIO(nested_r)) as rz:
        assert rz.testzip() is None,'Invalid nested R archive'
        rr=json.loads(rz.read(RROOT+'/results/report.json'))
        assert rr['version']=='v18.51r' and all(rr['checks'].values())
        rb=rz.read(RROOT+'/results/extended_paired_kappa_nonlinear_states.npz')
    with np.load(io.BytesIO(rb),allow_pickle=False) as obj:
        midpoint={seed:{'Y':obj[f'{seed}_midpoint_Y'].copy(),
                        'Q':obj[f'{seed}_midpoint_Q'].copy()} for seed in SEEDS}
        previous_xi=obj['xi'].copy()
    xi=a['xi']
    assert xi.shape==(513,) and np.max(np.abs(xi-previous_xi))<1e-14
    assert np.allclose(xi,np.linspace(0.,.2,513),rtol=0,atol=1e-14)
    for seed in SEEDS:
        assert a[seed+'_source_J'].shape==(257,10,10)
        assert a[seed+'_source_b'].shape==(257,10)
        assert a[seed+'_green_Y'].shape==(513,10,4)
        assert a[seed+'_forced_Y'].shape==(513,10,4)
        assert a[seed+'_radial2'].shape==(513,4)
        assert a[seed+'_azimuthal2'].shape==(513,4)
        assert midpoint[seed]['Y'].shape==(513,10)
        assert midpoint[seed]['Q'].shape==(513,5)
        assert np.all(np.isfinite(a[seed+'_source_J']))
        assert np.all(np.isfinite(a[seed+'_source_b']))
        assert np.max(np.abs(a[seed+'_source_xi']-xi[::2]))<1e-14
    return a,midpoint,report,{'parent_sha256':digest_bytes(TPARENT.read_bytes()),
           'source_hashes':len(src['files']),'output_hashes':len(out['file_sha256']),
           't_report_classification':report['classification']}


def make_B(b):
    """Four measured columns: v2, omega2, other, total (10D)."""
    v=np.zeros_like(b);v[5]=b[5]
    w=np.zeros_like(b);w[7]=b[7]
    other=b-v-w
    return np.column_stack((v,w,other,b))


def fundamental(xi,source_xi,Js,rtol=3e-11,atol=3e-13):
    Jx=CubicSpline(source_xi,Js,axis=0)
    def rhs(t,z):return (Jx(t)@z.reshape(10,10)).reshape(100)
    sol=solve_ivp(rhs,(float(xi[0]),float(xi[-1])),np.eye(10).reshape(100),
                  method='DOP853',t_eval=xi,rtol=rtol,atol=atol,max_step=.001)
    if not sol.success:raise RuntimeError('Fundamental IVP: '+sol.message)
    Phi=sol.y.T.reshape(len(xi),10,10)
    cond=np.linalg.cond(Phi)
    if not np.all(np.isfinite(cond)) or np.max(cond)>=1e11:
        raise RuntimeError('Fundamental representation ill-conditioned')
    return Phi,Jx,int(sol.nfev),float(np.max(cond))


def observation_covectors(Y,Q):
    """Rows of the 10D -> moving physical second-plane radial/azimuthal map.

    Construct by pushing forward the actual nonlinear Cartesian immersion:
    Q3=rho2*cos(theta2), Q4=rho2*sin(theta2). The radial covector
    is rho2*d(logrho2), the azimuthal covector rho2*d(theta2).
    """
    n=len(Y);L=np.zeros((n,2,10))
    for i,(y,q) in enumerate(zip(Y,Q)):
        th=float(y[6]);er=np.array([np.cos(th),np.sin(th)])
        ea=np.array([-np.sin(th),np.cos(th)])
        D=np.zeros((2,10))
        D[:,4]=q[2:4];D[:,6]=np.array([-q[3],q[2]])
        L[i,0]=er@D;L[i,1]=ea@D
    return L


def kernel_from_phi(xi,source_xi,bs,Phi,L):
    """K[i,u,p,c]= L_p(i) Phi(i) Phi(u)^(-1) b_c(u), 0 for u>i."""
    Bx=CubicSpline(source_xi,bs,axis=0)
    B=np.asarray([make_B(Bx(float(t))) for t in xi])
    # Interaction picture is more stable than materializing every inverse.
    H=np.linalg.solve(Phi,B)
    rows=np.einsum('nia,nab->nib',L,Phi)
    n=len(xi)
    K=np.zeros((n,n,2,4),dtype=np.float64)
    for i in range(n):
        K[i,:i+1,:,:]=np.einsum('pk,ukc->upc',rows[i],H[:i+1])
    return K,H,B


def simpson_weights(n,h):
    """Positive composite Simpson weights for 0..n; odd uses final 3/8."""
    w=np.zeros(n+1)
    if n==0:return w
    if n==1:return np.array([.5*h,.5*h])
    even_n=n if n%2==0 else n-3
    if even_n>=2:
        w[:even_n+1]+=h/3
        w[1:even_n:2]+=3*h/3
        w[2:even_n:2]+=h/3
    if n%2==1:
        w[n-3:n+1]+=3*h/8*np.array([1.,3.,3.,1.])
    return w


def integrated_kernel(K,xi):
    n=len(xi);h=float(xi[1]-xi[0]);out=np.zeros((n,2,4))
    for i in range(1,n):
        w=simpson_weights(i,h)
        out[i]=np.einsum('u,upc->pc',w,K[i,:i+1])
    return out


def signed_birth(xi,integrated,endpoint_kernel):
    """Operational onsets and history weights, not fundamental birth instants."""
    rv=integrated[:,0,0];rw=integrated[:,0,1];allr=integrated[:,0,3]
    summed=np.sum(np.abs(integrated[:,0,:3]),axis=1)
    frac=np.full(len(xi),np.nan)
    valid=summed>1e-14
    frac[valid]=1-np.abs(allr[valid])/summed[valid]
    def first(mask):
        jj=np.flatnonzero(mask)
        return float(xi[jj[0]]) if len(jj) else None
    onset={str(th):first(np.isfinite(frac)&(frac>=th)) for th in (.01,.05,.10)}
    meaningful=1e-5*max(np.max(np.abs(rw)),1e-14)
    onset['radial_opposition_above_1e-5_of_peak_omega2']=first((rv>0)&(rw< -meaningful))
    n=len(xi)-1;w=simpson_weights(n,float(xi[1]-xi[0]))
    ages={}
    for p,pr in enumerate(PROJECTIONS):
        for c,ch in enumerate(CHANNELS):
            weighted=w*np.abs(endpoint_kernel[:,p,c]);tot=float(weighted.sum())
            if tot<=1e-24:
                ages[pr+'_'+ch]={'median_source_u':None,'mean_source_lag':None,
                         'first_half_abs_mass_fraction':None,'positive_signed_mass':0.,'negative_signed_mass':0.}
                continue
            cdf=np.cumsum(weighted)
            medidx=int(np.searchsorted(cdf,.5*tot,side='left'))
            sig=w*endpoint_kernel[:,p,c]
            ages[pr+'_'+ch]={
              'median_source_u':float(xi[medidx]),
              'mean_source_lag':float(np.dot(weighted,xi[-1]-xi)/tot),
              'first_half_abs_mass_fraction':float(weighted[xi<=.1].sum()/tot),
              'positive_signed_mass':float(np.maximum(sig,0).sum()),
              'negative_signed_mass':float(np.minimum(sig,0).sum())}
    return frac,onset,ages


def window_records(xi,K,seed):
    """Fixed 0.025 source windows; exact summation of 64-step Simpson panels."""
    rows=[];h=float(xi[1]-xi[0]);span=64
    for obs in (128,256,384,512):
        for start in range(0,obs,span):
            end=start+span;w=simpson_weights(span,h)
            contrib=np.einsum('u,upc->pc',w,K[obs,start:end+1])
            for p,projection in enumerate(PROJECTIONS):
                for c,channel in enumerate(CHANNELS):
                    rows.append({'seed_family':seed,'xi_observation':float(xi[obs]),
                        'u_start':float(xi[start]),'u_end':float(xi[end]),
                        'projection':projection,'source_channel':channel,
                        'signed_window_response':float(contrib[p,c]),
                        'abs_signed_window_response':float(abs(contrib[p,c]))})
    return rows


def restart_impulse_checks(xi,Jx,source_xi,bs,L,K,restart_indices=(0,128,256,384)):
    """Independent short tangent IVPs from distinct source instants."""
    bx=CubicSpline(source_xi,bs,axis=0)
    errors=[]
    for j in restart_indices:
        t0=float(xi[j]);t1=float(xi[-1]);y0=make_B(bx(t0))[:,:3]
        sol=solve_ivp(lambda t,z:(Jx(t)@z.reshape(10,3)).reshape(30),
            (t0,t1),y0.ravel(),method='DOP853',rtol=2e-11,
            atol=2e-13,max_step=.001)
        if not sol.success:raise RuntimeError('Impulse restart IVP '+sol.message)
        z=sol.y[:,-1].reshape(10,3)
        actual=L[-1]@z
        expected=K[-1,j,:,:3]
        assert actual.shape==expected.shape==(2,3), (actual.shape,expected.shape)
        errors.append(float(np.max(np.abs(actual-expected))))
    return float(max(errors)),len(errors)


def analysis(seed,xi,a,mid):
    source_xi=a[seed+'_source_xi'];Js=a[seed+'_source_J'];bs=a[seed+'_source_b']
    Y=mid[seed]['Y'];Q=mid[seed]['Q'];L=observation_covectors(Y,Q)
    Phi,Jx,nfev,cond=fundamental(xi,source_xi,Js)
    K,H,B=kernel_from_phi(xi,source_xi,bs,Phi,L)
    integral=integrated_kernel(K,xi)
    parent=np.stack([a[seed+'_radial2'],a[seed+'_azimuthal2']],axis=1)
    parent_err=float(np.max(np.abs(integral-parent)))
    sum_err=float(np.max(np.abs(integral[:,:,:3].sum(axis=2)-integral[:,:,3])))
    future_err=float(np.max(np.abs(K[np.triu_indices(len(xi),k=1)])))
    diag_mag=float(np.max(np.abs(np.diagonal(K,axis1=0,axis2=1))))
    restart_err,restarts=restart_impulse_checks(xi,Jx,source_xi,bs,L,K)
    # Coarser source sampling for an independently integrated Phi.
    Phi_c,_,_,_=fundamental(xi,source_xi[::2],Js[::2],rtol=3e-11,atol=3e-13)
    Kc,_,_=kernel_from_phi(xi,source_xi[::2],bs[::2],Phi_c,L)
    coarse_int=integrated_kernel(Kc,xi)
    coarse_err=float(np.max(np.abs(coarse_int[-1]-integral[-1])))
    # Coarsen the quadrature nodes only (same Phi); this isolates quadrature.
    h=float(xi[1]-xi[0])
    half_weights=simpson_weights(256,2*h)
    half_integral=np.einsum('u,upc->pc',half_weights,K[-1,::2])
    quad_err=float(np.max(np.abs(half_integral-integral[-1])))
    # Source-impulse kernel and moving basis are on the observation orbit.
    frac,onset,ages=signed_birth(xi,integral,K[-1])
    out={'end_radial':integral[-1,0].tolist(),
         'end_azimuthal':integral[-1,1].tolist(),
         'end_cancellation_fraction_radial':float(frac[-1]),
         'compensation_onset_thresholds':onset,
         'source_age_at_final_observation':ages,
         'max_parent_channel_response_absolute_error':parent_err,
         'max_channel_additivity_error':sum_err,
         'max_future_source_kernel_abs':future_err,
         'max_same_time_positional_kernel_abs':diag_mag,
         'max_independent_impulse_restart_absolute_error':restart_err,
         'independent_restart_ivps':restarts,
         'source_grid_endpoint_max_abs_difference':coarse_err,
         'quadrature_endpoint_max_abs_difference':quad_err,
         'maximum_fundamental_condition_number':cond,
         'fundamental_ode_evals':nfev}
    return out,K,integral,frac,window_records(xi,K,seed)


def make_figures(xi,out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for seed,d in out.items():
        K=d['kernel'];I=d['integral'];frac=d['fraction']
        fig,axes=plt.subplots(2,2,figsize=(12,9),sharex=True,sharey=True)
        for ax,(p,c) in zip(axes.flat,[(0,0),(0,1),(1,0),(1,1)]):
            field=np.ma.masked_where(np.triu(np.ones((len(xi),len(xi)),bool),k=1),K[:,:,p,c])
            im=ax.pcolormesh(xi,xi,field,shading='auto',cmap='coolwarm')
            ax.plot(xi,xi,ls='--',linewidth=.7,c='black',alpha=.55)
            ax.set(title=f'{PROJECTIONS[p]} from {CHANNELS[c]}',
                   xlabel='Source time u',ylabel='Observation xi')
            fig.colorbar(im,ax=ax,shrink=.78,label='Signed Green kernel')
        fig.suptitle(f'Q550 v18.51u · time-resolved Green kernel · {seed}',fontsize=13)
        fig.tight_layout();fig.savefig(OUT/f'green_kernel_causal_map_{seed}.png',dpi=145);plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(12.8,5.2))
    for seed,d in out.items():
        I=d['integral'];frac=d['fraction'];axs[0].plot(xi,I[:,0,0],label=seed+' v2')
        axs[0].plot(xi,I[:,0,1],ls='--',label=seed+' omega2')
        axs[0].plot(xi,I[:,0,3],lw=1.3,label=seed+' total')
        axs[1].plot(xi,100*frac,label=seed)
    axs[0].set(title='Signed radial accumulated response',xlabel='Observation xi',ylabel='dQ_r / dκ₀')
    axs[1].set(title='Radial compensation across observations',xlabel='Observation xi',ylabel='Cancellation (%)')
    for ax in axs:
        ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.tight_layout();fig.savefig(OUT/'compensation_birth_and_channel_history.png',dpi=165);plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=('preflight','full'),default='full')
    args=parser.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    a,mid,treport,provenance=authenticate_parent()
    xi=a['xi']
    base_checks={
        'authenticated_complete_parent_t':True,
        'parent_t_numerical_audit_passed':treport['all_passed'],
        'two_seed_midpoint_physical_states':all(mid[s]['Y'].shape==(513,10) for s in SEEDS),
        'matched_stored_model_jacobian_and_sources':all(a[s+'_source_J'].shape==(257,10,10) for s in SEEDS),
        'causal_uniform_513_section_grid':len(xi)==513 and abs(xi[-1]-.2)<1e-14}
    (OUT/'LOCAL_PREFLIGHT.json').write_text(json.dumps({'version':'v18.51u','mode':'preflight',
        'checks':{k:bool(v) for k,v in base_checks.items()},'provenance':provenance},indent=2)+'\n')
    if args.mode=='preflight':
        print('Q550 v18.51u PREFLIGHT',sum(base_checks.values()),'/',len(base_checks),'PASS; no kernel calculated')
        return
    results={};arrays={'xi':xi};window=[];ledger=[]
    for seed in SEEDS:
        print('Time-resolved Green kernel:',seed,flush=True)
        summary,K,I,fraction,windows=analysis(seed,xi,a,mid)
        results[seed]={'summary':summary,'kernel':K,'integral':I,'fraction':fraction}
        arrays[seed+'_kernel']=K
        arrays[seed+'_integrated_radial_azimuthal']=I
        arrays[seed+'_cancellation_fraction']=fraction
        window.extend(windows)
        for i,t in enumerate(xi):
            ledger.append({'seed_family':seed,'xi_observation':float(t),
                'radial_v2':float(I[i,0,0]),'radial_omega2':float(I[i,0,1]),
                'radial_other':float(I[i,0,2]),'radial_total':float(I[i,0,3]),
                'azimuthal_v2':float(I[i,1,0]),'azimuthal_omega2':float(I[i,1,1]),
                'azimuthal_other':float(I[i,1,2]),'azimuthal_total':float(I[i,1,3]),
                'radial_signed_cancellation_fraction':float(fraction[i]) if np.isfinite(fraction[i]) else '',
                'signed_radial_recomposition_error':float(abs(I[i,0,:3].sum()-I[i,0,3]))})
        print('  terminal radial:',summary['end_radial'],flush=True)
        print('  birth markers:',summary['compensation_onset_thresholds'],flush=True)
    checks=dict(base_checks)
    checks.update({
      'causality_future_support_zero':all(r['summary']['max_future_source_kernel_abs']<1e-15 for r in results.values()),
      'source_channel_kernel_additivity':all(r['summary']['max_channel_additivity_error']<5e-11 for r in results.values()),
      'two_time_green_recovers_parent_t_response':all(r['summary']['max_parent_channel_response_absolute_error']<2e-7 for r in results.values()),
      'independent_impulse_restart_transport':all(r['summary']['max_independent_impulse_restart_absolute_error']<2e-7 for r in results.values()),
      'source_sampling_grid_refinement':all(r['summary']['source_grid_endpoint_max_abs_difference']<2e-5 for r in results.values()),
      'simpson_quadrature_refinement':all(r['summary']['quadrature_endpoint_max_abs_difference']<2e-5 for r in results.values()),
      'observed_source_lag_finite':all(r['summary']['source_age_at_final_observation']['radial_omega2']['mean_source_lag'] is not None for r in results.values()),
      'compensation_birth_marked_operationally':all(r['summary']['compensation_onset_thresholds']['0.01'] is not None for r in results.values()),
      'stable_fundamental_representation':all(r['summary']['maximum_fundamental_condition_number']<1e11 for r in results.values()),
      'zero_simultaneous_position_source_direct_feedthrough':all(r['summary']['max_same_time_positional_kernel_abs']<1e-9 for r in results.values()),
      'no_global_memory_or_material_return_claim':True})
    write_csv(OUT/'compensation_birth_section_ledger.csv',ledger)
    write_csv(OUT/'causal_source_window_contributions.csv',window)
    np.savez_compressed(OUT/'green_causal_two_time_kernels.npz',**arrays)
    make_figures(xi,results)
    public={seed:rec['summary'] for seed,rec in results.items()}
    report={'version':'v18.51u','title':'Time-Resolved Green Kernel & Compensation-Birth Audit',
        'mode':'full', 'classification': 'TWO_TIME_CAUSAL_GREEN_NUMERICAL_DIAGNOSTIC_NOT_MATERIAL_MEMORY_CERTIFICATE',
        'scientific_checks':{k:bool(v) for k,v in checks.items()},
        'all_passed':bool(all(checks.values())),'parent_t':provenance,
        'source_windows_per_seed':len(window)//2,'section_records':len(ledger),
        'seeds':public,'source_channels':list(CHANNELS),
        'observable_projections':list(PROJECTIONS),
        'limitations':[
        'The causal Green response is a linear on-trajectory first variation at the historical midpoint κ₀, not a physical independent actuator intervention.',
        'Source channels v2, omega2, and other and moving radial/azimuthal outputs are conditional on the declared 10D coordinate chart.',
        'The nonzero compensation onset is a threshold crossing of an accumulated response; it is not a microscopic material birth event.',
        'Past-window signs and source ages are response-weighted mathematical attributions, not intrinsic energy flow, entropy production, memory or delay clocks.',
        'All direct nonlinear field calculations are inherited from the verified t/s/r archives. This audit only repropagates their saved analytic Jacobians.',
        'Fundamental and impulse IVPs plus composite Simpson quadrature are high-precision numerical, not interval-certified.',
        'No return of the complete 10D material state or global two-kappa shadowing is established outside xi in [0,0.2].'],
        'environment':{'python':sys.version,'numpy':np.__version__,
             'scipy':scipy.__version__,'platform':platform.platform()}}
    (OUT/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    lines=['# Q550 v18.51u — server scientific report','',
        '**Status:** '+('PASS' if report['all_passed'] else 'FAIL'),
        '**Classification:** '+report['classification'],'',
        'This document is generated by the full server run, not pre-populated with predictions.','',
        '## Channel-resolved observations','']
    for seed,summary in public.items():
        lines.extend(['### '+seed,'',
          '- Terminal radial channel responses (v2, omega2, other, total): '+str(summary['end_radial']),
          '- Terminal azimuthal channel responses (v2, omega2, other, total): '+str(summary['end_azimuthal']),
          '- Terminal radial cancellation fraction: '+str(summary['end_cancellation_fraction_radial']),
          '- Numerical compensation onset markers: '+json.dumps(summary['compensation_onset_thresholds']),
          '- Source-age summary at final observation: '+json.dumps(summary['source_age_at_final_observation']),
          '- Maximum parent-response difference: '+str(summary['max_parent_channel_response_absolute_error']),
          '- Maximum independent impulse-restart difference: '+str(summary['max_independent_impulse_restart_absolute_error']),
          '- Source-grid and quadrature endpoint errors: '+str((summary['source_grid_endpoint_max_abs_difference'],summary['quadrature_endpoint_max_abs_difference'])),''])
    lines.extend(['## Interpretation limits','']+['- '+t for t in report['limitations']])
    (OUT/'SCIENTIFIC_REPORT_SERVER.md').write_text('\n'.join(lines)+'\n')
    print('Q550 v18.51u FULL',sum(checks.values()),'/',len(checks),
         'PASS' if report['all_passed'] else 'FAIL',flush=True)
    if not report['all_passed']:
        raise RuntimeError('Scientific checks failed: '+str([k for k,v in checks.items() if not v]))

if __name__=='__main__':main()
