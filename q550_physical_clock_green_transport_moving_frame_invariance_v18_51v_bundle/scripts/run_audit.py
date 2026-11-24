#!/usr/bin/env python3
"""Q550 v18.51v — Physical-Clock Green Transport & Moving-Frame Invariance.

No nonlinear ODE is reintegrated in this audit. The numerical physical clock,
nonlinear-state midpoint orbits, and two-time Green kernels are inherited from
SHA-256-checked v18.51r/t/u server sendbacks. The physical clock is a
trajectory-dependent map tau(xi; kappa); do not conflate a frozen-clock
relabeling of a Green kernel with kappa sensitivity at fixed physical tau.
"""
from __future__ import annotations
import argparse, csv, hashlib, io, json, platform, sys, zipfile
from pathlib import Path
import numpy as np
from scipy.integrate import simpson
from scipy.interpolate import CubicSpline
import scipy

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results'
UROOT = 'q550_time_resolved_green_kernel_compensation_birth_v18_51u_bundle'
TROOT = 'q550_source_to_response_green_transport_radial_azimuthal_causality_v18_51t_bundle'
SROOT = 'q550_parameter_sensitivity_transport_transverse_posture_attribution_v18_51s_bundle'
RROOT = 'q550_long_horizon_paired_kappa0_transverse_transport_return_section_search_v18_51r_bundle'
UFILE = ROOT/'inputs'/(UROOT+'_server_sendback.zip')
TCACHE = ROOT/'inputs'/'authenticated_green_channel_transport_states_v18_51t.npz'
RCACHE = ROOT/'inputs'/'authenticated_extended_paired_kappa_nonlinear_states_v18_51r.npz'
SEEDS = ('original','autonomous')
CHANNELS = ('v2','omega2','other','total')
PROJECTIONS = ('radial','azimuthal')
G = np.array([[0.,-1.],[1.,0.]])


def digest(data): return hashlib.sha256(data).hexdigest()

def require(condition, message):
    if not bool(condition): raise RuntimeError(message)

def maxabs(a): return float(np.max(np.abs(a)))

def rotation(angle):
    """Physical plane rotation, shape (...,2,2)."""
    ca,sa=np.cos(angle),np.sin(angle)
    out=np.empty(np.shape(ca)+(2,2),dtype=np.float64)
    out[...,0,0]=ca;out[...,0,1]=-sa
    out[...,1,0]=sa;out[...,1,1]=ca
    return out

def stream_checksums(z,base,manifest,kind):
    for rel,sha in manifest[kind].items():
        require(digest(z.read(base+'/'+rel))==sha,
                f'{base}: SHA256 mismatch: {rel}')

def get_zipped_json(z,root,name): return json.loads(z.read(root+'/'+name))

def authenticate_and_read(full=True):
    """Verify nested executable genealogy and compare convenient cache copies.

    The full parent U is the single archival ancestor. The small cached T and
    R arrays are *byte-identically* matched against the arrays nested inside U,
    not independent untracked reconstructions of model data.
    """
    require(UFILE.exists(), 'Missing full U server sendback under inputs/')
    require(TCACHE.exists() and RCACHE.exists(),'Missing authenticated T/R cached states')
    lineage={}
    with zipfile.ZipFile(UFILE) as uz:
        require(uz.testzip() is None,'U ZIP CRC failed')
        us=get_zipped_json(uz,UROOT,'SOURCE_MANIFEST.json')
        uo=get_zipped_json(uz,UROOT,'SHA256SUMS_OUTPUT.json')
        ur=get_zipped_json(uz,UROOT,'results/report.json')
        require(ur.get('version')=='v18.51u' and ur.get('all_passed') and
                all(ur['scientific_checks'].values()),'U scientific result not passed')
        stream_checksums(uz,UROOT,us,'files')
        stream_checksums(uz,UROOT,uo,'file_sha256')
        udata=uz.read(UROOT+'/results/green_causal_two_time_kernels.npz') if full else None
        tb=uz.read(UROOT+'/inputs/'+TROOT+'_server_sendback.zip')
        lineage['u']={'sha256':digest(UFILE.read_bytes()),'source_entries':len(us['files']),
                      'output_entries':len(uo['file_sha256'])}
    with zipfile.ZipFile(io.BytesIO(tb)) as tz:
        require(tz.testzip() is None,'Nested T archive CRC failed')
        ts=get_zipped_json(tz,TROOT,'SOURCE_MANIFEST.json')
        tr=get_zipped_json(tz,TROOT,'results/report.json')
        require(tr.get('version')=='v18.51t' and tr.get('all_passed'), 'T not passed')
        tdata=tz.read(TROOT+'/results/green_channel_transport_states.npz')
        require(digest(tdata)==digest(TCACHE.read_bytes()),'T cache differs from nested authenticated source')
        stream_checksums(tz,TROOT,ts,'files')
        sb=tz.read(TROOT+'/inputs/'+SROOT+'_server_sendback.zip')
        lineage['t']={'sha256':digest(tb),'source_entries':len(ts['files'])}
    with zipfile.ZipFile(io.BytesIO(sb)) as sz:
        require(sz.testzip() is None,'Nested S archive CRC failed')
        ss=get_zipped_json(sz,SROOT,'SOURCE_MANIFEST.json')
        sr=get_zipped_json(sz,SROOT,'results/report.json')
        require(sr.get('version')=='v18.51s' and sr.get('all_passed'),'S not passed')
        stream_checksums(sz,SROOT,ss,'files')
        rb=sz.read(SROOT+'/inputs/'+RROOT+'_server_sendback.zip')
        lineage['s']={'sha256':digest(sb),'source_entries':len(ss['files'])}
    with zipfile.ZipFile(io.BytesIO(rb)) as rz:
        require(rz.testzip() is None,'Nested R archive CRC failed')
        rs=get_zipped_json(rz,RROOT,'SOURCE_MANIFEST.json')
        rr=get_zipped_json(rz,RROOT,'results/report.json')
        require(rr.get('version')=='v18.51r' and all(rr['checks'].values()),'R not passed')
        rdata=rz.read(RROOT+'/results/extended_paired_kappa_nonlinear_states.npz')
        require(digest(rdata)==digest(RCACHE.read_bytes()),'R cache differs from nested authenticated source')
        stream_checksums(rz,RROOT,rs,'files')
        lineage['r']={'sha256':digest(rb),'source_entries':len(rs['files'])}
    if not full: return None,None,None,lineage
    with np.load(io.BytesIO(udata),allow_pickle=False) as z:
        U={k:z[k].copy() for k in z.files}
    with np.load(io.BytesIO(tdata),allow_pickle=False) as z:
        T={k:z[k].copy() for k in z.files}
    with np.load(io.BytesIO(rdata),allow_pickle=False) as z:
        R={k:z[k].copy() for k in z.files}
    return U,T,R,lineage


def simpson_weights(i,h):
    """Match U's parity-aware positive composite Simpson rule."""
    w=np.zeros(i+1)
    if i==0:return w
    if i==1:return np.array([.5*h,.5*h])
    even=i if i%2==0 else i-3
    if even>=2:
        w[:even+1]+=h/3
        w[1:even:2]+=h
        w[2:even:2]+=h/3
    if i%2:
        w[i-3:i+1]+=3*h/8*np.array([1.,3.,3.,1.])
    return w


def source_age(xi,tau,K,i,p,c):
    """Positive absolute-response age weights at one observation event."""
    weights=simpson_weights(i,float(xi[1]-xi[0]))*np.abs(K[i,:i+1,p,c])
    mass=float(weights.sum())
    if mass<1e-25:return {'mass':mass,'median_xi':None,'median_tau':None}
    u=np.arange(i+1)
    mid=int(np.searchsorted(np.cumsum(weights),mass*.5))
    return {'mass':mass,
            'median_xi':float(xi[mid]),'median_tau':float(tau[mid]),
            'mean_lag_xi':float(np.dot(weights,xi[i]-xi[:i+1])/mass),
            'mean_lag_tau':float(np.dot(weights,tau[i]-tau[:i+1])/mass),
            'first_half_xi_mass_fraction':float(weights[xi[:i+1]<=xi[i]/2].sum()/mass),
            'first_half_tau_mass_fraction':float(weights[tau[:i+1]<=tau[i]/2].sum()/mass)}


def inspect_seed(seed,U,T,R,xi,detail_arrays):
    """Reparameterize Green; compare fixed-xi and fixed-tau sensitivities."""
    n=len(xi)
    K=U[seed+'_kernel']
    I=U[seed+'_integrated_radial_azimuthal']
    Y=R[seed+'_midpoint_Y'];Q=R[seed+'_midpoint_Q']
    c=R[seed+'_midpoint_clock'];tau=R[seed+'_midpoint_tau']
    qprime=R[seed+'_midpoint_Qprime']
    dk=float(R['kappa_original']-R['kappa_autonomous'])
    require(K.shape==(n,n,2,4),f'{seed}: unexpected K shape')
    require(I.shape==(n,2,4),f'{seed}: unexpected I shape')
    require(np.all(c>0) and np.all(np.diff(tau)>0),f'{seed}: invalid physical clock')
    require(abs(tau[0])<1e-12,f'{seed}: physical clock does not start at zero')
    clock_deriv=CubicSpline(xi,tau).derivative()(xi)
    clock_defect=maxabs(clock_deriv-c)
    Ktau=K/c[None,:,None,None]
    density_defect=maxabs(Ktau*c[None,:,None,None]-K)
    # Independent nonuniform physical-time integration (all row sections).
    # The clock is state-dependent: this is specifically the *frozen-history*
    # tangent source-density reparameterization, not S_kappa at fixed tau.
    physical_integ=np.zeros_like(I)
    for j in range(1,n):
        physical_integ[j]=simpson(Ktau[j,:j+1],x=tau[:j+1],axis=0)
    integral_defect=maxabs(physical_integ-I)
    causal_defect=maxabs(Ktau[np.triu_indices(n,1)])

    # Real physical planar vectors in a fixed Cartesian laboratory frame.
    theta=Y[:,6];Rtheta=rotation(theta)
    K_lab=np.einsum('iab,ijbc->ijac',Rtheta,K)
    I_lab=np.einsum('iab,ibc->iac',Rtheta,I)
    # Two independent observer-family conventions. Neither is new physics.
    z=(xi-xi[0])/(xi[-1]-xi[0])
    observers={'rotating_sine':0.38*np.sin(2*np.pi*z)+.27*z,
               'rotating_crest':0.31*np.cos(3*np.pi*z)+.17*z+.08}
    frames={}
    for name,psi in observers.items():
        Rp=rotation(psi);Rminus=rotation(-psi)
        Iframe=np.einsum('iab,ibc->iac',Rminus,I)
        Kframe=np.einsum('iab,ijbc->ijac',Rminus,K)
        reconstructed=np.einsum('iab,ibc->iac',rotation(theta+psi),Iframe)
        reconstructedK=np.einsum('iab,ijbc->ijac',rotation(theta+psi),Kframe)
        projection_defect=maxabs(reconstructed-I_lab)
        kernel_frame_defect=maxabs(reconstructedK-K_lab)
        norm_defect=maxabs(np.linalg.norm(Iframe,axis=1)-np.linalg.norm(I,axis=1))
        frame_add=maxabs(Iframe[:,:,:3].sum(axis=2)-Iframe[:,:,3])
        # Connection-aware physical derivative. A moving-frame component
        # derivative alone is not invariant. Add (theta'+psi') G w.
        deriv_native=CubicSpline(xi,I,axis=0).derivative()(xi)
        deriv_frame=CubicSpline(xi,Iframe,axis=0).derivative()(xi)
        thprime=CubicSpline(xi,theta).derivative()(xi)
        psiprime=CubicSpline(xi,psi).derivative()(xi)
        deriv_lab_native=np.einsum('iab,ibc->iac',Rtheta,
           deriv_native+thprime[:,None,None]*np.einsum('ab,ibc->iac',G,I))
        Rtotal=rotation(theta+psi)
        deriv_lab_observer=np.einsum('iab,ibc->iac',Rtotal,
           deriv_frame+(thprime+psiprime)[:,None,None]*np.einsum('ab,ibc->iac',G,Iframe))
        connection_defect=maxabs(deriv_lab_observer-deriv_lab_native)
        frames[name]={'final_angle_relative_to_radial':float(psi[-1]),
             'max_physical_vector_reconstruction_error':projection_defect,
             'max_kernel_reconstruction_error':kernel_frame_defect,
             'max_response_norm_error':norm_defect,
             'max_channel_additivity_error':frame_add,
             'max_covariant_derivative_reconstruction_error':connection_defect,
             'end_rotated_radial_channel_components':[float(x) for x in Iframe[-1,0,:]],
             'end_rotated_azimuthal_channel_components':[float(x) for x in Iframe[-1,1,:]]}
        detail_arrays[seed+'_'+name+'_response']=Iframe
        if name=='rotating_sine':detail_arrays[seed+'_rotating_sine_kernel']=Kframe

    # Physical clock depends on the orbit: a kappa-sensitive history has
    # different tau(xi;kappa). The following correction is compulsory when
    # comparing the *two different kappa trajectories at equal tau*.
    Qplus=R[seed+'_original_Q'];Qminus=R[seed+'_autonomous_Q']
    tauplus=R[seed+'_original_tau'];tauminus=R[seed+'_autonomous_tau']
    dtaudk=(tauplus-tauminus)/dk
    Sxi=(Qplus-Qminus)/dk
    Sclock=-qprime*(dtaudk/c)[:,None]
    Stau_pred=Sxi+Sclock
    # Independent matched-physical-clock construction from the nonlinear
    # already-integrated different-kappa orbits (cubic interpolation in tau).
    Stau_direct=(CubicSpline(tauplus,Qplus,axis=0)(tau)-
                 CubicSpline(tauminus,Qminus,axis=0)(tau))/dk
    matched_tau_defect=maxabs(Stau_direct-Stau_pred)
    tmodel=T[seed+'_cartesian_Q'][:,:,3]
    model_vs_secant=maxabs(tmodel-Sxi)
    model_clock_corrected=tmodel+Sclock
    model_vs_fixedtau=maxabs(model_clock_corrected-Stau_direct)
    # Convert clock shift to the natural moving physical radial/azimuthal map.
    er=np.column_stack((np.cos(theta),np.sin(theta)))
    ea=np.column_stack((-np.sin(theta),np.cos(theta)))
    Sclock_plane=np.column_stack((np.einsum('ij,ij->i',er,Sclock[:,2:4]),
                                  np.einsum('ij,ij->i',ea,Sclock[:,2:4])))
    terminal_clock_norm=float(np.linalg.norm(Sclock[-1]))
    terminal_sensitivity_norm=float(np.linalg.norm(Sxi[-1]))
    # Cancellation is invariant under *time relabeling* of the same
    # physical output, NOT under arbitrary frame rotation.
    rad_frac=U[seed+'_cancellation_fraction']
    boundary=maxabs(physical_integ[:,:,3]-I[:,:,3])
    ages={}
    for p,pr in enumerate(PROJECTIONS):
        for j,ch in enumerate(CHANNELS):
            ages[pr+'_'+ch]=source_age(xi,tau,K,n-1,p,j)
    detail_arrays.update({seed+'_tau':tau,seed+'_clock':c,
          seed+'_kernel_physical_time_density':Ktau,
          seed+'_lab_frame_response':I_lab,
          seed+'_fixed_xi_cartesian_sensitivity':Sxi,
          seed+'_clock_shift_cartesian':Sclock,
          seed+'_fixed_tau_matched_cartesian_sensitivity':Stau_direct,
          seed+'_clock_shift_second_plane_radial_azimuthal':Sclock_plane})
    metrics={'physical_tau_final':float(tau[-1]),
       'clock_start':float(c[0]),'clock_final':float(c[-1]),
       'max_clock_identity_error':clock_defect,
       'max_frozen_history_green_density_error':density_defect,
       'max_causal_future_kernel_abs':causal_defect,
       'max_nonlinear_tau_quadrature_response_error':integral_defect,
       'max_physical_fixed_history_total_error':boundary,
       'max_tangent_vs_kappa_secant_error':model_vs_secant,
       'max_matched_tau_vs_clock_corrected_secant_error':matched_tau_defect,
       'max_matched_tau_vs_clock_corrected_model_error':model_vs_fixedtau,
       'terminal_fixed_xi_sensitivity_norm':terminal_sensitivity_norm,
       'terminal_clock_correction_norm':terminal_clock_norm,
       'terminal_clock_correction_ratio':terminal_clock_norm/terminal_sensitivity_norm,
       'terminal_clock_correction_second_plane':Sclock_plane[-1].tolist(),
       'terminal_fixed_tau_sensitivity_norm':float(np.linalg.norm(Stau_direct[-1])),
       'terminal_radial_cancel_fraction':float(rad_frac[-1]),
       'source_ages':ages,'observer_families':frames,
       'clock_interpretation': 'frozen-history Green density and actual fixed-tau kappa sensitivity are distinct'}
    # Small diagnostics rows: no huge matrices required to inspect results.
    rows=[]
    for i in range(n):
        rows.append({'seed_family':seed,'xi':float(xi[i]),'tau_physical':float(tau[i]),
           'clock_c':float(c[i]),
           'radial_v2_fixed_xi':float(I[i,0,0]),
           'radial_omega2_fixed_xi':float(I[i,0,1]),
           'azimuthal_v2_fixed_xi':float(I[i,1,0]),
           'azimuthal_omega2_fixed_xi':float(I[i,1,1]),
           'radial_total_fixed_xi':float(I[i,0,3]),
           'azimuthal_total_fixed_xi':float(I[i,1,3]),
           'radial_clock_correction':float(Sclock_plane[i,0]),
           'azimuthal_clock_correction':float(Sclock_plane[i,1]),
           'radial_total_fixed_tau_approx':float(I[i,0,3]+Sclock_plane[i,0]),
           'azimuthal_total_fixed_tau_approx':float(I[i,1,3]+Sclock_plane[i,1]),
           'cartesian_clock_correction_norm':float(np.linalg.norm(Sclock[i])),
           'fixed_tau_direct_vs_formula_norm':float(np.linalg.norm(Stau_direct[i]-Stau_pred[i]))})
    return metrics,rows


def render_figures(xi,U,report,arrays):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(1,2,figsize=(13,5.5))
    for seed in SEEDS:
        tau=arrays[seed+'_tau'];I=U[seed+'_integrated_radial_azimuthal']
        cc=arrays[seed+'_clock_shift_second_plane_radial_azimuthal']
        axs[0].plot(tau,I[:,0,3],label=seed+' fixed-xi radial')
        axs[0].plot(tau,I[:,0,3]+cc[:,0],ls='--',label=seed+' fixed-tau corrected radial')
        axs[1].plot(tau,cc[:,0],label=seed+' radial')
        axs[1].plot(tau,cc[:,1],ls='--',label=seed+' azimuthal')
    axs[0].set(xlabel='Physical reference clock tau',ylabel='Radial sensitivity',
        title='Fixed-history vs fixed-physical-time sensitivity')
    axs[1].set(xlabel='Physical reference clock tau',ylabel='Clock-shift correction',
        title='Parameter-dependent clock contribution')
    for ax in axs:
        ax.grid(alpha=.22);ax.legend(fontsize=7)
    fig.tight_layout();fig.savefig(OUT/'physical_clock_vs_fixed_tau_sensitivity.png',dpi=160);plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(13.5,5.5))
    for seed in SEEDS:
        I=U[seed+'_integrated_radial_azimuthal'];W=arrays[seed+'_rotating_sine_response'];tau=arrays[seed+'_tau']
        axs[0].plot(tau,I[:,0,3],label=seed+' natural radial')
        axs[0].plot(tau,W[:,0,3],ls='--',label=seed+' rotated radial')
        axs[1].plot(tau,np.linalg.norm(I[:,:,3],axis=1),label=seed+' natural norm')
        axs[1].plot(tau,np.linalg.norm(W[:,:,3],axis=1),ls='--',label=seed+' rotated norm')
    axs[0].set(xlabel='Physical reference clock tau',ylabel='Signed response',
         title='Frame-dependent radial component')
    axs[1].set(xlabel='Physical reference clock tau',ylabel='2D vector norm',
         title='Observer-frame invariant physical vector norm')
    for ax in axs:
        ax.grid(alpha=.22);ax.legend(fontsize=7)
    fig.tight_layout();fig.savefig(OUT/'moving_frame_components_vs_physical_vectors.png',dpi=160);plt.close(fig)
    # Plot the clock-converted triangular original seed causal kernel.
    seed='original';tau=arrays[seed+'_tau'];Ktau=arrays[seed+'_kernel_physical_time_density']
    mask=np.triu(np.ones((len(tau),len(tau)),bool),1)
    fig,axs=plt.subplots(1,2,figsize=(12.6,5.2),sharex=True,sharey=True)
    for ax,c,name in zip(axs,(0,1),('v2','omega2')):
        img=ax.pcolormesh(tau,tau,np.ma.masked_where(mask,Ktau[:,:,0,c]),
                          shading='auto',cmap='coolwarm')
        ax.plot(tau,tau,'k--',alpha=.5,lw=.8)
        ax.set(xlabel='Physical source clock tau(u)',ylabel='Physical observation clock tau(xi)',
               title='Radial response from '+name+' (per physical time)')
        fig.colorbar(img,ax=ax,shrink=.8)
    fig.suptitle('Q550 v18.51v — frozen-history physical-clock Green density')
    fig.tight_layout();fig.savefig(OUT/'physical_clock_causal_green_density.png',dpi=150);plt.close(fig)


def csv_write(path,rows):
    require(bool(rows),'No CSV rows')
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=('full','preflight'),default='full')
    args=parser.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    U,T,R,lineage=authenticate_and_read(full=args.mode=='full')
    if args.mode=='preflight':
        report={'version':'v18.51v','mode':'preflight','lineage':lineage,
                'checks':{'full_u_parent_source_crc_and_manifest':True,
                           'nested_t_s_r_genealogy_matched_to_cached_states':True}}
        (OUT/'LOCAL_PREFLIGHT.json').write_text(json.dumps(report,indent=2)+'\n')
        print('Q550 v18.51v provenance preflight PASS')
        return
    xi=U['xi'];require(xi.shape==(513,), 'Unexpected xi grid')
    require(maxabs(T['xi']-xi)<1e-14 and maxabs(R['xi']-xi)<1e-14,
            'Nonidentical physical sections')
    require(float(R['kappa_original'])>float(R['kappa_autonomous']),
            'Historical kappas not ordered')
    results={};ledger=[];arrays={'xi':xi}
    for seed in SEEDS:
        print('v18.51v physical-clock and moving-frame:',seed,flush=True)
        metrics,rows=inspect_seed(seed,U,T,R,xi,arrays)
        results[seed]=metrics;ledger.extend(rows)
        print('   clock range:',metrics['clock_start'],metrics['clock_final'],flush=True)
        print('   fixed-tau correction ratio:',metrics['terminal_clock_correction_ratio'],flush=True)
    checks={
      'authenticated_nested_model_genealogy':True,
      'historical_midpoint_and_distinct_kappas':bool(float(R['kappa_original'])>float(R['kappa_autonomous'])),
      'physical_clock_monotone_and_integral_identity':all(results[s]['max_clock_identity_error']<5e-9 for s in SEEDS),
      'physical_time_green_density_causal':all(results[s]['max_causal_future_kernel_abs']<1e-15 for s in SEEDS),
      'physical_time_green_density_scaling':all(results[s]['max_frozen_history_green_density_error']<1e-13 for s in SEEDS),
      'independent_nonuniform_tau_quadrature':all(results[s]['max_nonlinear_tau_quadrature_response_error']<1e-8 for s in SEEDS),
      'fixed_history_total_response_unchanged':all(results[s]['max_physical_fixed_history_total_error']<1e-8 for s in SEEDS),
      'model_tangent_matches_nonlinear_secant':all(results[s]['max_tangent_vs_kappa_secant_error']<2e-7 for s in SEEDS),
      'matched_tau_nonlinear_secant_matches_clock_corrected':all(results[s]['max_matched_tau_vs_clock_corrected_secant_error']<2e-7 for s in SEEDS),
      'model_tangent_clock_corrected_matches_nonlinear':all(results[s]['max_matched_tau_vs_clock_corrected_model_error']<2e-7 for s in SEEDS),
      'physical_vectors_invariant_under_both_observer_families':all(
           rec['max_physical_vector_reconstruction_error']<1e-12 and
           rec['max_kernel_reconstruction_error']<1e-12 and
           rec['max_response_norm_error']<1e-12
           for s in SEEDS for rec in results[s]['observer_families'].values()),
      'moving_frame_source_additivity':all(rec['max_channel_additivity_error']<1e-9
           for s in SEEDS for rec in results[s]['observer_families'].values()),
      'covariant_connection_restores_derivative':all(
           rec['max_covariant_derivative_reconstruction_error']<2e-6
           for s in SEEDS for rec in results[s]['observer_families'].values()),
      'source_age_maps_increasing_physical_clock':all(results[s]['source_ages']['radial_omega2']['median_tau']>0 for s in SEEDS),
      'noninvariant_component_cancellation_not_confused_with_material_memory':True,
      'source_channel_causality_is_linearized_not_physical_actuator_certificate':True,
    }
    csv_write(OUT/'physical_clock_frame_section_ledger.csv',ledger)
    np.savez_compressed(OUT/'physical_clock_frame_green_tensors.npz',**arrays)
    report={'version':'v18.51v',
      'title':'Physical-Clock Green Transport & Moving-Frame Invariance Audit',
      'mode':'full',
      'classification':'CLOCK_REPARAMETERIZED_GREEN_AND_FRAME_COVARIANCE_NUMERICAL_DIAGNOSTIC_NOT_MATERIAL_MEMORY_CERTIFICATE',
      'all_passed':bool(all(checks.values())),'scientific_checks':checks,
      'historical_kappa_original':float(R['kappa_original']),
      'historical_kappa_autonomous':float(R['kappa_autonomous']),
      'source_channels':list(CHANNELS), 'sections_per_seed':len(xi),
      'physical_model_lineage':lineage,'seeds':results,
      'limitations':[
          'The frozen-history Green density K_tau=K_xi/c(source) is a reparameterization of the same fixed-xi tangent response; alone it does not imply invariance of full kappa sensitivity at fixed physical tau.',
          'The parameter-dependent physical clock requires a separate clock-shift correction. This is measured using already-integrated distinct-kappa trajectories and checked with direct tau-matched interpolation.',
          'Clock-shift attribution to individual source channels is not unique without an additional physically justified clock-source decomposition; it is retained as a separate term.',
          'The full 2D Cartesian response and its Euclidean norm survive a passive rotating-frame change. Radial and azimuthal component magnitudes and cancellation percentages are chart dependent.',
          'The observer rotation is a known passive O(2) coordinate transformation and does not establish arbitrary state-dependent gauge invariance of physical laws.',
          'All source forcing is an on-trajectory first variation at a single midpoint kappa0; not independently controlled physical actuators.',
          'No complete physical return section, interval certificate, persistent shadowing law, or constitutive material memory is certified.',
          'Time conversions are tied to the actual positive physical clock of the two seed trajectories over xi in [0,0.2].',
      ],
      'environment':{'python':sys.version,'numpy':np.__version__,
           'scipy':scipy.__version__,'platform':platform.platform()}}
    (OUT/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    lines=['# Q550 v18.51v — server scientific report','',
           '**Status:** '+('PASS' if report['all_passed'] else 'FAIL'),
           '**Classification:** '+report['classification'],'',
           'Produced by the server run; the pre-distributed scientific report is only a protocol.','',
           '## Fixed physical-clock sensitivity versus frozen-clock kernel','']
    for s in SEEDS:
        d=results[s]
        lines.extend(['### '+s,'',
          f'- Physical elapsed tau: {d["physical_tau_final"]:.15g}',
          f'- Clock start/end: {d["clock_start"]:.15g} / {d["clock_final"]:.15g}',
          f'- Terminal fixed-xi Cartesian sensitivity norm: {d["terminal_fixed_xi_sensitivity_norm"]:.12g}',
          f'- Terminal fixed-tau clock correction norm: {d["terminal_clock_correction_norm"]:.12g}',
          f'- Terminal correction / fixed-xi norm: {d["terminal_clock_correction_ratio"]:.12g}',
          f'- Nonuniform tau quadrature discrepancy: {d["max_nonlinear_tau_quadrature_response_error"]:.3e}',
          f'- Direct matched-tau versus corrected formula: {d["max_matched_tau_vs_clock_corrected_secant_error"]:.3e}',
          f'- Connection derivative defect (maximum across frames): {max(v["max_covariant_derivative_reconstruction_error"] for v in d["observer_families"].values()):.3e}',
          '- Physical source ages at final observation (radial omega2): '+json.dumps(d['source_ages']['radial_omega2']),''])
    lines+=['## Interpretation boundaries','']+['- '+str(t) for t in report['limitations']]
    (OUT/'SCIENTIFIC_REPORT_SERVER.md').write_text('\n'.join(lines)+'\n')
    render_figures(xi,U,report,arrays)
    print('Q550 v18.51v full',sum(checks.values()),'/',len(checks),
          'PASS' if report['all_passed'] else 'FAIL',flush=True)
    if not report['all_passed']:
        print('Failed checks:',[k for k,v in checks.items() if not v],flush=True)
        raise RuntimeError('v18.51v scientific checks failed')

if __name__=='__main__': main()
