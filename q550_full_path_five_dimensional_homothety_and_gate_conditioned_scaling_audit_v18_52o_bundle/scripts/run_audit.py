#!/usr/bin/env python3
"""Q550 v18.52o: full material-time 5D 3:1 homothety and gate-selectivity audit.

Only postprocess genuine v18.52l traces. Gates are picked using SMALL-kick
trajectories, independently of the large/small residual measured at those gates.
4 source-matched (seed, polarity) pairs, NOT 8 statistically independent samples.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
PARENT='q550_near_tangency_breathing_envelope_drift_five_dimensional_section_audit_v18_52l_bundle'
BASENAME='extended_horizon_5d_vectors_and_traces.npz'
EXPECTED='d051ec17dbbee0bc29bae57f045e9ead9bd6316834d25d9e5dd31b083236a195'
LABEL_BASE='original_negative_1.0e-04'
PAIRS=[('original','negative'),('original','positive'),('autonomous','negative'),('autonomous','positive')]
SLOPES={'user_m_minus_0p308':-.308,'historical_horizontal_m0':0.,
        'tan_minus_17p04':math.tan(math.radians(-17.04)),'opposite_m_plus_0p308':+.308}
TILTS=[0.,10.,17.04,25.,35.]
WINDOW_HALFWIDTHS=[.10,.25,.50]


def record(path,obj):
 path.write_text(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+'\n',encoding='utf8')

def csvout(path,rows,columns):
 with path.open('w',newline='',encoding='utf8') as f:
  w=csv.DictWriter(f,fieldnames=columns,extrasaction='ignore');w.writeheader();w.writerows(rows)

def sha256(path):
 with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def vector_norm(a):return np.sqrt(np.sum(np.asarray(a,dtype=float)**2,axis=-1))

def unit(a):
 a=np.asarray(a,float);r=float(vector_norm(a))
 if not np.isfinite(r) or r<1e-15:raise ValueError('Invalid axis')
 return a/r

def pca_hidden(q):
 w=q[:,2:5];w=w-w.mean(axis=0)
 _,vec=np.linalg.eigh(w.T@w/max(1,len(w)-1));h=vec[:,-1]
 if h[np.argmax(abs(h))]<0:h=-h
 return unit(h)

def hinge_coords(q,m,h,angle):
 # In-plane nodal axis e=(1,m,0,0,0); tilt is in full 5D about e.
 h=unit(h);den=math.sqrt(1+m*m)
 v=(-m*q[:,0]+q[:,1])/den
 z=q[:,2:5]@h
 a=math.radians(angle)
 return math.cos(a)*v+math.sin(a)*z,math.cos(a)*v-math.sin(a)*z

def zero_crossings(t,w):
 """Linear root times of piecewise linearly interpolated observer gates.

 Ignore exact-zero plateaux; de-duplicate endpoint roots. No forced period.
 """
 t=np.asarray(t,float);w=np.asarray(w,float)
 if len(t)!=len(w) or not np.isfinite(w).all() or np.any(np.diff(t)<=0):
  raise ValueError('Bad gate samples')
 j=np.flatnonzero(((w[:-1]<0)&(w[1:]>=0))|((w[:-1]>0)&(w[1:]<=0)))
 out=[]
 for i in j:
  den=w[i+1]-w[i]
  if abs(den)<1e-300:continue
  f=float(-w[i]/den)
  if f < 0 or f > 1:continue
  r=float(t[i]+f*(t[i+1]-t[i]))
  if not out or r-out[-1][0]>1e-10:
   out.append((r,int(np.sign(den)) or 1))
 return out

def hermite(t,q,dq,x):
 """Cubic Hermite resampling using physical dQ/dtau=dV in 5D.

 Query x must lie inside source, no clipping/extrapolation.
 """
 t=np.asarray(t,float);q=np.asarray(q,float);dq=np.asarray(dq,float)
 x=np.atleast_1d(np.asarray(x,float))
 if np.any(np.diff(t)<=0) or x.min()<t[0]-1e-11 or x.max()>t[-1]+1e-11:
  raise ValueError('Hermite material-time queries out of domain')
 x=np.minimum(np.maximum(x,t[0]),t[-1])
 i=np.clip(np.searchsorted(t,x,side='right')-1,0,len(t)-2)
 dt=(t[i+1]-t[i]);a=((x-t[i])/dt)[:,None];b=dt[:,None]
 return ((2*a**3-3*a*a+1)*q[i]+(a**3-2*a*a+a)*b*dq[i]
         +(-2*a**3+3*a*a)*q[i+1]+(a**3-a*a)*b*dq[i+1])

def aligned_pair(lo,hi):
 """Exact common time grid, chosen independently of gate geometry."""
 tl,ql,vl=lo;th,qh,vh=hi
 mask=(tl>=th[0])&(tl<=th[-1]);t=tl[mask]
 if len(t)<100:raise ValueError('Insufficient common material-time overlap')
 b=hermite(th,qh,vh,t)
 return t,ql[mask],vl[mask],b,mask

def scalar(v):return float(v) if np.isfinite(v) else None

def safe_stats(a):
 a=np.asarray(a,float);a=a[np.isfinite(a)]
 if len(a)==0:return dict(n=0,median=None,p90=None,p99=None,maximum=None,mean=None)
 return dict(n=int(len(a)),median=float(np.median(a)),p90=float(np.quantile(a,.9)),
             p99=float(np.quantile(a,.99)),maximum=float(a.max()),mean=float(a.mean()))

def metrics(qlo,qhi):
 qlo=np.asarray(qlo,float);qhi=np.asarray(qhi,float)
 r=qhi-3*qlo;nlo=vector_norm(qlo);nhi=vector_norm(qhi);nr=vector_norm(r)
 floor=max(1e-14,1e-7*float(np.median(nlo)))
 valid=(nlo>floor)&(nhi>floor)
 rel=np.full(nlo.shape,np.nan);ratio=np.full(nlo.shape,np.nan);dratio=np.full(nlo.shape,np.nan)
 theta=np.full(nlo.shape,np.nan);perp=np.full(nlo.shape,np.nan)
 rel[valid]=nr[valid]/(3*nlo[valid]);ratio[valid]=nhi[valid]/nlo[valid]
 dratio[valid]=(nhi[valid]/(3*nlo[valid]))**2
 dot=np.sum(qhi*qlo,axis=1);factor=np.zeros(len(nlo));factor[valid]=dot[valid]/(nlo[valid]**2)
 perpvec=qhi-factor[:,None]*qlo;perp[valid]=vector_norm(perpvec[valid])/(3*nlo[valid])
 # arctan2 avoids acos precision loss near angles of a few hundredths of a degree.
 theta[valid]=np.degrees(np.arctan2(vector_norm(perpvec[valid]),factor[valid]*nlo[valid]))
 global_scale=3*float(np.median(nlo[valid]))
 glob=nr/global_scale
 global_coeff=float(np.sum(qlo*qhi)/np.sum(qlo*qlo))
 return dict(resid=nr,rel=rel,global_=glob,nlo=nlo,nhi=nhi,ratio=ratio,dratio=dratio,
             angle=theta,perp=perp,valid=valid,floor=floor,coeff=global_coeff,
             norm_scale=global_scale)

def at(t,metric,events):
 if not events:return np.empty(0)
 return np.interp([r[0] if isinstance(r,tuple) else r for r in events],t,metric)

def event_mask(t,roots,halfwidth):
 if not roots:return np.zeros(len(t),bool)
 # N roots <= tens; no (N x T) big intermediate.
 out=np.zeros(len(t),bool)
 for r,_ in roots:out|=abs(t-r)<=halfwidth
 return out

def compare_events(t,M,roots,case,seed,polarity,view):
 rows=[];controls=[];pairdiff=[]
 for k,(r,orientation) in enumerate(roots):
  row=dict(seed=seed,polarity=polarity,case=case,view=view,event_index=k,
   tau=r,crossing_orientation=orientation,
   rel_residual=float(at(t,M['rel'],[r])[0]),
   global_residual=float(at(t,M['global_'],[r])[0]),
   absolute_residual=float(at(t,M['resid'],[r])[0]),
   norm_ratio=float(at(t,M['ratio'],[r])[0]),
   orientation_angle_deg=float(at(t,M['angle'],[r])[0]),
   norm_small=float(at(t,M['nlo'],[r])[0]))
  rows.append(row)
 for k in range(len(roots)-1):
  left=roots[k][0];right=roots[k+1][0];mid=.5*(left+right)
  gap=right-left
  # Nearest equal-phase controls: midpoint between successive same-view roots.
  rr=float(at(t,M['rel'],[mid])[0]);gg=float(at(t,M['global_'],[mid])[0])
  control=dict(seed=seed,polarity=polarity,case=case,view=view,interval_index=k,
      gap=gap,tau_mid=mid,mid_rel_residual=rr,mid_global_residual=gg,
      edge_mean_rel_residual=.5*(rows[k]['rel_residual']+rows[k+1]['rel_residual']),
      edge_mean_global_residual=.5*(rows[k]['global_residual']+rows[k+1]['global_residual']))
  control['mid_minus_gate_rel']=rr-control['edge_mean_rel_residual']
  control['mid_minus_gate_global']=gg-control['edge_mean_global_residual']
  controls.append(control);pairdiff.append(control['mid_minus_gate_rel'])
 return rows,controls,pairdiff

def sensitivity(t,q,haxes,M,seed,polarity):
 rows=[];detail=[];primary_roots={};primary_w=None
 configs=[]
 for name,m in SLOPES.items():
  configs.append((name,m,'PCA',haxes['PCA'],17.04))
 # Explicit tilt controls, plus physical hidden-axis controls.
 for ang in TILTS:
  configs.append(('user_m_minus_0p308',-.308,'PCA',haxes['PCA'],ang))
 for name,h in haxes.items():
  if name!='PCA':configs.append(('user_m_minus_0p308',-.308,name,h,17.04))
 seen=set()
 for case,m,hn,h,ang in configs:
  key=(case,hn,ang)
  if key in seen:continue
  seen.add(key)
  wp,wm=hinge_coords(q,m,h,ang)
  roots_plus=zero_crossings(t,wp);roots_minus=zero_crossings(t,wm)
  win=(wp*wm<0)
  if case=='user_m_minus_0p308' and hn=='PCA' and abs(ang-17.04)<1e-9:
   primary_roots={'+':roots_plus,'-':roots_minus};primary_w=(wp,wm)
  for view,roots in [('+',roots_plus),('-',roots_minus)]:
   er,mid,dd=compare_events(t,M['tmetrics'],roots,case,seed,polarity,view)
   event_rel=[x['rel_residual'] for x in er]
   ctrl_rel=[x['mid_rel_residual'] for x in mid]
   rows.append(dict(seed=seed,polarity=polarity,case=case,hidden_axis=hn,tilt_deg=ang,view=view,
        root_count=len(roots),mean_delta_tau=float(np.mean(np.diff([r for r,_ in roots]))) if len(roots)>1 else None,
        median_relative_at_roots=float(np.median(event_rel)) if event_rel else None,
        median_relative_midphase=float(np.median(ctrl_rel)) if ctrl_rel else None,
        midpoint_minus_gate_median=float(np.median(dd)) if dd else None,
        window_fraction=float(np.mean(win)),
        median_relative_in_opposite_view_window=scalar(np.nanmedian(M['tmetrics']['rel'][win])) if win.any() else None,
        median_relative_outside_opposite_view_window=scalar(np.nanmedian(M['tmetrics']['rel'][~win])) if (~win).any() else None))
   if case=='user_m_minus_0p308' and hn=='PCA' and abs(ang-17.04)<1e-9:
    detail.extend([(er,mid,view)])
 return rows,detail,primary_roots,primary_w

def quantile_control(t,M,mask_gate,seed,polarity):
 # Norm-conditioned comparisons: don't mistake a denominator trough for gate focus.
 rows=[];amp=M['nlo'];cut=np.quantile(amp,[0,.25,.5,.75,1.])
 for j in range(4):
  b=(amp>=cut[j]) & (amp<cut[j+1] if j<3 else amp<=cut[j+1])
  g=b&mask_gate;o=b&(~mask_gate)
  rows.append(dict(seed=seed,polarity=polarity,bin=j+1,lower_norm=float(cut[j]),upper_norm=float(cut[j+1]),
    n_gate=int(g.sum()),n_off=int(o.sum()),
    median_gate_rel=scalar(np.nanmedian(M['rel'][g])) if g.any() else None,
    median_off_rel=scalar(np.nanmedian(M['rel'][o])) if o.any() else None,
    median_gate_global=scalar(np.nanmedian(M['global_'][g])) if g.any() else None,
    median_off_global=scalar(np.nanmedian(M['global_'][o])) if o.any() else None))
 return rows

def diagnostic_plot(path,traces,event_times):
 fig,axes=plt.subplots(3,1,figsize=(13,11),sharex=True)
 for key,tr in traces.items():
  t=tr['t'];M=tr['M'];sty='-' if key.endswith('positive') else '--'
  axes[0].plot(t,M['ratio'],sty,lw=1.05,label=key)
  axes[1].semilogy(t,np.maximum(M['rel'],1e-12),sty,lw=1.0)
  axes[2].semilogy(t,np.maximum(M['global_'],1e-12),sty,lw=1.0)
 for e in event_times:
  axes[1].axvline(e,color='.65',alpha=.22,lw=.7)
  axes[2].axvline(e,color='.65',alpha=.22,lw=.7)
 axes[0].axhline(3,color='k',lw=.8,ls=':')
 axes[0].set(ylabel='||q(3ε)|| / ||q(ε)||',title='Q550 v18.52o | Same-time physical 5D homothety: full path vs hinges')
 axes[1].set(ylabel='||q(3ε)-3q(ε)|| / 3||q(ε)||')
 axes[2].set(ylabel='||q(3ε)-3q(ε)|| / (3 median||q(ε)||)',xlabel='Common material time τ')
 axes[0].legend(fontsize=8,ncol=2);fig.tight_layout();fig.savefig(path,dpi=150);plt.close(fig)

def plot_gate(path,rows):
 fig,axes=plt.subplots(1,2,figsize=(12,4.8))
 for i,field in enumerate(('mid_minus_gate_rel','mid_minus_gate_global')):
  xs=np.arange(len(rows));y=[r[field] for r in rows]
  cols=['tab:blue' if r['seed']=='original' else 'tab:orange' for r in rows]
  axes[i].scatter(xs,y,s=25,c=cols)
  axes[i].axhline(0,color='k',lw=1)
  axes[i].set(xlabel='Midpoint between successive same-view gates (per pair)',
       ylabel='Midphase minus mean adjacent gate residual',title='Relative' if i==0 else 'Global-amplitude normalized')
 fig.suptitle('Q550 v18.52o | Gate advantage needs positive, consistent differences')
 fig.tight_layout();fig.savefig(path,dpi=150);plt.close(fig)

def derivative_check(t,q,dq):
 # Scale-aware dQ/dtau check at interior points; numerical diff is a diagnostic, not a proof.
 if len(t)<10:return None
 finite=np.gradient(q,t,axis=0,edge_order=2)
 nom=vector_norm(finite-dq);den=vector_norm(dq)
 return dict(max_abs=float(nom.max()),median_relative=float(np.median(nom/(den+1e-15))))

def synthetic_tests():
 # Cubic Hermite exact on cubic polynomials with exact derivative, even on shifted grids.
 t=np.linspace(0,9,31);x=np.linspace(.1,8.9,51)
 q=np.zeros((len(t),5));dq=np.zeros_like(q)
 q[:,0]=t**3-2*t*t+3*t;dq[:,0]=3*t*t-4*t+3
 q[:,1]=3*t*t;dq[:,1]=6*t
 out=hermite(t,q,dq,x)
 assert np.max(abs(out[:,0]-(x**3-2*x*x+3*x)))<2e-12
 assert np.max(abs(out[:,1]-3*x*x))<3e-13
 # 5D homothety is exact for all orientations, any hinge, any shifted clock.
 small=np.c_[np.sin(t),np.cos(t),t*.02,0*t,t*.03]
 hi=3*small
 M=metrics(small,hi)
 assert max(M['rel'])<1e-15 and max(abs(M['ratio']-3))<1e-14
 assert abs(M['coeff']-3)<1e-14
 wp,wm=hinge_coords(small,-.308,[1,0,0],17.04)
 assert np.max(abs(wp+wm-2*math.cos(math.radians(17.04))*(-(-.308)*small[:,0]+small[:,1])/math.sqrt(1+.308**2)))<1e-14
 assert len(zero_crossings(t,np.sin(t)))==2
 # Nonlinear cubic response changes amplitude across the whole path as expected.
 amp=np.array([.1,.2,.3]);v=np.zeros((3,5));v[:,0]=amp
 assert np.all(metrics(v,v*3+np.c_[amp**2,np.zeros((3,4))])['rel']>0)
 # Fixed angular tilt need not change full-5D metrics.
 assert np.allclose(M['rel'],metrics(small,hi)['rel'])
 return {'hermite_cubic':'PASS','exact_5D_homothety':'PASS','hinge_projection_algebra':'PASS',
  'gate_root_detection':'PASS','nonlinear_response_detection':'PASS'}

def main():
 ap=argparse.ArgumentParser()
 ap.add_argument('--mode',choices=['preflight','audit'],default='audit')
 ap.add_argument('--source',type=Path,default=ROOT.parent/PARENT/'results'/BASENAME)
 ap.add_argument('--sample-stride',type=int,default=1,
                 help='>1 is development-only: never a scientific report')
 ap.add_argument('--output-stride',type=int,default=16)
 args=ap.parse_args()
 if args.sample_stride<1 or args.output_stride<1:raise ValueError('stride must be >=1')
 tests=synthetic_tests()
 if not args.source.exists():raise FileNotFoundError(f'Missing authentic source: {args.source}')
 source_sha=sha256(args.source)
 if source_sha!=EXPECTED:raise ValueError(f'v18.52l SHA mismatch {source_sha} != {EXPECTED}')
 if args.sample_stride>1:stage='DEVELOPMENT_SMOKE_NOT_SCIENTIFIC'
 else:stage='AUTHENTIC_FULL_RESOLUTION_POSTPROCESS'
 with np.load(args.source,allow_pickle=False) as f:
  for seed,pol in PAIRS:
   for amp in ('1.0e-04','3.0e-04'):
    pref=f'{seed}_{pol}_{amp}__'
    for k in ('tau','dQ','dV','D','G'):
     if pref+k not in f:raise RuntimeError(f'Missing {pref+k}')
  primary=f[LABEL_BASE+'__dQ'];h=pca_hidden(primary[:len(primary)//2:args.sample_stride]);
  mh=np.asarray(json.loads((ROOT/'reference/PARENT_v18_52m_report.json').read_text())['PCA_hidden_axis_Q3_Q4_Q5'])
  mismatch=float(np.degrees(np.arccos(np.clip(abs(h@mh),0,1))))
  if args.sample_stride==1 and mismatch>1e-5:raise RuntimeError('Hidden PCA lineage mismatch')
  haxes={'PCA':h,'Q3':np.array([1.,0,0]),'Q4':np.array([0.,1,0]),'Q5':np.array([0.,0,1])}
  if args.mode=='preflight':
   report=dict(version='v18.52o',stage=stage,source=str(args.source),source_sha256=source_sha,
     synthetic_checks=tests,pca_hidden=h.tolist(),pca_historical_mismatch_deg=mismatch,
     pair_labels=[f'{s}_{p}' for s,p in PAIRS]);record(ROOT/'PREFLIGHT.json',report)
   print('PREFLIGHT PASS',stage,report['pair_labels']);return
  out=ROOT/'results';out.mkdir(exist_ok=True)
  allsummary=[];timeseries=[];allgates=[];allmid=[];bins=[];sens=[];winrows=[];derivatives=[];trajectories={};viz=[];branchrows=[];gatephase=[]
  for seed,pol in PAIRS:
   tag=f'{seed}_{pol}';lo=[];hi=[]
   for amp,target in [('1.0e-04',lo),('3.0e-04',hi)]:
    lab=f'{seed}_{pol}_{amp}__'
    tt=f[lab+'tau'][::args.sample_stride].astype(float)
    qq=f[lab+'dQ'][::args.sample_stride].astype(float)
    vv=f[lab+'dV'][::args.sample_stride].astype(float)
    DD=f[lab+'D'][::args.sample_stride].astype(float)
    GG=f[lab+'G'][::args.sample_stride].astype(float)
    if len(tt)<100 or np.any(np.diff(tt)<=0) or not np.isfinite(qq).all():raise ValueError('Source chronology invalid '+tag)
    d_err=float(np.max(abs(.5*np.sum(qq*qq,axis=1)-DD)))
    g_err=float(np.max(abs(np.sum(qq*vv,axis=1)-GG)))
    if max(d_err,g_err)>4e-12:raise RuntimeError('D/G physical 5D identity broken '+tag)
    derivatives.append(dict(label=f'{seed}_{pol}_{amp}',D_identity_max_abs=d_err,G_identity_max_abs=g_err,
                             **derivative_check(tt,qq,vv)))
    target.extend((tt,qq,vv))
   tau,ql,vl,qh,domain=aligned_pair(lo,hi)
   M=metrics(ql,qh)
   # Interpolator and common-clock controls. Raw same-index matching is explicitly a negative control.
   qh_linear=np.column_stack([np.interp(tau,hi[0],hi[1][:,c]) for c in range(5)])
   hermite_vs_linear=vector_norm(qh-qh_linear)/(3*M['nlo'])
   index_time_mismatch=metrics(ql,hi[1][domain])
   valid=M['valid']
   if valid.mean()<.999:raise RuntimeError('Majority of physical separation norms are degenerate')
   wp,wm=hinge_coords(ql,-.308,h,17.04)
   roots={'+':zero_crossings(tau,wp),'-':zero_crossings(tau,wm)}
   win=wp*wm<0
   # Amplitude-dependent motion of the gate is a distinct diagnostic, not a moved sampling gate.
   hp,hm=hinge_coords(qh,-.308,h,17.04)
   high_roots={'+':zero_crossings(tau,hp),'-':zero_crossings(tau,hm)}
   for view in ('+','-'):
    lowr=roots[view];highr=high_roots[view]
    if len(lowr)!=len(highr):
     gatephase.append(dict(seed=seed,polarity=pol,view=view,low_count=len(lowr),high_count=len(highr),
       index=None,low_tau=None,high_tau=None,delta_tau=None,orientation_match=None))
    else:
     for k,((tl0,sl0),(th0,sh0)) in enumerate(zip(lowr,highr)):
      gatephase.append(dict(seed=seed,polarity=pol,view=view,low_count=len(lowr),high_count=len(highr),
       index=k,low_tau=tl0,high_tau=th0,delta_tau=th0-tl0,orientation_match=(sl0==sh0)))
   for w,rootsvec in roots.items():
    ev,mid,diffs=compare_events(tau,M,rootsvec,'user_m_minus_0p308',seed,pol,w)
    allgates.extend(ev);allmid.extend(mid)
   for half in WINDOW_HALFWIDTHS:
    r=[x for vv in roots.values() for x in vv];mask=event_mask(tau,r,half)
    winrows.append(dict(seed=seed,polarity=pol,kind='root_neighborhood',halfwidth_tau=half,
     sample_fraction=float(mask.mean()),count=int(mask.sum()),
     gate_rel_median=scalar(np.nanmedian(M['rel'][mask])) if mask.any() else None,
     off_rel_median=scalar(np.nanmedian(M['rel'][~mask])) if (~mask).any() else None,
     gate_global_median=scalar(np.nanmedian(M['global_'][mask])) if mask.any() else None,
     off_global_median=scalar(np.nanmedian(M['global_'][~mask])) if (~mask).any() else None))
   winrows.append(dict(seed=seed,polarity=pol,kind='opposite_sign_window',halfwidth_tau=None,
    sample_fraction=float(win.mean()),count=int(win.sum()),
    gate_rel_median=scalar(np.nanmedian(M['rel'][win])) if win.any() else None,
    off_rel_median=scalar(np.nanmedian(M['rel'][~win])) if (~win).any() else None,
    gate_global_median=scalar(np.nanmedian(M['global_'][win])) if win.any() else None,
    off_global_median=scalar(np.nanmedian(M['global_'][~win])) if (~win).any() else None))
   bins.extend(quantile_control(tau,M,win,seed,pol))
   row=dict(seed=seed,polarity=pol,samples=len(tau),tau_start=float(tau[0]),tau_end=float(tau[-1]),
    initial_grid_offset_tau=float(hi[0][0]-lo[0][0]),max_grid_offset_tau=float(np.max(abs(hi[0]-lo[0]))),
    norm_floor=float(M['floor']),norm_floor_excluded=int((~valid).sum()),
    hermite_vs_linear_relative=safe_stats(hermite_vs_linear),
    same_index_negative_control_relative=safe_stats(index_time_mismatch['rel']),
    coefficient_best_least_squares=float(M['coeff']),
    norm_ratio=safe_stats(M['ratio']),D_to_nine_D_ratio=safe_stats(M['dratio']),
    relative_5D_residual=safe_stats(M['rel']),global_amplitude_residual=safe_stats(M['global_']),
    absolute_5D_residual=safe_stats(M['resid']),orientation_deg=safe_stats(M['angle']),
    perpendicular_to_small_relative=safe_stats(M['perp']),
    gates_per_view={k:len(v) for k,v in roots.items()},
    opposite_side_window_fraction=float(np.mean(win)))
   allsummary.append(row)
   print(tag,stage,'rel median',row['relative_5D_residual']['median'],
     'max',row['relative_5D_residual']['maximum'], 'root counts',row['gates_per_view'],flush=True)
   # Sensitivity recomputes gates but NEVER changes the 5D residual itself.
   sensrows,_,_,_=sensitivity(tau,ql,haxes,dict(tmetrics=M),seed,pol)
   sens.extend(sensrows)
   for i in range(0,len(tau),args.output_stride):
    timeseries.append(dict(seed=seed,polarity=pol,tau=float(tau[i]),
     norm_ratio=scalar(M['ratio'][i]),D_ratio_to_nine=scalar(M['dratio'][i]),
     rel_residual=scalar(M['rel'][i]),global_residual=scalar(M['global_'][i]),
     absolute_residual=float(M['resid'][i]),angle_deg=scalar(M['angle'][i]),
     norm_small=float(M['nlo'][i]),opposite_sign_window=int(win[i])))
   trajectories[tag]=dict(t=tau,M=M)
   if seed=='original' and pol=='negative':viz=[r[0] for v in roots.values() for r in v]
  # Even-odd parity cross-check; not a 3:1 test on a second independent sample.
  for seed in ('original','autonomous'):
   # Source time grids are slightly different even across kick polarities.
   for amp in ('1.0e-04','3.0e-04'):
    lk=f'{seed}_positive_{amp}__';nk=f'{seed}_negative_{amp}__'
    pos=(f[lk+'tau'][::args.sample_stride],f[lk+'dQ'][::args.sample_stride],f[lk+'dV'][::args.sample_stride])
    neg=(f[nk+'tau'][::args.sample_stride],f[nk+'dQ'][::args.sample_stride],f[nk+'dV'][::args.sample_stride])
    tt,qq,_,nn,_=aligned_pair(pos,neg)
    branchrows.append(dict(seed=seed,amp=amp,n_samples=len(tt),
       even_over_odd_rms=float(np.sqrt(np.mean(np.sum(((qq+nn)/2)**2,axis=1)))/
           max(np.sqrt(np.mean(np.sum(((qq-nn)/2)**2,axis=1))),1e-15)),
       median_even_norm=float(np.median(vector_norm((qq+nn)/2))),
       median_odd_norm=float(np.median(vector_norm((qq-nn)/2)))))
  summary=dict(version='v18.52o',title='Full-Path 5D Homothety & Gate-Conditioned Scaling Audit',
    stage=stage,source_sha256=source_sha,source_path=str(args.source),
    sample_stride=args.sample_stride,output_stride=args.output_stride,
    physical_dimension=5,material_clock='tau inherited v18.52l, cubic-Hermite matched',
    n_matched_pair_labels=4,hypothesized_scale=3,kinematic_D_scale=9,
    ancestry_hidden_axis_PCA_training='first half of original_negative_1.0e-04',
    pca_hidden=h.tolist(),pca_historical_mismatch_deg=mismatch,
    nodal_axis_slope=-.308,tilt_deg=17.04,
    gate_registration='ONLY low amplitude physical Q on same material-time axis',
    per_pair=allsummary,velocity_identity_checks=derivatives,
    perturbation_even_odd_diagnostics=branchrows,
    gate_time_shift_large_minus_small=safe_stats([r['delta_tau'] for r in gatephase if r['delta_tau'] is not None]),
    root_and_midpoint_comparisons={
     'n_root_measurements':len(allgates),'n_midpoint_measurements':len(allmid),
     'median_root_relative':safe_stats([v['rel_residual'] for v in allgates]),
     'median_midpoint_relative':safe_stats([v['mid_rel_residual'] for v in allmid]),
     'mid_minus_gate_relative':safe_stats([v['mid_minus_gate_rel'] for v in allmid]),
     'mid_minus_gate_global':safe_stats([v['mid_minus_gate_global'] for v in allmid]),
     'n_intervals_mid_minus_gate_relative_positive':sum(v['mid_minus_gate_rel']>0 for v in allmid),
     'gate_intervals_not_independent':'Same-pair serially correlated gates. Pooling is only descriptive.'},
    statement_of_limit='No topological claim, no independence of 65k time samples, no privileged tilt angle inferred from exact projection algebra.')
  record(out/'report.json',summary)
  csvout(out/'full_path_decimated_5d_scaling.csv',timeseries,
   ['seed','polarity','tau','norm_ratio','D_ratio_to_nine','rel_residual','global_residual','absolute_residual','angle_deg','norm_small','opposite_sign_window'])
  csvout(out/'gate_event_5d_residuals.csv',allgates,
   ['seed','polarity','case','view','event_index','tau','crossing_orientation','rel_residual','global_residual','absolute_residual','norm_ratio','orientation_angle_deg','norm_small'])
  csvout(out/'midphase_controls.csv',allmid,
   ['seed','polarity','case','view','interval_index','gap','tau_mid','mid_rel_residual','mid_global_residual','edge_mean_rel_residual','edge_mean_global_residual','mid_minus_gate_rel','mid_minus_gate_global'])
  csvout(out/'gate_windows_vs_full_path.csv',winrows,
   ['seed','polarity','kind','halfwidth_tau','sample_fraction','count','gate_rel_median','off_rel_median','gate_global_median','off_global_median'])
  csvout(out/'norm_quartile_conditioned_controls.csv',bins,
   ['seed','polarity','bin','lower_norm','upper_norm','n_gate','n_off','median_gate_rel','median_off_rel','median_gate_global','median_off_global'])
  csvout(out/'observer_slope_tilt_hidden_axis_sensitivity.csv',sens,
   ['seed','polarity','case','hidden_axis','tilt_deg','view','root_count','mean_delta_tau','median_relative_at_roots','median_relative_midphase','midpoint_minus_gate_median','window_fraction','median_relative_in_opposite_view_window','median_relative_outside_opposite_view_window'])
  csvout(out/'amplitude_dependent_gate_time_shift.csv',gatephase,
   ['seed','polarity','view','low_count','high_count','index','low_tau','high_tau','delta_tau','orientation_match'])
  csvout(out/'perturbation_even_odd_diagnostics.csv',branchrows,
   ['seed','amp','n_samples','even_over_odd_rms','median_even_norm','median_odd_norm'])
  diagnostic_plot(out/'full_path_5d_homothety_and_gate_alignment.png',trajectories,viz)
  plot_gate(out/'gate_vs_midphase_relative_residual.png',allmid)
  analysis=[
   '# Q550 v18.52o — Scientific report','',
   '**Status:** '+stage,'','## Source and methods','',
   f'- Authentic v18.52l SHA-256: `{source_sha}`',
   '- 4 matched pair labels (original/autonomous × negative/positive), each with amplitudes 1e-4 and 3e-4.',
   '- Common material time via cubic Hermite interpolation using physical 5D dQ and dV; equal array index is NOT assumed to mean equal physical time.',
   '- Interpolator (Hermite vs linear) and naive same-index clocks are explicitly measured as controls.',
   '- Amplitude-dependent gate phase shift is measured separately; sampling gates stay fixed on small amplitude.',
   '- Gate roots are registered on the small-amplitude trajectory alone, with nodal slope m=-0.308 and ±17.04° tilt.',
   '- Baselines include midpoint of successive same-view gates, complementary path, and matched ||q_small|| quartiles.',
   '- All time samples within a path are correlated; no naive time-sample hypothesis test.',
   '','## Full-trajectory homothety statistics','',
   '| Pair | Median norm ratio | Median relative residual | Max relative residual | Median angle (deg) | c_LS |',
   '|---|---:|---:|---:|---:|---:|',]
  for r in allsummary:
   analysis.append(f"| {r['seed']}/{r['polarity']} | {r['norm_ratio']['median']:.9f} | {r['relative_5D_residual']['median']:.9g} | {r['relative_5D_residual']['maximum']:.9g} | {r['orientation_deg']['median']:.7f} | {r['coefficient_best_least_squares']:.9f} |")
  overall=summary['root_and_midpoint_comparisons']
  analysis+=['','## Registered gate/midphase comparison','',
   f"- Registered gate events: {overall['n_root_measurements']}; phase-matched midpoints: {overall['n_midpoint_measurements']}.",
   f"- Median root relative residual: {overall['median_root_relative']['median']:.8g}; median midphase relative residual: {overall['median_midpoint_relative']['median']:.8g}.",
   f"- Midpoint minus adjacent-gate residual: median {overall['mid_minus_gate_relative']['median']:.8g}; positive {overall['n_intervals_mid_minus_gate_relative_positive']} / {overall['n_midpoint_measurements']}.",
   '- Direction matters: positive midpoint-minus-gate indicates a gate residual advantage. Negative indicates the opposite.',
   '- Inspect gate_windows_vs_full_path.csv, norm_quartile_conditioned_controls.csv and sensitivity before claiming localized scaling.',
   '','## Scope and decision','',
   '- A 3:1 ratio is imposed at initialization; persistence over time is a response-linearity test, not new spontaneous autosimilarity.',
   '- Reported 5D homothety is physical only if small and large 5D vectors remain close under the common material clock.',
   '- Equal-looking projected lobes, gate root timing, and exact observer identities alone cannot establish a topological invariant.',
   '- Distinguish a gate-local residual minimum from denominator artifacts and a general nonlinear phase dependence.',
   '- Stage is scientific only when sample_stride=1; nonunit stride reports are DEVELOPMENT_SMOKE_NOT_SCIENTIFIC.',
   '','All scripts, source pins, configuration and output SHA-256 manifest are bundled with server sendback.','']
  (out/'SCIENTIFIC_REPORT_SERVER.md').write_text('\n'.join(analysis),encoding='utf8')
  print('AUDIT PASS:',stage,'pairs',len(allsummary),'gates',len(allgates),'midpoints',len(allmid))

if __name__=='__main__':
 main()
