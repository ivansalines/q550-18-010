#!/usr/bin/env python3
"""Q550 v18.52n: four-point provenance and physical 5D paired-hinge observer audit.

This is only postprocessing of authenticated v18.52l 5D traces. A user-observed
four-point line is NOT identified from a screenshot unless the four coordinates
are supplied separately. The four historical Q2=0 crossings are negative controls.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,platform
from pathlib import Path
import numpy as np
import scipy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
PARENT='q550_near_tangency_breathing_envelope_drift_five_dimensional_section_audit_v18_52l_bundle'
SOURCE='extended_horizon_5d_vectors_and_traces.npz'
LABELS=(
 'original_negative_1.0e-04','original_positive_1.0e-04',
 'original_negative_3.0e-04','original_positive_3.0e-04',
 'autonomous_negative_1.0e-04','autonomous_positive_1.0e-04',
 'autonomous_negative_3.0e-04','autonomous_positive_3.0e-04')
PRIMARY=LABELS[0]
THETA_GRID=(0.,10.,17.04,25.,35.)


def dump(path,obj):
 path.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf8')

def csv_write(path,rows,columns=None):
 cols=columns or (list(rows[0]) if rows else ['label','tau'])
 with path.open('w',newline='',encoding='utf8') as f:
  w=csv.DictWriter(f,fieldnames=cols,extrasaction='ignore');w.writeheader();w.writerows(rows)

def norm(x):
 x=np.asarray(x,float);n=np.linalg.norm(x)
 if n<1e-15 or not np.isfinite(n): raise ValueError('Invalid basis norm')
 return x/n

def basis(m):
 """One specified physical in-plane line Q2=m*Q1+b, and its unit normal."""
 if not np.isfinite(m): raise ValueError('Invalid slope')
 e=np.array([1.,float(m),0,0,0]);n=np.array([-float(m),1.,0,0,0]);
 return norm(e),norm(n)

def coordinates(q,m,b):
 """Orthogonal planar coordinates relative to a hinge through (0,b,0,0,0)."""
 e,n=basis(m)
 x=np.asarray(q,float)
 # Avoid a full copy of 65537x5 per projection.
 u=x[...,0]*e[0]+(x[...,1]-b)*e[1]
 v=x[...,0]*n[0]+(x[...,1]-b)*n[1]
 return u,v

def projection(q,m,b,h,theta):
 u,v=coordinates(q,m,b)
 z=np.asarray(q)[:,2:5]@norm(h)
 a=math.radians(theta)
 wp=math.cos(a)*v+math.sin(a)*z
 wm=math.cos(a)*v-math.sin(a)*z
 return u,v,z,wp,wm

def normal_frames(m,h,theta):
 """Orthonormal observer axes in 5D with the *same sloped physical hinge*."""
 e,n=basis(m);hh=np.r_[np.zeros(2),norm(h)]
 a=math.radians(theta)
 return e,n,math.cos(a)*n+math.sin(a)*hh,math.cos(a)*n-math.sin(a)*hh

def pca_hidden(q):
 v=q[:,2:5];v=v-v.mean(axis=0)
 w,vec=np.linalg.eigh(v.T@v/max(1,len(v)-1))
 h=vec[:,-1];h*=1. if h[np.argmax(np.abs(h))]>=0 else -1.
 return norm(h),w[::-1].tolist()

def fit_four_points(x,y):
 """Both constrained through-origin and affine OLS, plus line orientation TLS.

 Four observed points need not lie on a unique line; no topological claim.
 """
 x=np.asarray(x,float);y=np.asarray(y,float)
 if x.shape!=(4,) or y.shape!=(4,) or not np.isfinite(x).all() or not np.isfinite(y).all():
  raise ValueError('Exactly four finite 2D coordinates required')
 denom=float(x@x)
 if denom<1e-22: raise ValueError('Four points lack useful Q1 for origin-constrained slope')
 m0=float((x@y)/denom)
 x0=x-x.mean();y0=y-y.mean();vv=float(x0@x0)
 if vv<1e-22: raise ValueError('Nearly vertical four-point set; cannot estimate Q2/Q1 finite slope')
 m=float((x0@y0)/vv);b=float(y.mean()-m*x.mean())
 resid=y-(m*x+b)
 resid0=y-m0*x
 scale=max(float(np.max(np.linalg.norm(np.c_[x,y],axis=1))),1e-15)
 lo=[]
 for omit in range(4):
  ids=np.arange(4)!=omit;z=x[ids];t=y[ids];zc=z-z.mean();tc=t-t.mean();d=float(zc@zc)
  lo.append(float((zc@tc)/d) if d>1e-22 else None)
 cov=np.cov(np.c_[x,y],rowvar=False)
 ev,evec=np.linalg.eigh(cov);d=evec[:,-1]
 if d[0]<0: d=-d
 tls_angle=float(np.degrees(math.atan2(d[1],d[0])))
 return dict(n=4,through_origin_slope=m0,through_origin_angle_deg=math.degrees(math.atan(m0)),
             affine_slope=m,affine_intercept=b,affine_angle_deg=math.degrees(math.atan(m)),
             tls_angle_deg=tls_angle,
             affine_rmse=float(np.sqrt(np.mean(resid**2))),
             affine_max_abs_residual=float(np.max(abs(resid))),
             through_origin_rmse=float(np.sqrt(np.mean(resid0**2))),
             normalized_affine_rmse=float(np.sqrt(np.mean(resid**2))/scale),
             per_point_affine_residuals=resid.tolist(),
             leave_one_out_affine_slopes=lo,
             x_span=float(np.ptp(x)),y_span=float(np.ptp(y)),
             angle_difference_ols_tls_deg=float(abs(math.degrees(math.atan(m))-tls_angle)))

def load_user_points(path):
 with path.open(newline='',encoding='utf8') as f:
  rows=list(csv.DictReader(f))
 if len(rows)!=4: raise ValueError('--four-points must be a CSV containing exactly four data rows')
 try:x=[float(r['Q1']) for r in rows];y=[float(r['Q2']) for r in rows]
 except (KeyError,ValueError) as e:raise ValueError('CSV must include finite numeric Q1,Q2 columns') from e
 names=[r.get('point_id',f'P{i+1}') for i,r in enumerate(rows)]
 return np.array(x),np.array(y),names

def crossing_events(t,q,dq,n,anchor,m):
 """Piecewise-linear roots of (q-anchor)·n=0, with 5D posture.

 n is unit normal to the tilted observer's transverse axis. Interpolation
 is an approximation, *not* an interval-certified ODE event root.
 """
 nn=norm(n);w=(q[:,0]-anchor[0])*nn[0]+(q[:,1]-anchor[1])*nn[1]+q[:,2:5]@nn[2:5]
 # Strict sign changes; exact zero nodes incorporated separately.
 ix=np.flatnonzero(w[:-1]*w[1:]<0)
 tol=max(1e-16,1e-13*float(np.max(np.abs(w))))
 ids=list(ix)
 zs=np.flatnonzero(abs(w)<tol)
 ev=[]
 for i in ids:
  f=float(-w[i]/(w[i+1]-w[i]));tt=float(t[i]+f*(t[i+1]-t[i]));qq=q[i]*(1-f)+q[i+1]*f;vv=dq[i]*(1-f)+dq[i+1]*f
  ev.append((tt,qq,vv))
 for i in zs:
  if i and i<len(t)-1:ev.append((float(t[i]),q[i],dq[i]))
 ev.sort(key=lambda z:z[0]);out=[]
 e,_=basis(m)
 for tt,qq,vv in ev:
  if out and tt-out[-1]['tau']<1e-8:continue
  a=qq-anchor;u=float(a@e);perp=float(np.linalg.norm(a-u*e))
  g=float(vv@nn);sp=float(np.linalg.norm(vv))
  out.append(dict(tau=tt,u=u,v_normal_derivative=g,
     transversality=abs(g)/sp if sp>1e-18 else None,
     crossing='up' if g>0 else ('down' if g<0 else 'grazing'),
     full_5d_radius=float(np.linalg.norm(qq)),D=float(qq@qq/2),
     hinge_perp_radius=perp,**{f'q{i+1}':float(qq[i]) for i in range(5)}))
 return out

def window_stats(t,wa,wb):
 """Opposite-sign dwell time of two piecewise-linear projections; exact for the interpolants.

 Gate endpoints are union of zero crossings (in continuous linear interpolation);
 no coarse grid-bin durations or expensive per-timepoint Python loops.
 """
 t=np.asarray(t,float);a=np.asarray(wa,float);b=np.asarray(wb,float)
 def roots(x):
  ix=np.flatnonzero(x[:-1]*x[1:]<0)
  return t[ix]+(t[ix+1]-t[ix])*(-x[ix]/(x[ix+1]-x[ix]))
 cuts=np.unique(np.r_[t[0],t[-1],roots(a),roots(b)])
 if len(cuts)<2:return dict(opposite_side_duration_fraction=0.,opposite_side_runs=0,windows=[])
 mid=.5*(cuts[:-1]+cuts[1:]);mask=(np.interp(mid,t,a)*np.interp(mid,t,b)<0)
 duration=float(np.sum(np.diff(cuts)[mask]));windows=[]
 idx=np.flatnonzero(mask)
 if len(idx):
  breaks=np.flatnonzero(np.diff(idx)>1)
  starts=np.r_[idx[0],idx[breaks+1]];ends=np.r_[idx[breaks],idx[-1]]
  for lo,hi in zip(starts,ends):
   s0=float(cuts[lo]);s1=float(cuts[hi+1]);windows.append(dict(start_tau=s0,end_tau=s1,duration=s1-s0))
 return dict(opposite_side_duration_fraction=duration/float(t[-1]-t[0]),opposite_side_runs=len(windows),windows=windows)

def shape_and_section(t,q,m,b,h,theta):
 u,v,z,wp,wm=projection(q,m,b,h,theta)
 def planar_ellipse(x,y):
  cc=np.cov(np.c_[x,y],rowvar=False);eig=np.linalg.eigvalsh(cc)
  return float(math.sqrt(max(eig[-1],1e-30)/max(eig[0],1e-30)))
 rperp=np.sqrt(v*v+np.sum(q[:,2:]**2,axis=1));r2=np.sqrt(v*v+z*z)
 au=np.unwrap(np.arctan2(z,v))
 rr=float(np.corrcoef(u,rperp)[0,1]) if np.std(u)>0 and np.std(rperp)>0 else 0.
 return dict(original_planar_axis_ratio=planar_ellipse(q[:,0],q[:,1]),
             sloped_untilted_axis_ratio=planar_ellipse(u,v),
             plus_ellipse_axis_ratio=planar_ellipse(u,wp),
             minus_ellipse_axis_ratio=planar_ellipse(u,wm),
             intrinsic_hinge_transverse_radius_cv=float(np.std(rperp)/max(np.mean(rperp),1e-30)),
             transverse_2d_radius_cv=float(np.std(r2)/max(np.mean(r2),1e-30)),
             cone_linear_radius_vs_axis_R2=rr*rr,
             axial_positive_step_fraction=float(np.mean(np.diff(u)>0)),
             transverse_reduced_net_phase_turns=float((au[-1]-au[0])/(2*math.pi)),
             phase_turns_topological_invariant=False)

def nearest_pairs(a,b):
 if not a or not b:return []
 aa=np.array([v['tau'] for v in a]);bb=np.array([v['tau'] for v in b]);out=[]
 for i,row in enumerate(a):
  j=int(np.argmin(abs(bb-row['tau'])))
  if int(np.argmin(abs(aa-bb[j])))!=i:continue
  peer=b[j];q1=np.array([row[f'q{k}'] for k in range(1,6)]);q2=np.array([peer[f'q{k}'] for k in range(1,6)])
  x=float(np.linalg.norm(q1));y=float(np.linalg.norm(q2));c=float(np.dot(q1,q2)/(x*y)) if x*y>1e-30 else None
  out.append(dict(tau_plus=row['tau'],tau_minus=peer['tau'],delta_tau=peer['tau']-row['tau'],
                  full_5d_angle_deg=float(np.degrees(np.arccos(np.clip(c,-1,1)))) if c is not None else None,
                  full_5d_distance=float(np.linalg.norm(q1-q2)),
                  D_ratio=peer['D']/max(row['D'],1e-30)))
 return out

def serial_pairs(events):
 res=[]
 for a,b in zip(events[:-1],events[1:]):
  v=np.array([a[f'q{i}'] for i in range(1,6)]);w=np.array([b[f'q{i}'] for i in range(1,6)])
  nv=np.linalg.norm(v);nw=np.linalg.norm(w);c=float(v@w/(nv*nw)) if nv*nw>1e-30 else None
  res.append(dict(tau0=a['tau'],tau1=b['tau'],delta_tau=b['tau']-a['tau'],
           D_ratio=b['D']/max(a['D'],1e-30),
           full_5d_angle_deg=float(np.degrees(np.arccos(np.clip(c,-1,1)))) if c is not None else None,
           same_posture_residual=float(np.linalg.norm(w-v)/max(nv,nw,1e-30))))
 return res

def slope_controls(m):
 return [('historical_Q1_axis',0.),('observed_cartesian_m',m),
         ('17p04_geometry_m',math.tan(math.radians(-17.04))),
         ('opposite_slope_control',-m)]

def load_source(src,hash_verify=True):
 if not src.exists():raise FileNotFoundError(f'Authentic v18.52l source absent: {src}')
 p=json.loads((ROOT/'reference/PARENT_v18_52l_OUTPUT_MANIFEST.json').read_text())
 expected=p['file_sha256']['results/'+SOURCE]
 if hash_verify:
  with src.open('rb') as fh: actual=hashlib.file_digest(fh,'sha256').hexdigest()
 else: actual='HASH_NOT_CHECKED'
 if hash_verify and actual!=expected:raise ValueError(f'Actual source SHA {actual} != pinned ancestral SHA {expected}')
 with np.load(src,allow_pickle=False) as f:
  for label in LABELS:
   for key in ('tau','dQ','dV','D','G'):
    if label+'__'+key not in f:raise ValueError('Missing '+label+'__'+key)
   n=len(f[label+'__tau'])
   if n<1000 or f[label+'__dQ'].shape!=(n,5):raise ValueError('Bad source shape')
  count=len(f[LABELS[0]+'__tau'])
 return dict(version='v18.52l',path=str(src),sha256=expected,observed_sha256=actual,
             trajectories=len(LABELS),samples_per_trajectory=count)

def synthetic_tests():
 # Basis: orthogonality and line origin / affine offset.
 for m,b in ((0.,0.),(-.308,0.),(-.308,.1),(1.,-.2)):
  e,n=basis(m)
  assert abs(e@n)<1e-13 and abs(np.linalg.norm(e)-1)<1e-13
  h=np.array([1.,0,0]);ee,nn,wp,wm=normal_frames(m,h,17.04)
  assert abs(ee@wp)<1e-13 and abs(ee@wm)<1e-13
  q=np.array([[0,b,0,0,0],[1,b+m,0,0,0.]])
  u,v=coordinates(q,m,b)
  assert max(abs(v))<1e-12
  assert np.isclose(u[1],np.sqrt(1+m*m))
  planar=np.tile(np.array([1.,2.,0.,0.,0.]),(10,1))
  _,_,_,a,c=projection(planar,m,b,h,17.04)
  assert max(abs(a-c))<1e-14
  nonplanar=planar.copy();nonplanar[:,2]=.2
  _,_,_,a,c=projection(nonplanar,m,b,h,17.04)
  assert max(abs(a-c))>.05
 # Four collinear observed points: complete affine / origin response.
 xx=np.array([-2.,-1.,1.,2.]); yy=-.308*xx+.002
 fit=fit_four_points(xx,yy)
 assert abs(fit['affine_slope']+.308)<1e-14 and abs(fit['affine_intercept']-.002)<1e-14
 assert fit['affine_rmse']<1e-13 and fit['through_origin_rmse']>.001
 assert abs(fit['through_origin_slope']+.308)<1e-14
 # Same 5D trajectory, tilted-root interpolation, absolute time order.
 t=np.linspace(0,5,1001)
 q=np.zeros((len(t),5));q[:,0]=.1*t;q[:,1]=np.sin(5*t);q[:,2]=.3*np.cos(5*t)
 dq=np.gradient(q,t,axis=0);ep,n,wp,wm=normal_frames(-.308,np.array([1.,0,0]),17.04)
 evp=crossing_events(t,q,dq,wp,np.zeros(5),-.308)
 evm=crossing_events(t,q,dq,wm,np.zeros(5),-.308)
 assert len(evp)>=6 and len(evm)>=6
 for row in evp:
  vec=np.array([row[f'q{i}'] for i in range(1,6)])
  assert abs(vec@wp)<1e-12
 u,v,z,a,c=projection(q,-.308,0,np.array([1.,0,0]),17.04)
 assert np.allclose(a+c,2*math.cos(math.radians(17.04))*v)
 assert np.allclose(a-c,2*math.sin(math.radians(17.04))*z)
 w=window_stats(t,a,c)
 assert 0<w['opposite_side_duration_fraction']<1 and w['opposite_side_runs']>0
 w0=window_stats(t,v,v)
 assert w0['opposite_side_duration_fraction']==0 and w0['opposite_side_runs']==0
 return dict(slope_hinge_basis='PASS',sloped_affine_hinge='PASS',four_point_fits='PASS',
             planar_null_and_nonplanar_control='PASS',single_path_two_views='PASS',
             five_d_crossing_interpolation='PASS',window_integration='PASS')

def plot_views(out,primary,records,m,b,theta,points=None):
 tt=records['tau'];q=records['q'];h=records['h'];step=max(1,len(tt)//3000)
 fig,axs=plt.subplots(2,2,figsize=(13.3,10.5),layout='constrained')
 scaler=1./max(float(np.max(np.linalg.norm(q,axis=1))),1e-20)
 for ax,mm,tag in ((axs[0,0],0.,'v18.52m baseline: Q1 horizontal hinge'),
                   (axs[0,1],m,'v18.52n: slope-constrained nodal hinge')):
  u,v,z,wp,wm=projection(q,mm,b if mm==m else 0,h,theta)
  aa=ax.scatter(u[::step]*scaler,wp[::step]*scaler,c=tt[::step],s=2.5,cmap='viridis',rasterized=True)
  ax.scatter(u[::step]*scaler,wm[::step]*scaler,c=tt[::step],s=2.5,cmap='plasma',marker='x',alpha=.40,rasterized=True)
  ax.axhline(0,c='gray',lw=.9,ls='--');ax.set_title(tag)
  ax.set(xlabel='Axial coordinate (normalized)',ylabel='Two +/- tilted shadows (normalized)')
  ax.set_aspect('equal',adjustable='datalim');ax.grid(alpha=.16)
 ax=axs[1,0]
 ax.plot(q[::step,0]*scaler,q[::step,1]*scaler,lw=.85,label='Single physical 5D path, projected to (Q1,Q2)')
 grid=np.array([np.min(q[:,0]),np.max(q[:,0])]);ax.plot(grid*scaler,(m*grid+b)*scaler,'--',label=f'Proposed line: Q2={m:.4f} Q1 + {b:.2e}')
 if points is not None:
  x,y=points;ax.scatter(x*scaler,y*scaler,s=35,c='red',zorder=5,label='Provided FOUR observed points')
 ax.set(xlabel='Physical relative Q1 / scale',ylabel='Physical relative Q2 / scale',title='Nodal cut: a geometric candidate, not a 5D invariant')
 ax.grid(alpha=.2);ax.legend(fontsize=7);ax.set_aspect('equal',adjustable='datalim')
 ax=axs[1,1];u,v,z,wp,wm=projection(q,m,b,h,theta)
 ax.plot(tt[::step],wp[::step]*scaler,label=f'+{theta:g} deg',lw=.95)
 ax.plot(tt[::step],wm[::step]*scaler,label=f'-{theta:g} deg',lw=.95)
 ax.axhline(0,c='gray',lw=.75,ls='--');ax.set(xlabel='Matched material time',ylabel='Tilted transverse / scale',title='One trajectory; gate openings are observer events')
 ax.grid(alpha=.2);ax.legend()
 fig.suptitle(f'Q550 v18.52n | {primary} | m={m:.6f} | theta={theta:.3f} deg')
 fig.savefig(out/'four_point_nodal_hinge_paired_observer_comparison.png',dpi=155)
 plt.close(fig)

def plot_sensitivity(out,rows,primary,m,theta):
 a=[r for r in rows if r['label']==primary and r['hidden_direction']=='PCA_train']
 fig,axs=plt.subplots(1,2,figsize=(12.5,4.8),layout='constrained')
 for tag,mm in slope_controls(m):
  v=sorted((r for r in a if r['slope_case']==tag),key=lambda r:r['theta_deg'])
  axs[0].plot([r['theta_deg'] for r in v],[r['opposite_side_duration_fraction'] for r in v],'-o',ms=3,label=tag)
  axs[1].plot([r['theta_deg'] for r in v],[r['plus_projected_axis_ratio'] for r in v],'-o',ms=3,label=tag)
 for ax in axs:ax.axvline(theta,color='gray',ls='--',lw=.8);ax.grid(alpha=.2);ax.legend(fontsize=7)
 axs[0].set(xlabel='Observer tilt magnitude (deg)',ylabel='Fraction of material time',title='Opposite-sign gate fraction')
 axs[1].set(xlabel='Observer tilt magnitude (deg)',ylabel='PCA semiaxis ratio (+ view)',title='Projection eccentricity diagnostics')
 fig.suptitle('Q550 v18.52n: distinct effects of hinge slope and observer tilt')
 fig.savefig(out/'nodal_slope_and_tilt_sensitivity.png',dpi=150);plt.close(fig)

def plot_four_points(out,plotdata,m,b):
 if plotdata is None:return
 x,y=plotdata;fig,ax=plt.subplots(figsize=(6.5,5.5),layout='constrained')
 ax.scatter(x,y,s=42,label='Supplied four observations',zorder=3)
 gx=np.linspace(min(x.min(),0),max(x.max(),0),100)
 fit=fit_four_points(x,y)
 ax.plot(gx,m*gx+b,ls='--',label='Used nodal hinge')
 ax.plot(gx,fit['through_origin_slope']*gx,label='OLS constrained to origin')
 ax.plot(gx,fit['affine_slope']*gx+fit['affine_intercept'],label='OLS affine')
 ax.set(xlabel='Observed Q1',ylabel='Observed Q2',title='Four-point line reconstruction: input coordinates')
 ax.legend(fontsize=8);ax.grid(alpha=.2)
 fig.savefig(out/'four_point_line_fit_and_origin_constraint.png',dpi=150);plt.close(fig)

def main():
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument('--mode',choices=['preflight','audit','synthetic'],default='audit')
 ap.add_argument('--source',type=Path,default=ROOT.parent/PARENT/'results'/SOURCE)
 ap.add_argument('--four-points',type=Path,default=None,help='External CSV with exactly four observed Q1,Q2 coordinates')
 ap.add_argument('--adopt-four-point-fit',action='store_true',help='Explicitly replace candidate slope AND intercept by affine OLS fit to supplied four points')
 ap.add_argument('--slope',type=float,default=-0.308)
 ap.add_argument('--intercept',type=float,default=0.,help='Nodal line intercept Q2 at Q1=0; default 0 is an assumption')
 ap.add_argument('--theta-deg',type=float,default=17.04,help='Tilt around hinge; independent of the slope angle')
 ap.add_argument('--primary',choices=LABELS,default=PRIMARY)
 ap.add_argument('--sample-stride',type=int,default=1,help='Development smoke only if >1')
 ap.add_argument('--output-dir',type=Path,default=ROOT/'results')
 args=ap.parse_args()
 if not 0<=args.theta_deg<80 or not np.isfinite(args.slope) or not np.isfinite(args.intercept):ap.error('Nonfinite or invalid geometric parameter')
 if args.sample_stride<1:ap.error('--sample-stride must be >=1')
 if args.adopt_four_point_fit and args.four_points is None:ap.error('--adopt-four-point-fit requires --four-points')
 out=args.output_dir;out.mkdir(parents=True,exist_ok=True)
 tests=synthetic_tests()
 if args.mode=='synthetic':print('SYNTHETIC TESTS PASS',tests);return
 src=load_source(args.source)
 points=None;four=dict(status='NO_USER_FOUR_POINT_COORDINATES_SUPPLIED',
   message='The historical four Q2=0 crossings are not the visually observed four off-axis points; slope -0.308 is a prespecified hypothesis, not a fitted result.')
 m=args.slope;b=args.intercept
 if args.four_points is not None:
  x,y,ids=load_user_points(args.four_points);points=(x,y)
  fit=fit_four_points(x,y)
  four=dict(status='USER_FOUR_POINT_COORDINATES_SUPPLIED',input_file=str(args.four_points),
            input_sha256=hashlib.sha256(args.four_points.read_bytes()).hexdigest(),point_ids=ids,
            points=[dict(id=ids[i],Q1=float(x[i]),Q2=float(y[i])) for i in range(4)],
            fit=fit,difference_from_user_slope=fit['affine_slope']-args.slope,
            used_for_geometry=bool(args.adopt_four_point_fit))
  if args.adopt_four_point_fit:m=fit['affine_slope'];b=fit['affine_intercept']
 if args.mode=='preflight':
  dump(out/'PREFLIGHT.json',dict(version='v18.52n',source=src,four_points=four,
    synthetic_checks=tests,theta_deg=args.theta_deg,slope=m,intercept=b))
  print('PREFLIGHT PASS',src['trajectories'],src['samples_per_trajectory'],four['status']);return
 # PCA hidden axis trained exclusively from first half of one historical trajectory.
 with np.load(args.source,allow_pickle=False) as data:
  primary_full=data[args.primary+'__dQ'];half=len(primary_full)//2
  ph,eig=pca_hidden(primary_full[:half:args.sample_stride]);test,_=pca_hidden(primary_full[half::args.sample_stride]);
  principal_ang=math.degrees(math.acos(float(np.clip(abs(ph@test),0,1))))
  m_report=json.loads((ROOT/'reference/PARENT_v18_52m_report.json').read_text())
  mh=np.array(m_report['PCA_hidden_axis_Q3_Q4_Q5']);baseline_angular_mismatch=math.degrees(math.acos(float(np.clip(abs(ph@mh),0,1))))
  if args.sample_stride==1 and baseline_angular_mismatch>1e-5:raise RuntimeError('PCA ancestry differs from genuine v18.52m')
  hidden_axes=dict(Q3=np.array([1.,0,0]),Q4=np.array([0.,1,0]),Q5=np.array([0.,0,1]),PCA_train=ph)
  rows=[];events=[];pairs=[];serial=[];windows=[];per_trace=[];negative=[];primary_record=None
  for label in LABELS:
   tt=np.asarray(data[label+'__tau'][::args.sample_stride],float)
   q=np.asarray(data[label+'__dQ'][::args.sample_stride],float)
   dq=np.asarray(data[label+'__dV'][::args.sample_stride],float)
   D=np.asarray(data[label+'__D'][::args.sample_stride],float)
   G=np.asarray(data[label+'__G'][::args.sample_stride],float)
   if len(q)<30 or np.any(np.diff(tt)<=0) or not np.isfinite(q).all():raise ValueError('Invalid historical 5D sequence '+label)
   chkD=float(np.max(np.abs(.5*np.sum(q*q,axis=1)-D)))
   chkG=float(np.max(np.abs(np.sum(q*dq,axis=1)-G)))
   if max(chkD,chkG)>4e-12:raise ValueError('5D D/G source identity failure '+label)
   # Historical Q2=0 crossings: four actual points, all by definition Q2~0.
   if label==args.primary:
    n0=np.array([0.,1.,0,0,0]);ev0=crossing_events(tt,q,dq,n0,np.zeros(5),0.)
    four_zeros=[(v['q1'],v['q2']) for v in ev0]
    negative=dict(status='HISTORICAL_ZERO_LINE_NEGATIVE_CONTROL',n_crossings=len(ev0),
      note='These roots are all Q2≈0; they cannot reconstruct nonzero user-observed slope. They are not the four user points.',
      events=[dict(tau=v['tau'],Q1=v['q1'],Q2=v['q2']) for v in ev0])
    if len(four_zeros)==4:
     fx=np.array([p[0] for p in four_zeros]);fy=np.array([p[1] for p in four_zeros]);negative['four_point_fit']=fit_four_points(fx,fy)
   for hidden_name,h in hidden_axes.items():
    z=q[:,2:5]@h
    for slope_case,sm in slope_controls(m):
     # Deliberate intercept: only the candidate geometry uses it; controls through zero.
     sb=b if slope_case=='observed_cartesian_m' else 0.
     u,v=coordinates(q,sm,sb)
     # Shape axes result from the observed projection itself, not hidden true metric.
     for theta in sorted(set(THETA_GRID+(args.theta_deg,))):
      if theta==0:wp=v;wm=v
      else:
       a=math.radians(theta);wp=math.cos(a)*v+math.sin(a)*z;wm=math.cos(a)*v-math.sin(a)*z
      w=window_stats(tt,wp,wm)
      # Shape/ellipse anisotropy needs covariance only; a control on apparent unrolling.
      cov=np.cov(np.stack((u,wp)));eigval=np.linalg.eigvalsh(cov)
      aniso=math.sqrt(max(float(eigval[-1]),1e-30)/max(float(eigval[0]),1e-30))
      row=dict(label=label,hidden_direction=hidden_name,slope_case=slope_case,slope=float(sm),intercept=float(sb),
        slope_angle_deg=math.degrees(math.atan(sm)),theta_deg=float(theta),
        opposite_side_duration_fraction=w['opposite_side_duration_fraction'],opposite_side_runs=w['opposite_side_runs'],
        plus_projected_axis_ratio=aniso)
      rows.append(row)
      if hidden_name=='PCA_train' and slope_case=='observed_cartesian_m' and abs(theta-args.theta_deg)<1e-9:
       e,n,npb,nmb=normal_frames(sm,h,theta)
       anchor=np.array([0.,sb,0.,0.,0.]);pe=crossing_events(tt,q,dq,npb,anchor,sm);me=crossing_events(tt,q,dq,nmb,anchor,sm)
       for nm,es in (('plus',pe),('minus',me)):
        for ev in es:events.append(dict(label=label,view=nm,theta_deg=(theta if nm=='plus' else -theta),**ev))
        for ev in serial_pairs(es):serial.append(dict(label=label,view=nm,**ev))
       for ev in nearest_pairs(pe,me):pairs.append(dict(label=label,**ev))
       for i,ww in enumerate(w['windows']):windows.append(dict(label=label,window_index=i+1,**ww))
       metrics=shape_and_section(tt,q,sm,sb,h,theta)
       rr=dict(label=label,plus_crossings=len(pe),minus_crossings=len(me),cross_view_pairs=len(nearest_pairs(pe,me)),
               opposite_side_duration_fraction=w['opposite_side_duration_fraction'],gate_windows=w['opposite_side_runs'],
               **metrics)
       per_trace.append(rr)
       if label==args.primary:primary_record=dict(tau=tt,q=q,h=h)
  csv_write(out/'nodal_hinge_slope_and_tilt_sensitivity.csv',rows)
  csv_write(out/'nodal_hinge_section_crossings_5d.csv',events)
  csv_write(out/'nodal_hinge_cross_view_posture_pairs.csv',pairs)
  csv_write(out/'nodal_hinge_serial_posture_transport.csv',serial)
  csv_write(out/'paired_observer_gate_windows.csv',windows)
  csv_write(out/'per_trajectory_geometry_diagnostics.csv',per_trace)
  dump(out/'four_point_reconstruction.json',four)
  dump(out/'historical_four_Q2_zero_negative_control.json',negative)
  plot_views(out,args.primary,primary_record,m,b,args.theta_deg,points)
  plot_sensitivity(out,rows,args.primary,m,args.theta_deg)
  plot_four_points(out,points,m,b)
  summary=dict(version='v18.52n',title='Four-Point Nodal-Axis Reconstruction & Paired-Hinge Observer Audit',
       stage='AUTHENTIC_FULL_POSTPROCESS' if args.sample_stride==1 else 'DEVELOPMENT_SMOKE_SUBSAMPLED_NOT_SCIENTIFIC_RESULT',
       source=src,theta_deg=args.theta_deg,observed_m_reference=args.slope,
       observed_m_angle_deg=math.degrees(math.atan(args.slope)),
       tan_of_tilt_deg=math.tan(math.radians(-args.theta_deg)),
       difference_tilt_vs_slope_angle_deg=math.degrees(math.atan(args.slope))+args.theta_deg,
       used_nodal_slope=m,used_nodal_intercept=b,
       used_nodal_slope_angle_deg=math.degrees(math.atan(m)),
       axis_source=('FOUR_POINT_AFFINE_OLS_EXPLICIT_USER_ADOPTION' if args.adopt_four_point_fit else 'USER_PRESPECIFIED_SLOPE_FIXED_INTERCEPT'),
       four_point_reconstruction=four,historical_four_zero_negative_control=negative,
       pca_hidden=ph.tolist(),pca_first_half_eigenvalues=eig,
       pca_heldout_angle_deg=principal_ang,
       lineage_v18_52m_pca_angle_mismatch_deg=baseline_angular_mismatch,
       primary=args.primary,trajectories=per_trace,
       controls=dict(slope_cases=[dict(name=name,m=s) for name,s in slope_controls(m)],
        theta_grid_deg=sorted(set(THETA_GRID+(args.theta_deg,))),
        hidden_directions={name:h.tolist() for name,h in hidden_axes.items()},
        basis_orthogonality_and_synthetic=tests,physical_D_G_identities=True,
        parent_l_source_sha_checked=True,m_server_pca_genealogy_checked=True),
       interpretive_limits=[
        'Four user-observed point coordinates are absent unless supplied with --four-points; m=-0.308 is not independently re-fitted.',
        'The four Q2=0 roots are NOT the same four visually observed points and cannot generate a nonzero Q2/Q1 slope.',
        'A line Q2=mQ1+b in the 2D projection is a one-dimensional *candidate hinge in R5*, not necessarily material or observer-invariant.',
        'Synchronous +/- observer projections algebraically generate sum/difference identities: complementary signs do not certify physical lobes.',
        'Apparent eccentricity/axis ratio of a projection is view-dependent; true in-plane rigid rotation leaves ellipse eigenvalues unchanged.',
        'No observer-derived single path winding in R5 is a topological invariant.',
        'Four-point line fit uses only four observations, thus is sensitive to how those points were chosen and to coordinate uncertainty.',
        'Cross-view event pairs occur at different times, so never establish exact recurrence; crossings are piecewise-linear interpolants.',
        'Cone/tube character would need independent geometry in full R5 and robust sections across observer changes.',
        'No new numerical integration, no physical two-filament model, and no topology certificate.'])
  dump(out/'report.json',summary)
  # compact scientific report with exact status, not a speculative promotion.
  avg=np.mean([x['opposite_side_duration_fraction'] for x in per_trace]);nr=np.mean([x['gate_windows'] for x in per_trace]);
  parent_m_frac=m_report['results'][args.primary]['opposite_side_duration_fraction']
  target=[r for r in rows if r['label']==args.primary and r['hidden_direction']=='PCA_train'
          and r['slope_case']=='observed_cartesian_m' and abs(r['theta_deg']-args.theta_deg)<1e-8][0]
  base=[r for r in rows if r['label']==args.primary and r['hidden_direction']=='PCA_train'
        and r['slope_case']=='historical_Q1_axis' and abs(r['theta_deg']-args.theta_deg)<1e-8][0]
  txt=[f'# Q550 v18.52n — Server scientific report','',
    '**Four-Point Nodal-Axis Reconstruction & Paired-Hinge Observer Audit**','',
    f'- Stage: `{summary["stage"]}`; no new ODE integration.',
    f'- Authentic v18.52l 5D trace SHA-256: `{src["sha256"]}`.',
    f'- Nodal user hypothesis: m = {args.slope:+.9f}, theta_slope = {math.degrees(math.atan(args.slope)):+.6f} deg.',
    f'- Separate observer tilt: +/-{args.theta_deg:.5f} deg; tan(-tilt)={math.tan(math.radians(-args.theta_deg)):+.9f}.',
    f'- Angle magnitude difference: {abs(summary["difference_tilt_vs_slope_angle_deg"]):.6f} deg.',
    f'- Axis source: `{summary["axis_source"]}`; line intercept {b:.5e}.',
    f'- Four-point coordinates status: `{four["status"]}`.',
    f'- Historical Q2=0 roots: {negative["n_crossings"]}; they are a negative control and not the user four.',
    f'- Ancestor v18.52m PCA-axis agreement: {baseline_angular_mismatch:.6g} deg.',
    f'- Hidden PCA train/heldout angle: {principal_ang:.6f} deg.','',
    '## Crucial one-trajectory observer comparison','',
    f'- Parent v18.52m primary Q1-axis gate duration fraction: {parent_m_frac:.7f}.',
    f'- Reproduced original Q1-axis gate fraction in current frame: {base["opposite_side_duration_fraction"]:.7f} (piecewise-linear exact gate intervals; v18.52m used sampled trapezoidal integration).',
    f'- Proposed sloped nodal axis gate fraction: **{target["opposite_side_duration_fraction"]:.7f}**.',
    f'- Proposed nodal hinge average gate fraction over eight trajectories: {avg:.7f}; mean gate windows: {nr:.3f}.',
    '','| Label | Opposite-side fraction | Windows | + crossings | - crossings | Cone R² (diagnostic) |',
    '|---|---:|---:|---:|---:|---:|']
  for r in per_trace:
   txt.append(f'| {r["label"]} | {r["opposite_side_duration_fraction"]:.7f} | {r["gate_windows"]} | {r["plus_crossings"]} | {r["minus_crossings"]} | {r["cone_linear_radius_vs_axis_R2"]:.4f} |')
  txt+=['','## Scientific decision','',
   'There is **no reconstructed four-point user line** without actual four coordinate pairs; the user-supplied slope remains a prior hypothesis.',
   'The nodal axis and the separate angular hinge tilt are now correctly distinguished and tested against the v18.52m Q1-axis control.',
   'Opposite-side gate fractions are observer diagnostics, not evidence of two physical filaments or invariant section topology.',
   'The intrinsic hinge-perpendicular radius can vary even where the projected curves look paired.','',
   '## Data lineage and limits','',
   'The entire analysis uses a **single physical 5D trajectory per kick** in matched material time, eight inherited kick labels.',
   'The ancestor v18.52l vector SHA and v18.52m PCA direction were checked. No high-cost nonlinear integration was performed.',
   'Four-point coordinates may be supplied by --four-points; explicit --adopt-four-point-fit replaces the candidate slope and intercept.',
   'Curvature, winding, ellipses, cones and tubes are not promoted to observer-independent dynamical/topological invariants.','']
  (out/'SCIENTIFIC_REPORT_SERVER.md').write_text('\n'.join(txt),encoding='utf8')
 print('v18.52n AUDIT PASS',len(LABELS),'trajectories;',len(tt),'samples each;',
       'm=',m,'b=',b,'angle=',args.theta_deg,'four_point=',four['status'],'output=',out)

if __name__=='__main__':main()
