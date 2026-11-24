#!/usr/bin/env python3
"""Q550 v18.53b: fixed-charge relative-equilibrium branch and spectral boundary.

The authentic archived 11D physical RHS is the authority whenever available.
The stand-alone reduced expression is independently checked, never promoted
from offline algebra into a certification of the genuine source. No IVPs.
"""
from __future__ import annotations
import argparse,csv,hashlib,importlib.util,json,math,sys,traceback
from pathlib import Path
import numpy as np
from scipy.optimize import root,least_squares
from scipy.optimize import linear_sum_assignment
from energy_math import Params,sigmoid,dU
from return_map import lift8,wrap,charge

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT.parent;OUT=ROOT/'results'
C=json.loads((ROOT/'inputs_config.json').read_text());P=Params(**C['model_params']);SEEDS=('original','autonomous')

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for chunk in iter(lambda:f.read(1<<20),b''):h.update(chunk)
 return h.hexdigest()

def importer(path,name):
 spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec)
 sys.modules[name]=m;spec.loader.exec_module(m);return m

def emit(name,content):
 path=OUT/name;path.parent.mkdir(parents=True,exist_ok=True)
 if isinstance(content,str):path.write_text(content)
 else:path.write_text(json.dumps(content,indent=2,ensure_ascii=False,allow_nan=False)+'\n')

def csvout(name,rows):
 if not rows:return
 with (OUT/name).open('w',newline='') as f:
  fields=list(dict.fromkeys(key for item in rows for key in item));writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows(rows)

def extract_params(axis,value,anchor):
 J,k,s=anchor
 return {'charge':(float(value),k,s),'kappa':(J,float(value),s),'signed_drive':(J,k,float(value))}[axis]

def u_from_a(a,J):
 l1,l2,w1=map(float,a);r1,r2=math.exp(l1),math.exp(l2)
 w2=(r1*r1*w1-J)/(r2*r2)
 return np.asarray([l1,0.,w1,l2,0.,0.,0.,math.pi],float),w2

def radial_residual(a,J,k,s,delta=math.pi):
 if np.any(~np.isfinite(a)) or max(abs(float(a[0])),abs(float(a[1])))>16:return np.asarray([1e6,1e6,1e6],float)
 l1,l2,w1=map(float,a);r1,r2=math.exp(l1),math.exp(l2)
 w2=(r1*r1*w1-J)/(r2*r2)
 if not math.isfinite(w2) or abs(w2)>1e8:return np.array([1e6]*3,float)
 x=.5;c=math.cos(delta)
 return np.array([w1*w1-dU(r1,P.e_reserve)/r1-k*x*r2/r1*c,
  w2*w2-dU(r2,P.e_reserve)/r2-k*x*r1/r2*c,
  w1+w2-s*x],float)

def solve(a0,J,k,s):
 """Follow delta=pi sheet continuously. A failed corrector is not a fold."""
 a0=np.asarray(a0,float)
 z=root(lambda a:radial_residual(a,J,k,s),a0,method='hybr',tol=1e-11)
 best=np.asarray(z.x,float);e=float(np.max(np.abs(radial_residual(best,J,k,s))))
 if not np.isfinite(e) or e>C['root_tol']:
  try:
   q=least_squares(lambda a:radial_residual(a,J,k,s),np.clip(a0,[-11.,-11.,-20.],[11.,11.,20.]),
     bounds=([-11.,-11.,-20.],[11.,11.,20.]),max_nfev=125,xtol=3e-13,ftol=3e-13,gtol=3e-13)
   test=float(np.max(np.abs(radial_residual(q.x,J,k,s))))
   if test<e:best=np.asarray(q.x,float);e=test
  except (ValueError,OverflowError,FloatingPointError):pass
 good=bool(e<C['root_tol'] and np.linalg.norm(best-a0)<1.3 and math.exp(best[0])>C['positive_radius_floor'] and math.exp(best[1])>C['positive_radius_floor'])
 return good,best,e

def reduced_field(u,J,k,s,omega):
 """Reconstructed xi chart; on-server crosschecked against original RHS."""
 l1,v1,w1,l2,v2,z,vz,D=np.asarray(u,float)
 r1,r2=math.exp(l1),math.exp(l2);x=float(sigmoid(z));g=x*(1-x)
 w2=(r1*r1*w1-J)/r2**2;sin,cos=math.sin(D),math.cos(D)
 dr1,dr2=r1*v1,r2*v2;dx=g*vz
 a1=(r1*w1*w1-dU(r1,P.e_reserve)-k*x*r2*cos-P.gamma_rho*dr1)/r1-v1*v1
 a2=(r2*w2*w2-dU(r2,P.e_reserve)-k*x*r1*cos-P.gamma_rho*dr2)/r2-v2*v2
 dw1=k*x*r2*sin/r1-2*v1*w1
 xdd=(-P.barrier_B*(2*x-1)/(x*(1-x))+P.mu_x*r1*r2*sin-P.gamma_x*dx)/P.m_x
 dvz=xdd/g-(1-2*x)*vz*vz
 return (omega/x)*np.array([v1,a1,dw1,v2,a2,vz,dvz,w1+w2-s*x],float)

def jacfd(fn,x,step):
 x=np.asarray(x,float);out=np.zeros((len(x),len(x)))
 for j in range(len(x)):
  h=step*max(1.,abs(x[j]));d=np.zeros(len(x));d[j]=h
  out[:,j]=(fn(x+d)-fn(x-d))/(2*h)
 return out

def native_field(u,J,k,s,eng):
 """Exact staged 11D source at arbitrary phase, reduced on fixed J leaf."""
 a=np.asarray(u,float).copy();a[7]=a[7]+float(eng.params.phi0)
 yy=lift8(a,J,0.)
 FF=np.asarray(eng.physical.rhs_from_original(np.r_[yy,0.],float(k),eng.params,eng.slope,float(s),eng.support,eng.audit),float)
 if FF.shape!=(11,):raise RuntimeError('AUTHENTIC_RHS_SHAPE_NOT_11')
 return np.asarray([FF[0],FF[1],FF[3],FF[4],FF[5],FF[8],FF[9],FF[2]+FF[6]-FF[10]],float),FF

def correct_and_audit(a,J,k,s,omega,eng=None):
 u,w2=u_from_a(a,J);F=reduced_field(u,J,k,s,omega)
 hbig,hsmall=C['jacobian_steps'];M=jacfd(lambda v:reduced_field(v,J,k,s,omega),u,hsmall)
 spectrum=np.linalg.eigvals(M);absc=float(np.max(spectrum.real))
 # Local 3x3 equilibrium solvability is separate from reduced 8D stability.
 D=jacfd(lambda aa:radial_residual(aa,J,k,s),a,1e-5)
 sv=np.linalg.svd(D,compute_uv=False)
 native={'available':False};spec=mval=specerr=None
 if eng is not None:
  fn=lambda v:native_field(v,J,k,s,eng)[0]
  Fauth,_=native_field(u,J,k,s,eng)
  Mbig=jacfd(fn,u,hbig);Msmall=jacfd(fn,u,hsmall)
  jac_rel=float(np.linalg.norm(Mbig-Msmall)/max(1.,np.linalg.norm(Msmall)))
  source_rel=float(np.linalg.norm(Fauth-F)/max(1.,np.linalg.norm(Fauth),np.linalg.norm(F)))
  matrix_rel=float(np.linalg.norm(Msmall-M)/max(1.,np.linalg.norm(Msmall),np.linalg.norm(M)))
  eigsrc=np.linalg.eigvals(Msmall)
  ii,jj=linear_sum_assignment(np.abs(eigsrc[:,None]-spectrum[None,:]))
  eigerr=float(np.max(np.abs(eigsrc[ii]-spectrum[jj])))
  native={'available':True,'native_equilibrium_rhs_max':float(np.max(abs(Fauth))),
   'source_reconstructed_field_relative':source_rel,'source_reconstructed_jacobian_relative':matrix_rel,
   'native_jacobian_step_relative':jac_rel,'native_local_eigen_max_abs_difference':eigerr,
   'source_match':bool(source_rel<C['source_rhs_tol'] and matrix_rel<C['native_jacobian_rel_tol'] and jac_rel<C['native_jacobian_rel_tol']),
   'native_spectral_abscissa':float(np.max(eigsrc.real))}
  spectrum=eigsrc;absc=float(np.max(eigsrc.real))
  # Only native eigvals can certify a genuine Q550 stability statement.
  vals=[{'re':float(v.real),'im':float(v.imag)} for v in eigsrc]
 else:vals=[{'re':float(v.real),'im':float(v.imag)} for v in spectrum]
 phase_turn=None if abs(s)<1e-10 else 2*math.pi/(omega*abs(s))
 return {'root_8d_max_abs':float(np.max(abs(F))),'singular_values_equilibrium_3x3':[float(x) for x in sv],
   'smallest_singular_equilibrium_3x3':float(sv[-1]),'spectral_abscissa':absc,
   'max_real_mode_imag_abs':float(max((abs(v.imag) for v in spectrum if abs(v.real-absc)<1e-6),default=0.)),
   'n_eigen_positive_real':int(sum(spectrum.real>1e-7)),
   'phase_turn_xi':phase_turn,'one_turn_spectral_radius':None if phase_turn is None else float(math.exp(min(300.,absc*phase_turn))),
   'eigenvalues':vals,'native':native,'charge_reconstruction_error':0.}

def make_row(axis,val,J,k,s,a,eps,qa,seed='common'):
 u,w2=u_from_a(a,J)
 return dict(axis=axis,parameter=float(val),J=float(J),kappa=float(k),signed_drive=float(s),
  r1=math.exp(float(a[0])),r2=math.exp(float(a[1])),omega1=float(a[2]),omega2=float(w2),
  x=0.5,delta=math.pi,root_residual=eps,smallest_equilibrium_singular=qa['smallest_singular_equilibrium_3x3'],
  spectral_abscissa=qa['spectral_abscissa'],max_real_mode_imag_abs=qa['max_real_mode_imag_abs'],
  n_positive_real_modes=qa['n_eigen_positive_real'],phase_turn_xi=qa['phase_turn_xi'],
  one_turn_spectral_radius=qa['one_turn_spectral_radius'],native_checked=qa['native']['available'],
  native_source_match=qa['native'].get('source_match',False),
  native_rhs_max=qa['native'].get('native_equilibrium_rhs_max'),
  native_jacobian_step_relative=qa['native'].get('native_jacobian_step_relative'))

def scan_axis(axis,anchor,a0,omega,eng):
 c=C['branch_scans'][axis]
 values=np.unique(np.r_[np.linspace(c['lo'],c['hi'],int(c['n'])),anchor[('charge','kappa','signed_drive').index(axis)]])
 refval=anchor[('charge','kappa','signed_drive').index(axis)]
 records=[];details=[];fails=[]
 # Two one-sided continuations avoid jumping across the anchor into a different sheet.
 for direction in [-1,1]:
  seq=sorted((v for v in values if (v-refval)*direction>1e-13),reverse=(direction==-1))
  prev=a0.copy();prevvalue=refval;prevprev=None
  for val in seq:
   J,k,s=extract_params(axis,val,anchor)
   if prevprev is None:estimate=prev.copy()
   else:
    old_a,old_value=prevprev
    estimate=prev+(prev-old_a)*(val-prevvalue)/(prevvalue-old_value)
    if np.linalg.norm(estimate-prev)>0.7:estimate=prev.copy()
   success,a,err=solve(estimate,J,k,s)
   if not success:
    success,a,err=solve(prev,J,k,s)
   if not success:
    fails.append({'axis':axis,'direction':direction,'value':float(val),'status':'CORRECTOR_FAILED_NOT_PROVEN_FOLD','last_accepted_value':float(prevvalue),'last_residual':err})
    break
   qa=correct_and_audit(a,J,k,s,omega,eng)
   if qa['root_8d_max_abs']>C['root_tol']*20:
    fails.append({'axis':axis,'direction':direction,'value':float(val),'status':'REDUCED_EQUILIBRIUM_RHS_NOT_ZERO','last_residual':qa['root_8d_max_abs']});break
   records.append(make_row(axis,val,J,k,s,a,err,qa))
   details.append({'axis':axis,'parameter':float(val),'root_a':a.tolist(),'analytic_diagnostics':qa})
   prevprev=(prev.copy(),prevvalue);prev=a.copy();prevvalue=val
 # anchor once
 J,k,s=anchor
 qa=correct_and_audit(a0,J,k,s,omega,eng)
 records.append(make_row(axis,refval,J,k,s,a0,float(np.max(abs(radial_residual(a0,J,k,s)))),qa))
 details.append({'axis':axis,'parameter':float(refval),'root_a':a0.tolist(),'analytic_diagnostics':qa})
 records.sort(key=lambda x:x['parameter']);details.sort(key=lambda x:x['parameter'])
 return records,details,fails

def pseudo_arclength_drive(ordinary,anchor,omega,eng,steps=60,h=0.02):
 """Follow the same root sheet beyond a signed-drive parameter turning point.

 4D augmented predictor-corrector: (three equilibrium eqns, one secant gauge).
 This is a numerical fold *candidate*, not a bifurcation certificate.
 """
 last=sorted(ordinary,key=lambda r:r['signed_drive'])[-2:]
 if len(last)<2:return [],[],[],{'status':'INSUFFICIENT_ORDINARY_BRANCH'}
 def make_z(r):return np.array([math.log(r['r1']),math.log(r['r2']),r['omega1'],r['signed_drive']],float)
 prev,curr=map(make_z,last)
 J,k,s0=anchor
 records=[];details=[];fail=[]
 for i in range(int(steps)):
  tan=curr-prev;norm=float(np.linalg.norm(tan))
  if norm<1e-14:fail.append({'step':i,'status':'ZERO_ARCLENGTH_SECANT'});break
  tan/=norm;pred=curr+h*tan
  def aug(z):return np.r_[radial_residual(z[:3],J,k,z[3]),np.dot(z-pred,tan)]
  try:sol=root(aug,pred,method='hybr',tol=1e-11)
  except (OverflowError,FloatingPointError,ValueError) as e:
   fail.append({'step':i,'status':'PSEUDOARCLENGTH_EXCEPTION','message':str(e)});break
  z=np.asarray(sol.x,float);err=float(np.max(abs(aug(z))))
  if not np.isfinite(err) or err>2e-8:
   fail.append({'step':i,'status':'PSEUDOARCLENGTH_CORRECTOR_NOT_CONVERGED','residual':err});break
  kk=float(z[3]);a=z[:3]
  qa=correct_and_audit(a,J,k,kk,omega,eng)
  if qa['root_8d_max_abs']>C['root_tol']*20:
   fail.append({'step':i,'status':'PSEUDOARCLENGTH_RHS_NOT_ZERO'});break
  row=make_row('signed_drive_pseudoarclength',kk,J,k,kk,a,err,qa)
  row['pseudo_step']=i;row['tangent_signed_drive']=float(tan[3]);row['arc_step']=h
  row['root_corrector_residual_augmented']=err
  records.append(row)
  details.append({'pseudo_step':i,'parameter':kk,'root_a':a.tolist(),'analytic_diagnostics':qa})
  prev,curr=curr,z
 summary={'status':'NO_TURN_OBSERVED','n_steps':len(records),'n_failed_steps':len(fail)}
 if records:
  ss=[last[-1]['signed_drive']]+[r['signed_drive'] for r in records]
  maxidx=int(np.argmax(ss));turn=0<maxidx<len(ss)-1
  summary.update({'max_signed_drive_including_tip':float(max(ss)),
    'max_signed_drive_position_in_sequence':maxidx,
    'parameter_turn_observed':bool(turn),
    'min_equilibrium_3x3_singular':float(min(r['smallest_equilibrium_singular'] for r in records)),
    'min_spectral_abscissa':float(min(r['spectral_abscissa'] for r in records)),
    'max_spectral_abscissa':float(max(r['spectral_abscissa'] for r in records)),
    'native_all_checked':bool(all(r['native_checked'] and r['native_source_match'] for r in records)),
    'status':('NUMERICAL_PSEUDOARCLENGTH_FOLD_AND_REAL_STABILITY_CROSSING_CANDIDATE' if turn and min(r['spectral_abscissa'] for r in records)<0<max(r['spectral_abscissa'] for r in records) else
      'PSEUDOARCLENGTH_BRANCH_CONTINUED_WITHOUT_FULL_FOLD_CROSSING')})
 return records,details,fail,summary

def refine_algebraic_fold(pseudo,anchor,omega,eng):
 """Solve equilibrium F=0 and det(D_a F)=0; quantify local saddle-node conditions.

 This is a high-accuracy *numerical* normal-form diagnostic, not interval proof.
 """
 if len(pseudo)<3:return {'status':'INSUFFICIENT_PSEUDO_POINTS'}
 near=min(pseudo,key=lambda r:r['smallest_equilibrium_singular'])
 J,k,_=anchor
 z0=np.array([math.log(near['r1']),math.log(near['r2']),near['omega1'],near['signed_drive']],float)
 def algebra(z):return radial_residual(z[:3],J,k,z[3])
 def a_jac(z,h):return jacfd(lambda a:radial_residual(a,J,k,z[3]),z[:3],h)
 def aug(z):return np.r_[algebra(z),np.linalg.det(a_jac(z,1e-5))]
 rr=root(aug,z0,method='hybr',tol=1e-11)
 z=np.asarray(rr.x,float);error=float(np.max(abs(aug(z))))
 if error>1e-6:return {'status':'REFINED_FOLD_CORRECTOR_NOT_CONVERGED','residual':error,'initial_z':z0.tolist()}
 U,sv,Vh=np.linalg.svd(a_jac(z,2e-6));w=U[:,-1];v=Vh[-1,:]
 eps=1e-5
 Fs=(radial_residual(z[:3],J,k,z[3]+eps)-radial_residual(z[:3],J,k,z[3]-eps))/(2*eps)
 trans=float(np.dot(w,Fs));quad=[]
 for h in [1e-3,3e-4,1e-4]:
  q=(radial_residual(z[:3]+h*v,J,k,z[3])+radial_residual(z[:3]-h*v,J,k,z[3])-2*algebra(z))/(h*h)
  quad.append(float(np.dot(w,q)))
 qa=correct_and_audit(z[:3],J,k,float(z[3]),omega,eng)
 eig=np.array([complex(t['re'],t['im']) for t in qa['eigenvalues']]);order=np.argsort(abs(eig))
 eig0=eig[order[0]];other=eig[order[1:]]
 authentic=qa['native'].get('available',False) and qa['native'].get('source_match',False)
 checks={'root_augmented_max_abs_lt_1e_7':error<1e-7,
  '3x3_null_singular_value_lt_1e_6':float(sv[-1])<1e-6,
  'parametric_transversality_nonzero':abs(trans)>1e-4,
  'quadratic_fold_curvature_nonzero':min(abs(q) for q in quad)>1e-4,
  'quadratic_stencil_consistent':abs(quad[1]-quad[2])/max(1.,abs(quad[2]))<1e-3,
  '8d_one_near_zero_real_mode':abs(eig0)<1e-5 and abs(eig0.imag)<1e-5,
  'other_7_modes_left_half_plane':bool(np.max(other.real)<-1e-4),
  'authentic_source_match':bool(authentic)}
 checks={key:bool(value) for key,value in checks.items()}
 return {'status':('NUMERICALLY_NONDEGENERATE_LOCAL_SADDLE_NODE_CANDIDATE__NOT_INTERVAL_CERTIFICATE' if all(v for kk,v in checks.items() if kk!='authentic_source_match') else 'INCOMPLETE_FOLD_NONDGENERACY'),
  'J':J,'kappa':k,'s_fold':float(z[3]),'core_a':z[:3].tolist(),
  'r1':float(math.exp(z[0])),'r2':float(math.exp(z[1])),
  'root_augmented_residual_max':error,'jacobian_equilibrium_singular_values':[float(x) for x in sv],
  'parametric_transversality_left_null_dot_Fs':trans,'quadratic_fold_curvature_left_null_dot_Faa_vv_by_stencil':quad,
  'smallest_8d_eigenvalue':{'real':float(eig0.real),'imag':float(eig0.imag)},
  'other_7_max_real':float(max(other.real)),'conditions':checks,
  'authentic_original_rhs_at_fold':qa['native'],
  'classification':'NUMERICAL_LOCAL_SADDLE_NODE_NONDGENERACY_DIAGNOSTIC_NOT_GLOBAL_PROOF'}

def find_candidates(axis,records):
 out=[]
 for a,b in zip(records[:-1],records[1:]):
  r1,r2=a['spectral_abscissa'],b['spectral_abscissa']
  if not (np.isfinite(r1) and np.isfinite(r2)):continue
  cross=bool((r1<-C['eigen_zero_real_tol'] and r2>C['eigen_zero_real_tol']) or
    (r2<-C['eigen_zero_real_tol'] and r1>C['eigen_zero_real_tol']))
  borderline=(min(abs(r1),abs(r2))<C['eigen_zero_real_tol'])
  fold=(min(a['smallest_equilibrium_singular'],b['smallest_equilibrium_singular'])<C['singular_value_warning'])
  if cross or borderline or fold:
   approx=(a['parameter']-r1*(b['parameter']-a['parameter'])/(r2-r1) if cross and r2!=r1 else .5*(a['parameter']+b['parameter']))
   out.append({'axis':axis,'from':a['parameter'],'to':b['parameter'],
    'spectral_abscissa_from':r1,'spectral_abscissa_to':r2,
    'interpolated_stability_boundary_parameter':float(approx),
    'spectral_sign_change':cross,'near_imag_axis':borderline,'small_equilibrium_jacobian_singular':fold,
    'classification':'NUMERICAL_LOCAL_BOUNDARY_CANDIDATE_NOT_BIFURCATION_THEOREM'})
 return out

def plot_points(rows,candidates):
 try:
  import matplotlib
  matplotlib.use('Agg')
  import matplotlib.pyplot as plt
 except ImportError:return
 fig,ax=plt.subplots(3,2,figsize=(12,10),constrained_layout=True)
 for i,axis in enumerate(('charge','kappa','signed_drive')):
  data=[r for r in rows if r['axis']==axis]
  if not data:continue
  x=[r['parameter'] for r in data]
  ax[i,0].plot(x,[r['r1'] for r in data],label='r1',lw=1.6)
  ax[i,0].plot(x,[r['r2'] for r in data],label='r2',lw=1.6)
  ax[i,0].set(title=axis+' — radial branch',xlabel=axis,ylabel='r',yscale='log');ax[i,0].legend()
  ax[i,1].plot(x,[r['spectral_abscissa'] for r in data],lw=1.5)
  ax[i,1].axhline(0,color='black',lw=.8,ls='--')
  ax[i,1].set(title=axis+' — 8D spectral abscissa',xlabel=axis,ylabel='max Re(lambda)')
  if axis=='signed_drive':
   ps=[r for r in rows if r['axis']=='signed_drive_pseudoarclength']
   if ps:
    px=[r['signed_drive'] for r in ps]
    ax[i,0].plot(px,[r['r1'] for r in ps],ls='--',lw=1.5,label='pseudo-arclength r1')
    ax[i,0].plot(px,[r['r2'] for r in ps],ls='--',lw=1.5,label='pseudo-arclength r2')
    ax[i,1].plot(px,[r['spectral_abscissa'] for r in ps],ls='--',lw=1.5,label='post-turn branch')
    ax[i,0].legend(fontsize=8);ax[i,1].legend(fontsize=8)
  for q in candidates:
   if q['axis']==axis:ax[i,1].axvspan(q['from'],q['to'],alpha=.15)
 fig.suptitle('Q550 v18.53b — relative equilibrium continuation (fixed J for kappa/drive)',fontsize=13)
 fig.savefig(OUT/'branch_radius_and_spectral_abscissa.png',dpi=150);plt.close(fig)

def main():
 ap=argparse.ArgumentParser()
 ap.add_argument('--parent-v1851s',type=Path)
 ap.add_argument('--require-parent',action='store_true')
 ap.add_argument('--offline-only',action='store_true',help='Algebraic screening; NOT genuine generator certification')
 args=ap.parse_args();OUT.mkdir(exist_ok=True,parents=True)
 ref=ROOT/'reference/v18_53a_authenticated_report.json'
 if sha(ref)!=C['v1853a_report_sha256']:raise RuntimeError('V1853A_AUTHENTICATED_REPORT_SHA_MISMATCH')
 report53=json.loads(ref.read_text())
 if report53['status']!=C['expected_v1853a_status']:raise RuntimeError('V1853A_STATUS_UNEXPECTED')
 A=report53['results']['original'];B=report53['results']['autonomous']
 if not report53['cross_seed']['equivalent_on_fixed_charge_and_phase_frame']:raise RuntimeError('BASE_SEEDS_NOT_EQUIVALENT')
 J0=float(A['J']);k0=float(A['kappa']);s0=-1.
 anch=tuple([J0,k0,s0]);omega=0.844357915384206
 a0=np.asarray(A['algebraic_branch_root']['core_root'],float)
 if float(np.max(abs(radial_residual(a0,*anch))))>1e-9:raise RuntimeError('ANCHOR_RADIAL_RESIDUAL_TOO_HIGH')
 if np.max(abs(np.asarray(A['equilibrium_reduced_chart'])-np.asarray(B['equilibrium_reduced_chart'])))>1e-8:
  raise RuntimeError('INTRINSIC_SEED_ANCHOR_MISMATCH')
 rep={'version':'v18.53b','scientific_decision':'OFFLINE_ANALYTIC_BRANCH_SCREENING_ONLY__AUTHENTIC_PARENT_REQUIRED',
  'authenticated_53a_anchor_sha256':sha(ref),'parent_v1851s_verified':False,
  'anchor':{'J':J0,'kappa':k0,'signed_drive':s0,'r1':math.exp(a0[0]),'r2':math.exp(a0[1]),'omega1':a0[2]},
  'full_historical_ramp_bifurcation_certified':False,'global_bifurcation_theorem_certified':False,
  'two_seed_branches_independent':False,'axes':{},'warnings':[]}
 eng=None;temp=None;source={}
 try:
  if not args.offline_only:
   stage=importer(ROOT/'reference/v18_52v_r1_run_audit.py','stage_q550_53b')
   parent=stage.find_parent(args.parent_v1851s)
   if parent is not None:
    if sha(parent)!=C['v1851s_sha256']:raise RuntimeError('GENEALOGICAL_PARENT_SHA_MISMATCH')
    temp,eng,_,ks,ctx,source=stage.stage_parent(parent,np.asarray([12.8,12.8001],float))
    params={n:float(getattr(eng.params,n)) for n in C['model_params']}
    if any(abs(params[n]-v)>1e-12 for n,v in C['model_params'].items()):raise RuntimeError('AUTHENTIC_PARAM_MISMATCH')
    omega=float(eng.omega)
    rep['parent_v1851s_verified']=True
    rep['source_provenance']={'parent_sha256':sha(parent),'parent_file':str(parent),
     'original_engine_omega':omega,'model_parameters':params,
     'signed_rhs_source':source['rhs_from_original'],
     'historical_clock_source':source['cubic_clock_sign']}
   elif args.require_parent:raise RuntimeError('AUTHENTIC_V1851S_PARENT_REQUIRED_FOR_BRANCH_CERTIFICATION')
  elif args.require_parent:raise RuntimeError('CONTRADICTORY_REQUIRE_PARENT_AND_OFFLINE_ONLY')
  rows=[];cands=[];fails=[];details=[]
  for axis in ('charge','kappa','signed_drive'):
   rr,dd,ff=scan_axis(axis,anch,a0,omega,eng)
   cc=find_candidates(axis,rr)
   rows+=rr;details+=dd;cands+=cc;fails+=ff
   rep['axes'][axis]={'n_points':len(rr),'n_failed_correctors':len(ff),
    'n_stable':sum(r['spectral_abscissa'] < -C['eigen_zero_real_tol'] for r in rr),
    'n_unstable':sum(r['spectral_abscissa'] > C['eigen_zero_real_tol'] for r in rr),
    'n_near_neutral':sum(abs(r['spectral_abscissa']) <= C['eigen_zero_real_tol'] for r in rr),
    'n_native_source_mismatch':sum(not r['native_source_match'] for r in rr if r['native_checked']),
    'radial_range_r1':[min(r['r1'] for r in rr),max(r['r1'] for r in rr)],
    'radial_range_r2':[min(r['r2'] for r in rr),max(r['r2'] for r in rr)],
    'stability_boundary_candidate_intervals':cc,
    'continued_range':[rr[0]['parameter'],rr[-1]['parameter']]}
   print(axis,'points',len(rr),'stability boundary brackets',len(cc),'source mismatches',rep['axes'][axis]['n_native_source_mismatch'],flush=True)
   emit('report.json',rep)
  # Continue through ordinary signed-drive parameter failure with pseudo-arclength.
  pseudo,pseudod,pseudof,psummary=pseudo_arclength_drive([r for r in rows if r['axis']=='signed_drive'],anch,omega,eng)
  rep['pseudoarclength_signed_drive']=psummary
  fold_refined=refine_algebraic_fold(pseudo,anch,omega,eng)
  rep['refined_local_fold']=fold_refined
  emit('refined_fold_nondegeneracy.json',fold_refined)
  rows+=pseudo;details+=pseudod;fails+=pseudof
  csvout('pseudo_arclength_drive.csv',pseudo)
  emit('pseudo_arclength_fold_candidates.json',{'summary':psummary,'correction_failures':pseudof})
  print('pseudo-arclength',psummary['status'],'max s',psummary.get('max_signed_drive_including_tip'),flush=True)
  csvout('branch_continuation.csv',rows)
  emit('branch_detailed_spectrum.json',details)
  emit('stability_boundary_candidates.json',cands)
  emit('corrector_stops_not_proven_folds.json',fails)
  plot_points(rows,cands)
  authentic=eng is not None
  checks=authentic and all(r['native_source_match'] and r['native_rhs_max']<C['source_rhs_tol'] for r in rows)
  rep['total_points']=len(rows);rep['candidate_boundary_count']=len(cands)
  rep['pseudo_points']=len(pseudo)
  rep['source_checked_all_equilibrium_points']=bool(checks)
  rep['source_checked_refined_fold']=bool(fold_refined.get('conditions',{}).get('authentic_source_match',False))
  checks=checks and rep['source_checked_refined_fold']
  rep['corrector_stops']=len(fails)
  rep['scientific_decision']=('AUTHENTIC_FIXED_CHARGE_RELATIVE_EQUILIBRIUM_BRANCHES_WITH_LOCAL_SPECTRAL_BOUNDARY_DIAGNOSTICS__NO_GLOBAL_BIFURCATION_THEOREM'
      if checks else ('AUTHENTIC_SOURCE_PRESENT_BUT_SOME_BRANCH_OR_JACOBIAN_CHECKS_FAILED' if authentic else
      'OFFLINE_ANALYTIC_BRANCH_SCREENING_ONLY__AUTHENTIC_PARENT_REQUIRED'))
  summary=[f'# Q550 v18.53b — {C["title"]}','','**'+rep['scientific_decision']+'**','',
   'The physical source has been checked at each root only if the authenticated v18.51s parent was available.','','## Branch scans']
  for ax,d in rep['axes'].items():
   summary+=['',f'### {ax}',f'- Continuation points: {d["n_points"]}',
    f'- Stable reduced 8D: {d["n_stable"]}; unstable: {d["n_unstable"]}',
    f'- Candidate stability-boundary intervals: {len(d["stability_boundary_candidate_intervals"])}',
    f'- Corrector stops: {d["n_failed_correctors"]} (not proved folds)']
  summary+=['','## Interpretation','',
   'Continuation in charge compares different invariant charge leaves, with stability always computed tangent to each leaf.',
   'Continuation in kappa or constant signed-drive holds charge fixed. A signed drive value is a counterfactual constant protocol parameter, not an actual reconstructed ramp trajectory.',
   'A local spectral sign crossing does not alone prove a Hopf, pitchfork, saddle-node or global bifurcation. A failed corrector is not evidence of a fold.',
   'No Poincare orbit was reintegrated here. The inherited v18.53a one-turn Floquet match applies only to its authenticated anchor.',
   'At s=0 the external phase does not turn: no one-phase-turn Floquet multiplier is defined.',
   'The two seed offsets are gauge/phase equivalent at the anchor and not separately continued as different physical branches.',
   'The signed-drive pseudo-arclength branch is followed through the apparent turn; its numerical parameter reversal and real eigenvalue crossing are candidates for a saddle-node, NOT a nondegenerate bifurcation certificate.',
   'Pseudo-arclength candidate: '+json.dumps(psummary),
   'Refined local fold diagnostic (not interval certificate): '+json.dumps(fold_refined)]
  emit('SCIENTIFIC_REPORT_SERVER.md','\n'.join(summary)+'\n');emit('report.json',rep)
  if args.require_parent and not checks:raise RuntimeError('STRICT_AUTHENTIC_BRANCH_SOURCE_CHECK_FAILED')
 finally:
  if temp is not None:temp.cleanup()

if __name__=='__main__':
 try:main()
 except Exception as exc:
  OUT.mkdir(exist_ok=True,parents=True);emit('audit_failure.txt',traceback.format_exc());print('V18_53B_FAILURE',repr(exc),flush=True);sys.exit(1)
