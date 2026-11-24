#!/usr/bin/env python3
"""Q550 v18.52u: audited orientation genealogy of external phase and 11D conjugacy.

No trajectory integration. External phase from authenticated v18.51s, native states
from authenticated v18.52l. Positive source clock comes from staged structural_jet.
"""
from __future__ import annotations
import argparse,csv,hashlib,importlib.util,inspect,io,json,math,re,sys,tempfile,traceback,zipfile
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from energy_math import Params,sigmoid
from clock_math import (field_plus,field_external_minus,relabelled_field,
    FLIP_THETA,FLIP_ALL_PHASES,general_conjugacy_residual,angular_accelerations,normres)

ROOT=Path(__file__).resolve().parents[1]; BASE=ROOT.parent; OUT=ROOT/'results'
C=json.loads((ROOT/'inputs_config.json').read_text());P=Params(**C['model_params'])
S51='q550_parameter_sensitivity_transport_transverse_posture_attribution_v18_51s_bundle'
R51='q550_long_horizon_paired_kappa0_transverse_transport_return_section_search_v18_51r_bundle'
L52='q550_near_tangency_breathing_envelope_drift_five_dimensional_section_audit_v18_52l_bundle'
DEFAULT_STATES=BASE/L52/'results/extended_horizon_10d_states.npz'
DEFAULT_PARENT=BASE/(S51+'_server_sendback.zip')
SEEDS=('original','autonomous')

def digest(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for chunk in iter(lambda:f.read(1<<20),b''):h.update(chunk)
 return h.hexdigest()

def dynamic_import(path,name):
 s=importlib.util.spec_from_file_location(name,path)
 if s is None or s.loader is None:raise RuntimeError('UNIMPORTABLE_SOURCE '+str(path))
 m=importlib.util.module_from_spec(s);sys.modules[name]=m;s.loader.exec_module(m);return m

def dump_json(path,obj):
 path.parent.mkdir(parents=True,exist_ok=True)
 path.write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False)+'\n')

def write_csv(path,rows):
 if not rows:return
 path.parent.mkdir(parents=True,exist_ok=True)
 with path.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def callable_metadata(f):
 info={'type':type(f).__module__+'.'+type(f).__qualname__,
       'module':getattr(f,'__module__',None),
       'name':getattr(f,'__qualname__',None)}
 if hasattr(f,'__code__'):
  code=f.__code__;info['source_file']=code.co_filename;info['first_line']=int(code.co_firstlineno)
  info['freevars']=list(code.co_freevars)
  closure=getattr(f,'__closure__',None) or []
  info['closure_types']=[type(c.cell_contents).__name__ for c in closure]
 try: info['source_excerpt']=inspect.getsource(f)[:18000]
 except (OSError,TypeError):info['source_excerpt']='UNAVAILABLE_FROM_INSPECT'
 return info

def locate_parent(source):
 search=([source] if source else [])+[DEFAULT_PARENT,ROOT/'inputs'/(S51+'_server_sendback.zip'),
           Path('/root/elab1/hysteresis')/(S51+'_server_sendback.zip'),
           Path('/root/Downloads')/(S51+'_server_sendback.zip')]
 return next((Path(p) for p in search if p is not None and Path(p).is_file()),None)

def source_excerpts(tree, limit=60):
 """Collect source, not fabricated narrative: exact lines from staged original files."""
 excerpt=[];index=[]
 pattern=re.compile(r'(?i)(def make_context|def .*theta|\btheta\s*=|\[.theta.\]|phase.clock|clock.direction|orientation|time.reverse|event\s*\+|CubicSpline|\btheta\(|\bTheta\b|structural_jet)')
 paths=sorted(tree.rglob('*.py'))
 ranked=[]
 for p in paths:
  try:
   data=p.read_text(encoding='utf-8',errors='replace');lines=data.splitlines()
  except OSError:continue
  hits=[i for i,line in enumerate(lines) if pattern.search(line)]
  if not hits:continue
  idx={'relative_path':str(p.relative_to(tree)),'sha256':digest(p),'bytes':p.stat().st_size,'matched_lines':len(hits)}
  index.append(idx)
  priority=0 if 'make_context' in data else 1 if 'theta' in data.lower() else 2
  ranked.append((priority,str(p),p,lines,hits))
 for _,_,p,lines,hits in sorted(ranked):
  if len(excerpt)>=limit:break
  # Present the most context-relevant relevant source excerpt, with line numbers.
  chosen=[]
  for i in hits:
   if any(abs(i-j)<=3 for j in chosen):continue
   chosen.append(i)
   if len(chosen)>=15:break
  selected=set()
  for i in chosen:
   selected.update(range(max(0,i-3),min(len(lines),i+4)))
  txt='\n'.join(f'{i+1:5d}: {lines[i]}' for i in sorted(selected))
  excerpt.append(f'### {p.relative_to(tree)} SHA256 {digest(p)}\n{txt}\n')
 return index,'\n'.join(excerpt)

def original_engine_genealogy(parent,xi,ref_states):
 if digest(parent)!=C['v1851s_sha256']:
  raise RuntimeError('PARENT_51S_SHA256_MISMATCH '+digest(parent))
 with zipfile.ZipFile(parent) as z:
  if z.testzip() is not None:raise RuntimeError('PARENT_51S_BAD_ZIP')
  nested=z.read(S51+'/inputs/'+R51+'_server_sendback.zip')
  srunsrc=z.read(S51+'/scripts/run_audit.py')
 with zipfile.ZipFile(io.BytesIO(nested)) as z:
  if z.testzip() is not None:raise RuntimeError('PARENT_51R_BAD_ZIP')
  a=z.read(R51+'/results/extended_paired_kappa_nonlinear_states.npz')
 with np.load(io.BytesIO(a),allow_pickle=False) as n:
  ks={seed:float(n['kappa_'+seed]) for seed in SEEDS}
 assert abs((ks['original']+ks['autonomous'])/2-C['v1852l_fixed_model_kappa'])<1e-12
 loader=dynamic_import(ROOT/'reference/parent_v18_51t_run_audit.py','q550_v1852u_staging_loader')
 with tempfile.TemporaryDirectory(prefix='q550_v1852u_original_') as temp:
  tmp=Path(temp)
  parts=loader.stage_engine(tmp,nested,srunsrc)[4]
  _,_,seed_api,observer,engine,unit,drpt,creport,_=parts
  got={name:float(getattr(engine.params,name)) for name in C['model_params']}
  diffs={name:got[name]-float(val) for name,val in C['model_params'].items()}
  if any(abs(d)>1e-12 for d in diffs.values()):
   dump_json(OUT/'parameter_mismatch.json',{'actual':got,'expected':C['model_params'],'diff':diffs})
   raise RuntimeError('AUTHENTIC_MODEL_PARAMETERS_DISAGREE see results/parameter_mismatch.json')
  stageindex,excerpts=source_excerpts(tmp)
  theta_source=callable_metadata(seed_api.make_context)
  (OUT/'parent_clock_source_excerpt.txt').write_text(
      '# Exact staged source matches (relative to temporary validated ancestor staging)\n'
      +excerpts+'\n\n# Exact make_context source excerpt\n'+theta_source.get('source_excerpt','')+'\n')
  phase={};metadata={};eng_rows=[]
  for seed in SEEDS:
   ctx=seed_api.make_context(observer,seed,engine,unit,drpt,creport)
   fn=ctx['theta'];event=float(ctx['event'])
   phase[seed]=np.asarray([float(fn(event+float(t))) for t in xi],float)
   md=callable_metadata(fn)
   metadata[seed]={'event':event,'theta_callable':md,
      'phase_at_first_xi':float(phase[seed][0]),'phase_at_last_xi':float(phase[seed][-1]),
      'context_keys':sorted(map(str,ctx.keys()))}
   # Independent original 11D generator at genuine native reference states,
   # plus arbitrary states off the orbit; original clock normalized to physical time.
   ks0=ks[seed]
   y=ref_states[seed]
   for i in np.unique(np.linspace(0,len(xi)-1,32,dtype=int)):
    z=np.r_[y[i,:10],phase[seed][i]]
    f=np.asarray(engine.jet.structural_jet(z,ks0,engine.params,engine.slope)[0],float)
    x=float(sigmoid(z[8]));f=f*x/(2*math.pi*float(engine.slope))
    archival=field_plus(z,ks0,P)
    eng_rows.append({'seed':seed,'kind':'native_reference','i':int(i),
      'theta_engine_normalized_rate':float(f[10]),'x':x,
      'physical_11d_relative_error':normres(f,archival),
      'original_source_clock_over_x':float(f[10]/x)})
   rng=np.random.default_rng(18525222+len(seed))
   for i in range(32):
    z=np.zeros(11);z[[0,4]]=rng.uniform(-1,1.5,2)
    z[[1,3,5,7,9]]=rng.uniform(-.8,.8,5)
    z[[2,6,10]]=rng.uniform(-5,5,3);z[8]=rng.uniform(-2,2)
    f=np.asarray(engine.jet.structural_jet(z,ks0,engine.params,engine.slope)[0],float)
    x=float(sigmoid(z[8]));f=f*x/(2*math.pi*float(engine.slope))
    eng_rows.append({'seed':seed,'kind':'off_orbit','i':i,
      'theta_engine_normalized_rate':float(f[10]),'x':x,
      'physical_11d_relative_error':normres(f,field_plus(z,ks0,P)),
      'original_source_clock_over_x':float(f[10]/x)})
  write_csv(OUT/'original_engine_11d_clock_probes.csv',eng_rows)
  dump_json(OUT/'source_genealogy.json',{'authenticated_parent':str(parent),
    'parent_sha256':digest(parent),'nested_v1851r_sha256':hashlib.sha256(nested).hexdigest(),
    'original_source_v1851s_script_sha256':hashlib.sha256(srunsrc).hexdigest(),
    'engine_slope':float(engine.slope),'original_engine_parameters':got,
    'staged_source_index':stageindex,'make_context':theta_source,'context_by_seed':metadata,
    'original_engine_probe_max_field_relative':max(r['physical_11d_relative_error'] for r in eng_rows),
    'original_engine_clock_ratio_range':[min(r['original_source_clock_over_x'] for r in eng_rows),max(r['original_source_clock_over_x'] for r in eng_rows)],
    'original_engine_probes':len(eng_rows)})
 return phase,ks,metadata,eng_rows

def evaluate_native(ref,phase,ks,xi,stride):
 force=[];trace=[];summary={}
 for seed in SEEDS:
  y=ref[seed];t=np.asarray(y[:,10],float);x=np.asarray(sigmoid(y[:,8]),float)
  ph=np.asarray(phase[seed],float)
  assert np.all(np.isfinite(ph)) and np.all(np.diff(t)>0)
  i=np.arange(0,len(t),max(1,stride),dtype=int)
  if i[-1]!=len(t)-1:i=np.r_[i,len(t)-1]
  yy=y[i];tt=t[i];pp=ph[i];xx=x[i]
  # The saved trajectory is evaluated only on its own physical time grid.
  domega=np.column_stack([np.gradient(yy[:,3],tt,edge_order=2),np.gradient(yy[:,7],tt,edge_order=2)])
  fit_minus=angular_accelerations(yy,ks[seed],P,pp,-1)
  fit_plus=angular_accelerations(yy,ks[seed],P,pp,+1)
  rms=lambda z:float(np.sqrt(np.mean(np.asarray(z,dtype=float)**2)))
  scale=max(1e-14,rms(domega))
  res_minus=rms(domega-fit_minus)/scale
  res_plus=rms(domega-fit_plus)/scale
  phase_integral=float(np.trapezoid(xx,tt))
  phase_delta=float(pp[-1]-pp[0])
  gradient=np.gradient(pp,tt,edge_order=2)
  ext_minus=rms(gradient+xx);ext_plus=rms(gradient-xx)
  summary[seed]={'n_native_samples':len(i),'theta_delta':phase_delta,
    'x_reference_integral':phase_integral,
    'signed_clock_endpoint_error':phase_delta+phase_integral,
    'external_clock_plus_x_rate_rms':ext_plus,'external_clock_minus_x_rate_rms':ext_minus,
    'minus_external_phase_angular_rms_relative':res_minus,
    'plus_external_phase_angular_rms_relative':res_plus,
    'reference_tau_start':float(tt[0]),'reference_tau_end':float(tt[-1])}
  force.append({'seed':seed,**summary[seed]})
  ss=np.unique(np.linspace(0,len(i)-1,min(301,len(i)),dtype=int))
  for j in ss:
   trace.append({'seed':seed,'sample_index':int(i[j]),'material_tau':float(tt[j]),
     'theta_external':float(pp[j]),'x_reference':float(xx[j]),
     'theta_rate_numerical':float(gradient[j]),
     'theta_rate_plus_x_error':float(gradient[j]-xx[j]),
     'theta_rate_minus_x_error':float(gradient[j]+xx[j])})
 write_csv(OUT/'clock_force_discrimination.csv',force)
 write_csv(OUT/'reference_clock_phase_traces.csv',trace)
 return summary

def full_field_maps(ks):
 rng=np.random.default_rng(18525223);rows=[]
 for seed in SEEDS:
  k=ks[seed]
  for i in range(128):
   y=np.zeros(11);y[[0,4]]=rng.uniform(-1,1.5,2)
   y[[1,3,5,7,9]]=rng.uniform(-.7,.7,5)
   y[[2,6,10]]=rng.uniform(-4,4,3);y[8]=rng.uniform(-2,2)
   for name,M in (('theta_sign_only',FLIP_THETA),('theta_local_and_omega_sign',FLIP_ALL_PHASES)):
    d,p,obs=general_conjugacy_residual(y,k,P,M)
    rows.append({'seed':seed,'test':'physical_full_field_conjugacy_'+name,'sample':i,
      'relative_residual':normres(p,obs),'max_component_abs':float(np.max(np.abs(d))),
      'dx_coordinate_residual':float(d[8]),'dvz_residual':float(d[9]),
      'omega1_residual':float(d[3]),'omega2_residual':float(d[7]),
      'theta_component_residual':float(d[10])})
   e=relabelled_field(y,k,P)
   rows.append({'seed':seed,'test':'passive_relabel_algebraic_identity','sample':i,
     'relative_residual':normres(e,FLIP_THETA*field_plus(FLIP_THETA*y,k,P)),
     'max_component_abs':0.,'dx_coordinate_residual':0.,'dvz_residual':0.,
     'omega1_residual':0.,'omega2_residual':0.,'theta_component_residual':0.})
 write_csv(OUT/'full_field_map_residuals.csv',rows)
 return {x:{'n':sum(r['test']==x for r in rows),
   'max_relative':max(r['relative_residual'] for r in rows if r['test']==x),
   'median_relative':float(np.median([r['relative_residual'] for r in rows if r['test']==x]))}
   for x in sorted(set(r['test'] for r in rows))}

def plot_reference_curves(ref,phase,ks,step=32):
 try:
  import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
 except ImportError:return []
 out=[]
 for seed in SEEDS:
  y=ref[seed];ix=np.arange(0,len(y),max(1,step));
  if ix[-1]!=len(y)-1:ix=np.r_[ix,len(y)-1]
  z=y[ix];t=z[:,10];ph=phase[seed][ix];x=sigmoid(z[:,8]);dph=np.gradient(ph,t,edge_order=2)
  d_omega=np.column_stack((np.gradient(z[:,3],t,edge_order=2),np.gradient(z[:,7],t,edge_order=2)))
  fm=angular_accelerations(z,ks[seed],P,ph,-1)
  fp=angular_accelerations(z,ks[seed],P,ph,+1)
  fig,axes=plt.subplots(3,1,figsize=(11,9),sharex=True,constrained_layout=True)
  axes[0].plot(t,dph,label='External clock rate')
  axes[0].plot(t,-x,'--',label='−x reference')
  axes[0].plot(t,x,':',label='+x reference');axes[0].legend(fontsize=8)
  axes[1].plot(t,d_omega[:,0],label='Native dω1/dτ')
  axes[1].plot(t,fm[:,0],'--',label='Force with −Θext')
  axes[1].plot(t,fp[:,0],':',label='Force with +Θext');axes[1].legend(fontsize=8)
  axes[2].plot(t,d_omega[:,1],label='Native dω2/dτ')
  axes[2].plot(t,fm[:,1],'--',label='Force with −Θext')
  axes[2].plot(t,fp[:,1],':',label='Force with +Θext');axes[2].legend(fontsize=8)
  for ax in axes:ax.grid(alpha=.2)
  axes[-1].set_xlabel('Material time τ')
  fig.suptitle(f'Q550 v18.52u: {seed} — phase genealogy, observed vs relabeled coupling')
  p=OUT/(seed+'_phase_orientation_force_comparison.png');fig.savefig(p,dpi=140);plt.close(fig)
  out.append(p.name)
 return out

def main():
 parser=argparse.ArgumentParser()
 parser.add_argument('--states',type=Path,default=DEFAULT_STATES)
 parser.add_argument('--parent-v1851s',type=Path)
 parser.add_argument('--stride',type=int,default=8,help='Native force finite-difference decimation stride')
 parser.add_argument('--require-parent',action='store_true')
 parser.add_argument('--require-native',action='store_true')
 args=parser.parse_args();OUT.mkdir(exist_ok=True)
 if args.stride<1:raise ValueError('--stride must be positive')
 parent=locate_parent(args.parent_v1851s)
 states=args.states;ref={}
 status={'version':'v18.52u','parent':'MISSING','native':'MISSING',
   'lineage':'NOT_STAGED','scientific_decision':'NOT_EVALUATED',
   'global_noncollapse_bound_certified':False,'full_nonlinear_conjugacy_certified':False}
 if states.is_file():
  h=digest(states)
  if h!=C['v1852l_10d_sha256']:raise RuntimeError('v18.52l NATIVE SHA256 MISMATCH '+h)
  with np.load(states,allow_pickle=False) as z:
   xi=np.asarray(z['common_xi'],float)
   assert len(xi)==65537 and np.all(np.diff(xi)>0)
   status['native_total_series']=len(z.files)-1
   for seed in SEEDS: ref[seed]=np.asarray(z[seed+'__reference__Y_tau'],float)
   status['native_kick_shared_phase_x_difference']={}
   for seed in SEEDS:
    r_x=sigmoid(ref[seed][:,8]);diff=[]
    for name in z.files:
     if name.startswith(seed) and '__kick__Y_tau' in name:
      X=np.asarray(z[name][:,8],float)
      diff.append(float(np.max(np.abs(sigmoid(X)-r_x))))
    status['native_kick_shared_phase_x_difference'][seed]=max(diff) if diff else None
  status['native']='AUTHENTIC_SHA256_VERIFIED'
  status['native_reference_xi_length']=len(xi)
 else:
  xi=np.linspace(0,1,17)
  if args.require_native:raise RuntimeError('AUTHENTIC_NATIVE_STATES_MISSING '+str(states))
 status['full_field_candidate_maps']=full_field_maps(C['kappa_seeds'])
 phase=None
 if parent is not None and parent.is_file():
  if not ref:
   status['parent']='AVAILABLE_NOT_STAGED_NATIVE_MISSING'
  else:
   phase,kappa,context,probes=original_engine_genealogy(parent,xi,ref)
   status['parent']='SHA256_AUTHENTICATED'
   status['lineage']='ORIGINAL_ENGINE_AND_THETA_FUNCTION_STAGED'
   status['authentic_engine_probes_max_relative']=max(row['physical_11d_relative_error'] for row in probes)
   status['original_engine_clock_ratio_range']=[min(row['original_source_clock_over_x'] for row in probes),max(row['original_source_clock_over_x'] for row in probes)]
   clock=evaluate_native(ref,phase,kappa,xi,args.stride)
   status['native_clock_force_discrimination']=clock
   status['plots']=plot_reference_curves(ref,phase,kappa)
 else:
  if args.require_parent:raise RuntimeError('AUTHENTIC_51S_PARENT_MISSING '+str(parent))
 summary=status.get('native_clock_force_discrimination')
 if summary:
  minus_clock=all(abs(s['signed_clock_endpoint_error'])<5e-7 and s['external_clock_minus_x_rate_rms']<2e-4 for s in summary.values())
  engine_plus=all(abs(x-1)<1e-10 for x in status['original_engine_clock_ratio_range'])
  mechanics_archive=all(s['minus_external_phase_angular_rms_relative']<1e-3 for s in summary.values())
  mechanics_relabel=all(s['plus_external_phase_angular_rms_relative']<1e-3 for s in summary.values())
  if minus_clock and engine_plus and mechanics_archive and not mechanics_relabel:
   decision='AUTHENTIC_GENEALOGY_WITH_EXTERNAL_DRIVE_DISTINCT_FROM_SIMPLE_COORDINATE_RELABEL'
  elif minus_clock and engine_plus and mechanics_relabel and not mechanics_archive:
   decision='AUTHENTIC_GENEALOGY_COMPATIBLE_WITH_PASSIVE_RELABEL'
  else:decision='AUTHENTIC_GENEALOGY_SIGN_DISCORDANCE_NOT_YET_CLASSIFIED'
  status['criteria']={'minus_clock':minus_clock,'engine_plus_clock':engine_plus,
    'archived_minus_phase_force_fit':mechanics_archive,'passive_plus_phase_force_fit':mechanics_relabel}
 elif parent is not None and parent.is_file():decision='AUTHENTIC_PARENT_BUT_NATIVE_TRACES_MISSING'
 else:decision='ANALYTIC_CANDIDATE_MAP_TESTS_ONLY_PARENT_MISSING'
 status['scientific_decision']=decision
 status['source_of_negative_external_clock_not_assumed']='Inspect results/parent_clock_source_excerpt.txt and results/source_genealogy.json before attributing the sign to a particular ancestor function.'
 status['no_claims']=['No proof of conjugacy by an arbitrary nonlinear transformation',
    'No automatic identification of external reference clock with per-kick x',
    'No topological integer winding or all-time global noncollapse theorem']
 dump_json(OUT/'report.json',status)
 lines=['# Q550 v18.52u — Phase-Clock Orientation Genealogy & Full-State Conjugacy Audit','',
    '## Scientific decision',decision,'',
    '## Exact passive coordinate identity','Theta_ext=-Theta_engine implies Theta_ext_dot=-x AND Delta=theta1+theta2-phi0+Theta_ext.',
    'This does not establish the actual archived external drive, which was evaluated using Delta=theta1+theta2-phi0-Theta_ext.','',
    '## Original ancestry','Original v18.51s must be SHA-256 checked and staged before its context theta is attributed to any source.',
    'Full-state map tests are limited to the explicitly listed linear transformations; arbitrary nonlinear conjugacy is untested.','',
    '## Results','```json',json.dumps(status,indent=2),'```','',
    '## Scientific restrictions','No new integrations. No global bound or topological charge concluded.',
    'The eleventh v18.52l column is material tau, never Theta.']
 (OUT/'SCIENTIFIC_REPORT_SERVER.md').write_text('\n'.join(lines)+'\n')
 print(json.dumps({'decision':decision,'parent':status['parent'],'native':status['native'],
   'criteria':status.get('criteria'),
   'native_clock_force_discrimination':status.get('native_clock_force_discrimination')},indent=2))
 if args.require_parent and status['parent']!='SHA256_AUTHENTICATED':raise SystemExit(3)
 if args.require_native and status['native']!='AUTHENTIC_SHA256_VERIFIED':raise SystemExit(3)

if __name__=='__main__':main()
