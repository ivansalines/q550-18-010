from pathlib import Path
import importlib.util,sys
import numpy as np
base=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(base/'scripts'))
import run_audit as n
from full_net_geometry import make_eval,local_values
x=np.linspace(12.8,56,2701)
C=np.zeros((len(x),10))
C[:,0]=.3*np.sin(0.63*x)
C[:,1]=.16*np.cos(.73*x)
C[:,2]=.01*(x-12.8)
C[:,3]=.05*np.sin(.17*x)
C[:,4]=.013*np.sin(.8*x)
C[:,5]=.008*(x-12.8)
C[:,6]=.07*np.sin(.37*x)
C[:,7]=.014*np.cos(.55*x)
C[:,8]=.001*x
C[:,9]=.013*np.sin(.41*x)
at,der=make_eval(x,C)
rows=[]
for a in [13.6,14.4,15.2]:
 for o,w in [(a,12.8),(a+.025,12.805)]:
  _,_,_,net=local_values((o,w),at,der)
  rows.append(dict(seed='original',kappa_rule='shared',solved_projected_root='True',anchor_origin_xi=a,origin_xi=o,width_xi=w,full_bivector_norm=float(np.linalg.norm(net))))
class ParentMock:
 def primitive_from_acc(self,acc):return None,C
results=dict(m=ParentMock(),grid=[],multistart=[],minima=[],best=[],parent_baseline=[])
B=np.gradient(C,x,axis=0)
primitive,best,base_rows=n.analyze_case('original','shared',x,None,B,np.ones(len(x)),rows,results)
assert len(results['grid']) == 53*43
assert results['best'] and len(results['parent_baseline']) == 6
assert np.isfinite(best['material_source_relative_error'])
assert np.isfinite(best['pchip_reoptimized_net_norm'])
assert np.isfinite(best['gain_over_parent_projected_root'])
print('SYNTHETIC_END_TO_END_ANALYZE_CASE PASS', 'grid',len(results['grid']), 'runs',len(results['multistart']),'best_norm',round(best['net_norm'],9),'source_error',round(best['material_source_relative_error'],6))
