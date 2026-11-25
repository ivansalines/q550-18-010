#!/usr/bin/env python3
"""Seven short post-processing checks on authenticated archived data; no model IVPs."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import numpy as np
from run_audit import authenticated_o
from valley_convergence_math import Primitive,optimum_width,eval_curve
x,rows,_,_,_,_=authenticated_o()
old=[r for r in rows if r['seed']=='original' and r['kappa_rule']=='shared']
old=sorted(old,key=lambda r:int(r['index']))
o0=float(old[70]['origin_xi']);w0=float(old[70]['width_xi'])
for mode in ['linear','pchip']:
 p=Primitive(x['common_xi'],x['original__shared__primitive'],mode)
 origins=np.linspace(o0-0.025,o0+0.025,7)
 w=np.array([optimum_width(p,o,w0-.05,w0+.05)[0] for o in origins])
 result,curve=eval_curve(p,origins,w)
 assert all(np.isfinite(v) for v in [result['arc_geodesic'],result['arc_chain_gradient'],result['center_norm']])
 assert np.max(np.abs(np.linalg.norm(curve['u'],axis=1)-1))<1e-12
print('V1853P_AUTHENTIC_SMALL_SAMPLE_POSTPROCESS_PASS_NO_IVPS')
