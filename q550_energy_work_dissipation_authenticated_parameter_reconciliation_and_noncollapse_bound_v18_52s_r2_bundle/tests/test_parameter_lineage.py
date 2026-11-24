#!/usr/bin/env python3
"""Fail closed if the pre-repair parameters leak back into the audit."""
import json,sys
from pathlib import Path
R=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(R/'scripts'))
from energy_math import Params,a_priori_bound
c=json.loads((R/'inputs_config.json').read_text())
p=Params(**c['model_params'])
assert vars(p)==c['model_params'],(vars(p),c['model_params'])
assert p.e_reserve==0.01 and p.phi0==0 and p.mu_x==0.05 and p.barrier_B==0.3
assert c['supersedes_wrong_parameter_config']['e_reserve']==0.26
b=a_priori_bound(c['v1852l_fixed_model_kappa'],p)
assert b['status']=='FINITE_TIME_GRONWALL_BOUND'
assert 0.69 < b['K'] < 0.71,b
print('PARAMETER_LINEAGE_PRECHECK_PASS',c['model_params'],'K',b['K'])
