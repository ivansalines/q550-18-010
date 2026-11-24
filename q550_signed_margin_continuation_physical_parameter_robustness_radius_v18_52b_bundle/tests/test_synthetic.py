#!/usr/bin/env python3
"""Cheap synthetic sanity checks; no Q550 physical IVP is run here."""
from pathlib import Path
import importlib.util,sys
import numpy as np
p=Path(__file__).resolve().parents[1]/'scripts/run_audit.py'
spec=importlib.util.spec_from_file_location('q550v1852b_tested',p)
m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m)
xx=np.linspace(0,.2,513)
def fake(margin,total,area,selff):
 return {'xi':xx,'rate_mixed':np.full(len(xx),margin-selff),
         'rate_self_ab':np.full(len(xx),selff),
         'rate_total':np.full(len(xx),total),'oriented_area':np.full(len(xx),area)}
c={
 'stable_signed_margin':m.classify({},fake(.2,.199,.1,-.3))=='SAMPLED_PASS',
 'signed_margin_loss':m.classify({},fake(-.02,.05,.1,-.3))=='MARGIN_SIGN_LOSS',
 'full_turning_loss':m.classify({},fake(.2,-.03,.1,-.3))=='TOTAL_TURNING_SIGN_LOSS',
 'area_sign_loss':m.classify({},fake(.2,.2,-.03,-.3))=='ORIENTED_AREA_SIGN_LOSS',
 'self_sign_loss':m.classify({},fake(.2,.2,.1,.03))=='SELF_CHANNEL_SIGN_LOSS',
 'near_zero_unknown':m.classify({},fake(1e-8,1e-7,.1,-.3))=='NEAR_ZERO_UNRESOLVED',
 'infinity_unknown':m.classify({},fake(float('inf'),.2,.1,-.3))=='NUMERICAL_UNRESOLVED',
 'negative_radius_not_used':m.REQUIRED_NUMERIC_RADIUS==0,
 'parent_root_named':m.AROOT.startswith('q550_dynamically_admissible_'),
 'two_physical_seeds':len(m.SEEDS)==2,
 'two_velocity_source_supports_in_protocol':'5 and 7' in (p.parents[1]/'MATHEMATICAL_PROTOCOL.md').read_text(),
 'dynamic_code_not_precomputation':callable(m.inspect_real),
}
print('Q550 v18.52b synthetic',sum(c.values()),'/',len(c),'PASS',flush=True)
for k,v in c.items():
 if not v:print('FAIL',k,flush=True)
if not all(c.values()):sys.exit(2)
