#!/usr/bin/env python3
"""Strict checks of the unbalanced-torque control and source symmetry."""
import importlib.util,sys
from pathlib import Path
root=Path(__file__).resolve().parents[1]
p=root/'scripts/run_audit.py'
spec=importlib.util.spec_from_file_location('v1852q',p);m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m)
print('Q550 v18.52q synthetic tests')
result=m.archived_probes(100)
assert result['passes']
assert abs(result['counterfactual_break_dC'])>1e-6
assert abs(result['counterfactual_break_dC']-result['counterfactual_predicted'])<1e-12
print('PASS: exact archived equation, drive-equivariance, and deliberate symmetry breaking',result)
