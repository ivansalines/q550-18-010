#!/usr/bin/env python3
"""Regression for the v18.52i h-seam orientation CSV name mismatch."""
from pathlib import Path
import importlib.util
import numpy as np
import sys
root=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('v1852i_ledger',root/'scripts/run_audit.py')
module=importlib.util.module_from_spec(spec)
sys.modules[spec.name]=module
spec.loader.exec_module(module)
sample={'angle_from_seam_deg': np.array([0.0,118.21]),
        'D':np.array([2.0,1.0]),'G':np.array([-0.7,0.5])}
assert module.ledger_value(sample,'angle_from_h_seam_deg',1)==118.21
assert module.ledger_value(sample,'D',0)==2.0
assert module.ledger_value(sample,'G',1)==0.5
assert 'angle_from_h_seam_deg' not in sample
print('LEDGER H-SEAM ANGLE ALIAS REGRESSION PASS')
