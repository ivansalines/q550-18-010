#!/usr/bin/env python3
"""Local deterministic unit tests of geometry, four-point regression and gates."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from run_audit import synthetic_tests
checks=synthetic_tests()
assert all(x=='PASS' for x in checks.values())
print('v18.52n SYNTHETIC TESTS PASS',checks)
