#!/usr/bin/env python3
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from run_audit import synthetic_tests
assert all(v=='PASS' for v in synthetic_tests().values())
print('SYNTHETIC TESTS PASS',synthetic_tests())
