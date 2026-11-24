#!/usr/bin/env python3
"""Check immutable inputs, including original v18.52c report/CRC."""
from pathlib import Path
import sys
p=Path(__file__).resolve().parents[1]/'scripts/run_audit.py'
spec=__import__('importlib.util',fromlist=['x']).spec_from_file_location('audit52d_verify',p)
m=__import__('importlib.util',fromlist=['x']).module_from_spec(spec)
sys.modules[spec.name]=m;spec.loader.exec_module(m)
m.verify_sources(False)
print('Q550 v18.52d immutable inputs and v18.52c prior report verified')
