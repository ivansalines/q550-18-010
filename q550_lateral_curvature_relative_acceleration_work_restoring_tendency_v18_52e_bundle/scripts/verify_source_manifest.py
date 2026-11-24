#!/usr/bin/env python3
"""Confirm every versioned source file before physical computation."""
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from run_audit import verify_sources
_,_,_,manifest=verify_sources()
print('Q550 v18.52e source verified',len(manifest['file_sha256']),'SHA-256')
