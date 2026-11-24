#!/usr/bin/env python3
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from run_audit import verify_source_manifest
_,rep,rows,_,hashes=verify_source_manifest()
print('Q550 v18.52c: source manifest verified; B physical rows',len(rows),'B exact result hashes',len(hashes),flush=True)
