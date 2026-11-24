#!/usr/bin/env python3
"""Independent exact algebra certificate using sympy when it is installed.

The proof itself is also spelled out in MATHEMATICAL_PROTOCOL.md and requires
no computer algebra package or numerical orbit data.
"""
import json
from pathlib import Path
p=Path(__file__).resolve().parents[1]/'results'/'symbolic_identity.json'
p.parent.mkdir(parents=True,exist_ok=True)
try:
 import sympy as S
except ImportError:
 p.write_text(json.dumps(dict(status='SYMPY_NOT_AVAILABLE',proof_document='MATHEMATICAL_PROTOCOL.md'),indent=2)+'\n')
 print('Sympy optional; formal derivation preserved in protocol')
else:
 r1,r2,v1,v2,w1,w2,k,x,Delta=S.symbols('r1 r2 v1 v2 omega1 omega2 k x Delta', real=True,positive=False)
 Fw1=k*x*r2/r1*S.sin(Delta)-2*v1*w1
 Fw2=k*x*r1/r2*S.sin(Delta)-2*v2*w2
 J1dot=S.expand(2*r1*r1*v1*w1+r1*r1*Fw1)
 J2dot=S.expand(2*r2*r2*v2*w2+r2*r2*Fw2)
 c=S.simplify(J1dot-J2dot)
 assert c==0 and S.simplify(J1dot-k*x*r1*r2*S.sin(Delta))==0
 eps=S.symbols('epsilon', real=True)
 delta_shift=(S.Symbol('theta1')+eps)+(S.Symbol('theta2')-eps)-S.Symbol('phi0')-S.Symbol('Theta')
 assert S.simplify(delta_shift-(S.Symbol('theta1')+S.Symbol('theta2')-S.Symbol('phi0')-S.Symbol('Theta')))==0
 result=dict(status='EXACT_SYMBOLIC_PASS',J1dot=str(J1dot),J2dot=str(J2dot),dC=str(c),counterrotation_delta_invariant=True,sympy=S.__version__)
 p.write_text(json.dumps(result,indent=2)+'\n')
 print(result)
