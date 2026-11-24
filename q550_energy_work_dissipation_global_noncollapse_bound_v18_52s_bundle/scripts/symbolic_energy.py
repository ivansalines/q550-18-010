#!/usr/bin/env python3
"""Exact symbolic energy ledger; optional SymPy, with documented analytic proof."""
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
report={'status':'SYMPY_NOT_INSTALLED','exact_model_identity':'see MATHEMATICAL_PROTOCOL.md'}
try:
 import sympy as S
 r1,r2,u1,u2,w1,w2,x,vx,d,k,mu,gr,gx,mx,e,B = S.symbols('r1 r2 u1 u2 w1 w2 x vx d k mu gr gx mx e B',nonzero=True,real=True)
 U=lambda r: r*r-S.Rational(4,5)*r**3+S.Rational(1,5)*r**4+e*e/(2*r*r)
 V=-B*S.log(x*(1-x));C=k*x*r1*r2*S.cos(d)
 E=S.Rational(1,2)*(u1*u1+u2*u2+r1*r1*w1*w1+r2*r2*w2*w2+mx*vx*vx)+U(r1)+U(r2)+V+C
 df={r1:u1,r2:u2,u1:r1*w1*w1-S.diff(U(r1),r1)-k*x*r2*S.cos(d)-gr*u1,
 u2:r2*w2*w2-S.diff(U(r2),r2)-k*x*r1*S.cos(d)-gr*u2,
 w1:k*x*r2*S.sin(d)/r1-2*u1*w1/r1,w2:k*x*r1*S.sin(d)/r2-2*u2*w2/r2,
 x:vx,vx:(-S.diff(V,x)+mu*r1*r2*S.sin(d)-gx*vx)/mx,
 d:w1+w2-x}
 derivative=S.expand(sum(S.diff(E,a)*b for a,b in df.items()))
 budget=-gr*(u1*u1+u2*u2)-gx*vx*vx+r1*r2*(mu*S.sin(d)+k*S.cos(d))*vx+k*x*x*r1*r2*S.sin(d)
 assert S.simplify(derivative-budget)==0
 Jsum=r1*r1*w1+r2*r2*w2
 djsum=S.simplify(sum(S.diff(Jsum,z)*df[z] for z in [r1,r2,w1,w2]))
 assert S.simplify(djsum-2*k*x*r1*r2*S.sin(d))==0
 H=E-x*Jsum/2
 dH=S.simplify(sum(S.diff(H,z)*df[z] for z in df))
 hbudget=-gr*(u1*u1+u2*u2)-gx*vx*vx+vx*(r1*r2*(mu*S.sin(d)+k*S.cos(d))-Jsum/2)
 assert S.simplify(dH-hbudget)==0
 report={'status':'SYMPY_EXACT_PASS','energy_residual':str(S.simplify(derivative-budget)),
   'rotating_frame_residual':str(S.simplify(dH-hbudget)),
   'individual_angular_moment_transfer':str(djsum)}
except ImportError:pass
(ROOT/'results').mkdir(exist_ok=True)
(ROOT/'results/symbolic_energy_identity.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
