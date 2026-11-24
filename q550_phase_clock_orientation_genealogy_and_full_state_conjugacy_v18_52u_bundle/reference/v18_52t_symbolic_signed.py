#!/usr/bin/env python3
"""Optional independent SymPy verification with arbitrary external phase speed."""
import json
from pathlib import Path
try:
 import sympy as S
except ImportError:
 result={'status':'SYMPY_NOT_INSTALLED','identity':'tested independently by algebraic and numeric checks'}
else:
 r1,r2,rd1,rd2,omega1,omega2,xd,x,ddx,rate,k,mu,gammar,gammax,mx=S.symbols('r1 r2 rd1 rd2 omega1 omega2 xdot x xddot rate k mu gammar gammax mx',real=True)
 U1,U2,Vprime,c,s=S.symbols('U1 U2 Vprime cosD sinD',real=True)
 Ddot=omega1+omega2-rate
 # Use independent kinetic/potential power plus native equations.
 dd1=r1*omega1**2-U1-k*x*r2*c-gammar*rd1
 dd2=r2*omega2**2-U2-k*x*r1*c-gammar*rd2
 dw1=k*x*r2*s/r1-2*rd1*omega1/r1
 dw2=k*x*r1*s/r2-2*rd2*omega2/r2
 ddx=(-Vprime+mu*r1*r2*s-gammax*xd)/mx
 Eder=rd1*dd1+rd2*dd2+r1*rd1*omega1**2+r1**2*omega1*dw1+r2*rd2*omega2**2+r2**2*omega2*dw2+mx*xd*ddx+U1*rd1+U2*rd2+Vprime*xd+k*(xd*r1*r2*c+x*(rd1*r2+r1*rd2)*c-x*r1*r2*s*Ddot)
 expected=-gammar*(rd1**2+rd2**2)-gammax*xd**2+xd*r1*r2*(mu*s+k*c)+k*x*r1*r2*s*rate
 rE=S.expand(Eder-expected)
 Jsum=r1**2*omega1+r2**2*omega2
 dJsum=S.expand(2*r1*rd1*omega1+r1**2*dw1+2*r2*rd2*omega2+r2**2*dw2)
 rJ=S.expand(dJsum-2*k*x*r1*r2*s)
 sigma=S.symbols('sigma',real=True)
 Hder=S.expand(Eder-sigma*(xd*Jsum+x*dJsum)/2)
 ideal=-gammar*(rd1**2+rd2**2)-gammax*xd**2+xd*(r1*r2*(mu*s+k*c)-sigma*Jsum/2)
 rH=S.expand((Hder-ideal).subs(rate,sigma*x))
 result={'status':'SYMPY_EXACT_PASS' if rE==rJ==rH==0 else 'SYMPY_IDENTITY_FAIL',
         'residual_dE':str(rE),'residual_dJsum':str(rJ),'residual_dH_signed':str(rH),
         'scope':'arbitrary external phase rate; signed frame requires phase_rate=sigma*x'}
 out=Path(__file__).resolve().parents[1]/'results/symbolic_signed_identity.json';out.parent.mkdir(exist_ok=True)
 out.write_text(json.dumps(result,indent=2)+'\n')
 print('SYMBOLIC',result['status'])
 if result['status']=='SYMPY_IDENTITY_FAIL':raise SystemExit(2)
